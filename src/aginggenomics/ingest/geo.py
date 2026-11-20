"""Bounded acquisition of pinned, public GEO artifacts."""

import json
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

from ..reporting import sha256


@dataclass(frozen=True)
class GeoFile:
    url: str
    digest: str
    max_bytes: int = 25 * 1024 * 1024

    @property
    def filename(self):
        return self.url.rsplit("/", 1)[-1]


def acquire(file: GeoFile, cache: Path, *, offline=False):
    """Return verified bytes and the recorded acquisition manifest."""
    cache.mkdir(parents=True, exist_ok=True)
    destination = cache / file.filename
    manifest = destination.with_suffix(".json")
    if destination.is_file() and manifest.is_file():
        try:
            recorded = json.loads(manifest.read_text(encoding="utf-8"))
            valid = (recorded["url"] == file.url and recorded["sha256"] == file.digest
                     and bool(recorded["downloaded_at_utc"])
                     and destination.stat().st_size == recorded["bytes"] <= file.max_bytes
                     and sha256(destination) == file.digest)
            if valid:
                return destination, recorded
        except (OSError, ValueError, TypeError, KeyError):
            pass
    if offline:
        raise OSError(f"No verified cache for {file.filename}. Run again with downloads enabled.")

    for attempt in range(3):
        try:
            with tempfile.TemporaryDirectory(dir=cache) as staging:
                pending = Path(staging) / file.filename
                size = 0
                deadline = time.monotonic() + 90
                with urlopen(file.url, timeout=15) as response, pending.open("wb") as output:
                    while chunk := response.read1(256 * 1024):
                        size += len(chunk)
                        if size > file.max_bytes or time.monotonic() > deadline:
                            raise OSError("Download exceeded the size or transfer deadline.")
                        output.write(chunk)
                if sha256(pending) != file.digest:
                    raise OSError("SHA-256 differs from the pinned GEO artifact.")
                recorded = {"url": file.url, "sha256": file.digest, "bytes": size,
                            "downloaded_at_utc": datetime.now(timezone.utc).isoformat()}
                pending_manifest = Path(staging) / manifest.name
                pending_manifest.write_text(json.dumps(recorded, indent=2), encoding="utf-8")
                pending.replace(destination)
                pending_manifest.replace(manifest)
            return destination, recorded
        except OSError as exc:
            if attempt == 2:
                raise OSError(f"Unable to download {file.filename}: {exc}") from exc
            time.sleep(2 ** attempt)
