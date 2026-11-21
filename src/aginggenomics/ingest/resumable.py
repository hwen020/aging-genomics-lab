"""Resumable, bounded HTTPS transfers for large, pinned GEO artifacts."""

import json
import os
import re
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from urllib.request import Request, urlopen

from ..reporting import sha256


@contextmanager
def transfer_lock(path):
    """Use a process-held lock that releases when an interrupted process exits."""
    with path.open("a+b") as handle:
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            lock = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            unlock = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            lock = lambda: fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            unlock = lambda: fcntl.flock(handle, fcntl.LOCK_UN)
        try:
            lock()
        except OSError as exc:
            raise OSError(f"Another process is acquiring {path.stem}. Retry after it finishes.") from exc
        try:
            yield
        finally:
            unlock()


def transfer_range(url, partial, size, source_size, *, seconds=1200):
    offset = partial.stat().st_size if partial.exists() else 0
    end = size - 1
    request = Request(url, headers={"Range": f"bytes={offset}-{end}", "Accept-Encoding": "identity"})
    deadline = time.monotonic() + seconds
    with urlopen(request, timeout=25) as response:
        content_range = response.headers.get("Content-Range", "")
        match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
        if response.status == 206:
            if not match or tuple(map(int, match.groups())) != (offset, end, source_size):
                raise OSError("Unexpected HTTP Content-Range. Partial bytes were preserved.")
        elif response.status == 200 and size == source_size:
            # Some proxies ignore Range. Restart the full transfer only.
            offset = 0
        else:
            raise OSError("Server did not honor the requested bounded byte range.")
        if response.headers.get("Content-Length") not in (None, str(size - offset)):
            raise OSError("Unexpected HTTP Content-Length.")
        with partial.open("ab" if offset else "wb") as output:
            remaining = size - offset
            while remaining:
                if time.monotonic() > deadline:
                    raise OSError("Transfer deadline exceeded; partial bytes were preserved.")
                chunk = response.read1(min(1024 * 1024, remaining))
                if not chunk:
                    raise OSError("Incomplete transfer; partial bytes were preserved.")
                output.write(chunk)
                remaining -= len(chunk)


def acquire_large(artifact, cache, size, *, source_size=None, offline=False):
    """Publish a pinned full artifact or a pinned prefix only after verification."""
    if not 0 < size <= artifact.max_bytes:
        raise ValueError("Artifact length exceeds the configured transfer limit.")
    source_size = source_size or size
    cache.mkdir(parents=True, exist_ok=True)
    destination = cache / artifact.filename
    manifest = destination.with_name(destination.name + ".json")
    partial = destination.with_name(destination.name + ".part")
    expected = {"url": artifact.url, "sha256": artifact.digest, "bytes": size,
                "source_bytes": source_size}
    if size != source_size:
        expected["http_range"] = f"bytes=0-{size - 1}"
    with transfer_lock(destination.with_name(destination.name + ".lock")):
        try:
            record = json.loads(manifest.read_text(encoding="utf-8"))
            if (all(record.get(key) == value for key, value in expected.items())
                    and record.get("downloaded_at_utc") and destination.stat().st_size == size
                    and sha256(destination) == artifact.digest):
                return destination, record
        except (OSError, ValueError, AttributeError):
            pass
        if offline:
            raise OSError(f"No verified cache for {artifact.filename}. Enable downloads to repair it.")
        for attempt in range(3):
            try:
                if partial.exists() and partial.stat().st_size > size:
                    partial.unlink()
                if not partial.exists() or partial.stat().st_size < size:
                    transfer_range(artifact.url, partial, size, source_size)
                if sha256(partial) != artifact.digest:
                    partial.unlink()
                    raise OSError("Download SHA-256 differs from the pinned artifact.")
                record = {**expected, "downloaded_at_utc": datetime.now(timezone.utc).isoformat()}
                pending = manifest.with_suffix(".tmp")
                pending.write_text(json.dumps(record, indent=2), encoding="utf-8")
                partial.replace(destination)
                pending.replace(manifest)
                return destination, record
            except OSError as exc:
                if attempt == 2:
                    raise OSError(f"{artifact.filename}: {exc}") from exc
                time.sleep(2 ** attempt)
