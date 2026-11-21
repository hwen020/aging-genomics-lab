# Portfolio Demo deployment

The public container uses `AGING_GENOMICS_MODE=portfolio`, starts on the platform's `PORT` (8501 locally) and serves `/_stcore/health`. The build context is allowlisted to the app, package and synthetic fixture. Public studies render bundled aggregate JSON. The full research workflow belongs on a local machine with access to the pinned GEO data.

```sh
docker build -t aging-genomics-lab:0.3.1 .
docker run --rm -p 8501:8501 aging-genomics-lab:0.3.1
curl -fsS http://localhost:8501/_stcore/health
```

`render.yaml` is a free-tier Blueprint for use after a repository is connected. As checked on 2026-11-20, Render lists 512 MB on the free web tier, spins down idle instances after 15 minutes, estimates about one minute to wake and uses ephemeral local storage. The local Windows Portfolio Demo server peaked at about 216 MiB in a single-browser walkthrough. The measured desktop EWAS previously peaked near 601 MiB, so the free tier is a candidate for the compact demo only. The local Docker engine was unavailable for validation. Linux container capacity, concurrent sessions and hosted cold starts still need measurement. [Render free services](https://render.com/docs/free) · [Render pricing](https://render.com/pricing) · [Streamlit container guide](https://docs.streamlit.io/deploy/tutorials/docker)
