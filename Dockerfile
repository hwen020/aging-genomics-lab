FROM python:3.12-slim
ENV AGING_GENOMICS_MODE=portfolio OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN python -m pip install --no-cache-dir .
COPY app.py ./
COPY .streamlit/config.toml ./.streamlit/config.toml
COPY data/demo ./data/demo
RUN useradd --create-home --uid 10001 portfolio
USER portfolio
EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s CMD python -c "import os, urllib.request; urllib.request.urlopen('http://127.0.0.1:' + os.getenv('PORT', '8501') + '/_stcore/health', timeout=3).close()"
CMD ["sh", "-c", "exec python -m streamlit run app.py --server.address=0.0.0.0 --server.port=${PORT:-8501}"]
