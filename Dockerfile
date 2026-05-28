# Container image for the Streamlit UI. Used by Railway / Render / any
# Docker host. Local CLI usage (python -m pipelines.runner) still works
# the same way outside the container.

FROM python:3.11-slim

# System deps RDKit needs at runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 libxext6 libsm6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python deps first so source changes don't bust the layer cache.
COPY requirements.txt .
# CPU-only torch from the PyTorch index (smaller image than the default
# CUDA build, which is ~2 GB). Other deps from PyPI.
RUN pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu \
    -r requirements.txt

COPY . .

# Railway / Render inject $PORT; default to 8501 for local docker run.
ENV PORT=8501
EXPOSE 8501

# --server.address=0.0.0.0 so the host can reach us; --server.headless=true
# so streamlit doesn't try to open a browser.
CMD streamlit run app.py \
    --server.port=${PORT} \
    --server.address=0.0.0.0 \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=true
