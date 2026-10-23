# Coherent MXene/metal interface builder + MACE relaxation on CPU
#   podman build -t coherent-interface-builder:2.1.0 .
#   podman run --rm -p 8000:8000 coherent-interface-builder:2.1.0
#   podman run --rm -p 8000:8000 -v "$PWD/data:/app/data:Z" coherent-interface-builder:2.1.0   # use your edited spreadsheet
FROM docker.io/library/python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    XDG_CACHE_HOME=/opt/cache

WORKDIR /app
RUN pip install --index-url https://download.pytorch.org/whl/cpu torch==2.14.1
COPY requirements.txt .
RUN pip install -r requirements.txt

# The model is part of the image, so relaxations never depend on a download.
RUN python -c "from mace.calculators.foundations_models import download_mace_mp_checkpoint as get; print(get('mh-1'))"

COPY . .

EXPOSE 8000
CMD ["gunicorn", "--workers", "1", "--threads", "4", "--timeout", "300", "--bind", "0.0.0.0:8000", "app:app"]
