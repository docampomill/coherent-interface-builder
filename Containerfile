# Coherent MXene/metal interface builder + MACE and Allegro relaxation on CPU
#   podman build -t mxbuilder:2.0.4 .
#   podman run --rm -p 8000:8000 mxbuilder:2.0.4
#   podman run --rm -p 8000:8000 -v "$PWD/data:/app/data:Z" mxbuilder:2.0.4   # use your edited spreadsheet
# The Allegro model is not in git: put the handoff in models/allegro/ before building (see README).
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
# Allegro in its own venv: nequip needs e3nn 0.6, and MACE's model loads only with 0.4.4. The venv reuses the image's
# torch and numpy (system site packages); relax.py talks to it through allegro_worker.py.
COPY requirements-allegro.txt .
RUN python -m venv --system-site-packages /opt/allegro && /opt/allegro/bin/pip install -r requirements-allegro.txt

# The model is part of the image, so relaxations never depend on a download.
RUN python -c "from mace.calculators.foundations_models import download_mace_mp_checkpoint as get; print(get('mh-1'))"

COPY . .

EXPOSE 8000
CMD ["gunicorn", "--workers", "1", "--threads", "4", "--timeout", "300", "--bind", "0.0.0.0:8000", "app:app"]
