# MXene/metal interface builder: Flask app + MACE-MH-1 on CPU.
#   podman build -t mxbuilder:0.1.0 .
#   podman run --rm -p 8000:8000 -v "$PWD/runs:/app/runs:Z" mxbuilder:0.1.0
FROM docker.io/library/python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    XDG_CACHE_HOME=/opt/cache \
    MXB_RUNS=/app/runs

WORKDIR /app

RUN pip install --index-url https://download.pytorch.org/whl/cpu torch==2.14.1
COPY requirements.txt .
RUN pip install -r requirements.txt && pip freeze > /app/pip-freeze.txt

# The model is part of the image, so runs never depend on a download.
RUN python -c "from mace.calculators.foundations_models import download_mace_mp_checkpoint as get; print(get('mh-1'))"

COPY mxbuilder mxbuilder
COPY webapp webapp
COPY tests tests
COPY examples examples

EXPOSE 8000
CMD ["gunicorn", "--workers", "1", "--threads", "4", "--timeout", "300", "--bind", "0.0.0.0:8000", "webapp.app:app"]
