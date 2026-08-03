# RepurposeIQ (CipherQ) — Flask app served by gunicorn.
# Python 3.12 to match the app's verified .venv312 build.
FROM python:3.12-slim

WORKDIR /app
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1

# RDKit / scientific wheels need a few shared libs at runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 libxext6 libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install -r requirements.txt gunicorn

COPY . .

EXPOSE 5002
CMD ["gunicorn", "-b", "0.0.0.0:5002", "-w", "1", "--threads", "4", \
     "--timeout", "180", "flask_app:app"]
