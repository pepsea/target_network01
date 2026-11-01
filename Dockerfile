# target_net: code + Python deps only. Databases, config and outputs are mounted (see docker-compose.yml).
FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends curl git && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
ENV TARGETNET_CONFIG=/config/target_net.yaml PYTHONUNBUFFERED=1
# default: CLI help; compose overrides the command
ENTRYPOINT ["python", "target_net.py"]
CMD ["--help"]
