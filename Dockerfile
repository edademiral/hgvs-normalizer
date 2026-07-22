FROM python:3.10-slim

WORKDIR /app

# psycopg2 (UTA baglantisi) icin derleme bagimliliklari
RUN apt-get update && apt-get install -y \
    libpq-dev \
    gcc \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# ENTRYPOINT sabit, CMD varsayilan argumanlar
ENTRYPOINT ["python", "hgvs_normalizer.py"]
CMD ["--input", "variants.tsv"]
