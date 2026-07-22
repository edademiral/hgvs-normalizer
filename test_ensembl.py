import requests

url = ("https://rest.ensembl.org/sequence/region/human/11:45917701..45917701"
       "?content-type=text/plain;coord_system_version=GRCh38")
try:
    r = requests.get(url, timeout=20)
    print("HTTP:", r.status_code)
    print("Referans baz:", r.text.strip())
    print("Bizim ref  : T")
    print("ESLESIYOR" if r.text.strip().upper() == "T" else "ESLESMIYOR")
except Exception as e:
    print("HATA:", type(e).__name__, str(e)[:150])
