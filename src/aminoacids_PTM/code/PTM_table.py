import csv
import io
import re
import sys
import tarfile
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from urllib.parse import urljoin
from urllib3.util.retry import Retry

BASE_URL = "https://biomics.lab.nycu.edu.tw/dbPTM/"
DOWNLOAD_PAGE = urljoin(BASE_URL, "download.php")
AMINO_ACIDS = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic acid",
    "E": "Glutamic acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine"
}


def make_session():
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"})
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.headers.update({"User-Agent": "dbPTM-human-matrix/1.0"})
    return session


def get_download_links(session):
    response = session.get(DOWNLOAD_PAGE, timeout=60)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    links = []
    for link in soup.select("a[href]"):
        url = urljoin(DOWNLOAD_PAGE, link["href"])
        if "/dbPTM/download/experiment/" in url and url.lower().endswith(".gz"):
            links.append(url)
    links = sorted(set(links))
    if not links:
        raise RuntimeError("Не найдены файлы экспериментальных PTM")
    return links


def process_archive(content, data):
    with tarfile.open(fileobj=io.BytesIO(content), mode="r:gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            with io.TextIOWrapper(extracted, encoding="utf-8-sig", errors="replace", newline="") as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if len(row) < 6 or not row[0].strip().endswith("_HUMAN"):
                        continue
                    uniprot_id = row[1].strip()
                    position = row[2].strip()
                    ptm = row[3].strip()
                    sequence = row[5].strip()
                    if not uniprot_id or not position or not ptm or len(sequence) <= 10:
                        continue
                    amino_acid = sequence[10].upper()
                    if amino_acid not in AMINO_ACIDS:
                        continue
                    pmids = {
                        value.strip()
                        for value in row[4].split(";")
                        if re.fullmatch(r"\d+", value.strip())
                    }
                    if not pmids:
                        continue
                    data[ptm]["pmids"].update(pmids)
                    data[ptm]["pair_pmids"][amino_acid].update(pmids)
                    data[ptm]["sites"][amino_acid].add((uniprot_id, position))


def build_tables(data):
    selected = sorted(ptm for ptm, values in data.items() if len(values["pmids"]) >= 2)
    matrix = pd.DataFrame(0, index=selected, columns=list(AMINO_ACIDS.values()), dtype=int)
    matrix.index.name = "PTM"
    summary_rows = []
    detail_rows = []
    for ptm in selected:
        values = data[ptm]
        summary_rows.append({
            "PTM": ptm,
            "Publication_count": len(values["pmids"]),
            "PMIDs": ";".join(sorted(values["pmids"], key=int))
        })
        for code, name in AMINO_ACIDS.items():
            sites = values["sites"][code]
            pair_pmids = values["pair_pmids"][code]
            matrix.loc[ptm, name] = len(sites)
            if sites:
                detail_rows.append({
                    "PTM": ptm,
                    "Amino_acid": name,
                    "Site_count": len(sites),
                    "Publication_count": len(pair_pmids),
                    "PMIDs": ";".join(sorted(pair_pmids, key=int))
                })
    summary = pd.DataFrame(summary_rows)
    details = pd.DataFrame(detail_rows)
    return matrix, summary, details


def main():
    output_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("dbptm_human_ptm_amino_acids.xlsx")
    if output_path.suffix.lower() != ".xlsx":
        raise ValueError("Имя выходного файла должно заканчиваться на .xlsx")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    session = make_session()
    links = get_download_links(session)
    data = defaultdict(lambda: {
        "pmids": set(),
        "pair_pmids": defaultdict(set),
        "sites": defaultdict(set)
    })
    for index, url in enumerate(links, start=1):
        print(f"[{index}/{len(links)}] {url.rsplit('/', 1)[-1]}", flush=True)
        response = session.get(url, timeout=(30, 300))
        response.raise_for_status()
        process_archive(response.content, data)
        time.sleep(0.2)
    matrix, summary, details = build_tables(data)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        matrix.to_excel(writer, sheet_name="PTM_by_amino_acid")
        summary.to_excel(writer, sheet_name="PTM_publications", index=False)
        details.to_excel(writer, sheet_name="Pair_details", index=False)
    print(f"Сохранено: {output_path.resolve()}")
    print(f"Модификаций с минимум двумя PMID: {len(matrix)}")


if __name__ == "__main__":
    main()