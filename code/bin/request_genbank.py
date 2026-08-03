#!/usr/bin/env python

import pandas as pd
import xml.etree.ElementTree as ET
import requests
import time
import argparse

parser = argparse.ArgumentParser(
    description="script to download taxonomy and yabby family from genbank"
)
parser.add_argument("uniprot_file", help="tsv file as output of uniprot DB"
)
args = parser.parse_args()

uniprot = args.uniprot_file

file = pd.read_csv(uniprot, sep="\t", header=None, names=["fasta", "accession", "score", "specie"])
data = ['AccessionVersion', 'TaxId', 'Title']
all_info = []

def get_summary(acc):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary?db=protein&id='+acc
    # creating HTTP response object from given url
    resp = requests.get(url)
    summary = ET.fromstring(resp.content)
    info = {}
    for item in list(summary)[0].findall('Item'):
        if item.attrib["Name"] in data:
            info[item.attrib["Name"]] = item.text
    if len(info) == 0:
        return
    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch?db=taxonomy&id="+str(info["TaxId"])
    resp2 = requests.get(url2)
    time.sleep(0.01)
    lineage2 = ET.fromstring(resp2.content)
    taxonomy = list(lineage2.iter("Lineage"))[0].text
    scientific_name = list(lineage2.iter("ScientificName"))[0].text
    info['lineage'] = taxonomy
    info['specie'] = scientific_name.replace(" ", "_")
    all_info.append([info[i] for i in ['AccessionVersion', 'specie','Title', 'lineage']])

for acc in file["accession"]:
    try:
        get_summary(acc)
    except IndexError:
        pass


lineage_df = pd.DataFrame(all_info, columns=["accession", "specie", "protein", "lineage"])

lineage_df.merge(file, on=["accession", "specie"]).drop_duplicates().to_csv("blast_genbank_lineage_species.uniq.tsv", sep="\t", index=None)
