#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="script to extract taxonomy and yabby family to use for tree annotation from genbank file"
)
parser.add_argument("uniprot_file", help="tsv file as output of uniprot DB"
)

parser.add_argument("taxonomy_file", help="file with list of taxonomic division to use"
)
args = parser.parse_args()

uniprot = args.uniprot_file
taxonomy = args.taxonomy_file

uniprot_DF = pd.read_csv(uniprot, sep="\t")

lineage_file = open(taxonomy, "r")
taxonomy = [i.strip() for i in lineage_file]


lineage = {spe:set() for spe in set(uniprot_DF["specie"])}
for tax in taxonomy:
    for specie, lin in zip(uniprot_DF["specie"], uniprot_DF["lineage"]):
        if tax.lower() in lin.lower():
            lineage[specie].add(tax)

lineage_list = []
for i in uniprot_DF.specie:
    lineage_list.append(list(lineage[i])[-1].lower())

uniprot_DF["classification"] = lineage_list

yabby_dict = {"FIL" : ["YABBY 1", "YABBY1", "YABBY 1-like", "YABBY 3", "YABBY3", "YABBY 3-like"],
 "YAB2" : ["YABBY 2", "YABBY2", "YABBY 2-like"],
 "INO" : ["YABBY 4", "YABBY4", "YABBY 4-like", "INNER NO OUTER"],
 "YAB5" : ["YABBY 5", "YABBY5", "YABBY 5-like"],
 "CRC" : ["DROOPING LEAF", "CRABS CLAW", "crabs claw"]
}

yabby = []
for i in uniprot_DF.protein:
    annot = set()
    for yab in yabby_dict:
        for name in yabby_dict[yab]:
            if name in i:
                annot.add(yab)
    yabby.append("".join(annot))

uniprot_DF["yabby"] = yabby

itol_df = uniprot_DF[["fasta", "specie", "yabby", "classification"]]
itol_df.to_csv("uniprot_annotation_lineage.itol.tsv", sep="\t", index=None)

annot_names = itol_df[itol_df.yabby.isin(list(yabby_dict.keys()))].fasta

with open("CRC_names.txt", 'w') as wf:
    wf.write("\n".join(annot_names))
    