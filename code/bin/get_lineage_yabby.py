#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="script to extract taxonomy and yabby family to use for tree annotation"
)
parser.add_argument("uniprot_file", help="tsv file as output of uniprot DB"
)

parser.add_argument("genbank_file", help="tsv file as output of genbank DB"
)

parser.add_argument("taxonomy_file", help="file with list of taxonomic division to use"
)
args = parser.parse_args()

uniprot = args.uniprot_file
genbank = args.genbank_file
taxonomy = args.taxonomy_file

uniprot_DF = pd.read_csv(uniprot, sep="\t", 
            header=None, names=["uniprot", "fasta", "score", "specie", "alphafold", "protein", "lineage", "pdb"])

genbank_DF = pd.read_csv(genbank, sep="\t")


lineage_file = open(taxonomy, "r")
taxonomy = [i.strip() for i in lineage_file]

yabby_dict = {"FIL" : ["YABBY 1", "YABBY1", "YABBY 1-like", "YABBY 3", "YABBY3", "YABBY 3-like"],
 "YAB2" : ["YABBY 2", "YABBY2", "YABBY 2-like"],
 "INO" : ["YABBY 4", "YABBY4", "YABBY 4-like", "INNER NO OUTER"],
 "YAB5" : ["YABBY 5", "YABBY5", "YABBY 5-like"],
 "CRC" : ["DROOPING LEAF", "CRABS CLAW", "crabs claw"]
}

### uniprot
lineage = {spe:set() for spe in set(uniprot_DF["specie"])}
for tax in taxonomy:
    for specie, lin in zip(uniprot_DF["specie"], uniprot_DF["lineage"]):
        if tax.lower() in lin.lower():
            lineage[specie].add(tax)

lineage_list = []
for i in uniprot_DF.specie:
    lineage_list.append(list(lineage[i])[-1].lower())

yabby = []
for i in uniprot_DF.protein:
    annot = set()
    for yab in yabby_dict:
        for name in yabby_dict[yab]:
            if name in i:
                annot.add(yab)
    yabby.append("".join(annot))

uniprot_DF["classification"] = lineage_list
uniprot_DF["yabby"] = yabby

### genbank
lineage = {spe:set() for spe in set(genbank_DF["specie"])}
for tax in taxonomy:
    for specie, lin in zip(genbank_DF["specie"], genbank_DF["lineage"]):
        if tax.lower() in lin.lower():
            lineage[specie].add(tax)

lineage_list = []
for i in genbank_DF.specie:
    lineage_list.append(list(lineage[i])[-1].lower())

yabby = []
for i in genbank_DF.protein:
    annot = set()
    for yab in yabby_dict:
        for name in yabby_dict[yab]:
            if name in i:
                annot.add(yab)
    yabby.append("".join(annot))

genbank_DF["classification"] = lineage_list
genbank_DF["yabby"] = yabby



subset_uniprot_df = uniprot_DF[["fasta", "specie", "yabby", "classification"]]
subset_genbank_df = genbank_DF[["fasta", "specie", "yabby", "classification"]]

for f in subset_genbank_df[subset_genbank_df.yabby.isna()].fasta: 
    if len(subset_uniprot_df.loc[subset_uniprot_df["fasta"] == f, "yabby"].values) != 0:
        subset_genbank_df.loc[subset_genbank_df["fasta"] == f, "yabby"] = subset_uniprot_df.loc[subset_uniprot_df["fasta"] == f, "yabby"].values[0]

itol_df = pd.concat([subset_genbank_df, subset_uniprot_df[~subset_uniprot_df.fasta.isin(subset_genbank_df.fasta)]])

itol_df.to_csv("annotation_lineage.itol.tsv", sep="\t", index=None)

phylo_names = itol_df[(itol_df.yabby.isin(list(yabby_dict.keys()))) | (itol_df.classification.isin(["acrogymnospermae", "amborella", "aristolochia", "nymphaea"]))].fasta

with open("phylo_names.txt", 'w') as wf:
    wf.write("\n".join(phylo_names))
    