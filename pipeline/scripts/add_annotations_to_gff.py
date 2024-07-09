#!/usr/bin/env python

# script to add annotations from DRAM to gff file
# 2022-11-23

import sys

annot_fl = sys.argv[1]
gff_fl = sys.argv[2]
output = sys.argv[3]

# first collect annotation info in a dictionary

# Annotations file doesn't contain barnap hits because they are rRNA hits and are ignored
# This is causing the error in the comparison stage between the annot_dict and the gff_fl file 
# as there are different numbers of lines (I think this is the proximate cause but it may just be
# a symptom.) I don't think that the rRNA is actually wanted considering we are doing viral discovery here
# So easiest solution is to just build out the append lines section to include a grep for rRNA?
annot_dict = {}
target_column_name = "pfam_hits"  # Specify the desired column name here
column_index = None

with open(annot_fl) as fl:
    header = next(fl).strip().split("\t")
    if target_column_name in header:
        column_index = header.index(target_column_name)
    else:
        print("Column name not found: " + target_column_name)
        column_index = -1
    
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        gene_id = cols[0]
        product = cols[column_index].split("[")[0].lstrip("HerpesvirusABC").lstrip("AlphaherpesvirusABC").strip()
        if product == "":
            product = "Hypothetical protein"
        annot_dict[gene_id] = product

# now read through gff file and add this information

output_lines = []

with open(gff_fl) as fl:
    for line in open(gff_fl):
        if line.startswith("#"):
            output_lines.append(line)
            continue
        if line.find("ribosomal RNA") != -1:
            continue
        line = line.strip()
        cols = line.split("\t")
        IDs = cols[8].split(";")
        gene = IDs[0].lstrip("ID=")
        # make my new strings to add to gff file
        prod = "product=" + annot_dict[gene]
        name = "Name=" + annot_dict[gene] + " CDS"
        IDs.insert(0, name)
        IDs.append(prod)
        IDs_str = ";".join(IDs).replace(";;", ";")
        # now put everything back together
        first_cols = cols[0:8]
        first_cols.append(IDs_str)
        first_cols.append("\n")
        CDS_line = "\t".join(first_cols).replace("\t\n", "\n")
        gene_line = CDS_line.replace("CDS", "gene")
        output_lines.append(gene_line)
        output_lines.append(CDS_line)
        

# go through the output_lines and write to new file

with open(output, "w") as out:
    for ln in output_lines:
        out.write(ln)




