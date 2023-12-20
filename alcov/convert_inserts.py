# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 2023

@author: Jenn Knapp
email: jknapp@uwaterloo.ca
"""
"""
Purpose: Turns the insert.bed file containing all primer pairs and positions covered into an amplicon.py file that is read by the amplicon_coverage function of alcov.

Reguires:
SARS-CoV-2.insert.bed downloaded from https://github.com/BioWilko/primer-schemes/blob/master/nCoV-2019/V4.1/SARS-CoV-2.insert.bed of the appropriate version.
"""

input_file = "SARS-CoV-2.insert.bed"
output_file = "amplicons.py"

inserts = []

with open(input_file, 'r') as f:
    for line in f:
        data = line.strip().split('\t')
        insert = [data[0], data[1], data[2], data[4], data[5], data[3][-1]]
        inserts.append(insert)

with open(output_file, 'w') as f:
    f.write(f"inserts = {str(inserts)}\n")
