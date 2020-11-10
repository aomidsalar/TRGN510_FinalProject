%%file mutationcalling.py
#!/usr/bin/python
import sys
from optparse import OptionParser
import fileinput
import re
import numpy as np
import pandas as pd
import csv

parser = OptionParser()
parser.add_option("--impact=", dest="impact")
parser.add_option("--quality", type="int", dest="quality")
parser.add_option("--denovo", action="store_true", dest="denovo", default=False)
parser.add_option("-f", "--file", dest="filename")
(options, args) = parser.parse_args()
impact=options.impact
quality=options.quality
denovo=options.denovo
vcfinput=options.filename

def read_comments(csv_file):
    for row in csv_file:
        if row[0] == '#':
            yield row.split('#')[1].strip()

def get_last_commented_line(filename):
    with open(filename, 'r', newline='') as f:
        decommented_lines = [line for line in csv.reader(read_comments(f))]
        header = decommented_lines[-1]
        header=header[0].split("\t")
        skiprows = len(decommented_lines)
        return header, skiprows

header, skiprows = get_last_commented_line(vcfinput)
vcf=pd.read_csv(vcfinput, sep="\t", names=header, skiprows=skiprows)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('expand_frame_repr', False)
if impact:
    if impact=="HIGH":
        vcf=vcf[vcf.iloc[:,7].astype(str).str.contains('HIGH')]
    if impact=="MODERATE":
        vcf=vcf[vcf.iloc[:,7].astype(str).str.contains('MODERATE')]
    if impact=="LOW":
        vcf=vcf[vcf.iloc[:,7].astype(str).str.contains('LOW')]
if quality:
    qualityscore=int(quality)
    vcf = vcf[vcf.iloc[:,5].ge(qualityscore)]
if denovo==True:
    vcf=vcf[vcf.iloc[:,9].astype(str).str.contains('0/1' or '1/1') & vcf.iloc[:,10].astype(str).str.contains('0/0') & vcf.iloc[:,11].astype(str).str.contains('0/0')]
vcf.to_csv(sys.stdout, sep="\t", quoting=csv.QUOTE_NONE, index=False)
