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

vcf = pd.read_csv(vcfinput, comment='#', sep='\t', header=0)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('expand_frame_repr', False)
if impact:
    if impact=="HIGH":
        filtered_highimpact=vcf[vcf.iloc[:,7].astype(str).str.contains('HIGH')]
        print(filtered_highimpact)
    if impact=="MODERATE":
        filtered_moderateimpact=vcf[vcf.iloc[:,7].astype(str).str.contains('MODERATE')]
        print(filtered_moderateimpact)
    if impact=="LOW":
        filtered_lowimpact=vcf[vcf.iloc[:,7].astype(str).str.contains('LOW')]
        print(filtered_lowimpact)
    else:
        print("Invalid impact")
if quality:
    qualityscore=int(quality)
    filtered_quality = vcf[vcf.iloc[:,5].ge(qualityscore)]
    print(filtered_quality)
if denovo==True:
    denovo_rows=vcf[vcf.iloc[:,9].astype(str).str.contains('0/1' or '1/1') & vcf.iloc[:,10].astype(str).str.contains('0/0') & vcf.iloc[:,11].astype(str).str.contains('0/0')]
    print(denovo_rows)
        
