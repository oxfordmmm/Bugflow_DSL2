#!/usr/bin/env python
#tsv2html
# Copyright (c) 2019, Gaussian Solutions LLC. All rights reserved
# 
"""
Converts a given tsv file to html table file
"""
import csv
import os
import sys
import pandas as pd
def tsv2html(tsv_file_name, html_file_name):
    df = pd.read_csv(tsv_file_name,sep='\t', header=0)
    old_width = pd.get_option('display.max_colwidth')
    #pd.set_option('display.max_colwidth', -1)    
    with open(html_file_name,'w') as html_file:
        html_file.write(df.to_html(index=False))
    pd.set_option('display.max_colwidth', old_width)
if __name__ == "__main__":
   tsv2html(sys.argv[1], sys.argv[2])