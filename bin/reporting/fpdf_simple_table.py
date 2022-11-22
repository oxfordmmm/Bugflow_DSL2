# simple_table.py
 
from fpdf import FPDF
import json
#import pandas as pd
import numpy as np


def simple_table(spacing=1):

#Sample name
    with open('sample_name.tsv') as h:
        sample = h.readlines()
    name = ("").join(sample)
    print(name)

#Read in QC stats in TSV
    with open('QC_summary_table.tsv') as g:
        lines = g.readlines()
    tsv = ("").join(lines)
    print(tsv)
   
#Read in QC stats in JSON
    with open('amr.json') as f:
        data = json.load(f) 
    print(data)

#Append logos
    WIDTH = 210
    #HEIGHT = 297
    pdf = FPDF()
    pdf.add_page()
    pdf.image('oxford.png', x=0, y=0, w=WIDTH/2-10)
    pdf.image('leeds.jpg', x=120, y=WIDTH/9.5, w=80)
    pdf.ln(30)  # move 30 down

#Show poject title
    pdf.set_font("Times", "B", size=17)
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Sequence Analysis Report", align = "L", ln=2)
    pdf.ln(5)

#Show sample details
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,255) #blue
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Sample details", align = "L", ln=2)
    pdf.ln(10)

    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = name)
    pdf.ln(10)

#Append QC stats
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,255) #blue
    pdf.cell(w=0, h=5, txt="Sequencing Quality Stats", align = "L", ln=2)
    
    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = tsv)
    pdf.ln(10)

#Append AMR Profile
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,255) #blue
    pdf.cell(w=0, h=5, txt="Antimicrobial Resistance Profile", align = "L", ln=2)
    pdf.set_font("Times", size=14)
    pdf.ln(10)

    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    col_width = pdf.w
    row_height = pdf.font_size
    for key, value in data.items():
        pdf.multi_cell(col_width, row_height*spacing,
                txt="{}\t\t\t{}".format(key, value))
    pdf.ln(10)
    
#Save to PDF
    
    pdf.output('simple_table.pdf')
 
if __name__ == '__main__':
    simple_table()

