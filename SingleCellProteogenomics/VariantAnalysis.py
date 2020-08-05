# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 21:50:29 2020

@author: antho
"""

import re
anncomma = re.compile(b',(?!\()')
def read_variants(filepath):
    gene_transcript = {}
    transcript_variant = {}
    lines=[]
    with gzip.open(filepath, 'rb') as file:
        for line in file:
            if line.startswith(b"#"): continue
            if b'MODERATE' in line: lines.append(line)
            annotations = anncomma.split(line.split(b'\t')[7].split(b'ANN=')[1])
            for ann in annotations:
                fields = ann.split(b'|')
                eff = fields[2]
                gene = fields[4]
                site = 
                transcript = fields[6]
                if eff in [b"MODERATE",b"HIGH"] and gene.startswith(b'ENSG'):
                    if gene in gene_variant: gene_variant[gene].append(ann)
                    else: gene_variant[gene] = []
    gene_variant = {}
