# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 21:50:29 2020

@author: antho
"""

import re, gzip

anncomma = re.compile(b',(?!\()')

# def read_variants(filepath):
    
transcript_variant = {}
with gzip.open(filepath, 'rb') as file:
    for line in file:
        if line.startswith(b"#"): continue
        infos = line.split(b'\t')[7].split(b';')
        annotation = [info for info in infos if info.startswith(b"ANN=")]
        if len(annotation) != 1: continue
        annotations = anncomma.split(annotation[0].split(b'ANN=')[1])
        for ann in annotations:
            allele, efftypes, putative_impact, geneName, geneId, featureType, featureId, biotype, exonIntronRank, hgvsDna, hgvsProtein, cdnaPosition, cdsPosition, protPos, distToFeature, warnings = ann.split(b'|')
            if putative_impact in [b"MODERATE", b"HIGH"] and featureId.startswith(b'transcript:ENST') and hgvsProtein: # skip splice acceptor/donator variations by requiring hgvsProtein
                transcriptId = featureId.strip(b'transcript:')
                if transcriptId in transcript_variant: transcript_variant[transcriptId].append(ann)
                else: transcript_variant[transcriptId] = [ann]

effs=[]
for ann in np.concatenate([v for v in transcript_variant.values()]):
    allele, efftypes, putative_impact, geneName, geneId, featureType, featureId, biotype, exonIntronRank, hgvsDna, hgvsProtein, cdnaPosition, cdsPosition, protPos, distToFeature, warnings = ann.split(b'|')
    effs.append(efftypes)
print(np.unique(effs))