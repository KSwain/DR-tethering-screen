#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 12:57:34 2019

@author: davidnoble
"""

import pandas as pd
import pybedtools 

domains = pd.read_table('559292.tsv', header=3, sep='\s',
                        names=['seq id', 'alignment start', 'alignment end', 'envelope start',
                               'envelope end', 'hmm acc', 'hmm name', 'type', 'hmm start', 'hmm end',
                               'hmm length', 'bit score', 'E-value', 'clan'])

annotations_bed = pybedtools.BedTool('saccharomyces_cerevisiae.bed')
annotations_df = annotations_bed.to_dataframe()

#Using UniProt data to get ordered locus gene name from seq id
uniprot = pd.read_table('uniprot-filtered-proteome_UP000002311+AND+organism__Saccharomyces+cerevisi--.tab')

domains = pd.merge(domains, uniprot, left_on='seq id', right_on='Entry')
domains = domains.drop('Entry', axis=1)
new_cols = domains.columns.values
new_cols[14] = 'Gene locus'
domains.columns = new_cols

domains_merged = pd.merge(domains, annotations_df, left_on='Gene locus', right_on='name')
domain_start = (domains_merged['alignment start'] * 3) + domains_merged['start']
domain_end = (domains_merged['alignment end'] * 3) + domains_merged['start']

domains_merged['domain start'] = domain_start
domains_merged['domain end'] = domain_end

domains_table = pd.DataFrame()
domains_table['chrom'] = domains_merged['chrom']
domains_table['chromStart'] = domains_merged['domain start']
domains_table['chromEnd'] = domains_merged['domain end']
domains_table['name'] = domains_merged['seq id']
domains_table['strand'] = domains_merged['strand']

domains_genomic = pybedtools.BedTool.from_dataframe(domains_table).saveas('domains_genomic.bed')


x = .75
sequences = sequences = pybedtools.BedTool('pacbio-190731-facs-assign.bed')

intersections = sequences.intersect(domains_genomic, f=x, wo=True, nonamecheck=True).to_dataframe()
intersections.columns = ['chrom', 'start', 'end', 'barcode', 'score', 'strand', 'domainChrom', 'domainStart', 'domainEnd', 
                         'domainID', 'strand', 'overlap']



fragment_domains = pd.merge(intersections, domains_merged, left_on=['domainID', 'domainStart'], 
                            right_on=['seq id', 'domain start'])
fragment_domains.to_csv(r'fragment_domains.csv')
