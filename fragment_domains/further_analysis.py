# David 6/17/2020

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pybedtools
import re
from scipy import stats
import seaborn as sns

# Loading the ANNOTATED GENOME dataset as a dataframe
# I added an empty column to better format it when I convert the DataFrame into a BedTool, it contains no information
annotations = pybedtools.BedTool('data/saccharomyces_cerevisiae.bed').to_dataframe()
annotations = annotations.drop(['score', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'], axis=1)
annotations['empty'] = np.zeros(annotations['chrom'].size)
column_titles = ['chrom','start','end','name','empty','strand']
annotations = annotations.reindex(columns=column_titles)
annotations_bed = pybedtools.BedTool.from_dataframe(annotations).sort()
annotations.head()

# Loading the FRAGMENT dataset as a BedTool
fragments = pybedtools.BedTool('../primary_analysis/work-cached/pacbio-190731-facs-assign.bed').sort()
fragments_df = fragments.to_dataframe()
fragments_df.head()

# Separating the annotations dataset into two, one for each strand
genes_plus = annotations['strand'] == '+'
genes_minus = annotations['strand'] == '-'
annotations_plus = annotations[genes_plus]
annotations_minus = annotations[genes_minus]
annotations_plus_bed = pybedtools.BedTool.from_dataframe(annotations_plus)
annotations_minus_bed = pybedtools.BedTool.from_dataframe(annotations_minus)

# Seperating the fragments dataset into two, one for each strand
frags_plus = fragments_df['strand'] == '+'
frags_minus = fragments_df['strand'] == '-'
fragments_plus = fragments_df[frags_plus]
fragments_minus = fragments_df[frags_minus]
fragments_plus_bed = pybedtools.BedTool.from_dataframe(fragments_plus)
fragments_minus_bed = pybedtools.BedTool.from_dataframe(fragments_minus)

# Finding the coverage of the positive strand fragments over positive strand genes
plus_covg = annotations_plus_bed.coverage(fragments_plus_bed, s=True, nonamecheck=True).to_dataframe()
plus_covg.columns = ['chrom', 'start', 'end', 'name', 'empty', 'strand', 
                    'frag coverage count', 'bases covered count', 'length', 'frag coverage fraction']

np.count_nonzero(plus_covg['frag coverage count']) / 3387

# Finding the coverage of the minus strand fragments over minus strand genes
minus_covg = annotations_minus_bed.coverage(fragments_minus_bed, s=True, nonamecheck=True).to_dataframe()
minus_covg.columns = ['chrom', 'start', 'end', 'name', 'empty', 'strand', 
                    'frag coverage count', 'bases covered count', 'length', 'frag coverage fraction']
np.count_nonzero(minus_covg['frag coverage count']) / 3309

# Finding the coverage of the pacbio dataset on the yeast proteome. Each row represents a gene.
frag_covg = annotations_bed.coverage(fragments, s=True, nonamecheck=True).to_dataframe()
frag_covg.columns = ['chrom', 'start', 'end', 'name', 'empty', 'strand', 
                    'frag coverage count', 'bases covered count', 'length', 'frag coverage fraction']
#frag_covg.dropna(axis=0, inplace=True)
print(frag_covg.shape)
frag_covg.head()

# Making a dataframe of the genes that are covered by at least one fragment, removing entries with 0 coverage
# to make a clearer scatter plot below.
at_least_one = frag_covg['frag coverage count'] != 0
genes_covered = frag_covg[at_least_one]
genes_covered.head()

# The fraction of genes represented by at least one fragment
np.count_nonzero(frag_covg['frag coverage count']) / frag_covg['frag coverage count'].size 

# Saving the fragment coverage dataset as a csv
frag_covg.to_csv('output/fragment_coverage.csv')

# A scatter plot showing gene length vs. fragment coverage
plt.scatter(genes_covered['length'], genes_covered['frag coverage fraction'], s=0.5)

plt.hist(frag_covg['frag coverage count'], bins=np.arange(0, 100, 2))

# Calculating the z scores on fragment coverage count for each gene to find outliers that are represented in
# the screen unusually highly. I used a threshold of a z-score of 3.
covg_z = stats.zscore(frag_covg['frag coverage count'])
frag_covg['z'] = covg_z
outliers = frag_covg[frag_covg['z'] > 3]
print('number of genes overrepresented:', outliers.shape[0])
outliers.head()

### Fragments in-frame

# Intersecting the fragments dataset with the annotated genome 
frag_annotations = fragments.intersect(annotations_bed, wb=True, s=True, nonamecheck=True).to_dataframe().drop(['thickStart','blockSizes'], axis=1)
frag_annotations.columns = ['chrom', 'fragStart', 'fragEnd', 'barcode', 'score', 'strand', 'geneStart', 'geneEnd', 'name', 'geneStrand']
print(frag_annotations.shape)
frag_annotations.head()

frag_annotations_num_rows = frag_annotations.shape[0]

frag_ant_plus = frag_annotations[frag_annotations['strand'] == '+']
frag_ant_minus = frag_annotations[frag_annotations['strand'] == '-']

frag_ant_plus['reading frame'] = np.abs((frag_ant_plus['fragStart'] - frag_ant_plus['geneStart']) % 3)
print(frag_ant_plus.shape)
frag_ant_plus.head()

frag_ant_minus['reading frame'] = np.abs(((frag_ant_minus['fragEnd']) - frag_ant_minus['geneEnd']) % 3)
print(frag_ant_minus.shape)
frag_ant_minus.head()

frag_plus_inframe = frag_ant_plus[frag_ant_plus['fragEnd'] < frag_ant_plus['geneEnd']]
print(frag_plus_inframe.shape)
frag_plus_inframe.head()

frag_plus_inframe = frag_plus_inframe[frag_plus_inframe['reading frame'] == 2]
frag_plus_inframe.shape[0] / frag_ant_plus.shape[0]

frag_minus_inframe = frag_ant_minus[frag_ant_minus['fragStart'] > frag_ant_minus['geneStart']]
print(frag_minus_inframe.shape)
frag_minus_inframe.head()

frag_minus_inframe = frag_minus_inframe[frag_minus_inframe['reading frame'] == 1]
frag_minus_inframe.shape[0] / frag_ant_minus.shape[0]

# A function that takes in rows from the frag_annotations dataset and returns the coordinate for the protein start of the fragment
def findProteinStart(row):
    if row['strand'] == '+':
        end = (row['fragEnd'] - row['geneStart'] + 1) // 3
        dist = (row['fragEnd'] - row['fragStart']) // 3
        return end - dist
    elif row['strand'] == '-':
        dist = (row['geneEnd'] - row['fragStart'] + 1) // 3
        tot_dist = (row['geneEnd'] - row['geneStart']) // 3
        return tot_dist - dist
    
# A function that takes in rows from the frag_annotations dataset after finding the proteinStarts and returns the coordinate for the 
# protein end of the fragment
def findProteinEnd(row):
    return row['proteinStart'] + ((row['fragEnd'] - row['fragStart']) // 3)

# Applying the above functions to append the proteinStart and proteinEnd indeces to each row of frag_annotations
frag_annotations['proteinStart'] = frag_annotations.apply(findProteinStart, axis=1)
frag_annotations['proteinEnd'] = frag_annotations.apply(findProteinEnd, axis=1)
frag_annotations.head()

# Finding the proportion of fragments that are in-frame.
1 - (np.count_nonzero(frag_annotations['geneStart'] + 3 * frag_annotations['proteinStart'] == frag_annotations['fragStart']) / frag_annotations_num_rows)

1 - (np.count_nonzero((frag_annotations['fragStart'] - frag_annotations['geneStart']) % 3 == 0) / frag_annotations['fragStart'].size)

frag_plus_inframe.to_csv('output/plus_inframe_frags.csv')
frag_minus_inframe.to_csv('output/minus_inframe_frags.csv')


from Bio.SeqIO.FastaIO import SimpleFastaParser

### Fragment-domains, domain presence in each fragment

# Loading the DOMAIN dataset into a dataframe
domains = pd.read_csv('data/domains.csv', names=['seq id', 'alignment start', 'alignment end', 'envelope start',
                               'envelope end', 'hmm acc', 'hmm name', 'type', 'hmm start', 'hmm end',
                               'hmm length', 'bit score', 'E-value', 'clan'])

# Loading table with SEQ ID and corresponding GENE NAMES (ORDERED LOCUS)
domain_loci = pd.read_csv('data/uniprot-filtered-proteome_UP000002311+AND+organism__Saccharomyces+cerevisi--.tab', sep='\t')
# Merging the domain dataset with domain_loci to add gene locus names onto the domain dataset; dropping uninformative columns
domains = pd.merge(domains, domain_loci, left_on='seq id', right_on='Entry').drop(['Entry', 'type', 'clan', 'envelope start', 'envelope end',
                                                                                  'hmm start', 'hmm end', 'hmm length', 'bit score', 'E-value'], axis=1)
new_cols = domains.columns.values
new_cols[5] = 'Gene locus'
domains.columns = new_cols
domains.head()

# Merging the domains dataset with the annotations dataset to align domains to their genomic coordinates
domains_genomic = pd.merge(domains, annotations, left_on='Gene locus', right_on='name').drop('name', axis=1)
domain_start = (domains_genomic['alignment start'] * 3) + domains_genomic['start']
domain_end = (domains_genomic['alignment end'] * 3) + domains_genomic['start']
domains_genomic['domain start'] = domain_start
domains_genomic['domain end'] = domain_end
domains_genomic.head()

# Making a dataset in BED format of domains and their genomic coordinates, to use for interescting with fragments later
domain_bed = pd.DataFrame()
domain_bed['chrom'] = domains_genomic['chrom']
domain_bed['domStart'] = domains_genomic['domain start']
domain_bed['domEnd'] = domains_genomic['domain end']
domain_bed = pybedtools.BedTool.from_dataframe(domain_bed)
domain_bed.head()

# Intersecting the fragments and domain datasets to find the presence of domains in each fragment
frag_domains = fragments.intersect(domain_bed, wo=True, nonamecheck=True).to_dataframe().drop(['thickStart'], axis=1)
frag_domains.columns = ['chrom','start','end','barcode','score','strand','domStart','domEnd','overlap']
print(frag_domains.shape)
frag_domains.head()

# Merging the fragment_domains dataframe with the domains dataset to add more information about the domains
frag_domain_presence = pd.merge(frag_domains, domains_genomic, left_on=['chrom', 'domStart', 'domEnd'], right_on=['chrom', 'domain start', 'domain end'])
frag_domain_presence = frag_domain_presence.drop(['start_y','end_y','strand_y', 'domain start', 'domain end','Gene locus'], axis=1)

frag_domain_presence.columns = ['chrom','fragStart','fragEnd','barcode','score','strand','domStart','domEnd','overlap',
                                'seq id','alignment start','alignment end','hmm acc','hmm name', '']
print(frag_domain_presence.shape)
frag_domain_presence.head()

domain_presence = frag_domain_presence['overlap'] / (frag_domain_presence['domEnd'] - frag_domain_presence['domStart'])
frag_domain_presence['domain presence'] = domain_presence
frag_domain_presence.head()



frag_domain_presence.to_csv('output/frag-domain-presence.csv')

