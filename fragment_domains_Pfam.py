
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
import pybedtools
import re
from scipy import stats


# In[6]:


# Loading the DOMAIN dataset into a dataframe
domains = pd.read_csv('domains.csv', names=['seq id', 'alignment start', 'alignment end', 'envelope start',
                               'envelope end', 'hmm acc', 'hmm name', 'type', 'hmm start', 'hmm end',
                               'hmm length', 'bit score', 'E-value', 'clan'])

# Loading table with SEQ ID and corresponding GENE NAMES (ORDERED LOCUS)
domain_loci = pd.read_csv('uniprot-filtered-proteome_UP000002311+AND+organism__Saccharomyces+cerevisi--.tab', sep='\t')
# Merging the domain dataset with domain_loci to add gene locus names onto the domain dataset; dropping uninformative columns
domains = pd.merge(domains, domain_loci, left_on='seq id', right_on='Entry').drop(['Entry', 'type', 'clan', 'envelope start', 'envelope end',
                                                                                  'hmm start', 'hmm end', 'hmm length', 'bit score', 'E-value'], axis=1)
new_cols = domains.columns.values
new_cols[5] = 'Gene locus'
domains.columns = new_cols
domains.head()


# In[7]:


# Loading the ANNOTATED GENOME dataset as a dataframe
annotations = pybedtools.BedTool('saccharomyces_cerevisiae.bed').to_dataframe()
annotations = annotations.drop(['score', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'], axis=1)
annotations.head()


# In[8]:


# Merging the domains dataset with the annotations dataset to align domains to their genomic coordinates
domains_genomic = pd.merge(domains, annotations, left_on='Gene locus', right_on='name').drop('name', axis=1)
domain_start = (domains_genomic['alignment start'] * 3) + domains_genomic['start']
domain_end = (domains_genomic['alignment end'] * 3) + domains_genomic['start']
domains_genomic['domain start'] = domain_start
domains_genomic['domain end'] = domain_end
domains_genomic.head()


# In[9]:


# Making a dataset in BED format of domains and their genomic coordinates, to use for interescting with fragments later
domain_bed = pd.DataFrame()
domain_bed['chrom'] = domains_genomic['chrom']
domain_bed['domStart'] = domains_genomic['domain start']
domain_bed['domEnd'] = domains_genomic['domain end']
domain_bed = pybedtools.BedTool.from_dataframe(domain_bed)
domain_bed.head()


# In[10]:


# Loading the FRAGMENT dataset as a BedTool
fragments = pybedtools.BedTool('pacbio-190731-facs-assign.bed')
fragments.head()


# In[11]:


# Let x be an adjustable parameter to determine the amount of overlap we require a fragment to have with a domain to consider them as a match
x = .75


# In[12]:


# Using the intersect tool to find overlap between fragments and domains, returning the result as a dataframe
fragment_domains = fragments.intersect(domain_bed, f=x, wo=True, nonamecheck=True).to_dataframe().drop('thickStart', axis=1)
fragment_domains.columns = ['chrom', 'fragStart', 'fragEnd', 'barcode', 'score', 'strand', 'domStart', 'domEnd', 'overlap']
fragment_domains.head()


# In[13]:


# Merging the fragment_domains dataframe with the domains dataset to add more information about the domains
fragment_domains = pd.merge(fragment_domains, domains_genomic, left_on=['chrom', 'domStart', 'domEnd'], right_on=['chrom', 'domain start', 'domain end'])
fragment_domains = fragment_domains.drop(['strand_y', 'domain start', 'domain end'], axis=1)
cols = fragment_domains.columns.values
cols[5] = 'strand'
cols[15] = 'geneStart'
cols[16] = 'geneEnd'
fragment_domains.columns = cols
fragment_domains.head()


# In[14]:


# Save the fragment_domains dataset as a CSV file
fragment_domains.to_csv('fragment_domains.csv')


# In[15]:


# A function that takes in rows from the fragment_domain dataset and returns the coordinate for the protein start of the fragment
def findProteinStart(row):
    if row['strand'] == '+':
        end = (row['fragEnd'] - row['geneStart'] + 1) // 3
        dist = (row['fragEnd'] - row['fragStart']) // 3
        return end - dist
    elif row['strand'] == '-':
        dist = (row['geneEnd'] - row['fragStart'] + 1) // 3
        tot_dist = (row['geneEnd'] - row['geneStart']) // 3
        return tot_dist - dist
    
# A function that takes in rows from the fragment_domain dataset after finding the proteinStarts and returns the coordinate for the 
# protein end of the fragment
def findProteinEnd(row):
    return row['proteinStart'] + ((row['fragEnd'] - row['fragStart']) // 3)


# In[16]:


# Applying the functions above to find the protein coordinates of each fragment
fragment_domains['proteinStart'] = fragment_domains.apply(findProteinStart, axis=1)
fragment_domains['proteinEnd'] = fragment_domains.apply(findProteinEnd, axis=1)
fragment_domains = fragment_domains[['chrom','fragStart','fragEnd','proteinStart','proteinEnd','barcode','score','strand','domStart','domEnd','overlap',
                                    'seq id','alignment start','alignment end','hmm acc','hmm name','Gene locus','geneStart','geneEnd']]
fragment_domains.head()


# In[17]:


# Sorting the fragments BED
frags_sorted = fragments.sort()
frags_sorted.head()


# In[18]:


# Picking the most sequenced fragments out of identical fragments
frags_sorted = frags_sorted.to_dataframe()
frags_sorted.columns = ['chrom','start','stop','barcode','score','strand']

frags_sorted = frags_sorted.sort_values('score', ascending=False).drop_duplicates(['chrom','start','stop']).sort_index()
frags_sorted = pybedtools.BedTool.from_dataframe(frags_sorted)
frags_sorted.head()


# In[19]:


# Intersecting the fragments with themselves to find 90% similar fragments
frags_intersected = frags_sorted.intersect(frags_sorted, wao=True, f=.9, r=True, s=True, sorted=True, nonamecheck=True)
frags_intersected.head()


# In[20]:


# Selecting 
frags_intersected_df = frags_intersected.to_dataframe(names=['chrom A','start A','end A','barcode A','score A','strand A','chrom B','start B','end B','barcode B','score B','strand B','overlap'])
frags_intersected_collapsed = frags_intersected_df.sort_values('score B',ascending=False).drop_duplicates(['chrom A', 'start A','end A']).drop_duplicates(['chrom B','start B','end B']).sort_index()
frags_intersected_collapsed.head()


# In[21]:


frags_collapsed = frags_intersected_collapsed.iloc[:,np.arange(6,12)]
frags_collapsed.columns = ['chrom','start','end','barcode','score','strand']
frags_collapsed.head()


# In[22]:


# Merging our frags_clustered dataset with the fragment_domains dataset to map domain information onto our clustered fragments
frag_domain_collapsed = pd.merge(frags_collapsed, fragment_domains, left_on=['barcode'],
                                 right_on=['barcode'], how='left')
frag_domain_collapsed = frag_domain_collapsed.drop(['strand_y','score_y','chrom_y','fragStart','fragEnd'], axis=1)
frag_domain_collapsed


# In[23]:


#This seemed to add some rows; checking that the additional rows are because some fragments contained more than one domain
pd.concat(g for _, g in frag_domain_collapsed.groupby("barcode") if len(g) > 1).head()


# In[24]:


# Loading the joint fragment-mle-peak dataset to map mle peak data onto our frag_domain_collapsed dataset 
mlep = pd.read_csv('joint-frag-mle-peak.csv')
mlep.head()


# In[25]:


# Functions that process fragment data from the mlep dataset and returns its individual components (chromosome, start, end, strand)
def fragChrom(row):
    return re.search("chr\d*", row['frag']).group()

def fragStart(row):
    return int(re.search(":\d*", row['frag']).group()[1:])

def fragEnd(row):
    return int(re.search("-\d*", row['frag']).group()[1:])

def fragStrand(row):
    return re.search("\(.\)", row['frag']).group()


# In[26]:


mlep['chrom'] = mlep.apply(fragChrom, axis=1)
mlep['start'] = mlep.apply(fragStart, axis=1)
mlep['end'] = mlep.apply(fragEnd, axis=1)
mlep['strand'] = mlep.apply(fragStrand, axis=1)
mlep = mlep.drop('frag', axis=1)
mlep


# In[27]:


# Merging the frag_domain_collapsed dataset with our mlep dataset to map on mlep data to our fragments
frag_domain_mlep = pd.merge(frag_domain_collapsed, mlep, left_on=['chrom_x','start','end'], right_on=['chrom','start','end'], how='left')
frag_domain_mlep = frag_domain_mlep.drop('strand', axis=1).dropna(subset=['mlePeak'])
frag_domain_mlep


# In[28]:


frag_domain_mlep.to_csv('frag_domain_mlep.csv')


# In[29]:


# Finding the occurence of domains in the frag_domain_mlep dataset
domain_occurence_groups = frag_domain_mlep.groupby(['hmm name'])
domain_occurence = pd.DataFrame()
domain_occurence['domain'] = domain_occurence_groups.groups
domain_occurence['count'] = domain_occurence_groups.size()
domain_occurence = domain_occurence.reset_index().drop('index', axis=1)
domain_occurence.sort_values('count', ascending=False).head()


# In[30]:


def mannwhitney(row):
    domain = row['domain']
    with_domain = frag_domain_mlep.loc[frag_domain_mlep['hmm name'] == domain]
    wo_domain = frag_domain_mlep.loc[frag_domain_mlep['hmm name'] != domain]
    return stats.mannwhitneyu(with_domain['mlePeak'], wo_domain['mlePeak'], alternative='two-sided')


# In[31]:


domains_totest = domain_occurence.loc[domain_occurence['count'] >= 2].reset_index()
tests = pd.DataFrame(domains_totest.apply(mannwhitney, axis=1).tolist())
domain_tests = domains_totest.join(tests).drop('index', axis=1)
domain_tests


# In[34]:


def get_mlep(row):
    domain = row['domain']
    frag = frag_domain_mlep.loc[frag_domain_mlep['hmm name'] == domain]
    return np.mean(frag['mlePeak'])


# In[35]:


valids = domain_tests.loc[domain_tests['pvalue'] <= .05]
mleps = valids.apply(get_mlep, axis=1)
valids['mlePeak'] = mleps
valids = valids.sort_values('mlePeak', axis=0)
valids


# In[36]:


valids.reset_index().to_csv('domain_mleps.csv')


# In[37]:


# Dropping duplicated fragments based on their start, end, and strand
frags_full = fragments.to_dataframe()
frags_full = frags_full.drop_duplicates(['start', 'end', 'strand'])
frags_full.shape


# In[38]:


# Turning the fragments dataframe into a BedTool for intersecting below
frags_full_bed = pybedtools.BedTool.from_dataframe(frags_full)
frags_full_bed.head()


# In[39]:


# Loading the ANNOTATED GENOME dataset as a dataframe
annotations = pybedtools.BedTool('saccharomyces_cerevisiae.bed').to_dataframe()
annotations = annotations.drop(['score', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'], axis=1)
annotations_bed = pybedtools.BedTool.from_dataframe(annotations)
annotations.head()


# In[40]:


# Using the intersect tool to attach gene coordinates onto the fragment dataset
frags_coordinates = frags_full_bed.intersect(annotations_bed, wo=True, nonamecheck=True).to_dataframe()
frags_coordinates.columns = ['chrom', 'fragStart', 'fragEnd', 'barcode', 'score', 'strand', 'chrom2', 'geneStart', 'geneEnd', 'geneName', 'strand2', 'overlap']
frags_coordinates = frags_coordinates.drop(['chrom2', 'strand2'], axis=1)
print(frags_coordinates.shape)
frags_coordinates.head()


# In[41]:


# Applying findProteinStart and findProteinEnd functions to attach protein coordinates onto frags_collapsed_coordinates
frags_coordinates['proteinStart'] = frags_coordinates.apply(findProteinStart, axis=1)
frags_coordinates['proteinEnd'] = frags_coordinates.apply(findProteinEnd, axis=1)
frags_coordinates.head()


# In[42]:


frags_coordinates.to_csv('fragment-protein-coordinates_updatedALL.csv')


# In[43]:


mlep.head()


# In[54]:


# Merging the frag_domain_collapsed dataset with our mlep dataset to map on mlep data to our fragments
frags_coord_mlep = pd.merge(frags_coordinates, mlep, left_on=['chrom','fragStart','fragEnd'], right_on=['chrom','start','end'], how='left')
#frags_coord_mlep_NA = frags_coord_mlep
frags_coord_mlep = frags_coord_mlep.dropna(subset=['mlePeak'])
frags_coord_mlep.shape


# In[55]:


frags_coord_mlep.to_csv('frags_coordinates_mlep.csv')
frags_coord_mlep_NA.to_csv('frags_coordinates_NAmlep.csv')

