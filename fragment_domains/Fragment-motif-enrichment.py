#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pybedtools
import re
import os
import glob
from Bio import SeqIO


# The following cells run code that take fragments we found to be in-frame from previous analysis, adjusts their coordinates to neglect partial codons, and create FASTA files of their sequences. These files can then be translated into protein sequences to feed into the MEME suite and look for motif enrichment.

# In[2]:


# Loading the datasets of in-frame fragments from the plus and minus strands
plus_inframe_df = pd.read_csv('output/plus_inframe_frags.csv').iloc[:,1:7]
minus_inframe_df = pd.read_csv('output/minus_inframe_frags.csv').iloc[:,1:7]


# In[3]:


# Dropping duplicates


# In[4]:


plus_inframe_df.drop_duplicates(['chrom','fragStart','fragEnd'], inplace=True)
print(plus_inframe_df.shape[0], 'rows')
plus_inframe_df.head()


# In[5]:


minus_inframe_df.drop_duplicates(['chrom','fragStart','fragEnd'], inplace=True)
print(minus_inframe_df.shape[0], 'rows')
minus_inframe_df.head()


# We only want to look at motif enrichment in fragments with an absolute MLE peak score of greater than 1. The following cells will load the dataset with MLE peak data for each fragment, so we can filter out irrelevant fragments and intersect with the in-frame fragments.

# In[6]:


# Loading the joint-frag-mle-peak dataset
frag_mlepeak_df = pd.read_csv('../primary_analysis/work-cached/joint-frag-mle-peak.csv')[['frag','mlePeak']]
frag_mlepeak_df.dropna(axis=0, inplace=True)
frag_mlepeak_df


# In[7]:


# Functions that process fragment data from the mlep dataset and returns its individual components (chromosome, start, end, strand)
def fragChrom(row):
    return re.search("chr\d*", row['frag']).group()

def fragStart(row):
    return int(re.search(":\d*", row['frag']).group()[1:])

def fragEnd(row):
    return int(re.search("-\d*", row['frag']).group()[1:])

def fragStrand(row):
    return re.search("\(.\)", row['frag']).group()[1:2]


# In[8]:


# Applying the above functions to the dataframe and rearranging the columns
frag_mlepeak_df['chrom'] = frag_mlepeak_df.apply(fragChrom, axis=1)
frag_mlepeak_df['start'] = frag_mlepeak_df.apply(fragStart, axis=1)
frag_mlepeak_df['end'] = frag_mlepeak_df.apply(fragEnd, axis=1)
frag_mlepeak_df['strand'] = frag_mlepeak_df.apply(fragStrand, axis=1)

frag_mlepeak_df.drop('frag', axis=1, inplace=True)

frag_mlepeak_df = frag_mlepeak_df[['chrom','start','end','strand','mlePeak']]

print(frag_mlepeak_df.shape[0], 'rows')
frag_mlepeak_df.head()


# In[9]:


# Merging the inframe fragment dataset with the mlepeak dataset to append mlepeak data onto the inframe fragments.


# In[10]:


plus_inframe_mlepeak = plus_inframe_df.merge(frag_mlepeak_df, left_on=['chrom','fragStart','fragEnd','strand'],
                                             right_on=['chrom','start','end','strand'], how='inner')
plus_inframe_mlepeak = plus_inframe_mlepeak[['chrom','start','end','barcode','mlePeak','strand']]
print(plus_inframe_mlepeak.shape[0], 'rows')
plus_inframe_mlepeak.head()


# In[11]:


minus_inframe_mlepeak = minus_inframe_df.merge(frag_mlepeak_df, left_on=['chrom','fragStart','fragEnd','strand'],
                                             right_on=['chrom','start','end','strand'], how='inner')
minus_inframe_mlepeak = minus_inframe_mlepeak[['chrom','start','end','barcode','mlePeak','strand']]
print(minus_inframe_mlepeak.shape[0], 'rows')
minus_inframe_mlepeak.head()


# In[12]:


#****************************************************************************


# In[13]:


# We don't care for the partial codons, so to make it easier to translate later we will adjust the 
# start and end coordinates of the fragments.
# For plus strand fragments, we will add 1 to the start coordinate.
# For minus strand fragments, we will subtract 1 from the end coordinate.
plus_inframe_mlepeak_adjusted = plus_inframe_mlepeak
plus_inframe_mlepeak_adjusted['start'] = plus_inframe_mlepeak['start'] + 1
minus_inframe_mlepeak_adjusted = minus_inframe_mlepeak
minus_inframe_mlepeak_adjusted['end'] = minus_inframe_mlepeak['end'] - 1


# In[14]:


#****************************************************************************


# In[15]:


# Concatenating the plus and minus dataframes into one.
inframe_mlepeak_adjusted = pd.concat([plus_inframe_mlepeak_adjusted,minus_inframe_mlepeak_adjusted])
inframe_mlepeak_adjusted


# In[16]:


# Creating two new dataframes, one for in-frame fragments with an mlePeak greater than 1, another for in-frame
# fragments less than -1.
pos_mlep_df = inframe_mlepeak_adjusted[inframe_mlepeak_adjusted['mlePeak'] >= 1]
neg_mlep_df = inframe_mlepeak_adjusted[inframe_mlepeak_adjusted['mlePeak'] <= -1]


# In[17]:


print(pos_mlep_df.shape[0], 'rows')
pos_mlep_df.head()


# In[18]:


print(neg_mlep_df.shape[0], 'rows')
neg_mlep_df.sort_values(['chrom','start','end'])


# In[19]:


# A function that takes in two fragments as Series objects and returns their minimum percentage overlap. 
# Assumes A comes before B. If A and share a start or stop index, returns 1.
#This is a helper function for COLLAPSE below.
def overlap(A, B):
    if A['chrom'] != B['chrom']:
        return 0
    elif (A['end'] < B['start']) or (A['start'] > B['end']):
        return 0
    elif (A['start'] == B['start']) or (A['end'] > B['end']):
        return 1
    else:
        lenA = A['end'] - A['start']
        lenB = B['end'] - B['start']
        overlap = A['end'] - B['start']
        AoverB = overlap / lenA
        BoverA = overlap / lenB
        return min(AoverB, BoverA)


# In[20]:


# Another helper function for collapse below. Takes in the CLUSTER list, iterates through them and returns 
# the index of the entry with the greatest absolute mlePeak score.
def maxindex(cluster):
    i = 0
    k = 0
    maxm = 0
    for entry in cluster:
        if np.abs(entry['mlePeak']) > maxm:
            maxm = np.abs(entry['mlePeak'])
            k = i
        i += 1
    return k


# Explanation for my algorithm below:
# 
# The function begins with two rows (A,B) of the dataframe, and creates an empty list called CLUSTER. It adds the first row A to CLUSTER, and if the second row B is similar enough given by the THRESHOLD, it will add B to CLUSTER. The function iterates on by assigning B to the next row, going down the dataframe until CLUSTER is filled with fragments similar to A by the THRESHOLD. It then picks the fragment with the greatest absolute mlePeak score and appends it to RESULT. In the next iteration, the function will look at the row after B and assign that as A. The iteration stops when there are no more rows to look at, and returns the RESULT as a DataFrame.

# In[21]:


# A function that takes in a dataframe DF of fragments and returns a collapsed copy. The collapsed copy 
# omits fragments that are similar to one another given by the threshhold TH between 0 and 1.
# Picks the fragment with the higher absolute MLE Peak score.
def collapse(df, th):
    sortd = df.sort_values(['chrom','start','end'])
    result = []
    indexA = 0
    indexB = 1
    while indexB < sortd.shape[0]:
        cluster = []
        A = sortd.iloc[indexA]
        B = sortd.iloc[indexB]
        cluster.append(A)
        while (overlap(A, B) >= th) and (indexB < sortd.shape[0] - 1):
            cluster.append(B)
            indexB += 1
            B = sortd.iloc[indexB]
        k = maxindex(cluster)
        result.append(cluster[k])
        indexA = indexB + 1
        indexB = indexA + 1
    return pd.DataFrame(result)
            


# In[22]:


# Creating collapsed copies of the dataframes
pos_frags_clpsd = collapse(pos_mlep_df, 0.5)
neg_frags_clpsd = collapse(neg_mlep_df, 0.5)


# In[23]:


print(pos_frags_clpsd.shape[0], 'rows')
pos_frags_clpsd.head()


# In[24]:


print(neg_frags_clpsd.shape[0], 'rows')
neg_frags_clpsd.head()


# In[25]:


# Making BedTools from the DataFrames and saving them as .bed files


# In[26]:


pos_frags_bed = pybedtools.BedTool.from_dataframe(pos_frags_clpsd)
neg_frags_bed = pybedtools.BedTool.from_dataframe(neg_frags_clpsd)


# In[27]:


pos_frags_bed.saveas('frag-seqs/pos_mlep_frags.bed')
neg_frags_bed.saveas('frag-seqs/neg_mlep_frags.bed')


# In[28]:


# Terminal commands using bedtools' getfasta command. Looks at the bed files we saved above, the sequences
# of the chromosomes, and returns a file of the sequences of the fragments.


# In[29]:


os.system('bedtools getfasta -fi data/saccharomyces_cerevisiae.fa -bed frag-seqs/pos_mlep_frags.bed -s -fo frag-seqs/pos_mlep_frag_seqs.fasta')


# In[30]:


os.system('bedtools getfasta -fi data/saccharomyces_cerevisiae.fa -bed frag-seqs/neg_mlep_frags.bed -s -fo frag-seqs/neg_mlep_frag_seqs.fasta')


# In[ ]:





# In[31]:


## fixing the headers on the fasta files of the chromosome's sequences to match with my bed files


# In[32]:


# chr_dir = 'yeastgenome/'


# # In[33]:


# with os.scandir(chr_dir) as chr_files:
#     for chr_entry in chr_files:
#         chrom = chr_entry.name[:5]
#         fixed = chrom + 'fixed.fsa'
#         with open(chr_entry) as original, open(fixed, 'w+') as fixed:
#             records = SeqIO.parse(original, 'fasta')
#             for record in records:
#                 record.id = chrom
#                 record.description = chrom
#                 SeqIO.write(record, fixed, 'fasta')


# In[34]:


# get_ipython().system('cat yeastgenomefixed/*.fsa > yeastgenomefixed/genomefixed.fsa')


# In[ ]:





# In[ ]:





# In[35]:


# ran EMBOSS transeq to translate the fragment sequences into their amino acid sequences,
# spot checking their validity on the UCSC Genome Browser


# In[36]:


# transeq discarded chromosome information for each of the fragments, fixing that below


# In[37]:


pos_seqs = 'frag-seqs/pos_mlep_frag_seqs.fasta'
neg_seqs = 'frag-seqs/neg_mlep_frag_seqs.fasta'
pos_aa_seqs = 'frag-seqs/pos_mlep_frag_aa_seqs.faa'
neg_aa_seqs = 'frag-seqs/neg_mlep_frag_aa_seqs.faa'


# In[38]:


with open(pos_seqs) as pos_nuc, open(pos_aa_seqs) as pos_aa, open('frag-seqs/pos_aa_seqs_fixed.faa', 'w') as fixed:
    nuc_records = SeqIO.parse(pos_nuc, 'fasta')
    chr_list = []
    for nuc_record in nuc_records:
        chr_list.append(nuc_record.id[0:5])
    aa_records = SeqIO.parse(pos_aa, 'fasta')
    index = 0
    for aa_record in aa_records:
        aa_record.id = aa_record.id[:-1] + chr_list[index]
        aa_record.description = aa_record.description[:-1] + chr_list[index]
        index += 1
        SeqIO.write(aa_record, fixed, 'fasta')


# In[39]:


with open(neg_seqs) as neg_nuc, open(neg_aa_seqs) as neg_aa, open('frag-seqs/neg_aa_seqs_fixed.faa', 'w') as fixed:
    nuc_records = SeqIO.parse(neg_nuc, 'fasta')
    chr_list = []
    for nuc_record in nuc_records:
        chr_list.append(nuc_record.id[0:5])
    aa_records = SeqIO.parse(neg_aa, 'fasta')
    index = 0
    for aa_record in aa_records:
        aa_record.id = aa_record.id[:-1] + chr_list[index]
        aa_record.description = aa_record.description[:-1] + chr_list[index]
        index += 1
        SeqIO.write(aa_record, fixed, 'fasta')


# In[ ]:

os.system('meme frag-seqs/neg_aa_seqs_fixed.faa -protein -oc frag-seqs/neg_mlep_meme -nostatus -time 18000 -mod anr -nmotifs 8 -minw 6 -maxw 50 -objfun classic -markov_order 0')
os.system('meme frag-seqs/pos_aa_seqs_fixed.faa -protein -oc frag-seqs/pos_mlep_meme -nostatus -time 18000 -mod anr -nmotifs 8 -minw 6 -maxw 50 -objfun classic -markov_order 0')



# In[40]:


#################################################################################
# Plugged motifs found from MEME into FIMO to collect genes those motifs are found in.
# Goal is to see which motifs are legitimate motifs that belong to a variety of genes, or if they are just
# repetitive sequences belonging to telomeric elements or families of similar genes.


# In[ ]:





# In[41]:


neg_tsv_path = r'FIMO/neg-mlep'
pos_tsv_path = r'FIMO/pos-mlep'

neg_tsv = glob.glob(os.path.join(neg_tsv_path, "*.tsv"))
pos_tsv = glob.glob(os.path.join(pos_tsv_path, "*.tsv"))

neg_motif_genes = pd.concat((pd.read_csv(f, sep='\t') for f in neg_tsv))
neg_motif_genes = neg_motif_genes[neg_motif_genes.motif_id.str.isnumeric()].drop('motif_alt_id', axis=1)
pos_motif_genes = pd.concat((pd.read_csv(f, sep='\t') for f in pos_tsv))
pos_motif_genes = pos_motif_genes[pos_motif_genes.motif_id.str.isnumeric()].drop('motif_alt_id', axis=1)


# In[51]:


neg_motif_genes.drop_duplicates(['sequence_name','motif_id'], inplace=True)
pos_motif_genes.drop_duplicates(['sequence_name','motif_id'], inplace=True)


# In[52]:


print(neg_motif_genes.shape[0], 'rows')
neg_motif_genes.head()


# In[53]:


print(pos_motif_genes.shape[0], 'rows')
pos_motif_genes.head()


# In[54]:


gene_dict = {}
genes_fasta = 'orf_trans.fasta'
with open(genes_fasta) as genes:
    records = SeqIO.parse(genes, 'fasta')
    for record in records:
        desc = re.search(r'"(.*?)"', record.description).group(0)
        gene_dict[record.id] = desc[1:-1]


# In[55]:


neg_description = []
for index, row in neg_motif_genes.iterrows():
    neg_description.append(gene_dict.get(row.sequence_name))
neg_motif_genes['gene_description'] = neg_description
neg_motif_genes.head()


# In[56]:


pos_description = []
for index, row in pos_motif_genes.iterrows():
    pos_description.append(gene_dict.get(row.sequence_name))
pos_motif_genes['gene_description'] = pos_description
pos_motif_genes.head()


# In[58]:


neg_motif_genes.to_csv('FIMO/neg_mlep_motif_genes.csv')
pos_motif_genes.to_csv('FIMO/pos_mlep_motif_genes.csv')


# In[ ]:




