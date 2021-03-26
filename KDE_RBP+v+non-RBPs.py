
# coding: utf-8

# In[1]:


import pandas as pd
import scipy as sp
import seaborn as sns;
import matplotlib.pyplot as plt


# In[2]:


#reading datasets as csv files into panda dataframes
RNA_binding = pd.read_csv('priorities.csv')
fragment_peaks = pd.read_csv('joint-frag-mle-peak.csv')


# In[3]:


#merging datasets bv systematic gene name, separating by screen 
merged = pd.merge(RNA_binding, fragment_peaks, on='yorf', how='inner')
merged.drop(['beckmannHi','beckmannLo','hogan','tsvetanova','mitchell','essential','desc_x','desc_y','mlePeak','gene_y'], axis=1, inplace=True)
screen_15 = merged.drop(['mlePeak.18','nread.18','nbc.18'], axis=1)
screen_18 = merged.drop(['mlePeak.15','nread.15','nbc.15'], axis=1)


# In[10]:


#taking the absolute value of the MLEP scores 
screen_15['mlePeak.15'] = screen_15['mlePeak.15'].apply(abs)   
screen_18['mlePeak.18'] = screen_18['mlePeak.18'].apply(abs)  

#MLEP from into RNA-binding (score >= 2) and non-binding (score <= 1)
dist_15_binding = screen_15.loc[screen_15['score'] >= 2]
dist_15_nonbinding = screen_15.loc[screen_15['score'] <= 1]
dist_18_binding = screen_18.loc[screen_18['score'] >= 2]
dist_18_nonbinding = screen_18.loc[screen_18['score'] <= 1]


# In[5]:


dist_15_binding.head()


# In[16]:


sns.set_style("ticks")

#KDE plots for screen: 
ax = sns.kdeplot(dist_18_binding['mlePeak.18'], color='salmon', shade=True)
ax = sns.kdeplot(screen_18['mlePeak.18'], color='cornflowerblue', linestyle=":")
ax = sns.kdeplot(dist_18_nonbinding['mlePeak.18'], color='c')

plt.legend(title=None, loc='upper right', labels=['RNA-Binding', 'All proteins', 'Non RNA-binding'],
           fontsize=16)

# plt.legend(title=None, loc='upper right', labels=['All proteins', 'Non RNA-binding'],
#            fontsize=16)



plt.ylim(0,2)
#plt.xlim(-0,2)
plt.xlabel("Absolute activity score", fontsize=20)
plt.ylabel("Density", fontsize=20)
plt.yticks([0,1,2], fontsize=16)
plt.xticks([0,0.5,1,1.5,2], fontsize=16)
sns.despine()
plt.show()
#plt.savefig("screen_allvnonbinding_kde.png", format='png', dpi=300, bbox_inches='tight')


# In[15]:


plt.show()


# In[8]:


#Statisical tests of binding vs. nonbinding (a)
MW_stat_a, MW_p_a = sp.stats.mannwhitneyu(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'], alternative='two-sided')
KS_stat_a, KS_p_a = sp.stats.ks_2samp(dist_15_binding['mlePeak.15'], dist_15_nonbinding['mlePeak.15'])

#KS_p_a = 3.168474975481327e-27
#KS_stat_a = 0.3146807109940751
#MW_p_a = 3.1892154629370733e-23
#MW_stat_a = 394413.0

#Statistical tests of binding vs. all known (b)
MW_stat_b, MW_p_b = sp.stats.mannwhitneyu(dist_15_binding['mlePeak.15'], screen_15['mlePeak.15'], alternative='two-sided')
KS_stat_b, KS_p_b = sp.stats.ks_2samp(dist_15_binding['mlePeak.15'], screen_15['mlePeak.15'])

#KS_p_b = 4.210325923173611e-18
#KS_stat_b = 0.250130821559393
#MW_p_b = 7.69044204590104e-16
#MW_stat_b = 471207.0

