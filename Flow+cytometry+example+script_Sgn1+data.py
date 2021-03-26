
# coding: utf-8

# In[1]:


import numpy as np
import FlowCytometryTools
from FlowCytometryTools import FCMeasurement, PolyGate, ThresholdGate
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


## load in control 1 and check channel names
halo1 = FCMeasurement(ID='halo1', 
                       datafile= '/Users/KendraSwain/Documents/Lab Files/Experiments/flow cytometry/180125 Sgn1 FC data/Halo only 1.fcs')
halo1.channel_names                      


# In[3]:


## load in the rest of the samples
halo2 = FCMeasurement(ID='halo2', 
                       datafile= '/Users/KendraSwain/Documents/Lab Files/Experiments/flow cytometry/180125 Sgn1 FC data/Halo only 2.fcs')
sgn1 = FCMeasurement(ID='sgn1', 
                       datafile= '/Users/KendraSwain/Documents/Lab Files/Experiments/flow cytometry/180125 Sgn1 FC data/old SGN1 1.fcs')
sgn2 = FCMeasurement(ID='sgn2', 
                       datafile= '/Users/KendraSwain/Documents/Lab Files/Experiments/flow cytometry/180125 Sgn1 FC data/old SGN1 2.fcs')
sgn1.channel_names


# In[5]:


## plot halo control 1
halo1.plot(['FSC-A', 'SSC-A'])
plt.show()


# In[6]:


## draw gate around healthiest cells
gate1= PolyGate([(50000,0),(80000,0), (80000,45000), (50000,45000)],
               ('FSC-A', 'SSC-A'), region='in')
halo1.plot(['FSC-A', 'SSC-A'], gates=gate1)
plt.show()


# In[7]:


##see population that falls inside the gate
halo1_gated=halo1.gate(gate1)
halo1_gated.plot(['FSC-A','SSC-A'], gates=gate1)
plt.show()


# In[8]:


sgn1.plot(['FSC-A', 'SSC-A'])
plt.show()


# In[9]:


## gate Sgn1 healthiest cells 
gate2= PolyGate([(75000,50000),(125000,50000), (185000,125000), (150000,125000), (75000,75000)],
               ('FSC-A', 'SSC-A'), region='in')
sgn1.plot(['FSC-A', 'SSC-A'], gates=gate2)
plt.show()


# In[10]:


## gate the rest samples:
halo2_gated=halo2.gate(gate1)
sgn1_gated=sgn1.gate(gate2)
sgn2_gated=sgn2.gate(gate2)


# In[11]:


## visualize YFP / FITC channel for all samples: 
sns.kdeplot(np.log10(halo1_gated['FITC-A']), color='lightgrey',shade=True)
sns.kdeplot(np.log10(halo2_gated['FITC-A']), color='silver',shade=True)
sns.kdeplot(np.log10(sgn1_gated['FITC-A']), color='gold', shade=True)
sns.kdeplot(np.log10(sgn2_gated['FITC-A']), color='orange', shade=True)

plt.legend(title=None, loc='upper left', labels=['Halo control 1', 'Halo control 2', 'Sgn1 replicate 1', 'Sgn1 replicate 2'],
     fontsize=14, framealpha=0.0)


plt.xlabel(" ")
plt.ylabel(" ")
plt.yticks([1], "")
plt.xticks([1,2,3,4,5], fontsize=18)
sns.despine()

plt.savefig("SGN1_YFPplot.png", format='png',dpi=300, bbox_inches='tight')
plt.show()


# In[12]:


## visualize RFP / TexR channel for all samples:

sns.kdeplot(np.log10(halo1_gated['PE-Tx-Red-YG-A']), color='silver', shade=True)
sns.kdeplot(np.log10(halo2_gated['PE-Tx-Red-YG-A']), color='lightgrey', shade=True)
sns.kdeplot(np.log10(sgn1_gated['PE-Tx-Red-YG-A']), color='indianred', shade=True)
sns.kdeplot(np.log10(sgn2_gated['PE-Tx-Red-YG-A']), color='tomato', shade=True)

plt.legend(title=None, loc='upper left', labels=['Halo control 1', 'Halo control 2', 'Sgn1 replicate 1', 'Sgn1 replicate 2'],
     fontsize=14, framealpha=0.0)

plt.xlabel(" ")
plt.ylabel(" ")
plt.yticks([1], "")
plt.xticks([1,2,3,4,5], fontsize=18)
sns.despine()

plt.savefig("SGN1_RFPplot.png", format='png',dpi=300, bbox_inches='tight')
plt.show()


# In[41]:


#Plot forward scatter for all cells per sample: 
sns.kdeplot((halo1['FSC-A']), color='silver', shade=True)
sns.kdeplot((halo2['FSC-A']), color='lightgrey', shade=True)

sns.kdeplot((sgn1['FSC-A']), color='cyan', shade=True)
sns.kdeplot((sgn2['FSC-A']), color='mediumturquoise', shade=True)

plt.legend(title=None, loc='upper right', labels=['Halo control 1', 'Halo control 2', 'Sgn1 replicate 1', 'Sgn1 replicate 2'],
     fontsize=14, framealpha=0.0)


plt.xlabel(" ")
plt.ylabel(" ")
plt.yticks([1], "")
plt.xlim(35000, 300000)
plt.xticks([50000, 100000, 150000, 200000, 250000], ("50", "100", "150", "200", "250"), fontsize=18)
sns.despine()



plt.savefig("SGN1_FSCplot.png", format='png',dpi=300, bbox_inches='tight')
plt.show()


# In[13]:


## Define YFP / RFP ratios and normalize to Halo control ratio ( 0.86)
halo1_gated_ratio = (halo1_gated['FITC-A']/halo1_gated['PE-Tx-Red-YG-A'])/0.86
halo1_gated_ratio = halo1_gated_ratio.apply(np.log2)
halo2_gated_ratio = (halo2_gated['FITC-A']/halo2_gated['PE-Tx-Red-YG-A'])/0.86
halo2_gated_ratio = halo2_gated_ratio.apply(np.log2)

sgn1_gated_ratio = (sgn1_gated['FITC-A']/sgn1_gated['PE-Tx-Red-YG-A'])/0.86
sgn1_gated_ratio = sgn1_gated_ratio.apply(np.log2)
sgn2_gated_ratio = (sgn2_gated['FITC-A']/sgn2_gated['PE-Tx-Red-YG-A'])/0.86
sgn2_gated_ratio = sgn2_gated_ratio.apply(np.log2)


# In[16]:


##plot Sgn1 and halo control data for YFP/RFP ratio
sns.kdeplot(halo1_gated_ratio, color='silver', shade=True, legend=False)
sns.kdeplot(halo2_gated_ratio, color='lightgrey', shade = True, legend=False)
sns.kdeplot(sgn1_gated_ratio, color='b', shade=True, legend=False)
sns.kdeplot(sgn2_gated_ratio, color='dodgerblue', shade=True, legend=False)


plt.legend(title=None, loc='upper left', labels=['Halo control 1', 'Halo control 2', 'Sgn1 replicate 1', 'Sgn1 replicate 2'], 
           fontsize=14, framealpha=0.0)

plt.ylim(0,2.5)
plt.xlim(-2,4)
plt.xlabel(" ")
plt.ylabel(" ")
plt.yticks([1], "")
plt.xticks([-1,0,1,2,3], ("0.5", "1", "2", "4", "8"), fontsize=18)
sns.despine()
plt.axvline(2.68, 0,4.5, linestyle='--', color='dodgerblue')


#plt.savefig("SGN1 activity plot_Hires.png", format='png',dpi=300, bbox_inches='tight')
plt.show()


# In[17]:


sgn1_logY = (sgn1_gated['FITC-A'])
sgn1_logY = sgn1_logY.apply(np.log10)
sgn1_logR = (sgn1_gated['PE-Tx-Red-YG-A'])
sgn1_logR = sgn1_logR.apply(np.log10)


# In[22]:


## define tighter Halo distribution for contour plot
gate3= PolyGate([(50000,7000),(75000,15000), (75000,25000), (50000,25000)],
               ('FSC-A', 'SSC-A'), region='in')
halo1_tight=halo1.gate(gate3)
halo1_tight.shape


# In[23]:


halo1_logY = (halo1_tight['FITC-A'])
halo1_logY = halo1_logY.apply(np.log10)
halo1_logR = (halo1_tight['PE-Tx-Red-YG-A'])
halo1_logR = halo1_logR.apply(np.log10)


# In[25]:


sns.kdeplot(sgn1_logY, sgn1_logR, cmap='Blues_r')
plt.show()  


# In[24]:


sns.kdeplot(halo1_logY, halo1_logR, cmap='gist_gray')
plt.show()


# In[26]:


sns.kdeplot(sgn1_logY, sgn1_logR, cmap='Blues_r', legend=True)
sns.kdeplot(halo1_logY, halo1_logR, cmap='gist_gray', legend=True)

plt.legend(title=None, loc='upper left', labels=['Halo control', 'Sgn1'],
           fontsize=14, framealpha=0.0)

plt.ylim(3,6)
plt.xlim(3,6)
plt.xlabel(" ")
plt.ylabel(" ")
plt.yticks([3,4,5,6], fontsize=18)
plt.xticks([3,4,5,6], fontsize=18)
sns.despine(trim=True)

#plt.savefig("SGN1vHalo_RvY_HiRes.png", format='png', dpi=300, bbox_inches='tight')
plt.show()

