# Author: Isaac Filella, Riccardo Pellarin -- isaac.filella-merce@pasteur.fr
# https://www.linkedin.com/in/isaacfilella/
# 2021-07-07 16:42:30 (UTC+0100)

#Scprit prepared to recover the results show un the manuscript

import fastaf
import closeness
import pylab
import os
import pickle
import ast

# Data loading and parsing
aliB=fastaf.Alignment(aliname="data/alignments/tssb.aln")
aliG=fastaf.Alignment(aliname="data/alignments/tssg.aln")
aliGf1=fastaf.Alignment(aliname="data/alignments/tssgf1.aln")
aliGf2=fastaf.Alignment(aliname="data/alignments/tssgf2.aln")
aliK=fastaf.Alignment(aliname="data/alignments/tssk.aln")
aliKntd=fastaf.Alignment(aliname="data/alignments/tsskntd.aln")
alis=[aliK,aliG,aliB,aliKntd,aliGf1,aliGf2]
keylegend=["TssK","TssG","TssB","TssK","TssG","TssG"]
legend=["TssK","TssG","TssB","TssK-NTD","TssG-foot1","TssG-foot2"]
colors = ["#F8EE9A", "#EE238E","#455676","#FF8300","#67FF79","#67EFFF"]
delimiters=["-","-","-","-","-","-"]
subtypes=["i1","i2","i3","i4b","i5"]
renaming = ["data/renaming/"+l+".rename.txt" for l in legend]
drenaming=[]
for r in renaming:
    f=open(r,"r")
    contents = f.read()
    drenaming.append(ast.literal_eval(contents))
    f.close()

seqidss=[]
for i,ali in enumerate(alis):
    seqids=[]
    for aliseq in ali.ali.keys():
        if keylegend[i].split("-")[0] in aliseq:
            seqids.append(aliseq)
    seqidss.append(seqids)
pathogens=[]
for i,seqids in enumerate(seqidss):
    for seqid in seqids:
        pat=seqid.split(delimiters[i])[0]
        if ";" in pat:
            pat=pat.split(";")
            for p in pat:
                if p not in pathogens: pathogens.append(p)
        else:
            if pat not in pathogens: pathogens.append(pat)
seqidspats={}
for pat in pathogens:
    seqidspat=[]
    for ali in alis:
        aux=False
        for aliseq in ali.ali.keys():
            if (pat+"-" in aliseq) or (pat+";" in aliseq):
                seqidspat.append(aliseq)
                aux=True
                break
        if not aux: seqidspat.append(None)
    seqidspats[pat]=seqidspat

# Obtain the clusters of variants per each protein
"""
for i,l in enumerate(legend):
    print("Clustering of "+l)
    closeness.get_clustered_bins(seqidss[i],alis[i],rename=drenaming[i],outname='results/'+l+'.clusters')
"""
clusterings = ['results/'+l+'.clusters.pkl' for l in legend]
"""
# Obtain the EAEC TssK-TssG interface conservation (Figure 2C)
EAEC3dic = {'EAEC3' : seqidspats["EAEC3"]}
closeness.plot_closeness_barplots(seqidss=EAEC3dic,alis=alis,legend=legend,colors=colors,out="results/EAEC3.TssKGBbarplots",clusterings = clusterings,annot=True)
"""
# Obtain TssK-TssG interface for all pathogens (Supplementary Figure 4)
#closeness.plot_closeness_barplots(seqidss=seqidspats,alis=alis,legend=legend,colors=colors,out="results/ALLpathogens.TssKGBbarplots",clusterings = clusterings,annot=True)

# Obtain heatmaps for all TssK-TssG proteins and protein domains as well as for TssB (Figure 6 and Supplementary Figure 6B)

for i,l in enumerate(legend):
    closeness.plot_closeness_heatmap(seqids=seqidss[i],ali=alis[i],clustering=clusterings[i],pout="results/heatmap.%s"%legend[i],ddout="results/heatmap.%s"%legend[i],rename=drenaming[i],subtypes=subtypes)
    exit()
