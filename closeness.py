# Author: Isaac Filella, Riccardo Pellarin -- isaac.filella-merce@pasteur.fr
# https://www.linkedin.com/in/isaacfilella/
# 2019-07-18 11:42:30 (UTC+0100)

import cPickle as pickle
import numpy as np
import pylab
import pandas as pd
import seaborn as sns
from Bio.SubsMat import MatrixInfo
import matplotlib
import matplotlib.pyplot as plt
import random
from matplotlib import colors as mcolors
import re
import math

matrix = MatrixInfo.blosum62
gap_s=-5
gap_e=-1
rex=re.compile('[A-Z]-')
aalist = list('ACDEFGHIKLMNPQRSTVWXY|-')
clim=5
mixvec=True

#enmax_palette = ["#F8EE9A", "#EE238E", "#455676","#FF8300","#67FF79","#67EFFF","#E8C7DE","#BCAF9C"]
#color_codes_wanted = ["TssK","TssG","TssB","TssK-NTD","TssG foot1","TssG foot2","NT foot1","NT foot2"]
#ccolor = lambda x: enmax_palette[color_codes_wanted.index(x)]

def _blosum_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def _seqs2vec(seqs):
    """
    Get the normalized vector representing a cluster of aligned sequences
    - seqs: list of strings. Each string is an aligned sequence
    """
    for n,seq in enumerate(seqs):
        if n==0: vec=_seq2vec(seq)
        else: vec+=_seq2vec(seq)
    vec /= vec.sum(axis=-1)[...,None]
    return vec

def _seq2vec(seq):
    """
    Transform a sequence to a vector
    - seq: string. Aligned sequence
    """
    seq=_substitute_opening_gap_char(seq)
    mapper = dict([(r, i) for i, r in enumerate(aalist)])
    naa = len(seq)
    naa_types = len(aalist)
    vec = np.zeros((naa, naa_types))
    for i, res in enumerate(list(seq)):
        ind = mapper[res]
        vec[i, ind] = 1.
    return vec

def _substitute_opening_gap_char(seq):
    """
    Add an specific character for gap opening. "|" by default
    - seq: string. Aligned sequence
    """
    newseq=list(seq)
    iterator=rex.finditer(seq)
    for match in iterator:
        try:
            newseq[match.span()[1]-1]="|"
        except:
            continue
    return "".join(newseq)

def _rename(name,rename):
    """
    Transform name by renaming it with rename dictionary
    - name: string
    - rename: dictionary
    """
    for k in rename.keys():
        if k==name:
            name=rename[k]
    return name

def get_subscore(seq1, seq2, matrix=matrix, gap_s=gap_s, gap_e=gap_e):
    """
    Get the substitution matrix score between two aligned sequences
    - seq1: string. Aligned sequence
    - seq2: string. Aligned sequence
    - matrix: substitution matrix (blosum62 by default)
    - gap_s: integer. gap opening penalty
    - gap_e: integer. gap extension penalty
    """
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += _blosum_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += _blosum_match(pair, matrix)
            else:
                score += gap_e
    return score

def get_subscore_mixvec(vec1, vec2, matrixi=matrix, gap_s=gap_s, gap_e=gap_e):
    """
    Get the score between two clusters of aligned sequences by using a normalized 
    vector per cluster.
    - vec1: numpy array. normalized vector of the first cluster
    - vec2: numpy array. normalized vector of the second cluster
    - matrix: substitution matrix (blosum62 by default)
    - gap_s: integer. gap opening penalty
    - gap_e: integer. gap extension penalty
    """
    score = 0
    for i in range(len(vec1)):
        n1s=np.nonzero(vec1[i][:-2])[0]
        n2s=np.nonzero(vec2[i][:-2])[0]
        for n1 in n1s:
            k1=vec1[i][n1]
            for n2 in n2s:
                k2=vec2[i][n2]
                pair=(aalist[n1],aalist[n2])
                score += _blosum_match(pair, matrix)*k1*k2
        score += gap_s*max(vec1[i][-2],vec2[i][-2])
        score += gap_e*max(vec1[i][-1],vec2[i][-1])
    return score

def get_subscore_pairwise(seqs1,seqs2, matrix=matrix, gap_s=gap_s, gap_e=gap_e):
    """
    Get the score between two clusters of aligned sequences by pairwisely computing its 
    substitution score and then returning the mean of all pairs.
    - seqs1: list of strings. Cluster of aligned sequences
    - seqs2: list of strings. Cluster of aligned sequences
    - matrix: substitution matrix (blosum62 by default)
    - gap_s: integer. gap opening penalty
    - gap_e: integer. gap extension penalty
    """
    print type(seqs1)
    score=0
    terms=0
    for seq1 in seqs1:
        for seq2 in seqs2:
            score+=get_subscore(seq1, seq2, matrix, gap_s, gap_e)
            terms+=1.0

    return score/terms

def get_identical_score(bin1,bin2=None):
    """
    Get the identical score (maximal score) for a set or between two sets of
    aligned sequences by computing the substitution matrix score between all sequences
    and themself.
    - bin1: tuple of strings. Cluster of aligned sequences
    - bin2: tuple of strinfs. Cluster of aligned sequences 
    """
    if bin2==None: bin2=[]
    tmpscore=0.0
    norm=0
    for ali1 in bin1:
        tmpscore+=get_subscore(ali1,ali1)
        norm+=1
    for ali2 in bin2:
        tmpscore+=get_subscore(ali2,ali2)
        norm+=1
    return tmpscore/norm

def get_random_score(bin1,bin2,nsample,bootstrap=1.0):
    """
    Get the random score (minimum score) betweem two sets of aligned sequences
    by randomizing all their sequences, pairwisely computing its substitution matrix score
    and computing the mean.
    - bin1: tuple of strings. Cluster of aligned sequences
    - bin2: tuple of strinfs. Cluster of aligned sequences
    - nsample: integer. Compute the random score nsample times to obtain an average random score
    - bootstrap: float. Select a different percentage of the aligned sequences each nsample time
    """
    totscore=0.0
    isdiagonal=False
    if bin1==bin2:
        isdiagonal=True
    bin1=random.sample(bin1,int(len(bin1)*bootstrap))
    bin2=random.sample(bin2,int(len(bin2)*bootstrap))
    for n in range(0,nsample):
        rbin1=[''.join(random.sample(ali1,len(ali1))) for ali1 in bin1]
        if isdiagonal:
            rbin2=rbin1 #if the two bins are identical, the randomization should also be
        else:
            rbin2=[''.join(random.sample(ali2,len(ali2))) for ali2 in bin2]
        if mixvec:
            vec1=_seqs2vec(rbin1)
            vec2=_seqs2vec(rbin2)
            score=get_subscore_mixvec(vec1,vec2)
        else:
            score=get_subscore_pairwise(rbin1,rbin2)
        totscore+=score
    return totscore/nsample

def get_clustered_bins(seqids,ali,delimiter=None,rename=None,outname=None):
    """
    Get clusters of aligned sequences by assigning them to one of the query sequences
    according to the substitution matrix score. Clusters are reduced up to 0.75 best 
    scoring.
    - seqids: list of strings. Aligned query sequences
    - ali: class fastaf.Alignment. Alignment of sequences
    - delimiter: string. Used to trim the query sequence titles by a specific delimiter 
    - rename: dictionary. If rename seqids are renamed according a given dictionary
    - outname: string. If outname the resuting clustering is saved into a pickle object
    """
    if delimiter!=None and rename!=None:
        raise ValueError('ERROR: both options not compatible either delimiter or rename') 
    seq_dict={}
    pats=[None]*len(seqids)
    for i,c in enumerate(seqids):
        if rename!=None: pats[i]=_rename(c,rename)
        elif delimiter!=None: pats[i]=c.split(delimiter)[0]
        else: pats[i]=c
    for seq in ali.ali:
        maxscore=None
        assigned=None
        for seqid in seqids:
            if rename!=None: pat=_rename(seqid,rename)
            elif delimiter!=None: pat=seqid.split(delimiter)[0]
            else: pat=seqid
            if seqid == seq:
                assigned=pat
                maxscore=get_subscore(ali.ali[seqid], ali.ali[seq])
                break
            if not maxscore:
                maxscore=get_subscore(ali.ali[seqid], ali.ali[seq])
                assigned=pat
            else:
                score=get_subscore(ali.ali[seqid], ali.ali[seq])
                if score>maxscore:
                    maxscore=score
                    assigned=pat
        if assigned not in seq_dict:
            seq_dict[assigned]=[(seq,ali.ali[seq],maxscore)]
        else:
            seq_dict[assigned].append((seq,ali.ali[seq],maxscore))
    for key in seq_dict:
        seq_dict[key].sort(key=lambda x: x[2],reverse=True)
        totlen=len(seq_dict[key])
        seq_dict[key]=seq_dict[key][0:int(0.75*totlen+1)]
        print key, "Cluster with %d total sequences and with %d when 0.75 best"%(totlen,len(seq_dict[key]))
    if outname!=None:
        with open(outname+".pkl","w") as f:
            pickle.dump(seq_dict,f)
    return pats,seq_dict

def get_individual_bins(seqids,ali,n=100,delimiter=None,rename=None,outname=None):
    """
    Get clusters of identical aligned sequences where each cluster contain n copies of a single query.
    - seqids: list of strings. Aligned query sequences
    - ali: class fastaf.Alignment. Alignment of sequences
    - n: integer. Number of query copies per query cluster
    - delimiter: string. Used to trim the query sequence titles by a specific delimiter
    - rename: dictionary. If rename seqids are renamed according a given dictionary
    - outname: string. If outname the resuting clustering is saved into a pickle object
    """ 
    if delimiter!=None and rename!=None:
        raise ValueError('ERROR: both options not compatible either delimiter or rename')
    seq_dict={}
    pats=[None]*len(seqids)
    for i,c in enumerate(seqids):
        if rename!=None: pats[i]=_rename(c,rename)
        elif delimiter!=None: pats[i]=c.split(delimiter)[0]
        else: pats[i]=c
    for seq in seqids:
        if rename!=None: pat=_rename(seq,rename)
        elif delimiter!=None: pat=seq.split(delimiter)[0]
        else: pat=seq
        seq_dict[pat]=[(seq,ali.ali[seq],0)]*n
        print pat, "Cluster with %d identical sequences"%(len(seq_dict[pat]))
    if outname!=None:
        with open(outname+".pkl","w") as f:
            pickle.dump(seq_dict,f)
    return pats,seq_dict

def write_seq_dict(seq_dict,clusters_dict):
    '''
    Write the sequences of single or combination of clusters of aligned sequence into multiple files. 
    The clusters and/or combinations of clusters are given by a diccionary where its keys are
    filenames and its values a list of the bins to be added per file.
    - seq_dict: dictionary. Clusters of aligned sequences per query sequence
    - clusters_dic: dictionary. Keys are filenames and values are lists of clusters 
                    e.g clusters_dic={"clusterEAEC3_PA1.fa":["EAEC3-i2","P.a1-i3"],"clusterVC.fa":["V.c-i1"]}
    '''
    for filename in clusters_dict:
        f=open(filename,"w")
        added_seqs=[]
        for pat in clusters_dict[filename]:
            seqs=seq_dict[pat]
            for seq in seqs:
                if seq  not in added_seqs:
                   f.write(">"+seq[0]+"\n")
                   f.write(seq[1]+"\n")
                   added_seqs.append(seq)
                else: print seq
        f.close()

def get_closeness(pats,seq_dict,isdiagonal=False,log=False):
    """
    Get the closeness metric beteween all possible pairs of clusters of aligned sequences
    - pats: list. Aligned query sequences representing the clusters
    - seq_dict: dictionary. Clusters of aligned sequences per query sequence
    - isdiagonal: Boolean. If True computes the closeness only between the clusters on itsel
    - log: Boolean. If True it computes -log(1-similarity)
    """
    similarities= np.zeros((len(pats),len(pats)))
    intersectCount = np.zeros((len(pats),len(pats)))
    for i,patI in enumerate(pats):
        seqsI=list(zip(*seq_dict[patI]))[1]
        for j,patJ in enumerate(pats):
            seqsJ=list(zip(*seq_dict[patJ]))[1]
            if isdiagonal and i!=j: continue
            if j>i: continue
            if i==j and len(seqsI) < clim:
                similarities[i][j]=0
                continue
            keys1=set(list(zip(*seq_dict[pats[i]]))[0])
            keys2=set(list(zip(*seq_dict[pats[j]]))[0])
            intersectCount[i][j] = len(list(keys1 & keys2))
            random_score=get_random_score(seqsI,seqsJ,nsample=10)
            identical_score=get_identical_score(seqsI,seqsJ)
            score=0.0
            norm=0
            if mixvec:
                vecI=_seqs2vec(seqsI)
                vecJ=_seqs2vec(seqsJ)
                score=get_subscore_mixvec(vecI,vecJ)
            else:
                score=get_subscore_pairwise(list(seqsI),list(seqsJ)) 
            print "idscore (max):", round(identical_score,4), "randscore (min):", round(random_score,4), "score:", round(score,4)
            if log: similarity = -math.log(1-((float(score)-float(random_score))/(float(identical_score)-float(random_score))))
            else: similarity = ((float(score)-float(random_score))/(float(identical_score)-float(random_score)))
            similarities[i][j]= similarity
            similarities[j][i]= similarity
            print patI,patJ,similarity
    dfDists=pd.DataFrame(similarities,columns=pats,index=pats)
    dfCount=pd.DataFrame(intersectCount,columns=pats,index=pats)
    return dfDists, dfCount

def plot_closeness_heatmap(seqids,ali,delimiter=None,rename=None,out=None,clustering=None,subtypes=None,log=False):
    """
    Plot the closeness hetamap. Each cell is the closeness between two clusters of aligned sequences.
    - seqids: list of strings. Aligned query sequences
    - ali: class fastaf.Alignment. Alignment of sequences
    - delimiter: string. Used to trim the query sequence titles by a specific delimiter
    - rename: dictionary. If rename seqids are renamed according a given dictionary
    - clustering: string. If 'clustered' compute clusters of aligned sequences using seqids
                          If 'individual' compute clusters of multiple copies of seqids
                          It can also be  pickable object containing a previously computed clustering
    - rename: dictionary. If rename seqids are renamed according a given dictionary
    - log: 
    """
    if clustering=="clustered":
        pats,seq_dict=get_clustered_bins(seqids,ali,delimiter=delimiter,rename=rename)
    elif clustering=="individual":
        pats,seq_dict=get_individual_bins(seqids,ali,delimiter=delimiter,rename=rename)
    else:
        f=open(clustering,"r")
        seq_dict=pickle.load(f)
        pats=seq_dict.keys()
    dfDists, dfCount = get_closeness(pats,seq_dict,log=log)
    if subtypes==None:
        cg = sns.clustermap(dfDists,vmin=0,vmax=1,cmap="RdBu_r",linewidths = 0.30,metric='cityblock')
    else:
        subtypes_pal = sns.color_palette("Set1", n_colors=len(subtypes), desat=.99)
        subtypes_lut = dict(zip(map(str, subtypes), subtypes_pal))
        columnsNames = dfDists.columns.values
        rowsNames = dfDists.index.values
        colscolor=[]
        rowscolor=[]
        for i,name in enumerate(zip(columnsNames,rowsNames)):
            colsubtype=name[0].split("-")[1]
            rowsubtype=name[1].split("-")[1]
            colscolor.append(subtypes_lut[colsubtype])
            rowscolor.append(subtypes_lut[rowsubtype])
        dfcolcolors=pd.DataFrame({'subtype':colscolor},index=columnsNames)
        dfrowcolors=pd.DataFrame({'subtype':rowscolor},index=rowsNames)
        cg = sns.clustermap(dfDists,vmin=0,vmax=1,cmap="RdBu_r",linewidths = 0.30,metric='cityblock',col_colors=dfcolcolors, row_colors=dfrowcolors)
        for label in subtypes:
            cg.ax_col_dendrogram.bar(0, 0, color=subtypes_lut[label],label=label, linewidth=0)
            cg.ax_col_dendrogram.legend(loc="best", bbox_to_anchor=(0, 1.2) ,ncol=1)
    if out!=None:
        if log:
            cg.savefig(out+".log.png")
        else:
            cg.savefig(out+".png")
    if log:
        idxr=cg.dendrogram_row.reordered_ind
        idxc=cg.dendrogram_col.reordered_ind
        dfDists, dfCount = get_closeness(pats,seq_dict)
        columnsNames = dfDists.columns.values
        rowsNames = dfDists.index.values
        columnsNames = [columnsNames[i] for i in idxc]
        rowsNames = [rowsNames[i] for i in idxr]
        dfDists=dfDists.reindex(columns=columnsNames,index=rowsNames)
        cg = sns.clustermap(dfDists,vmin=0,vmax=1,cmap="RdBu_r",linewidths = 0.30,metric='cityblock',col_colors=dfcolcolors, row_colors=dfrowcolors, row_cluster=False,col_cluster=False)
        cg.savefig(out+".png")
    return dfDists, dfCount

def get_pair_closeness(seqid1,seqid2,clustering,delimiter=None,n=10):
    """
    Get the closeness metric between two clusters of aligned sequences
    - seqid1:
    - seqid2:
    - clustering:
    - delimiter: 
    - n: 
    """
    f=open(clustering,"r")
    seq_dict=pickle.load(f)
    pats=[seqid1.split(".")[0].split(delimiter)[0], seqid2.split(".")[0].split(delimiter)[0]]
    closeness=[]
    seqs1=list(zip(*seq_dict[pats[0]]))[1]
    seqs2=list(zip(*seq_dict[pats[1]]))[1]
    for i in range(0,n):
        random_score=get_random_score(seqs1,seqs2,nsample=10)
        identical_score=get_identical_score(seqs1,seqs2)
        if mixvec:
            vec1=_seqs2vec(seqs1)
            vec2=_seqs2vec(seqs2)
            score=get_subscore_mixvec(vec1,vec2)
        else:
            score=get_subscore_pairwise(list(seqs1),list(seqs2))
        close = (float(score)-float(random_score))/(float(identical_score)-float(random_score))
        print close
        closeness.append(close)
    closeness=np.asarray(closeness)
    print "cl(" + pats[0] + "," + pats[1] + ")=" + str(np.mean(closeness)) + " +- " + str(np.std(closeness))

def plot_closeness_barplot(seqidss,aligs,legend,delimiters,limsy=[0,1],out="",clusterings=None,annot=True):
    """
    Plot a barplot with the closeness of multiple bins (seqidss) on themselfs for several alignments (aligs)
    """
    newpats=[]
    for i,seqids in enumerate(seqidss):
        for seqid in seqids:
            pat=seqid.split(delimiters[i])[0]
            if ";" in pat:
                pat=pat.split(";")
                for p in pat:
                    if p.split(".")[0] not in newpats: newpats.append(p.split(".")[0])
            else:
                if pat.split(".")[0] not in newpats: newpats.append(pat.split(".")[0])
    #All barplots in a single figure
    fig, ax = plt.subplots()
    rects= [None] * len(seqidss)
    labels=[]
    total_width=0.8
    single_width=1.
    n_bars = len(seqidss)
    bar_width = total_width / n_bars
    r = np.array(range(len(newpats)))+1
    x_offsets = [(j - n_bars / 2) * bar_width + bar_width / 2 for j in range(len(seqidss))]
    for i,(seqids,alig,delimiter) in enumerate(zip(seqidss,aligs,delimiters)):
        print legend[i]
        bardata=[]
        if clusterings==None:
            pats,seq_dict=get_clustered_bins(seqids,alig,delimiter,out=legend[i]+"."+str(len(newpats))+"."+".seq_dict",re=re)
        else:
            f=open(clusterings[i],"r")
            seq_dict=pickle.load(f)
            pats= seq_dict.keys()
        dfDists, dfCount = get_closeness(pats,seq_dict,isdiagonal=True)
        for npat in newpats:
            aux=False
            org = ''.join([j for j in npat if not j.isdigit()])
            npatre = npat.replace(org,rename[org])
            for col in dfDists.columns:
                if npatre in col:
                    npat= npatre+ " - " + subtype[npat]
                    if dfCount[col][col] < clim:
                        bardata.append((npat, 0, dfCount[col][col]))
                    else:
                        bardata.append((npat, dfDists[col][col], dfCount[col][col]))
                    aux=True
            if aux==False:
                npat= npatre+ " - " + subtype[npat]
                bardata.append((npat, 0, 0))
        #All barplots in a single figure
        print bardata
        color=[ccolor(legend[i])]*len(newpats)
        rects[i]=ax.bar(x_offsets[i]+r,list(zip(*bardata))[1],width=bar_width * single_width,label=legend[i],color=color)
        labels.append(list(zip(*bardata))[2])
        print "------------------------------------------------------------------" 
    ax.set_ylabel('Conservation level')
    ax.set_title('Conservation of pathogenic T6SS gene clusters')
    ax.set_xlabel('Pathogenic T6SS')
    ax.set_xticks(r)
    ax.set_xticklabels(list(zip(*bardata))[0])#,rotation='vertical')
    ax.set_ylim(bottom=limsy[0], top=limsy[1])
    ax.legend()
    #plt.axhline(y=1,color='r', linestyle='-')
    if annot==True:
        for i in range(len(labels)):
            for rect,label in zip(rects[i],labels[i]):
                height = rect.get_height()
                ax.annotate('{}'.format(int(label)),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')
    fig.tight_layout()
    fig.savefig("plot/barplot.%s.png"%(out))

def plot_closeness_barplots(seqidss,aligs,legend,delimiters,limsy=[0,1],out="",clusterings=None,annot=True):
    newpats=[]
    data=[]
    color=[]
    for prot in legend:
        color.append(ccolor(prot))
    for i,seqids in enumerate(seqidss):
        for seqid in seqids:
            pat=seqid.split(delimiters[i])[0]
            if ";" in pat:
                pat=pat.split(";")
                for p in pat:
                    if p.split(".")[0] not in newpats: newpats.append(p.split(".")[0])
            else:
                if pat.split(".")[0] not in newpats: newpats.append(pat.split(".")[0])
    for i,(seqids,alig,delimiter) in enumerate(zip(seqidss,aligs,delimiters)):
        print legend[i]
        bardata={}
        if clusterings==None:
            pats,seq_dict=get_cluster_bins(seqids,alig,delimiter,out=legend[i]+"."+str(len(newpats))+"."+".seq_dict",re=re)
        else:
            f=open(clusterings[i],"r")
            seq_dict=pickle.load(f)
            pats= seq_dict.keys()
        dfDists, dfCount = get_closeness(pats,seq_dict,isdiagonal=True)
        for npat in newpats:
            aux=False
            org = ''.join([j for j in npat if not j.isdigit()])
            npatre = npat.replace(org,rename[org])
            for col in dfDists.columns:
                if npatre in col:
                    npat= npatre+ " - " + subtype[npat]
                    if dfCount[col][col] < clim:
                        bardata[npat]=(0, dfCount[col][col])
                    else:
                        bardata[npat]=(dfDists[col][col], dfCount[col][col])
                    aux=True
            if aux==False:
                npat= npatre+ " - " + subtype[npat]
                bardata[npat]=(0, 0)
        data.append(bardata)
    if len(newpats) % 6 == 0: rows=len(newpats)/6
    else: rows=len(newpats)/6+1
    fig = plt.figure(figsize=(5*6, 4*rows))
    for j,pat in enumerate(data[0].keys()):
        values=[]
        for i,prot in enumerate(data):
            values.append(prot[pat][0])
        print pat,zip(legend,values)
        ax = fig.add_subplot(rows, 6, j+1)
        ax.bar(legend,values,color=color)
        ax.set_ylabel('Conservation level',fontsize=17)
        ax.set_title('%s'%pat,fontsize=20)
        ax.set_ylim(0, 1)
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)
        ax.set_aspect(abs(x1-x0)/abs(y1-y0))
    fig.savefig("plot/barplot.%s.png"%(out),bbox_inches='tight')