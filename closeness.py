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
enmax_palette = ["#F8EE9A", "#EE238E", "#455676","#FF8300","#67FF79","#67EFFF","#E8C7DE","#BCAF9C"]
color_codes_wanted = ["TssK","TssG","TssB","TssK-NTD","TssG foot1","TssG foot2","NT foot1","NT foot2"]
ccolor = lambda x: enmax_palette[color_codes_wanted.index(x)]
ncluster=150
clim=5
mixvec=True

used_subtypes = ["i1","i2","i3","i4b","i5"]
#subtypes_pal = sns.cubehelix_palette(len(used_subtypes),
#                                    light=.9, dark=.1, reverse=True,
#                                    start=1, rot=-2)
#subtypes_pal = sns.color_palette("hls", 5)
#subtypes_pal = sns.hls_palette(5, l=.3, s=.8)
#subtypes_pal = sns.color_palette("Set2",5)
#subtypes_pal = sns.color_palette("cubehelix", 5)
#subtypes_pal = sns.diverging_palette(290, 130, l=80, n=5, center="dark")
subtypes_pal = sns.color_palette("Set1", n_colors=8, desat=.99)
subtypes_lut = dict(zip(map(str, used_subtypes), subtypes_pal))

"""
clusters_dic={"clusterEAEC1.fa":["EAEC1-i1"],"clusterPA1.fa":["PA1-i3"],"clusterEAEC3.fa":["EAEC3-i2"],
              "clusterACB.fa":["ACB-i4b"],"clusterPA2.fa":["PA2-i1"],"clusterPA3.fa":["PA3-i4b"],
              "clusterVC.fa":["VC-i1"],"clusterEAEC2.fa":["EAEC2-i4b"],"clusterAT.fa":["AT-i5"],
              "clusterBM3BPM4.fa":["BM3;BPM4-i1"],"clusterAH.fa":["AH-i1"],"clusterBC2.fa":["BC2-i4b"],
              "clusterBC1.fa":["BC1-i2"],"clusterBPM1.fa":["BPM1-i3"],#"clusterBM1BPM5.fa":["BM1;BPM5-i3"],
              "clusterBPM6.fa":["BPM6-i4b"],"clusterYP5YPE1.fa":["YP5;YPE1-i3"],"clusterBM2BPM2.fa":["BM2;BPM2-i3"],
              "clusterET.fa":["ET-i4b"],"clusterYP1YPE4.fa":["YP1;YPE4-i1"],#"clusterYP2.fa":["YP2-i2"],
              "clusterSR2.fa":["SR2-i3"],"clusterSR1.fa":["SR1-i3"],"clusterKP1.fa":["KP1-i2"],
              #"clusterYEYP3YPE3.fa":["YE;YP3;YPE3-i3"],
              "clusterCJ.fa":["CJ-i1"],"clusterSEST.fa":["SEST-i3"]}
"""
clusters_dic={"clusterEAEC3.fa":["EAEC 3-i2"]}
#new
subtype={'EAEC1':'i1','EAEC2':'i4b','EAEC3':'i2','ACB':'i4b','FTN':'ii','BF':'iii','VC':'i1','PA1':'i3',
         'PA2':'i1','PA3':'i4b','SR1':'i3','SR2':'i3','SEST':'i3','KP1':'i2',
         'KP2':'i2','AH':'i1','CJ':'i1','YE':'i3','YPS1':'i1','YPS2':'i2','YPS3':'i2','YPS4':'i3','YPS5':'i2',
         'YPS6':'i3','YPE1':'i3','YPE2':'i2','YPE3':'i2','YPE4':'i3','YPE5':'i3','YPE6':'i1',
         'ET':'i4b','AT':'i5','BC1':'i2','BC2':'i4b','BM1':'i3','BM2':'i3','BM3':'i4b',
         'BM4':'i3','BM5':'i1','BPS1':'i3','BPS2':'i3','BPS3':'i3','BPS4':'i1','BPS5':'i3','BPS6':'i4b'}
#new
rename={'EAEC':'EAEC ','ACB':'A.b','VC':'V.c','PA':'P.a ','SR':'S.m ','SEST':'S.Tm','KP':'K.p ','AH':'A.h','CJ':'C.j','YE':'Y.e','YPE':'Y.p ','YPS':'Y.ps ','ET':'E.t','AT':'A.t','BC':'B.c ','BM':'B.m ','BPS':'B.ps '}
#old
#subtype={'EAEC1':'i1','EAEC2':'i4b','EAEC3':'i2','ACB':'i4b','FTN':'ii','BF':'iii','VC':'i1','PA1':'i3',
#         'PA2':'i1','PA3':'i4b','SR1':'i3','SR2':'i3','SEST':'i3','KP1':'i2',
#         'KP2':'i2','AH':'i1','CJ':'i1','YE':'i3','YP1':'i1','YP2':'i2','YP3':'i3','YP4':'i2','YP5':'i3',
#         'YPE1':'i1','YPE2':'i2','YPE3':'i3','YPE4':'i2','ET':'i4b','AT':'i5','BC1':'i2',
#         'BC2':'i4b','BM1':'i3','BM2':'i3','BM3':'i1','BPM1':'i3','BPM2':'i3','BPM3':'i3',
#         'BPM4':'i1','BPM5':'i3','BPM6':'i4b'}
#rename={'EAEC':'EAEC ','ACB':'A.b','VC':'V.c','PA':'P.a ','SR':'S.m ','SEST':'S.Tm','KP':'K.p ','AH':'A.h','CJ':'C.j','YE':'Y.e','YPE':'Y.p ','YP':'Y.ps ','ET':'E.t','AT':'A.t','BC':'B.c ','BPM':'B.ps ','BM':'B.m '}

def _blosum_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def get_sb62(seq1, seq2, matrix=matrix, gap_s=gap_s, gap_e=gap_e):
    """
    Blosum score between two sequences
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

def get_sb62_mixvec(vec1, vec2, matrixi=matrix, gap_s=gap_s, gap_e=gap_e):
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

def get_sb62_pairwise(seqs1,seqs2, matrix=matrix, gap_s=gap_s, gap_e=gap_e):
    score=0
    terms=0
    for seq1 in seqs1:
        for seq2 in seqs2:
            score+=get_sb62(seq1, seq2, matrix, gap_s, gap_e)
            terms+=1.0

    return score/terms

def _seqs2vec(sequences):
    """
    Get the normalized vectore of a cluster of sequences
    """
    for n,seq in enumerate(sequences):
        if n==0: vec=_seq2vec(seq)
        else: vec+=_seq2vec(seq)
    vec /= vec.sum(axis=-1)[...,None]
    return vec

def _seq2vec(sequence):
    """
    Transform a sequence to a vector
    - sequence: string
    """
    sequence=_substitute_opening_gap_char(sequence)
    mapper = dict([(r, i) for i, r in enumerate(aalist)])
    naa = len(sequence)
    naa_types = len(aalist)
    vec = np.zeros((naa, naa_types))
    for i, res in enumerate(list(sequence)):
        ind = mapper[res]
        vec[i, ind] = 1.
    return vec

def _substitute_opening_gap_char(seq):
    """
    Add a character for gap opening
    """
    newseq=list(seq)
    iterator=rex.finditer(seq)
    for match in iterator:
        try:
            newseq[match.span()[1]-1]="|"
        except:
            continue
    return "".join(newseq)

def _add_subtype(pre,suf):
    try:
        aux=pre+"-"+subtype[suf.split(";")[0]] 
        return aux
    except:
        return pre

def _rename(name):
    for k in rename.keys():
        if k in name: name=name.replace(k,rename[k])
    return name

def get_scores_seq(seqid,alig):
    """
    Given a sequence (seqid) of an alignment (alig) compute its scaled b62 score (range 0-1)
    against all other sequences in the alignment.
    identity_score as the maximum b62 score of a given alignment
    random_score as the maximum b62 score of a given alignment
    """
    scores=[]
    seqs=alig.ali.values()
    random_score=get_random_score(seqs,seqs,10,bootstrap=0.05)
    identity_score=get_identical_score([alig.ali[seqid]])
    for seq in alig.ali:
        score=get_sb62(alig.ali[seqid],alig.ali[seq])
        nscore=(score-random_score)/(identity_score-random_score)
        scores.append((seq,alig.ali[seq],nscore))
    return scores

def get_identical_score(bin1,bin2=None):
    """
    Compute the sb62 of all the aligned sequences of two bins (bin1, bin2) against themselves.
    Then return the mean of all computed sb62 (identical_score)
    """
    if bin2==None: bin2=[]
    tmpscore=0.0
    norm=0
    for ali1 in bin1:
        tmpscore+=get_sb62(ali1,ali1)
        norm+=1
    for ali2 in bin2:
        tmpscore+=get_sb62(ali2,ali2)
        norm+=1
    return tmpscore/norm

def get_random_score(bin1,bin2,nsample,bootstrap=1.0):
    """
    Randomize all the sequences of two sets and pairwisely compute its sb62.
    Then, compute the mean of all this sb62 repeat the whole process nsamples
    times and return the average of the obtained means (random_score)
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
            score=get_sb62_mixvec(vec1,vec2)
        else:
            score=get_sb62_pairwise(rbin1,rbin2)
        totscore+=score
    return totscore/nsample

def get_scores_seq_filtered(seqid,alig,dthr):
    """
    Given a sequence (seqid) of an alignmt (alig) compute its scaled b62 score (range 0-1) 
    against all other sequences in the alignment and return the ones below a given threshold (dthr)
    """
    scores=get_scores_seq(seqid,alig)
    filtscores=[]
    for i in scores:
        if i[2] >= dthr:
            filtscores.append(i)
    return filtscores

def get_sb62_bins(seqids,alig,dthr,thr,delimiter,filename=None):
    """
    Compute the bins of aligned sequences based on a sb62 threshold (dthr) or the % of
    sequences with higher sb62 (thr) per seqid (different seqids might share homologous seqs)
    """
    seq_dict={}
    pats=[]
    for seqid in seqids:
        pat=seqid.split(".")[0].split(delimiter)[0]
        pat=_add_subtype(pat,pat)
        if dthr:
            didsfilt=get_scores_seq_filtered(seqid,alig,dthr)
            print pat, len(didsfilt)
        elif thr:
            unsort_didsfilt=get_scores_seq(seqid,alig)
            unsort_didsfilt.sort(key=lambda x: x[2],reverse=True)
            totlen=len(unsort_didsfilt)
            didsfilt=unsort_didsfilt[0:int(thr*totlen)]
            print pat, len(didsfilt)
        else:
            print "You must provide a threshold"
            return 0
        if filename:
            fl=open(filename+"."+pat+".aln","w")
            for item in didsfilt:
                fl.write(">"+item[0]+"\n")
                fl.write(item[1]+"\n")	
        seq_dict[pat]=didsfilt
        pats.append(pat)
    return pats,seq_dict

def get_cluster_bins(seqids,alig,delimiter,out=None,re=False):
    """
    Compute the bins of aligned sequences assigning to each of them to a single seqid 
    (the one with highest sb62). Then select the ncluster*0.75 (global variable) highest
    sb62 aligned sequences
    """
    seq_dict={}
    pats=[seqid.split(delimiter)[0] for seqid in seqids]
    for i,c in enumerate(pats):
        pats[i]=_add_subtype(c,c)
        if re: pats[i]=_rename(pats[i])
    for seq in alig.ali:
        maxscore=None
        assigned=None
        for seqid in seqids:
            pat=seqid.split(delimiter)[0]
            pat=_add_subtype(pat,pat)
            if re: pat=_rename(pat)
            if seqid == seq:
                assigned=pat
                maxscore=get_sb62(alig.ali[seqid], alig.ali[seq])
                if assigned not in seq_dict:
                    seq_dict[assigned]=[(seq,alig.ali[seq],maxscore)]
                else:
                    seq_dict[assigned].append((seq,alig.ali[seq],maxscore))
                break
            if not maxscore:
                maxscore=get_sb62(alig.ali[seqid], alig.ali[seq])
                assigned=pat
            else:
                score=get_sb62(alig.ali[seqid], alig.ali[seq])
                if score>maxscore:
                    maxscore=score
                    assigned=pat
        if assigned not in seq_dict:
            seq_dict[assigned]=[(seq,alig.ali[seq],maxscore)]
        else:
            seq_dict[assigned].append((seq,alig.ali[seq],maxscore))
    for key in seq_dict:
        print "Before: ", key, len(seq_dict[key])
        seq_dict[key].sort(key=lambda x: x[2],reverse=True)
        totlen=len(seq_dict[key])
        seq_dict[key]=seq_dict[key][0:int(0.75*totlen+1)]
        if len(seq_dict[key])>ncluster:
            seq_dict[key]=random.sample(seq_dict[key],ncluster)
        print "After: ", key, len(seq_dict[key])
    if out!=None:
        with open(out+".pkl","w") as f:
            pickle.dump(seq_dict,f)
    return pats,seq_dict

def get_individual_bins(seqids,alig,delimiter,n=100):
    """
    Compute a bin per seqid generated with n copies of itself
    """ 
    seq_dict={}
    pats=[seqid.split(".")[0].split(delimiter)[0] for seqid in seqids]
    for i,c in enumerate(pats):
        pats[i]=_add_subtype(c,c)
    for seq in seqids:
        pat=seq.split(".")[0].split(delimiter)[0]
        pat=_add_subtype(pat,pat)
        seq_dict[pat]=[(seq,alig.ali[seq],0)]*n
    return pats,seq_dict

def get_closeness(pats,seq_dict,isdiagonal=False,log=False):
    """
    Compute the closeness metric beteween bins of aligned sequences
    """
    similarities= np.zeros((len(pats),len(pats)))
    intersectCount = np.zeros((len(pats),len(pats)))
    for i,patI in enumerate(pats):
        seqsI=list(zip(*seq_dict[patI]))[1]
        for j,patJ in enumerate(pats):
            if isdiagonal and i!=j: continue
            if j>i: continue
            if i==j and len(seqsI) < clim:
                similarities[i][j]=0
                continue
            keys1=set(list(zip(*seq_dict[pats[i]]))[0])
            keys2=set(list(zip(*seq_dict[pats[j]]))[0])
            intersectCount[i][j] = len(list(keys1 & keys2))
            seqsJ=list(zip(*seq_dict[patJ]))[1]
            random_score=get_random_score(seqsI,seqsJ,nsample=10)
            identical_score=get_identical_score(seqsI,seqsJ)
            score=0.0
            norm=0
            if mixvec:
                vecI=_seqs2vec(seqsI)
                vecJ=_seqs2vec(seqsJ)
                score=get_sb62_mixvec(vecI,vecJ)
            else:
                score=get_sb62_pairwise(list(seqsI),list(seqsJ)) 
            print "idscore (max):", round(identical_score,4), "randscore (min):", round(random_score,4), "score:", round(score,4)
            if log: similarity = -math.log(1-((float(score)-float(random_score))/(float(identical_score)-float(random_score))))
            else: similarity = ((float(score)-float(random_score))/(float(identical_score)-float(random_score)))
            similarities[i][j]= similarity
            similarities[j][i]= similarity
            print patI,patJ,similarity
    dfDists=pd.DataFrame(similarities,columns=pats,index=pats)
    dfCount=pd.DataFrame(intersectCount,columns=pats,index=pats)
    return dfDists, dfCount

def write_seq_dict(seq_dict,pathogens):
    '''
    Write the bins or the combinations of them into separated files. The bins or its combinations
    are given by a diccionary where its keys are filenames and its values a list of the bins to be
    added per file(pathogens).
    '''
    for filename in pathogens:
        f=open(filename,"w")
        added_seqs=[]
        for pat in pathogens[filename]:
            seqs=seq_dict[pat]
            for seq in seqs:
                if seq  not in added_seqs:
                   f.write(">"+seq[0]+"\n")
                   f.write(seq[1]+"\n")
                   added_seqs.append(seq)
        f.close()

def plot_scores_comulatives(seqids,alig):
    """
    Plot the b62 score comulative distributions for several sequences (seqids)
    against a given alignment (alig)
    """
    pylab.figure()
    for seqid in seqids:
        leg=seqid.split(".")[0]
        dids=get_scores_seq(seqid,alig)
        distances=list(zip(*dids))[2]
        ids=list(zip(*dids))[0]
        xaxis=np.linspace(0,max(distances),1000,endpoint=False)
        xaxis=sorted(xaxis,reverse=True)
        yaxis=[]
        for x in xaxis:
            y = distances < x
            yval = np.sum(y)
            yaxis.append(yval)
        pylab.plot(xaxis,yaxis,label=leg)
        pylab.legend()
    pylab.title("Comulative plots: b62 scaled score")
    pylab.xlabel("b62 scaled score")
    pylab.ylabel("Number of homologous sequences")
    pylab.savefig("plot/all.scaledb62.comulative.png")

def plot_closeness_heatmap(seqids,alig,dthr=None,thr=None,delimiter="-",out="",clustering=None,re=False,log=False):
    """
    Plot the pairwise closeness between bins of aligned sequences in a heatmap
    """
    if clustering==None:
        if dthr!=None or thr!=None:
            pats,seq_dict=get_sb62_bins(seqids,alig,dthr,thr,delimiter)
        else:
            pats,seq_dict=get_cluster_bins(seqids,alig,delimiter,out=str(ncluster)+".seq_dict",re=re)
            #pats,seq_dict=get_individual_bins(seqids,alig,delimiter) 
    else:
        f=open(clustering,"r")
        seq_dict=pickle.load(f)
        pats=[seqid.split(".")[0].split(delimiter)[0] for seqid in seqids]
    #write_seq_dict(seq_dict,clusters_dic)
    dfDists, dfCount = get_closeness(pats,seq_dict,log=log)
    columnsNames = dfDists.columns.values
    rowsNames = dfDists.index.values
    colscolor=[]
    rowscolor=[]
    #Color rows and columns by subtype
    for i,name in enumerate(zip(columnsNames,rowsNames)):
        colsubtype=name[0].split("-")[1]
        rowsubtype=name[1].split("-")[1]
        colscolor.append(subtypes_lut[colsubtype])
        rowscolor.append(subtypes_lut[rowsubtype])
    dfcolcolors=pd.DataFrame({'subtype':colscolor},index=columnsNames)
    dfrowcolors=pd.DataFrame({'subtype':rowscolor},index=rowsNames)
    #Plot heatmap
    cg = sns.clustermap(dfDists,vmin=0,vmax=1,cmap="RdBu_r",linewidths = 0.30,metric='cityblock',col_colors=dfcolcolors, row_colors=dfrowcolors)
    #Add subtype legend
    for label in used_subtypes:
        cg.ax_col_dendrogram.bar(0, 0, color=subtypes_lut[label],
                            label=label, linewidth=0)
        cg.ax_col_dendrogram.legend(loc="best", bbox_to_anchor=(0, 1.2) ,ncol=1)
    cg.savefig("plot/distanceHeatmaplog.%s.png"%(out))
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
        cg.savefig("plot/distanceHeatmap.%s.png"%(out))
    return dfDists, dfCount

def get_pair_closeness(seqid1,seqid2,delimiter,clustering,n=10):
    """
    compute the closeness metric in between two bins of aligned sequences
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
            score=get_sb62_mixvec(vec1,vec2)
        else:
            score=get_sb62_pairwise(list(seqs1),list(seqs2))
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
            pats,seq_dict=get_cluster_bins(seqids,alig,delimiter,out=legend[i]+"."+str(len(newpats))+"."+str(ncluster)+".seq_dict",re=re)
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
            pats,seq_dict=get_cluster_bins(seqids,alig,delimiter,out=legend[i]+"."+str(len(newpats))+"."+str(ncluster)+".seq_dict",re=re)
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
