import fastaf
import closeness
import pylab
import os


print "\n-------------------- Testing Fastaf class --------------------\n"
hbl90=fastaf.Fastaf(fastaname="test/data/test.blast.fa",searchtool="blast")
if len(hbl90.homolseqs) != 5770: raise ValueError("ERROR: while reading fastaf")
hbl90.do_filtbyquerylength_fasta()
if len(hbl90.homolseqs) != 5302: raise ValueError("ERROR: in do_filtbyquery")
hbl90.do_filt_fasta(ethr=1.e-30)
if len(hbl90.homolseqs) != 3660: raise ValueError("ERROR: in do_filt_fasta")
hbl90.do_filt_fasta(lthr=[350,490])
if len(hbl90.homolseqs) != 3418: raise ValueError("ERROR: in do_filt_fasta")
hbl90.do_filtbyquartile_fasta()
if len(hbl90.homolseqs) != 2563: raise ValueError("ERROR: in do_filtbyquartile")
hbl90.do_shuffle_homolseqs()
if len(hbl90.homolseqs) != 2563: raise ValueError("ERROR: in do_shuffle_homolseqs")
hbl90.write_fasta(output="test/output/out.blast.fa",evalue=True)
aux = os.popen('grep ">" test/output/out.blast.fa | wc').read()
if int(aux.split()[0]) != 2563: raise ValueError("ERROR: in write_fasta")
print "\n"
hja90=fastaf.Fastaf(fastaname="test/data/test.jack.fa",searchtool="jackhmmer",table="test/data/test.jack.table")
if len(hja90.homolseqs) != 8365: raise ValueError("ERROR: while reading fastaf")
hja90.do_filtbyquerylength_fasta()
if len(hja90.homolseqs) != 7813: raise ValueError("ERROR: in do_filtbyquery")
hja90.do_filt_fasta(ethr=1.e-30)
if len(hja90.homolseqs) != 7475: raise ValueError("ERROR: in do_filt_fasta")
hja90.do_filt_fasta(lthr=[350,490])
if len(hja90.homolseqs) != 6519: raise ValueError("ERROR: in do_filt_fasta")
hja90.do_filtbyquartile_fasta()
if len(hja90.homolseqs) != 4889: raise ValueError("ERROR: in do_filtbyquartile_fasta")
hja90.do_shuffle_homolseqs()
if len(hja90.homolseqs) != 4889: raise ValueError("ERROR: in do_shuffle_homolseqs")
hja90.write_fasta(output="test/output/out.jack.fa")
aux = os.popen('grep ">" test/output/out.jack.fa | wc').read()
if int(aux.split()[0]) != 4889: raise ValueError("ERROR: in write_fasta")
print "\n-------------------- Fastaf class testing COMPLETE--------------------\n"


print "\n-------------------- Testing Alignment class --------------------\n"
ali=fastaf.Alignment(aliname="test/data/test.alignment.aln")
if len(ali.ali.keys())!= 2536: raise ValueError("ERROR: while reading Alignment")
ali.write_ali(output="test/output/out.alignment.aln")
aux = os.popen('grep ">" test/output/out.alignment.aln | wc').read()
if int(aux.split()[0])!= 2536: raise ValueError("ERROR: in write_ali")
ali.write_fasta(output="test/output/out.alignment.fa")
aux = os.popen('grep ">" test/output/out.alignment.fa | wc').read()
if int(aux.split()[0])!= 2536: raise ValueError("ERROR: in write_fasta")
ali1to20=ali.get_trimmed_ali(selection=[1,8])
if ali1to20.ali["EAEC3-TssK|eck:EC55989_3338"]!="MKIYRPLW": raise ValueError("ERROR: in get_trimmed_ali")
ali50to70=ali.get_trimmed_ali(selection=[15,25])
if ali50to70.ali["EAEC3-TssK|eck:EC55989_3338"]!="MPQQFQQQAAW": raise ValueError("ERROR: in get_trimmed_ali")
aliconc1=ali1to20.get_concatenated_horizontally(Alignment2=ali50to70)
if aliconc1.ali["EAEC3-TssK|eck:EC55989_3338"]!="MKIYRPLWMPQQFQQQAAW": raise ValueError("ERROR: in get_concatenated_horizontally")
aliconc1.write_ali(output="test/output/out.alignment.conc1.aln")
aux = os.popen('grep ">" test/output/out.alignment.conc1.aln | wc').read()
if int(aux.split()[0])!= 2536: raise ValueError("ERROR: in write_ali or get_concatenated_horizontally")
aliconc2=ali1to20.get_concatenated_horizontally(Alignment2=ali50to70,spacer="WCKWCKWCK")
if aliconc2.ali["EAEC3-TssK|eck:EC55989_3338"]!="MKIYRPLWWCKWCKWCKMPQQFQQQAAW": raise ValueError("ERROR: in get_concatenated_horizontally")
aliconc2.write_ali(output="test/output/out.alignment.conc2.aln")
aux = os.popen('grep ">" test/output/out.alignment.conc1.aln | wc').read()
if int(aux.split()[0])!= 2536: raise ValueError("ERROR: in write_ali or get_concatenated_horizontally")
aliclean=ali.get_clean_EvalMSA_ali(Evalout="test/data/test.EvalMSA.out")
if len(aliclean.ali.keys())!=2273: raise ValueError("ERROR: in get_clean_EvalMSA_ali")
try:
    aliclean.ali["A0A318MPE3"]
    raise ValueError("ERROR: in get_clean_EvalMSA_ali")
except:
    pass
aliclean.write_ali(output="test/output/out.alignment.clean.aln")
aliid=fastaf.Alignment(aliname="test/data/test.alignment.aln",uniqueid=True)
if len(aliid.ali.keys())!=2536 or len(aliid.ali.keys()[0].split("."))<=1: raise  ValueError("ERROR: in reading Alignment class with uniqueid")
aliid.write_ali(output="test/output/out.alignmentid.aln")
aliid.write_fasta(output="test/output/out.alignmentid.fa")
aliupd=ali.__add__(Alignment2=aliid)
if len(aliupd.ali.keys())!=5072: raise ValueError("ERROR: in __add__ inside Alignment class")
print "\n-------------------- Alignment class testing COMPLETE--------------------\n"


print "\n-------------------- Testing Fastafcont class --------------------\n"
cont=fastaf.Fastafcont(name="test")
cont.add_fastaf(hbl90)
cont.add_fastaf(hja90)
cont.write_fasta(output="test/output/out.cont.fa")
print "\n-------------------- Fastafcont class testing COMPLETE--------------------\n"


print "\n-------------------- Testing closeness --------------------\n"
ali1=fastaf.Alignment(aliname="test/data/test.alignment1.aln")
ali2=fastaf.Alignment(aliname="test/data/test.alignment2.aln")
ali3=fastaf.Alignment(aliname="test/data/test.alignment3.aln")
ali4=fastaf.Alignment(aliname="test/data/test.alignment4.aln")
ali5=fastaf.Alignment(aliname="test/data/test.alignment5.aln")
ali6=fastaf.Alignment(aliname="test/data/test.alignment6.aln")
alis=[ali1,ali2,ali3,ali4,ali5,ali6]
legend=["TssK-NTD","TssK","TssG","TssB","TssG-foot1","TssG-foot2"]
colors = ["#FF8300", "#F8EE9A", "#EE238E","#455676","#67FF79","#67EFFF"]
delimiters=["-","-","-","-","-","-"]
rename={'EAEC2-TssK|eck:EC55989_3287':'EAEC2-i4b', 'ET-TssK|etr:ETAE_2441':'E.t-i4b', 'BM5;BPS4-TssK|bmal:DM55_4711':'B.m5;B.ps4-i1', 'ACB-TssK|abaz:P795_10845':'A.b-i4b', 'BPS6-TssK|bps:BPSL3110':'B.ps6-i4b', 'BC2-TssK|bceo:I35_0329':'B.c2-i4b', 'BC1-TssK|bceo:I35_4156':'B.c1-i2', 'KP1-TssK|kpu:KP1_2395':'K.p1-i2', 'YPE2;YPE3;YPS2;YPS5-TssK|ype:YPO0976':'Y.p2;Y.p3;Y.ps2;Y.ps5-i2', 'YPS3-TssK|ypi:YpsIP31758_0779':'Y.ps3-i2', 'EAEC3-TssK|eck:EC55989_3338':'EAEC3-i2', 'KP2-TssK|kpu:KP1_3387':'K.p2-i2', 'PA3-TssK|paei:N296_2437':'P.a3-i4b', 'CJ-TssK|cjl:PJ17_05100':'C.j-i1', 'SR1-TssK|serf:L085_12910':'S.m1-i3', 'YPE1;YPS6-TssK|ype:YPO0513':'Y.p1;Y.ps6-i3', 'BPS3-TssK|bps:BPSS0530':'B.ps3-i3', 'AH-TssK|ahj:V469_13135':'A.h-i1', 'VC-TssK|vcm:VCM66_A0112':'V.c-i1', 'EAEC1-TssK|eck:EC55989_0225':'EAEC1-i1', 'YPE6;YPS1-TssK|ype:YPO3597':'Y.p6;Y.ps1-i1', 'PA1-TssK|paei:N296_85':'P.a1-i3', 'AT-TssK|atf:Ach5_45010':'A.t-i5', 'PA2-TssK|paei:N296_1713':'P.a2-i1', 'BPS1-TssK|bps:BPSS0101':'B.ps1-i3', 'YE;YPE5;YPS4-TssK|yet:CH48_3197':'Y.e;Y.p5;Y.ps4-i3', 'BM2;BPS2-TssK|bmal:DM55_3873':'B.m2;B.ps2-i3', 'SEST-TssK|setc:CFSAN001921_16000':'S.Tm-i3', 'SR2-TssK|serf:L085_13805':'S.m2-i3', 'BM1;BPS5-TssK|bmal:DM55_3351':'B.m1;B.ps5-i3'}
clusters_dic={"test/output/clusterEAEC3_PA1.fa":["EAEC3-i2","P.a1-i3"],"test/output/clusterVC.fa":["V.c-i1"]}
subtypes=["i1","i2","i3","i4b","i5"]
seqidss=[]
for i,ali in enumerate(alis):
    seqids=[]
    for aliseq in ali.ali.keys():
        if legend[i].split("-")[0] in aliseq:
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
print pathogens
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

print "\nClustering: \n"
pats,seq_dict=closeness.get_clustered_bins(seqids=seqidss[0],ali=alis[0],outname="test/output/out.clusterdict")
if len(pats)!=30 or len(seq_dict["EAEC3-TssK|eck:EC55989_3338"])!=46: raise ValueError("ERROR: in get_clustered_bins")
pats,seq_dict=closeness.get_individual_bins(seqids=seqidss[0],ali=alis[0],outname="test/output/out.indvdict")
if len(pats)!=30 or len(seq_dict["EAEC3-TssK|eck:EC55989_3338"])!=100: raise ValueError("ERROR: in get_clustered_bins")

print "\nClustering delimiter: \n"
pats,seq_dict=closeness.get_clustered_bins(seqids=seqidss[0],ali=alis[0],delimiter=delimiters[0],outname="test/output/out.clusterdict.delimited")
if len(pats)!=30 or len(seq_dict["EAEC3"])!=46: raise ValueError("ERROR: in get_individual_bins")
pats,seq_dict=closeness.get_individual_bins(seqids=seqidss[0],ali=alis[0],delimiter=delimiters[0],outname="test/output/out.indvdict.delimited")
if len(pats)!=30 or len(seq_dict["EAEC3"])!=100: raise ValueError("ERROR: in get_individual_bins")

print "\nClustering renamed: \n"
pats,seq_dict=closeness.get_clustered_bins(seqids=seqidss[0],ali=alis[0],rename=rename,outname="test/output/out.clusterdict.renamed")
if len(pats)!=30 or len(seq_dict["EAEC3-i2"])!=46: raise ValueError("ERROR: in get_individual_bins")
closeness.write_seq_dict(seq_dict,clusters_dict=clusters_dic)
aux1 = os.popen('grep ">" test/output/clusterEAEC3_PA1.fa | wc').read()
aux2 = os.popen('grep ">" test/output/clusterVC.fa | wc').read()
if int(aux1.split()[0])!=593 or int(aux2.split()[0])!=64: raise ValueError("ERROR: in write_seq_dict")
dfDists, dfCount = closeness.get_closeness(pats,seq_dict)
columnsNames = dfDists.columns.values
rowsNames = dfDists.index.values
for c in columnsNames:
    for r in rowsNames:
        print dfDists[c][r]
        if dfDists[c][r]<0 or dfDists[c][r]>1: raise ValueError("ERROR: in get_closenesss")
pats,seq_dict=closeness.get_individual_bins(seqids=seqidss[0],ali=alis[0],rename=rename,outname="test/output/out.indvdict.renamed")
if len(pats)!=30 or len(seq_dict["EAEC3-i2"])!=100: raise ValueError("ERROR: in get_individual_bins")
pair=closeness.get_pair_closeness(seqid1='EAEC2-TssK|eck:EC55989_3287',seqid2='ACB-TssK|abaz:P795_10845',clustering="test/data/test.seq_dict.pkl",rename=rename,n=10)
if pair[0]<0 or pair[0]>1: raise ValueError("ERROR: in get_pair_closeness")


print "\nPlots: \n"
closeness.plot_closeness_heatmap(seqids=seqidss[0],ali=alis[0],rename=rename,clustering="clustered")
closeness.plot_closeness_heatmap(seqids=seqidss[0],ali=alis[0],rename=rename,clustering="individual")
closeness.plot_closeness_heatmap(seqids=seqidss[0],ali=alis[0],rename=rename,out="test/output/closeness.heatmap.renamed.clustering",clustering="test/data/test.seq_dict.pkl")
closeness.plot_closeness_heatmap(seqids=seqidss[0],ali=alis[0],rename=rename,out="test/output/closeness.heatmap.renamed.clustering",clustering="test/data/test.seq_dict.pkl",subtypes=subtypes)
closeness.plot_closeness_heatmap(seqids=seqidss[0],ali=alis[0],rename=rename,out="test/output/closeness.heatmap.renamed.clustering",clustering="test/data/test.seq_dict.pkl",subtypes=subtypes,log=True)
closeness.plot_closeness_barplot(seqidss=seqidspats,alis=alis,legend=legend,colors=colors,out="barplot")
closeness.plot_closeness_barplots(seqidss=seqidspats,alis=alis,legend=legend,colors=colors,out="barplots")
pylab.show()
print "\n-------------------- Closeness testing COMPLETE --------------------\n"

