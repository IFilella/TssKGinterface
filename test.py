import fastaf


print "\n-------------------- Testing Fastaf class --------------------\n"
hbl90=fastaf.Fastaf(fastaname="test/data/test.blast.uniref90.fa",searchtool="blast")
hbl90.do_filtbyquerylength_fasta()
hbl90.do_filt_fasta(ethr=1.e-30)
hbl90.do_filt_fasta(lthr=[350,490])
hbl90.do_filtbyquartile_fasta()
hbl90.do_shuffle_homolseqs()
hbl90.write_fasta(output="test/output/out.blast.uniref90.fa",evalue=True)
print "\n"
hblSGB=fastaf.Fastaf(fastaname="test/data/test.blast.SGB.fa",searchtool="blast")
hblSGB.do_filtbyquerylength_fasta()
hblSGB.do_filt_fasta(ethr=1.e-30)
hblSGB.do_filt_fasta(lthr=[350,490])
hblSGB.do_filtbyquartile_fasta()
hblSGB.do_shuffle_homolseqs()
hblSGB.write_fasta(output="test/output/out.blast.SGB.fa")
print "\n"
hmulti=fastaf.Fastaf(fastaname="test/data/test.mutiquery.fa",searchtool="blast")
hmulti.do_filtbyquerylength_fasta()
hmulti.do_filt_fasta(ethr=1.e-30)
hmulti.do_filt_fasta(lthr=[350,490])
hmulti.do_filtbyquartile_fasta()
hmulti.do_shuffle_homolseqs()
hmulti.write_fasta(output="test/output/out.multiquey.fa")
print "\n"
hja90=fastaf.Fastaf(fastaname="test/data/test.jack.uniref90.fa",searchtool="jackhmmer",table="test/data/test.jack.uniref90.table")
hja90.do_filtbyquerylength_fasta()
hja90.do_filt_fasta(ethr=1.e-30)
hja90.do_filt_fasta(lthr=[350,490])
hja90.do_filtbyquartile_fasta()
hja90.write_fasta(output="test/output/out.jack.uniref90.fa")
hja90.do_shuffle_homolseqs()
print "\n"
hjaSGB=fastaf.Fastaf(fastaname="test/data/test.jack.SGB.fa",searchtool="jackhmmer",table="test/data/test.jack.SGB.table")
hjaSGB.do_filtbyquerylength_fasta()
hjaSGB.do_filt_fasta(ethr=1.e-30)
hjaSGB.do_filt_fasta(lthr=[350,490])
hjaSGB.do_filtbyquartile_fasta()
hjaSGB.write_fasta(output="test/output/out.jack.SGB.fa")
hjaSGB.do_shuffle_homolseqs()
print "\n"
print "\n-------------------- Fastaf class testing COMPLETE--------------------\n"


print "\n-------------------- Testing Alignment class --------------------\n"
ali=fastaf.Alignment(aliname="test/data/test.alignment.aln")
ali.write_ali(output="test/output/out.alignment.aln")
ali.write_fasta(output="test/output/out.alignment.fa")
ali1to20=ali.get_trimmed_ali(selection=[1,20])
ali50to70=ali.get_trimmed_ali(selection=[50,70])
aliconc1=ali1to20.get_concatenated_horizontally(Alignment2=ali50to70)
aliconc1.write_ali(output="test/output/out.alignment.conc1.aln")
aliconc2=ali1to20.get_concatenated_horizontally(Alignment2=ali50to70,spacer="WCKWCKWCK")
aliconc2.write_ali(output="test/output/out.alignment.conc2.aln")
aliclean=ali.get_clean_EvalMSA_ali(Evalout="test/data/test.EvalMSA.out")
aliclean.write_ali(output="test/output/out.alignment.clean.aln")
print "\n"
aliid=fastaf.Alignment(aliname="test/data/test.alignment.aln",uniqueid=True)
aliid.write_ali(output="test/output/out.alignmentid.aln")
aliid.write_fasta(output="test/output/out.alignmentid.fa")
print "\n"
aliupd=ali.__add__(Alignment2=aliid)
print len(aliupd.ali.keys())
print "\n"
print "\n-------------------- Alignment class testing COMPLETE--------------------\n"


print "\n-------------------- Testing Fastafcont class --------------------\n"
cont=fastaf.Fastafcont(name="test")
cont.add_fastaf(hbl90)
cont.add_fastaf(hblSGB)
cont.add_fastaf(hja90)
cont.add_fastaf(hjaSGB)
cont.write_fasta(output="test/output/out.cont.fa")
print "\n-------------------- Fastafcont class testing COMPLETE--------------------\n"
