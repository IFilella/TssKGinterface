# Author: Isaac Filella, Riccardo Pellarin -- isaac.filella-merce@pasteur.fr
# https://www.linkedin.com/in/isaacfilella/
# 2019-07-18 11:42:30 (UTC+0100)

import os
import copy
from collections import OrderedDict
import hashlib
import random
import warnings
import math

class Homolseq(object):
    """
    Init the object homologous sequence from a fasta file
    """
    def __init__(self,seq,title,searchtool,table):
        self.seq=seq
        self.title=title
        self.length=len(seq)
        self.database=None
        self.key=None
        self.score=None
        self.evalue=None
        self.query=None
        while "/" in self.title:
            self.title=os.path.dirname(self.title)
        if "UniRef90" in title:
            self.key=self.title.split("UniRef90_")[1].split()[0]
            self.database="uniref90"
            self.query=False
        elif "SGB" in title:
            self.key=self.title.split(":")[0].split()[0]
            self.database="SGB"
            self.query=False
        else:
            self.key=title.split()[0]
            self.database="query"
            self.query=True
        if searchtool=="jackhmmer":
            if table!=None and self.query:
                self.evalue=0
                self.score=0
            elif table!=None and not self.query:
                command="grep %s %s"%(self.key,table)
                out=os.popen(command).read()
                self.evalue=float(out.split()[4])
                self.score=float(out.split()[5])
        elif searchtool=="blast":
            try:
                self.evalue=float(self.title.split()[1])
            except:
                self.evalue=0
            try:
                self.score=float(self.title.split()[2])
            except:
                self.evalue=0

class Alignment(object):
    """
    Init the object alignment 
    """
    def __init__(self,aliname=None,dictionary=None,uniqueid=False):
        """
        aliname: name of the file containing the original alignment
        dictionary: dicctionary containing the alignment (used to trimm or filter an existing alignment)
        uniqueid: to assign a uniqueid to each sequence of the alignment
        """
        if dictionary is not None:
            self.ali = dictionary
        else:
            self.ali = OrderedDict()
        if aliname is not None:
            self.aliname=aliname
            self.read_ali(self.aliname,uniqueid)
        else:
            self.aliname=None

    def __add__(self,Alignment2):
        alicopy=copy.copy(self.ali)
        alicopy.update(Alignment2.ali)
        return Alignment(aliname=None,dictionary=alicopy)

    def read_ali(self,aliname,uniqueid=False):
        f=open(aliname,"r")
        sequence=""
        for l in f:
            line=l.replace("\n","")
            if line[0] == ">":
                if len(sequence) > 0:
                    if not uniqueid:
                       self.ali[key]=sequence
                    else:
                       sequence_no_ali=sequence.replace("-","")
                       seqid = hashlib.md5(sequence_no_ali).hexdigest()
                       self.ali[key+"."+seqid]=sequence
                line=line.replace(">","")
                if "/" in line:
                    line=os.path.dirname(line)
                if "UniRef90" in line:
                    line=line.split("UniRef90_")[1].split()[0]
                elif "SGB" in line:
                    line=line.split(":")[0].split()[0]
                else:
                    line=line.split()[0]
                key = line
                sequence=""
            else:
                sequence+=line
        if not uniqueid:
           self.ali[key]=sequence
        else:
           sequence_no_ali=sequence.replace("-","")
           seqid = hashlib.md5(sequence_no_ali).hexdigest()
           self.ali[key+"."+seqid]=sequence
        print "Alignment with %d aligned sequences"%len(self.ali.keys()) 

    def write_ali(self,output):
        f=open(output,"w")
        for k in self.ali:
            f.write(">"+k+"\n")
            f.write(self.ali[k]+"\n")

    def write_fasta(self,output):
        f=open(output,"w")
        for k in self.ali:
            f.write(">"+k+"\n")
            f.write(self.ali[k].replace("-","")+"\n")
    
    def get_trimmed_ali(self,selection,aliname=None):
        """
        Get a trimmed alignment from one ali position to another.
        Selection must be a list with two integers giving the ali positions 
        for the cut (first ali position is 1)."
        """
        trimmedAli=copy.copy(self.ali)
        for seq in trimmedAli:
            trimmedAli[seq]=trimmedAli[seq][selection[0]-1:selection[1]]
        return Alignment(aliname=aliname,dictionary=trimmedAli)
    
    def get_concatenated_horizontally(self,Alignment2,spacer=""):
        """
        Get a new MSA by concatenating two MSA with same sequence ids.
        A spacer can be added as a separtor in the resulting MSA
        caveat: same size same ids order
        """
        alicopy=copy.copy(self.ali)
        for k in alicopy:
            seq2=Alignment2.ali[k]
            seq1=alicopy[k]
            alicopy[k]=seq1+spacer+seq2
        return Alignment(aliname=None,dictionary=alicopy)

    def get_clean_EvalMSA_ali(self,Evalout,aliname=None):
        """
        Get a new alignment without EvalMSA outlaiers 'Evalout'.
        EvalMSA is an independent software dedicated to search outliers in a given MSA.
        """
        
        cleanAli=copy.copy(self.ali)
        fclean=open(Evalout,"r")
        outliers=[]
        aux=False
        for l in fclean:
            if "ALIGNMENT ANALYSIS" in l:
                aux=False
            elif aux:
                line=l.split("(")[0]
                line=line.replace(" ","")
                line=line.replace("\n","")
                if line!="":
                    outliers.append(line)
            elif "Outliers:" in l:
                aux=True
                line=l.split()[1].split("(")[0]
                line=line.replace(" ","")
                line=line.replace("\n","")
                if line!="":
                    outliers.append(line)
        outliers=set(outliers)
        for o in outliers:
            if "NODE" in o:
                o=o.split(":")[0].split("/")[0]
            elif "UniRef90" in o:
                o=o.split("_")[1].split("/")[0]
            del cleanAli[o]
        return Alignment(aliname=aliname,dictionary=cleanAli)
 

class Fastafcont(object):
      """
      Init the object containing more than one fastaf that might come from different databases
      and diferent searchtools
      """
      def __init__(self,name):
          self.name=name
          self.searchtools=set()
          self.fastafs=[]
      
      def add_fastaf(self,fastaf):
          self.fastafs.append(fastaf)
          self.searchtools.add(fastaf.searchtool)
      
      def write_fasta(self,output):
          f=open(output,"w")
          if len(self.fastafs)!=0:
              for fa in self.fastafs:
                  for se in fa.homolseqs:
                      f.write(">%s\n"%se.title)
                      f.write("%s\n"%se.seq)

class Fastaf(object): 
    """
    Init the object using a fasta file coming from and homologous search done with
    jackhmmer or blast on uniref90 or SGB databases. Groups together all homologous sequences.
    """
    def __init__(self,fastaname,searchtool,table=None):
        """
        fastaname: fasta file 
        searchtool: homolog engine used, either jackhmmer or blast
        table: if jackhammer, we can add the e-value and score per homolog through 'table'
        """
        self.fastaname=fastaname
        self.searchtool=searchtool
        if searchtool!="jackhmmer" and searchtool!="blast":
            raise ValueError('ERROR: the searchtool must be either jackhmmer or blastp')
        self.homolseqs=[]
        self.table=table
        self.read_fasta()
       
    def read_fasta(self):
        homolseqs=[]
        f=open(self.fastaname,"r")
        count=0
        for l in f:
            l=l.replace("\n","")
            if l[0]==">":
                if count!=0:
                    homolseqs.append(Homolseq(seq,title,self.searchtool,self.table))
                title=l.replace(">","")
                seq=""
                count+=1
            else:
                seq=seq+l
        homolseqs.append(Homolseq(seq,title,self.searchtool,self.table))
        self.homolseqs=homolseqs
        queries=self._get_queries_title()
        print "Searchtool: " + self.searchtool
        print "#Homologs: " + str(len(homolseqs)-len(queries))
        if len(queries)>1:
            warnings.warn("WARNING: Multiple queries found")
        print "Queries: "
        for query in queries:
            print query
    
    def write_fasta(self,output,key=True,evalue=False):
        f=open(output,"w")
        for homol in self.homolseqs:
            if key:
                if evalue: f.write(">%s %.2e\n"%(homol.key,homol.evalue))
                else: f.write(">%s\n"%homol.key)
            else:
                if evalue: f.write(">%s %.2e\n"%(homol.title,homol.evalue))
                else: f.write(">%s\n"%homol.title)
            f.write("%s\n"%homol.seq)
        
    def _get_queries_title(self):
        queries = []
        for homolseq in self.homolseqs:
            if homolseq.query:
               queries.append(homolseq.title)
        return queries

    def _get_queries_seq(self):
        queries = []
        for homolseq in self.homolseqs:
            if homolseq.query:
               queries.append(homolseq.seq)
        return queries

    def _get_minseq_length(self):
        queries=self._get_queries_seq()
        if len(queries)==1:
            return int(len(queries[0])/2)
        else:
            mlen=float(sum(map(len,queries))) / len(queries)
            return math.ceil(mlen/2)

    def _get_maxseq_length(self):
        queries=self._get_queries_seq()
        if len(queries)==1:
            return int(len(queries[0])*3/2)
        else:
            mlen=float(sum(map(len,queries))) / len(queries)
            return math.ceil(mlen*3/2)
    
    def do_filtbyquerylength_fasta(self):
        filtindexes=[]
        queries=self._get_queries_seq()
        print "#Homologs BEFORE filt by query/ies lenght: " + str(len(self.homolseqs)-len(queries))
        for i,homol in enumerate(self.homolseqs):
            if self._get_minseq_length() > homol.length or homol.length > self._get_maxseq_length():
                filtindexes.append(i)
        filtindexes=set(filtindexes)
        filtindexes=list(filtindexes)
        filtindexes=sorted(filtindexes,reverse=True)
        for i in filtindexes:
            del self.homolseqs[i]
        print "#Homologs AFTER filt by query/ies lenght: " + str(len(self.homolseqs)-len(queries))

    def do_filt_fasta(self,ethr=None,lthr=[None,None]):
        filtindexes=[]
        queries=self._get_queries_seq()
        print "#Homologs BEFORE filt by lenght/e-val: " + str(len(self.homolseqs)-len(queries))
        for i,homol in enumerate(self.homolseqs):
            if ethr:
                if homol.evalue > ethr and not homol.query:
                    filtindexes.append(i)
            if lthr!=[None,None]:
                if lthr[0] > homol.length or homol.length > lthr[1]:
                    filtindexes.append(i)
        filtindexes=set(filtindexes)
        filtindexes=list(filtindexes)
        filtindexes=sorted(filtindexes,reverse=True)
        for i in filtindexes:
            del self.homolseqs[i]
        queries=self._get_queries_seq()
        print "#Homologs AFTER filt by lenght/e-val: " + str(len(self.homolseqs)-len(queries))

    def do_filtbyquartile_fasta(self,quart=0.75):
        evalues=[]
        queries=self._get_queries_seq()
        print "#Homologs BEFORE filt by %.2f best e-val: "%quart + str(len(self.homolseqs)-len(queries))
        for i,homol in enumerate(self.homolseqs):
            evalues.append((i,homol.evalue))
        evalues=sorted(evalues, key=lambda tup: tup[1])
        totlen=len(evalues)
        evaluesfilt=evalues[int(quart*totlen):]
        filtindex=list(zip(*evaluesfilt))[0]
        filtindex=sorted(filtindex,reverse=True)
        for i in filtindex:
            del self.homolseqs[i]
        print "#Homologs AFTER filt by %.2f best e-val: "%quart + str(len(self.homolseqs)-len(queries))
    
    def get_homol_keys(self):
        keys=[]
        for homol in self.homolseqs:
            keys.append(homol.key)
        return keys
    
    def get_homol_seq(self,key):
        for homol in self.homolseqs:
            if homol.key == key:
                return homol.seq

    def get_homol_title(self,key):
        for homol in self.homolseqs:
            if homol.key == key:
                return homol.title

    def do_shuffle_homolseqs(self):
        random.shuffle(self.homolseqs)
