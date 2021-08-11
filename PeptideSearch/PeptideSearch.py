#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 14:09:59 2021

@author: ekhaledian
"""


import pandas as pd
import numpy as np
import re
from itertools import chain
from more_itertools import locate
regex = re.compile('[^a-zA-Z]')
import sys


class PeptideSearch():
    
    def __init__(self, Peptides , Fasta_File):
        self.Peptides = Peptides
        self.Fasta_File = Fasta_File
        
    class Fasta:
         def __init__(self, sequence , ID):
             self.sequence = sequence
             self.ID = ID
    
    #read sequences from a fasta file and return them as a list of Fasta(sequence , ID)    
    def read_fasta(self):
        file1 = open(self.Fasta_File, 'r')
        Lines = file1.readlines()
        List=[]
        seq=""
        ID=""
        for line in Lines:
            if line[0]=='>':
                List.append(self.Fasta(seq, ID))
                if line[-1:]=="\n":
                    ID=line[:-1]
                else:
                    ID=line
                seq=""
            else:
                if line[-1:]=="\n":
                    seq=seq+line[:-1]
                elif line[-1:] in ['A', 'C', 'G', 'T','U']:
                     seq=seq+line
        List.append(self.Fasta(seq, ID))
        List=List[1:]
        return List
    
    #read sequences and IDs from a fasta and return a dataframe    
    def FastatoCsv(self):
        df2=pd.DataFrame()
        records = self.read_fasta()
        geneid= [str(records[j].ID) for j in range(len(records))]
        sequence=[str(records[j].sequence) for j in range(len(records))]
        df2["geneid"]=geneid
        df2[ "sequence"]= sequence
        return df2
    
    #read sequences from a .txt file or from a list 
    def read_Peptides(self):
        if type(self.Peptides)==list:
            return self.Peptides
        file1 = open(self.Peptides, 'r')
        Lines = file1.readlines()
        Lines= [i.strip() for i in Lines]
        return list(Lines)

     #find exact matches for a list of peptides in the Fasta_File
    def ExactMatch(self, Peptides=[]):
        if len(Peptides)==0:
            Peptides= self.read_Peptides()
        db=self.FastatoCsv()

        DBString=db['sequence']
        L=list()
        pept=list()
        start=list()
        count=0
        sys.stdout.write("[progress: ")
        for j in range(len(Peptides)):
            if j%800==0:
                sys.stdout.write(" %"+str(round(j/len(Peptides)*100)))
                sys.stdout.flush()
            i=Peptides[j]
            string_ind = list(locate(DBString, lambda a: i in a))
            if len(string_ind)>0:
                count=count+1
                for j in string_ind:
                    pept.append(i)
                    start.append(db['sequence'][j].find(i)) 
            L.append(string_ind)
        peptide = pept
        string_index = list(chain.from_iterable(L))
        start_index=list((start))
          
        sys.stdout.write("step1]\n")
        a = np.array(peptide)
        _, idx = np.unique(a, return_index=True)
        lista= a[np.sort(idx)]
        
        count=list()
        for i in range(len(lista)):
            count.append(0)
            count[i]= peptide.count(lista[i])
        
        
        len(count)      
    
        pList= [a for a in Peptides if a not in lista]
    
        df = pd.DataFrame( columns = ['Peptide','geneid/s', 'start/s', 'ExactMatch/OneMismatch'])
        j=0
        for i in range(len(lista)):
            if i%200==0:
                sys.stdout.write(":")
                sys.stdout.flush()
            start_ends=''
            geneid=''
            Occurances= count[i]
            while(Occurances>0):
                b= string_index[j]
                Occurances = Occurances-1
                geneid= db['geneid'][b][1:]
                start_ends = str(start_index[j])+ '--'+ str(start_index[j]+len(lista[i]))
                j=j+1
                df = df.append({'Peptide':lista[i],'geneid/s':geneid, 'start/s':start_ends, 'ExactMatch/OneMismatch':'ExactMatch'}, ignore_index=True)
        sys.stdout.write(" Done!]\n")
        return df, pList

    #find matches with one mismatch for a list of peptides in the Fasta_File
    def OneMismatch(self, peptides=[]):
        if len(peptides)==0:
            peptides= self.read_Peptides()

        db=self.FastatoCsv()
        DBString=db['sequence']
        df=pd.DataFrame()
        NoMatch=set()
        sys.stdout.write("This might take a while for large number of peptides, use ctrl+c if you want to interupt the script \n [Progress:")
        for i in range(len(peptides)):
            if i%300==0:
                sys.stdout.write(" %"+str(round(i/len(peptides)*100)))
                sys.stdout.flush()
            strp= peptides[i]            
            for q in range(len(strp)):
                str1= strp[0:q]
                str2= strp[q+1:]
                dfrelay1= pd.DataFrame( columns = ['Peptide', 'geneid/s', 'start/s', 'ExactMatch/OneMismatch'])
                c1 = list(locate(DBString, lambda a: str1 in a))
                c2 = list(locate(DBString, lambda a: str2 in a))
    
                common = set(c1) & set(c2)
                
                if(len(c1)>0 and len(c2)>0):
                    #check for the same sequence
                    if(len(common)>0):
                        idx= common.pop()
                        geneid = db['geneid'][idx][1:]
                        s1= db['sequence'][idx].find(str1)
                        s2= db['sequence'][idx].find(str2)
                        if s2==(s1+len(str1)+1) or s1==0 or s2==0:
                            start_ends= s1
                            dfrelay1 = dfrelay1.append({'Peptide':strp, 'geneid/s':geneid, 'start/s':start_ends, 'ExactMatch/OneMismatch':'OneMismatch'}, ignore_index=True)
                            df=df.append(dfrelay1,ignore_index=True)
                    else:
                        NoMatch.add(strp)
        sys.stdout.write(" Done!]\n")
                            


        return df, list(NoMatch)
    
    #Combine Exact match and one mismatch functions
    def MatchFinder(self):
        df, pList= self.ExactMatch()
        df2, N= self.OneMismatch(pList)
        dffinal= df.append(df2, ignore_index=True)
        dffinal.to_csv(str(self.Fasta_File)+"Matches.csv")
        return dffinal, N
    
    #search for a peptide in a sequence
    def SequenceSearch(Peptide, Sequence):
        Result, Index= -1, -1 
        x= Sequence.find(Peptide)
        if x>=0:
            Index=x
            Result=0
            return Result, Index
        
        for q in range(len(Peptide)):
            str1= Peptide[0:q]
            str2= Peptide[q+1:]
            s1= Sequence.find(str1)
            s2= Sequence.find(str2)

            if len(str1)==0 and s2>=0:
                Result=1
                if (s2>0):
                    Index=s2-1
                else:
                    Index=0
            elif len(str2)==0 and s1>=0:
                Result=1
                Index=s1
            elif s2==(s1+len(str1)+1):
                Result=1
                Index= s1
        return Result, Index   
                    

#List=["ACH", "K"]
#p= PeptideSearch("Peptides.txt", "Mtb.H37Rv.faFrame1Direct.txt")
#df, pList= p.OneMismatch(List)
# df2, N= p.OneMismatch(List)
# d, s= PeptideSearch.SequenceSearchs("sl", "abcsm")


