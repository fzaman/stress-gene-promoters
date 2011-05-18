__doc__="""
MotifStats

A Class to Implement necessary Functions for the Module of Motifs-Statistics Exploration in a Given Sequence Set

File/Module:        MotifStats.py                                          
Title/Description:  Functional Implementaion for Motifs-Statistics Exploration in a Sequence Set
Version: 	    v.0.0.1                                                      
Last-Modified:      17-May-2011
Author: 	    Saddam Hossain
Email:              saddam_raj@yahoo.com 
Institution:        Bio-Bio-1
Email:              bio-bio-1@googlegroups.com
Status: 	    Development                                                 
Type: 	            Class Object Definition & Functional Implementation
Created: 	    06-May-2011         
Python-Version:     2.7.1                                                       
Change History:     17-May-2011 (v.0.0.1) (successful use of re (regular expression) package, has reduced the search time tremendously)
                    16-May-2011 (v.0.0.1) (Bio.Motif has been incorporated more robustly)
                    15-May-2011 (v.0.0.1)
                    12-May-2011 (v.0.0.1)
                    06-May-2011 (v.0.0.1)
http://bio-bio-1.wikispaces.com/
"""

__version__= "MotifStats v.0.0.1"

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Motif
from Bio.Alphabet import IUPAC
import re


class MotifStats:
    def __init__(self):
        pass

    def FindMotifAllOccurances(self,pSeqRecord,pMotif):
        forwardOccurances=[]
        reverseOccurances=[]
        allOccurances=[]
        Occurances=[]
        forwardOccurances = self.FindForwardOccurances(pSeqRecord,pMotif)
        if(self.IsSelfComplementaryDNAMotif(pMotif)==False):
            reverseOccurances = self.FindReverseOccurances(pSeqRecord,pMotif)
        allOccurances = forwardOccurances + reverseOccurances
        allOccurances = list(set(allOccurances))
        forwardOccurances.sort()
        forwardOccurances.reverse()
        reverseOccurances.sort()
        reverseOccurances.reverse()
        allOccurances.sort()
        allOccurances.reverse()
        Occurances=[forwardOccurances,reverseOccurances,allOccurances]
        return Occurances
        pass


    def FindForwardOccurances(self,pSeqRecord,pMotif):
        occurance=[]
        listindex=self.FindMotif(pSeqRecord,pMotif)
        #Reversing the indexing
        for i in listindex:
            occurance.append(len(pSeqRecord.seq.tostring())-i)
        return occurance            
        pass


    def FindReverseOccurances(self,pSeqRecord,pMotif):
        occurance=[]
        sequence=(pSeqRecord.seq)[::-1]
        motif=pMotif.complement()
        tseqrecord = SeqRecord(Seq((pSeqRecord.seq.tostring())[::-1]),id=pSeqRecord.id,name=pSeqRecord.name,description=pSeqRecord.description)
        tmotif=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
        tmotif.add_instance(Seq(str(pMotif.complement()),tmotif.alphabet))    
        listindex=self.FindMotif(tseqrecord,(tmotif.instances)[0])
        #Converting 0 to 1 indexing
        for i in listindex:
            occurance.append(i+1)

        return occurance       
        pass


    #This function should run in a for- with a Look-Up table
    #will be implemented later
    def GetMotifRegExp(self,motif):
        motifRegExp = str(motif).replace("\n","")
        motifRegExp = motifRegExp.replace("W","[AT]")
        motifRegExp = motifRegExp.replace("Y","[CT]")
        motifRegExp = motifRegExp.replace("B","[TCG]")
        motifRegExp = motifRegExp.replace("D","[ATG]")
        motifRegExp = motifRegExp.replace("H","[ATC]")
        motifRegExp = motifRegExp.replace("K","[GT]")
        motifRegExp = motifRegExp.replace("M","[CA]")
        motifRegExp = motifRegExp.replace("R","[AG]")
        motifRegExp = motifRegExp.replace("S","[CG]")
        motifRegExp = motifRegExp.replace("V","[ACG]")
        motifRegExp="("+motifRegExp+")"
        return motifRegExp
        pass


    def FindMotif(self,seqrecord,motif):
        motifRegExp = self.GetMotifRegExp(motif)
        regexp = re.compile(motifRegExp)
        results = regexp.finditer(seqrecord.seq.tostring())
        matchIndexes=[]
        for result in results:
            matchIndexes.append(result.start())
        return matchIndexes
        pass



    def IsSelfComplementaryDNAMotif(self,pMotif):
        motif=str(pMotif)
        if(motif==((pMotif.complement())[::-1])):
            return True
        return False
        pass


    def CalculateMotifFrequencyMatrix(self,pSeqRecords,pMotifs):
        motifFrequencyMatrix=[]
        motifFrequencyMatrixRow=[]
        totaloccurance = []
        motifFrequencyMatrixRow.append(None)
        for s in pSeqRecords:
            motifFrequencyMatrixRow.append(s)
        motifFrequencyMatrix.append(motifFrequencyMatrixRow)
        motifFrequencyMatrixRow = []
        i=0
        for slMotifs in pMotifs:
            for m in slMotifs.instances:
                motifFrequencyMatrixRow.append(m)
                for s in pSeqRecords:
                    totaloccurance=(self.FindMotifAllOccurances(s,m))[2]
                    motifFrequencyMatrixRow.append(len(totaloccurance))
                motifFrequencyMatrix.append(motifFrequencyMatrixRow)
                motifFrequencyMatrixRow = []
        return motifFrequencyMatrix
        pass


    def CalculateTopNMotifFrequencyMatrix(self,pSeqRecords,pMotifs,pTopN):
        motifFrequencyMatrix=[]
        motifTopNFrequencyMatrix=[]
        motifFreqSum=[]
        motifFreqSumRow=[]
        motifSum=[]
        motifFrequencyMatrix = self.CalculateMotifFrequencyMatrix(pSeqRecords,pMotifs)
        for r in motifFrequencyMatrix[1:]:
            motifFreqSumRow.append(r[0])
            motifFreqSumRow.append(sum(r[1:]))
            motifSum.append(sum(r[1:]))
            motifFreqSum.append(motifFreqSumRow)
            motifFreqSumRow=[]
            
        motifSum.sort()
        motifSum.reverse()
        motifTopNSum=motifSum[:pTopN]
        motifTopN=[]
        for r in motifFreqSum:
            if(r[1] in motifTopNSum):
                motifTopN.append(str(r[0]))        

        motifTopNFrequencyMatrix.append(motifFrequencyMatrix[0])
        for r in motifFrequencyMatrix[1:]:
            if(str(r[0]) in motifTopN):
                motifTopNFrequencyMatrix.append(r)
            
        return motifTopNFrequencyMatrix
        pass



    def CalculatePositionalDistributionMatrix(self,pSeqRecords,pMotif):
        positionalDistributionMatrix=[]
        occuranceList=[]
        seqLengths=[]
        for s in pSeqRecords:
            occuranceList.append((self.FindMotifAllOccurances(s,pMotif))[2])
            seqLengths.append(len(s.seq))    

        positionalDistributionMatrixRow=[]
        for i in range(1,max(seqLengths)+1):
            positionalDistributionMatrixRow.append(i)
            freq=0
            for o in occuranceList:
                freq+=o.count(i)
            positionalDistributionMatrixRow.append(freq)
            positionalDistributionMatrix.append(positionalDistributionMatrixRow)
            positionalDistributionMatrixRow=[]

        return positionalDistributionMatrix
        pass









