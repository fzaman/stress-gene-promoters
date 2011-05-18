__doc__="""
FlankingsAlign

A Class to do Multiple Sequence Alignment for Flanking Regions of Length N of a given Motif in a Given Set of Sequences

File/Module:        FlankingsAlign.py                                          
Title/Description:  Functional Implementation for Flanking Regions 
Version: 	    v.0.0.1                                                      
Last-Modified:      15-May-2011                                                 
Author: 	    Saddam Hossain
Email:              saddam_raj@yahoo.com 
Institution:        Bio-Bio-1
Email:              bio-bio-1@googlegroups.com
Status: 	    Development                                                 
Type: 	            Class Object Definition & Function Implementation                                     
Created: 	    06-May-2011                                                 
Python-Version:     2.7.1                                                       
Change History:     15-May-2011 (v.0.0.1)
                    12-May-2011 (v.0.0.1)
                    06-May-2011 (v.0.0.1)
http://bio-bio-1.wikispaces.com/
"""

__version__= "FlankingsAlign v.0.0.1"

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Motif
from Bio.Alphabet import IUPAC

from MotifStats import *

class FlankingsAlign:
    def __init__(self):
        pass

    def CalculateFlankingsAlignment(self,pSeqRecords,pMotif,pFlankingLen):
        alignedFlankings = []
        
        for s in pSeqRecords:
            MotifStats1 = MotifStats()
            occurances = MotifStats1.FindMotifAllOccurances(s,pMotif)
            for foc in occurances[0]:
                alignedFlankings.append(self.GetFlankingRegion(s,pMotif,foc,1,pFlankingLen))
            for roc in occurances[1]:
                alignedFlankings.append(self.GetFlankingRegion(s,pMotif,roc,2,pFlankingLen))                
            
        return alignedFlankings    
        pass


    def GetFlankingRegion(self,pSeqRecord, pMotif,pMotifPosition,pScanningDirection,pFlankingLen):
        vSequence = str(pSeqRecord.seq)

        vFlnkRegID = "dummy"
        vFlnkRegName="dummy"
        vFlnkRegDesc="dummy"
        vFlnkRegConsensus = ""

        if(pScanningDirection==1):
            leftFlank = ((vSequence[:len(vSequence)-pMotifPosition])[-pFlankingLen:]).rjust(pFlankingLen,'-')
            motif = str(pMotif)
            tflank=(vSequence[((len(vSequence)-pMotifPosition)+len(motif)):]).ljust(pFlankingLen,'-')
            rightFlank=tflank[:pFlankingLen]
            vFlnkRegConsensus = leftFlank+motif+rightFlank            
            pass
        else:
            vSequence=vSequence[::-1]
            leftFlank = (vSequence[:pMotifPosition-1])[-pFlankingLen:]
            leftFlank = str((Seq(leftFlank)).complement())
            leftFlank = leftFlank.rjust(pFlankingLen,'-')
            motif = str(pMotif)
            tflank=(vSequence[pMotifPosition+len(motif)-1:])[:pFlankingLen]
            tflank = str((Seq(tflank)).complement())
            tflank = tflank.rjust(pFlankingLen,'-')
            rightFlank=tflank
            vFlnkRegConsensus = leftFlank+motif+rightFlank            
            pass
        
        flankingRegion = SeqRecord(Seq(vFlnkRegConsensus),id=vFlnkRegID,name=vFlnkRegName,description=vFlnkRegDesc)
        

        return flankingRegion
        pass



