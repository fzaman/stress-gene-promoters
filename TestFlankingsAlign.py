__doc__="""
Class for Testing Functionality of FlankingsAlign-Class & Functions inside
"""
__version__= "TestFlankingsAlign v.0.0.1"


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Motif
from Bio.Alphabet import IUPAC

import FlankingsAlign
from FlankingsAlign import *

class TestFlankingsAlign:
    def __init__(self):
        pass

    def PrintFlankingsAlignment(self,pAlignments):
        for s in pAlignments:
            print("\t\t"+s.seq)
        pass


    def testUseCase001MotifFlankingRegionAlignment(self,pSeqRecords,pMotif,pFlankLen):
        flankingsAlign = FlankingsAlign()
        print("\n\n\n")
        print("Test UseCase001: Multiple Sequence Alignment for Flanking Regions of a Given Motif")
        print("""USER INPUTS:
    1. fasta file (Containing a set of Promoter Sequences)
    2. Motif - whose flankng regions need to be aligned
    3. Flanking Region Length""")
        

        print("""PROGRAM OUTPUTS:
    1. Well Formated text file for Motif-Logo input 
    2. Graphical Presentation of Motif-Logo""")
        print("Sample Simulation-")
        print("User Input 01: a Set of Promoter Sequences-")
        for iSeqRecord in pSeqRecords:
            print("\t\t"+iSeqRecord.seq)            

        print("***The user may select a Motif from the existing data base or input a customize or new motif to search")
        print("Target Motif-")
        print("\t"+pMotif)

        print("Flanking Region Length: "+str(pFlankLen))

        print("Module Output:")
        print("\tTARGET MOTIF: "+pMotif)
        print("\tFLANKING REGIONS ALIGNMENT")
        print("\n")

        alignments = flankingsAlign.CalculateFlankingsAlignment(pSeqRecords,pMotif,pFlankLen)
        self.PrintFlankingsAlignment(alignments)
        print("\n\n\n")
        
        pass


    def testFlankingsAlignUseCases(self,pUseCaseID):
        if(pUseCaseID == 1):
            #A Sequence will be treated as a Bio.SeqRecord object instance
            #And a Collection of Sequences (Sequence Set) will be treated as list of Bio.SeqRecord
            SeqRecord01 = SeqRecord(Seq("TTTCAATGCATTAACCTGATGATTGGC"),id="S01",name="Seq01",description="Dummy Sequence 01")
            SeqRecord02 = SeqRecord(Seq("ATCGATCGATCGTATATATACGCGTATATACGCGGCGCGATCG"),id="S02",name="Seq02",description="Dummy Sequence 02")
            SeqRecord03 = SeqRecord(Seq("CGCGCGCGCGGCGCGCGATCGATCGATCGTATATACGCGTATA"),id="S03",name="Seq03",description="Dummy Sequence 03")
            SeqRecords=[SeqRecord01,SeqRecord02,SeqRecord03]                
            #A single Motif will be represented as a single instance of Bio.Motif Class
            #And a colletion of Motifs will be a collection of Bio.Motif objects instances
            motifs=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motifs.add_instance(Seq("TATA",motifs.alphabet))
            #Take input of Flanking Region Length
            print("Please Enter the Flankng Region Length:")
            flankLen=int(input())            
            self.testUseCase001MotifFlankingRegionAlignment(SeqRecords,(motifs.instances)[0],flankLen)
        else:
            print("Sorry Use-Case ID not recognized!")        
        pass



if __name__=='__main__':
    testFlankingsAlign1=TestFlankingsAlign()
    #Usecase 1 = Multiple Sequence Alignment for Flanking Regions Alignment of a Given Motif
    testFlankingsAlign1.testFlankingsAlignUseCases(1)
    pass


