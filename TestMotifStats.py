__doc__="""
Class for Testing Functionality of MotifStats-Class & Functions inside
"""
__version__= "TestMotifStats v.0.0.1"


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Motif
from Bio.Alphabet import IUPAC

import MotifStats
from MotifStats import *


class TestMotifStats:
    def __init__(self):
        pass

    def PrintMotifFrequencyMatrix(self,pMotifFrequencyMatrix):
        printRow=""
        for i in pMotifFrequencyMatrix:
            for j in i:
                if("Bio.Seq.Seq" in str(type(j))):
                    printRow+=((str(j)).center(5))
                elif("Bio.SeqRecord.SeqRecord" in str(type(j))):
                    printRow+=((str(j.id)).center(5))
                elif("None" in str(type(j))):
                    printRow+=(("").center(5))                    
                else:
                    printRow+=((str(j)).center(5))
            print("\t\t"+str(printRow))
            printRow=""
        pass

    def PrintMotifPositionalDistributionMatrix(self,vMotifPositionalDistributionMatrix):
        printRow=""
        for i in vMotifPositionalDistributionMatrix:
            for j in i:
                printRow+=((str(j)).center(5))
            print("\t\t"+str(printRow))
            printRow=""
                
        pass




    def testUseCase001FrequencyMatrix(self,pSeqRecords,pMotifs):
        motifStats = MotifStats()
        print("\n\n\n")
        print("Test UseCase001: Find Frequency Matrix for Motifs (from Motif Database) into a Input Set of Sequences")
        print("""USER INPUTS:
    1. fasta file (Containing a set of Promoter Sequences)""")

        print("""PROGRAM OUTPUTS:
    1. Well Formated text file for Frequency Distribution
    2. Graphical Presentation (2D Histograms)""")
        print("\n")
        print("Sample Simulation-")
        print("User Input 01: a Set of Promoter Sequences-")
        for iSeqRecord in pSeqRecords:
            print("\t\t"+iSeqRecord.seq)            
    
        print("\n")
        print("***The Motifs which will be considered for Frequency Search in Input Sequence]kjn s")
        print("Motifs-")
        for slMotifs in pMotifs:
            for iMotif in slMotifs.instances:
                print("\t\t"+iMotif)
        
        print("\n")
        print("Module Output:")
        print("\tFREQUENCY MATRIX")

        motifFrequencyMatrix=motifStats.CalculateMotifFrequencyMatrix(pSeqRecords,pMotifs)
        self.PrintMotifFrequencyMatrix(motifFrequencyMatrix)
        print("\n\n\n")
        
        pass
        
    def testUseCase002TopNFrequencyMatrix(self,pSeqRecords,pMotifs,pTopN):
        motifStats = MotifStats()        
        print("\n\n\n")
        print("Test UseCase002: Find Frequency Matrix for Top-N Motifs (from Motif Database) into a Input Set of Sequences")
        print("""USER INPUTS:
    1. fasta file (Containing a set of Promoter Sequences)
    2. Value of N for Top-N""")

        print("""PROGRAM OUTPUTS:
    1. Well Formated text file for Frequency Distribution
    2. Graphical Presentation (2D Histograms)""")
        print("Sample Simulation-")
        print("User Input 01: a Set of Promoter Sequences-")
        for iSeqRecord in pSeqRecords:
            print("\t\t"+iSeqRecord.seq)            

        print("N (Top-N)="+str(pTopN))
        print("***The Motifs which will be considered for Frequency Search in Input Sequences")
        print("Motifs-")
        for slMotifs in pMotifs:
            for iMotif in slMotifs.instances:
                print("\t\t"+iMotif)

        print()
        print("Module Output:")
        print("\tFREQUENCY MATRIX (Top-N)")

        motifTopNFrequencyMatrix=motifStats.CalculateTopNMotifFrequencyMatrix(pSeqRecords,pMotifs,pTopN)
        self.PrintMotifFrequencyMatrix(motifTopNFrequencyMatrix)
        print("\n\n\n")
        
        pass


    def testUseCase003PositionDistributionMatrix(self,pSeqRecords,pMotif):
        motifStats = MotifStats()        
        print("\n\n\n")
        print("Test UseCase003: Find Position Distribution Matrix for a Specific Motif into a Input Set of Sequences")
        print("""USER INPUTS:
    1. fasta file (Containing a set of Promoter Sequences)
    2. Motif - which position distribution needs to Find """)

        print("""PROGRAM OUTPUTS:
    1. Well Formated text file for Position Distribution 
    2. Graphical Presentation (2D Histograms)""")
        print("Sample Simulation-")
        print("User Input 01: a Set of Promoter Sequences-")
        for iSeqRecord in pSeqRecords:
            print("\t\t"+iSeqRecord.seq)            

        print("***The user may select a Motif from the existing data base or input a customize or new motif to search")
        print("Target Motif-")
        print("\t"+pMotif)

        print("Module Output:")
        print("\tTARGET MOTIF: "+pMotif)
        print("\tPOSITION DISTRIBUTION MATRIX")

        motifPositionalDistributionMatrix=motifStats.CalculatePositionalDistributionMatrix(pSeqRecords,pMotif)
        self.PrintMotifPositionalDistributionMatrix(motifPositionalDistributionMatrix)
        print("\n\n\n")
        
        pass

    def testMotifStatsUseCases(self,pUseCaseID):
        if(pUseCaseID == 1):
            #A Sequence will be treated as a Bio.SeqRecord object instance
            #And a Collection of Sequences (Sequence Set) will be treated as list of Bio.SeqRecord
            SeqRecord01 = SeqRecord(Seq("TTTCAATGCATTAACCTGATGATTGGC"),id="S01",name="Seq01",description="Dummy Sequence 01")
            SeqRecord02 = SeqRecord(Seq("ATCGATCGATCGTATATATACGCGTATATACGCGGCGCGATCG"),id="S02",name="Seq02",description="Dummy Sequence 02")
            SeqRecord03 = SeqRecord(Seq("CGCGCGCGCGGCGCGCGATCGATCGATCGTATATACGCGTATA"),id="S03",name="Seq03",description="Dummy Sequence 03")
            SeqRecords=[SeqRecord01,SeqRecord02,SeqRecord03]                
            #A single Motif will be represented as a single instance of Bio.Motif Class
            #And a colletion of Motifs will be a collection of Bio.Motif objects instances
            #Things to consider using Bio.Motif for Motif objects is that same-lenth(sl) motifs can only be
            #...added to a set of Motifs
            motif1=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif1.add_instance(Seq("TATA",motif1.alphabet))
            motif2=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif2.add_instance(Seq("ATCG",motif2.alphabet))
            motif3=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif3.add_instance(Seq("CG",motif3.alphabet))
            motif4=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif4.add_instance(Seq("CAAT",motif4.alphabet))
            motifs=[motif1,motif2,motif3,motif4]
            self.testUseCase001FrequencyMatrix(SeqRecords,motifs)
            pass
        elif(pUseCaseID == 2):
            print("Please Enter the N value for Top-N Analysis:")
            topN=int(input())
            #A Sequence will be treated as a Bio.SeqRecord object instance
            #And a Collection of Sequences (Sequence Set) will be treated as list of Bio.SeqRecord
            SeqRecord01 = SeqRecord(Seq("TTTCAATGCATTAACCTGATGATTGGC"),id="S01",name="Seq01",description="Dummy Sequence 01")
            SeqRecord02 = SeqRecord(Seq("ATCGATCGATCGTATATATACGCGTATATACGCGGCGCGATCG"),id="S02",name="Seq02",description="Dummy Sequence 02")
            SeqRecord03 = SeqRecord(Seq("CGCGCGCGCGGCGCGCGATCGATCGATCGTATATACGCGTATA"),id="S03",name="Seq03",description="Dummy Sequence 03")
            SeqRecords=[SeqRecord01,SeqRecord02,SeqRecord03]                
            #A single Motif will be represented as a single instance of Bio.Motif Class
            #And a colletion of Motifs will be a collection of Bio.Motif objects instances
            #Things to consider using Bio.Motif for Motif objects is that same-lenth(sl) motifs can only be
            #...added to a set of Motifs
            motif1=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif1.add_instance(Seq("TATA",motif1.alphabet))
            motif2=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif2.add_instance(Seq("ATCG",motif2.alphabet))
            motif3=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif3.add_instance(Seq("CGCG",motif3.alphabet))
            motif4=Motif.Motif(alphabet=IUPAC.unambiguous_dna)
            motif4.add_instance(Seq("CAAT",motif4.alphabet))
            motifs=[motif1,motif2,motif3,motif4]
            self.testUseCase002TopNFrequencyMatrix(SeqRecords,motifs,topN)
            pass
        elif(pUseCaseID == 3):
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

            self.testUseCase003PositionDistributionMatrix(SeqRecords,(motifs.instances)[0])            
            pass
        else:
            print("Sorry Use-Case ID not recognized!")
        
        pass


if __name__=="__main__":
    testMotifStats = TestMotifStats()
    #To Test the UseCases just change the parameter within (1,2,3)
    # 1: MotifStats UseCase001  : Motif Frequency Distribution
    # 2: MotifStats UseCase002  : Top-N Motif Frequency Distribution
    # 3: MotifStats UseCase003  : Position Distribution Matrix for a Specific Motif     
    testMotifStats.testMotifStatsUseCases(1)
    testMotifStats.testMotifStatsUseCases(2)
    testMotifStats.testMotifStatsUseCases(3)

    pass
