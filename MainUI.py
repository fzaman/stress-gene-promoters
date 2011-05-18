import wx
from wxPython.wx import *
from OccurenceFinder import *
from MotifStats import *

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import Motif 


import re
import csv

class MainFrame(wxFrame):

    motifList = []
    promoterList = []

    def __init__(self, parent, ID, title):
        wxFrame.__init__(self, parent, ID, title, wxDefaultPosition, wxSize(1000, 600))
        motifLoadButton = wx.Button(self, label = 'Browse Motif', pos = (20, 5), size = (100, 25))
        promoterLoadButton = wx.Button(self, label = 'Browse Promoter', pos = (515, 5), size = (120, 25))
        motifLoadButton.Bind(wx.EVT_BUTTON, self.loadMotif)
        promoterLoadButton.Bind(wx.EVT_BUTTON, self.loadPromoter)

        motifLabel = wx.StaticText(self, -1, "Selected Motifs" , wx.Point(20, 40))
        promoterLabel = wx.StaticText(self, -1, "Selected Promoter" , wx.Point(515, 40))
        self.motifListBox = wx.ListBox(self, -1, (20, 60), (400, 300), [], wx.LB_SINGLE)
        self.promoterListBox = wx.ListBox(self, -1, (515, 60), (400, 300), [], wx.LB_SINGLE)
    
        saveMotifFrequencyButton =wx.Button(self, label = 'Save Motif Occurances in File', pos = (25, 450), size = (220, 25))
        saveMotifFrequencyButton.Bind(wx.EVT_BUTTON, self.savePromoterFrequency)

        noOfMotifsLabel = wx.StaticText(self, -1, "Enter No of Motifs" , wx.Point(25, 375))
        self.noOfMotifsTextCtrl = wx.TextCtrl(self, pos=(125,370), size=(40,25))
          
        fileNameLabel = wx.StaticText(self, -1, "Enter File Name" , wx.Point(25, 400))
        self.fileNameTextCtrl = wx.TextCtrl(self, pos=(125,405), size=(500,25))

        pass
    
    def loadMotif(self, event):
        dialog = wxFileDialog ( None, style = wxOPEN)
        if dialog.ShowModal() == wxID_OK:
            file_name = dialog.GetPath()
	    self.motifList = list(SeqIO.parse(file_name, "fasta"))            
            self.bio_motifs = []
	    for a_seq_record in self.motifList:
                a_bio_motif = Motif.Motif(alphabet = IUPAC.unambiguous_dna)
		a_bio_motif.add_instance(Seq(str(a_seq_record.seq), a_bio_motif.alphabet))
		self.bio_motifs.append(a_bio_motif)
		#self.motifListBox.Append(aMotif.id)
            #self.noOfMotifsTextCtrl.SetValue(str(len(self.motifList)))
        pass

    def loadPromoter(self, event):
        dialog = wxFileDialog ( None, style = wxOPEN)
        if dialog.ShowModal() == wxID_OK:
            file_name = dialog.GetPath()
	    self.promoterList = list(SeqIO.parse(file_name, "fasta"))   
	    # self.promoter_seq_records = SeqIO.parse(file_name, "fasta")
            #for aPromoter in self.promoterList:
            #    self.promoterListBox.Append(aPromoter.id)
        pass

    def savePromoterFrequency(self, event):
        fileNameLength = len(str(self.fileNameTextCtrl.GetValue()))
        
	if ( fileNameLength > 0) :
            output_file  = open(self.fileNameTextCtrl.GetValue() + '.csv', "wb")
            writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            basicInfoList = [""]

	    aMotifStats = MotifStats()
	    freq = aMotifStats.CalculateMotifFrequencyMatrix(self.promoterList, self.bio_motifs)
	    print len(freq)

           # for row in self.promoterList:
           #     basicInfoList.append(row.id)
           # writer.writerow(basicInfoList)

            #for aMotif in self.motifList:
            #    motifOccurancesList = [aMotif.id]
            #    for aPromoter in self.promoterList:
            #        anOccurenceFinder = OccuranceFinder(str(aMotif.seq).strip(), str(aPromoter.seq).strip())
            #        noOfOccurances = anOccurenceFinder.getNumberOfMotifs()
            #        motifOccurancesList.append(noOfOccurances)
            #    writer.writerow(motifOccurancesList)
            #output_file.close()
        else:
            dlg = wxMessageDialog(self, "Please, enter the file name first", "File name missing", wxOK | wxICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()

class MyApp(wxApp):
    def OnInit(self):
        frame = MainFrame(NULL, -1, "Motif Stats")
        frame.Show(true)
        self.SetTopWindow(frame)
        return true
        pass
    
app = MyApp(0)
app.MainLoop()
