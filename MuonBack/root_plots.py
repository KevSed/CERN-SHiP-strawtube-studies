import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt

# Supress display output
if not '-showCanvas' in sys.argv:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

##########################################################
#                     INPUT                              #
##########################################################

# Load data from rootfiles
print('Loading datafiles...')
dataFile = ROOT.TFile.Open('ship.conical.MuonBack-TGeant4-100000.root', 'read')
dataTree = dataFile.Get('cbmsim')
print('Finished.')

print('Entries in data file: {}'.format(str(dataTree.GetEntries())))
hist1D = ROOT.TH1F('hist','muonPoint.fPz',100,-0.11,0.09)

# Load data from reco rootfiles
print('Loading recofiles...')
dataFile_rec = ROOT.TFile.Open('ship.conical.MuonBack-TGeant4-100000_rec.root', 'read')
dataTree_rec = dataFile_rec.Get('cbmsim;3')
print('Finished.')

print('Entries in reconstructed file: {}'.format(str(dataTree_rec.GetEntries())))


# Create a canvas
canvas = ROOT.TCanvas('can','can',650,650)

Pdg_val = np.zeros(1, dtype=float)
dataTree.SetBranchAddress("MCTrack.fPdgCode", Pdg_val)
nentries = 1514723 #dataTree.GetEntries()

Pdg = np.zeros(nentries, dtype=int)
for i in range(nentries):
    dataTree.GetEntry(i)
    Pdg[i] = Pdg_val
dataFile.Close()

print('Drawing the Histogram.')
plt.hist(Pdg, bins=20, label='MCTrack.fPdgCode')
plt.xlim(-10, 10)
plt.legend(loc='best')
plt.grid()
name = 'pdgcode_hist.pdf'
plt.savefig('pdgcode_hist.pdf')
print('Finished. Histogram saved in \'{}\''.format(name))
plt.show()
