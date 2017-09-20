from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
import math
#from decorators import *
import sys

if len(sys.argv)==1:
    sys.exit("Please choose one of the flags -new -old -all and repeat.")

##################################################################################
##################################################################################


# Print iterations progress
def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=50):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '#' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

#################################################################################
##################################################################################

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)


h00 = r.TH2D('h00','x:y hits in 1st plane of T1 B-field=0.0 (mc/weights)',160,-600,600,160,-600,600)
h01 = r.TH2D('h01','x:y hits in 1st plane of T1 B-field=0.1 (mc/weights)',160,-600,600,160,-600,600)
h02 = r.TH2D('h02','x:y hits in 1st plane of T1 B-field=0.2 (mc/weights)',160,-600,600,160,-600,600)
h04 = r.TH2D('h04','x:y hits in 1st plane of T1 B-field=0.4 (mc/weights)',160,-600,600,160,-600,600)
h06 = r.TH2D('h06','x:y hits in 1st plane of T1 B-field=0.6 (mc/weights)',160,-600,600,160,-600,600)
h08 = r.TH2D('h08','x:y hits in 1st plane of T1 B-field=0.8 (mc/weights)',160,-600,600,160,-600,600)
h10 = r.TH2D('h10','x:y hits in 1st plane of T1 B-field=1.0 (mc/weights)',160,-600,600,160,-600,600)
h12 = r.TH2D('h12','x:y hits in 1st plane of T1 B-field=1.2 (mc/weights)',160,-600,600,160,-600,600)
h14 = r.TH2D('h14','x:y hits in 1st plane of T1 B-field=1.4 (mc/weights)',160,-600,600,160,-600,600)
h16 = r.TH2D('h16','x:y hits in 1st plane of T1 B-field=1.6 (mc/weights)',160,-600,600,160,-600,600)
h18 = r.TH2D('h18','x:y hits in 1st plane of T1 B-field=1.8 (mc/weights)',160,-600,600,160,-600,600)

z00 = r.TH1D('z00','strawtubesPoint z B-field=0.0 T1',250,2580,2598)
z00.SetStats(False)

PDG = r.TH1D('PDG','particle distribution at T1',60,-30,30)
PDG.SetStats(False)

if '-new' in sys.argv:
    f00 = r.TFile("/eos/user/k/ksedlacz/files/new/1M-00-MuonBack-new_rec.root", "READ")
    f02 = r.TFile("/eos/user/k/ksedlacz/files/new/02-100000-MuonBack-new_rec.root", "READ")
    f04 = r.TFile("/eos/user/k/ksedlacz/files/new/04-100000-MuonBack-new_rec.root", "READ")
    f06 = r.TFile("/eos/user/k/ksedlacz/files/new/06-100000-MuonBack-new_rec.root", "READ")
    f08 = r.TFile("/eos/user/k/ksedlacz/files/new/08-100000-MuonBack-new_rec.root", "READ")
    f10 = r.TFile("/eos/user/k/ksedlacz/files/new/10-100000-MuonBack-new_rec.root", "READ")
    f14 = r.TFile("/eos/user/k/ksedlacz/files/new/14-100000-MuonBack-new_rec.root", "READ")
    f18 = r.TFile("/eos/user/k/ksedlacz/files/new/18-100000-MuonBack-new_rec.root", "READ")
    directory='new'
    t00 = f00.cbmsim
    t02 = f02.cbmsim
    t04 = f04.cbmsim
    t06 = f06.cbmsim
    t08 = f08.cbmsim
    t10 = f10.cbmsim
    t14 = f14.cbmsim
    t18 = f18.cbmsim
    trees = [t00,t02,t04,t06,t08,t10,t14,t18]
    hists = [h00,h02,h04,h06,h08,h10,h14,h18]
    fields = [0.0,0.2,0.4,0.6,0.8,1.0,1.4,1.8]
    texts = [0,0,0,0,0,0,0,0]
    xprojections = [0,0,0,0,0,0,0,0]
    yprojections = [0,0,0,0,0,0,0,0]
    xtexts = [0,0,0,0,0,0,0,0]
    ytexts = [0,0,0,0,0,0,0,0]

if '-old' in sys.argv:
    f00 = r.TFile("/eos/user/k/ksedlacz/files/old/00-100000-MuonBack_rec.root", "READ")
    f01 = r.TFile("/eos/user/k/ksedlacz/files/old/01-100000-MuonBack_rec.root", "READ")
    f02 = r.TFile("/eos/user/k/ksedlacz/files/old/02-100000-MuonBack_rec.root", "READ")
    f04 = r.TFile("/eos/user/k/ksedlacz/files/old/04-100000-MuonBack_rec.root", "READ")
    f06 = r.TFile("/eos/user/k/ksedlacz/files/old/06-100000-MuonBack_rec.root", "READ")
    f10 = r.TFile("/eos/user/k/ksedlacz/files/old/10-100000-MuonBack_rec.root", "READ")
    f12 = r.TFile("/eos/user/k/ksedlacz/files/old/12-100000-MuonBack_rec.root", "READ")
    directory='old'
    t00 = f00.cbmsim
    t01 = f01.cbmsim
    t02 = f02.cbmsim
    t04 = f04.cbmsim
    t06 = f06.cbmsim
    #t08 = f08.cbmsim
    t10 = f10.cbmsim
    t12 = f12.cbmsim
    trees = [t00,t01,t02,t04,t06,t10,t12]
    hists = [h00,h01,h02,h04,h06,h10,h12]
    fields = [0.0,0.1,0.2,0.4,0.6,1.0,1.2]
    texts = [0,0,0,0,0,0,0]
    xprojections = [0,0,0,0,0,0,0]
    yprojections = [0,0,0,0,0,0,0]
    xtexts = [0,0,0,0,0,0,0]
    ytexts = [0,0,0,0,0,0,0]

if '-all' in sys.argv:
    f00 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/00-100000-MuonBack_rec.root", "READ")
    f02 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/02-100000-MuonBack_rec.root", "READ")
    f04 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/04-100000-MuonBack_rec.root", "READ")
    f06 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/06-100000-MuonBack_rec.root", "READ")
    f08 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/08-100000-MuonBack_rec.root", "READ")
    f10 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/10-100000-MuonBack_rec.root", "READ")
    f12 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/12-100000-MuonBack_rec.root", "READ")
    f14 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/14-100000-MuonBack_rec.root", "READ")
    f16 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/16-100000-MuonBack_rec.root", "READ")
    f18 = r.TFile("/eos/user/k/ksedlacz/files/all-lep/18-100000-MuonBack_rec.root", "READ")
    directory='all'
    t00 = f00.cbmsim
    t02 = f02.cbmsim
    t04 = f04.cbmsim
    t06 = f06.cbmsim
    t08 = f08.cbmsim
    t10 = f10.cbmsim
    t12 = f12.cbmsim
    t14 = f14.cbmsim
    t16 = f16.cbmsim
    t18 = f18.cbmsim
    trees = [t00,t02,t04,t06,t08,t10,t12,t14,t16,t18]
    hists = [h00,h02,h04,h06,h08,h10,h12,h14,h16,h18]
    fields = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
    texts = [0,0,0,0,0,0,0,0,0,0]
    xprojections = [0,0,0,0,0,0,0,0,0,0]
    yprojections = [0,0,0,0,0,0,0,0,0,0]
    xtexts = [0,0,0,0,0,0,0,0,0,0]
    ytexts = [0,0,0,0,0,0,0,0,0,0]


for hist in hists:
    hist.SetStats(False)

# initialising histograms
n=0
for tree in trees:
    counter=0
    items = list(range(0, tree.GetEntries()))
    l = len(items)
    print '\n' "Analysing tree ", n+1, "    "
    muon = 0
    amuon = 0
    electron = 0
    aelectron = 0
    j=0
    for event in tree:
        j=j+1
        for hit in event.strawtubesPoint:
            x = hit.GetX()
            y = hit.GetY()
            z = hit.GetZ()
            pz = hit.GetPz()
            for track in event.MCTrack:
                w = track.GetWeight()
                pdg = track.GetPdgCode()
                if z<2581.5 and z>2580:
                    PDG.Fill(pdg,w)
                    if pdg ==13:
                        muon=muon+1
                    if pdg ==-13:
                        amuon=amuon+1
                    if pdg ==11:
                        electron=electron+1
                    if pdg ==-11:
                        aelectron=aelectron+1
                    hists[n].Fill(x,y,w*17786274/tree.GetEntries())
                if n==0:
                    if z<2598 and z>2580:
                        z00.Fill(z,w*17786274/tree.GetEntries())
        print_progress(j + 1, l, prefix = 'Progress:', suffix = 'Complete')
    texts[n] = r.TLatex(.13,.85,"entries: {:.2g}, #mu-: {}, #mu+: {}, e-:{}, e+:{}".format(hists[n].Integral(),muon,amuon,electron,aelectron))
    texts[n].SetTextSize(0.035)
    texts[n].SetNDC(r.kTRUE)
    xprojections[n] = hists[n].ProjectionX()
    xprojections[n].SetTitle('T1 x-projection')
    xprojections[n].SetStats(False)
    meanx = xprojections[n].GetMean()
    xtexts[n] = r.TLatex(.13, .85,"X for B = {} mean: {:.4g}".format(fields[n],meanx))
    xtexts[n].SetNDC(r.kTRUE)
    yprojections[n] = hists[n].ProjectionY()
    yprojections[n].SetTitle('T1 y-projection')
    yprojections[n].SetStats(False)
    meany = yprojections[n].GetMean()
    ytexts[n] = r.TLatex(.13, .85,"Y for B = {} mean: {:.4g}".format(fields[n],meany))
    ytexts[n].SetNDC(r.kTRUE)
    n=n+1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Combined
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
num = np.rint(len(trees)/2.)
nom = num*650

h = 0
for hist in hists:
    cans = r.TCanvas('cans', 'cans', 950,650)
    cans.cd(1)
    hist.Draw('COLZ')
    hist.GetXaxis().SetTitle('x-Point in strawtubes')
    hist.GetYaxis().SetTitle('y-Point')
    texts[h].Draw()
    cans.Draw()
    cans.SaveAs('hists/strawtubes/{}/xy_{}.pdf'.format(directory,int(fields[h]*10)))
    cans.Close()
    h = h+1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           z distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c_z = r.TCanvas('c_z', 'c_z', 950,650)

c_z.cd(1)
z00.Draw('HIST BAR')
z00.GetXaxis().SetTitle('z positions around T1')
z00.GetYaxis().SetTitle('hits/s')
c_z.Draw()
c_z.SaveAs('hists/strawtubes/{}/z_T1.pdf'.format(directory))
c_z.Close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           particle distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

t_mu = r.TLatex(13, 100e6,"#mu^{-}")
t_mu.SetTextSize(0.04)
t_mu.SetTextColor(r.kBlack)

t_amu = r.TLatex(-13, 100e6,"#mu^{+}")
t_amu.SetTextSize(0.04)
t_amu.SetTextColor(r.kBlack)

t_el = r.TLatex(9.5, 500e6,"e^{-}")
t_el.SetTextSize(0.04)
t_el.SetTextColor(r.kBlack)

t_ael = r.TLatex(-11, 500e6,"e^{+}")
t_ael.SetTextSize(0.04)
t_ael.SetTextColor(r.kBlack)

t_gam = r.TLatex(24, 500e6,"#gamma")
t_gam.SetTextSize(0.04)
t_gam.SetTextColor(r.kBlack)

c_pdg = r.TCanvas('c_pdg', 'c_pdg', 950,650)

c_pdg.cd(1)
PDG.Draw('HIST BAR')
PDG.GetXaxis().SetTitle('particle ID at T1')
t_mu.Draw()
t_amu.Draw()
t_el.Draw()
t_ael.Draw()
t_gam.Draw()
c_pdg.Draw()
c_pdg.SaveAs('hists/strawtubes/{}/pdg_T1.pdf'.format(directory))
c_pdg.Close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         X/Y projections
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for h in range(len(xprojections)):
    cproj = r.TCanvas('cproj', 'cproj', 2*950,650)
    cproj.Divide(2,1)
    cproj.cd(1)
    xprojections[h].Draw('HIST BAR')
    xprojections[h].GetXaxis().SetTitle('x')
    #xtexts[h].Draw()
    cproj.cd(2)
    yprojections[h].Draw('HIST BAR')
    yprojections[h].GetXaxis().SetTitle('y')
    #ytexts[h].Draw()
    cproj.Draw()
    cproj.SaveAs('hists/strawtubes/{}/proj_{}.pdf'.format(directory,int(fields[h]*10)))
    cproj.Close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analysis of magnetic field
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcan = r.TCanvas('bcan', 'bcan', 950,650)

text = r.TLatex(0.2,0.2,"Reconstructed from 100.000 events samples")
text.SetNDC(r.kTRUE)
text.SetTextSize(0.03)

bcan.cd(1)
bfield = r.TGraph()
for n in range(len(fields)):
    bfield.SetPoint(n,fields[n],hists[n].Integral()+1)
bfield.SetTitle('Number of total hits in first plane of T1 for different B-fields')
bfield.SetMarkerStyle(5)
bfield.SetMarkerSize(3)
bfield.SetMarkerColorAlpha(r.kGreen+1,1.0)
r.gPad.SetLogy()
bfield.SetMaximum(4*10e10)
bfield.SetMinimum(1)
bfield.GetXaxis().SetTitle('Magnetic field B [T]')
bfield.GetYaxis().SetTitle('hits/s')
bfield.Draw('AP')
#text.Draw()
bcan.Draw()
bcan.SaveAs('hists/strawtubes/{}/hits_bfield.pdf'.format(directory))
bcan.Close()
