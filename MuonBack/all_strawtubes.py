# -*- coding: utf-8 -*-
from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys


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


# read the rootfile
f00 = r.TFile("files/00-100000-MuonBack.root", "READ")
f01 = r.TFile("files/01-100000-MuonBack.root", "READ")
f02 = r.TFile("files/02-100000-MuonBack.root", "READ")
f04 = r.TFile("files/04-100000-MuonBack.root", "READ")
f06 = r.TFile("files/06-100000-MuonBack.root", "READ")
f10 = r.TFile("files/10-100000-MuonBack.root", "READ")
t00 = f00.cbmsim
t01 = f01.cbmsim
t02 = f02.cbmsim
t04 = f04.cbmsim
t06 = f06.cbmsim
t10 = f10.cbmsim

trees = [t00,t01,t02,t04,t06,t10]
fields = [0.0,0.1,0.2,0.4,0.6,1.0]

h00 = r.TH2D('h00','x:y hits in T1 B-field=0.0 (weights)',160,-600,600,160,-600,600)
h01 = r.TH2D('h01','x:y hits in T1 B-field=0.1 (weights)',160,-600,600,160,-600,600)
h02 = r.TH2D('h02','x:y hits in T1 B-field=0.2 (weights)',160,-600,600,160,-600,600)
h04 = r.TH2D('h04','x:y hits in T1 B-field=0.4 (weights)',160,-600,600,160,-600,600)
h06 = r.TH2D('h06','x:y hits in T1 B-field=0.6 (weights)',160,-600,600,160,-600,600)
h10 = r.TH2D('h10','x:y hits in T1 B-field=1.0 (weights)',160,-600,600,160,-600,600)

hists = [h00,h01,h02,h04,h06,h10]
texts = [0,0,0,0,0,0]
xprojections = [0,0,0,0,0,0]
yprojections = [0,0,0,0,0,0]
xtexts = [0,0,0,0,0,0]
ytexts = [0,0,0,0,0,0]

for hist in hists:
    hist.SetStats(False)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           b-field = 0.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                mu = track.GetPdgCode()
                if z<2600:
                    if mu==13:
                        muon=muon+1
                    if mu==-13:
                        amuon=amuon+1
                    if mu==11:
                        electron=electron+1
                    if mu==-11:
                        aelectron=aelectron+1
                    hists[n].Fill(x,y,w)
        print_progress(j + 1, l, prefix = 'Progress:', suffix = 'Complete')
    texts[n] = r.TLatex(.13,.85,"entries: {}, #mu-: {}, #mu+: {}, e-:{}, e+:{}".format(hists[n].GetEntries(),muon,amuon,electron,aelectron))
    texts[n].SetTextSize(0.04)
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
#           XY PLOTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lc1 = r.TCanvas('can', 'can', 2*950,3*650)
lc1.Divide(2,3)

h = 1
for hist in hists:
    lc1.cd(h)
    hist.Draw('COLZ')
    hist.GetXaxis().SetTitle('x-Point in strawtubes')
    hist.GetYaxis().SetTitle('y-Point')
    texts[h-1].Draw()
    h = h+1
lc1.Draw()
lc1.SaveAs('hists/strawtubesPoint_xy_all6.pdf')
lc1.Close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           PROJECTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cproj = r.TCanvas('cproj', 'cproj', 2*950,6*650)
cproj.Divide(2,6)

h = 0
for xhist in xprojections:
    cproj.cd(2*h+1)
    xhist.Draw('HIST BAR')
    xhist.GetXaxis().SetTitle('x')
    xtexts[h].Draw()
    h = h+1
h = 1
for yhist in yprojections:
    cproj.cd(2*h)
    yhist.Draw('HIST BAR')
    yhist.GetXaxis().SetTitle('y')
    ytexts[h-1].Draw()
    h = h+1
cproj.Draw()
cproj.SaveAs('hists/strawtubesPoint_proj_all6.pdf')
cproj.Close()
