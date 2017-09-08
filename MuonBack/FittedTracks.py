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
f00 = r.TFile("/eos/user/k/ksedlacz/files/no-mag_rec.root", "READ")
f02 = r.TFile("/eos/user/k/ksedlacz/files/02-100000-MuonBack_rec.root", "READ")
f04 = r.TFile("/eos/user/k/ksedlacz/files/04-100000-MuonBack_rec.root", "READ")
f06 = r.TFile("/eos/user/k/ksedlacz/files//06-100000-MuonBack_rec.root", "READ")
f10 = r.TFile("/eos/user/k/ksedlacz/files/10-100000-MuonBack_rec.root", "READ")
f12 = r.TFile("/eos/user/k/ksedlacz/files/12-100000-MuonBack_rec.root", "READ")
t00 = f00.cbmsim
t02 = f02.cbmsim
t04 = f04.cbmsim
t06 = f06.cbmsim
t10 = f10.cbmsim
t12 = f12.cbmsim
trees = [t00,t02,t04,t06,t10,t12]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# initialising histograms
h00 = r.TH2D('h00','x:y B-field=0.0 (reco+weights) T1',80,-300,300,160,-600,600)
h00.SetStats(False)
h02 = r.TH2D('h02','x:y B-field=0.2 (reco+weights) T1',80,-300,300,160,-600,600)
h02.SetStats(False)
h04 = r.TH2D('h04','x:y B-field=0.4 (reco+weights) T1',80,-300,300,160,-600,600)
h04.SetStats(False)
h06 = r.TH2D('h06','x:y B-field=0.6 (reco+weights) T1',80,-300,300,160,-600,600)
h06.SetStats(False)
h10 = r.TH2D('h10','x:y B-field=1.0 (reco+weights) T1',80,-300,300,160,-600,600)
h10.SetStats(False)
h12 = r.TH2D('h12','x:y B-field=1.2 (reco+weights) T1',80,-300,300,160,-600,600)
h12.SetStats(False)

z00 = r.TH1D('z00','z B-field=0.0 (reco)',250,2580.5,2584)
z00.SetStats(False)
z02 = r.TH1D('z02','z B-field=0.2 (reco)',200,-2200,3800)
z02.SetStats(False)
z04 = r.TH1D('z04','z B-field=0.4 (reco)',200,-2200,3800)
z04.SetStats(False)
z06 = r.TH1D('z06','z B-field=0.6 (reco)',200,-2200,3800)
z06.SetStats(False)
z10 = r.TH1D('z10','z B-field=1.0 (reco)',200,-2200,3800)
z10.SetStats(False)
z12 = r.TH1D('z12','z B-field=1.2 (reco)',200,-2200,3800)
z12.SetStats(False)

pz00 = r.TH1D('pz00','p_z B-field=0.0 (reco)',100,0,300)
pz00.SetStats(False)
pz02 = r.TH1D('pz02','p_z B-field=0.2 (reco)',100,0,300)
pz02.SetStats(False)
pz04 = r.TH1D('pz04','p_z B-field=0.4 (reco)',100,0,300)
pz04.SetStats(False)
pz06 = r.TH1D('pz06','p_z B-field=0.6 (reco)',100,0,300)
pz06.SetStats(False)
pz10 = r.TH1D('pz10','p_z B-field=1.0 (reco)',100,0,300)
pz10.SetStats(False)
pz12 = r.TH1D('pz12','p_z B-field=1.2 (reco)',100,0,300)
pz12.SetStats(False)

target_xy = r.TH2D('traget_xy','x:y target (reco/only muon shield off)',80,-300,300,160,-600,600)
target_xy.SetStats(False)
target_dist = r.TH1D('target_dist','dist to target (reco/only muon shield off)',200,0,400)
target_dist.SetStats(False)

weight = r.TH1D('weights','mc weights',200,2000,5500)
#weight = r.TH1D('weights','event weights',100,-10,10)

loweight = r.TH1D('loweight','p_{t} for mc weights low',80,0,10)
hiweight = r.TH1D('hiweight','p_{t} for mc weights high',80,0,10)
loweight.SetStats(False)
hiweight.SetStats(False)

hists = [h00,h02,h04,h06,h10,h12]
zhists = [z00,z02,z04,z06,z10,z12]
pzhists = [pz00,pz02,pz04,pz06,pz10,pz12]
texts = [0,0,0,0,0,0]

n=0
target_z = -7067.0
for tree in trees:
    counter=0
    items = list(range(0, tree.GetEntries()))
    l = len(items)
    print '\n' "Analysing tree ", n+1, "   "
    muon=0
    amuon=0
    j=0
    for event in tree:
        j=j+1
        w1 = event.MCTrack[1].GetWeight()
        evw = event.GetWeight()
        for track in event.FitTracks:
            for mctrack in event.MCTrack:
                w = mctrack.GetWeight()
            state = track.getFittedState()
            mom = state.getMom()
            gmom = np.sqrt(mom[0]**2+mom[1]**2)
            pos = state.getPos()
            pdg = state.getPDG()
            if pdg==13:
                muon=muon+1
            elif pdg==-13:
                amuon=amuon+1
            x = pos[0]
            y = pos[1]
            z = pos[2]
            pzhists[n].Fill(mom[2],w)
            zhists[n].Fill(z,w)
            if n==1:
                frac = (z-target_z)/mom[2]
                xtarget = x-frac*mom[0]
                ytarget = y-frac*mom[1]
                dist = np.sqrt(xtarget**2+ytarget**2)
                target_xy.Fill(xtarget,ytarget)
                target_dist.Fill(dist)
            if 2579.<= z <=2584.:
                hists[n].Fill(x,y,w)
            weight.Fill(w1)
            if w<3000: loweight.Fill(gmom)
            if w>3001: hiweight.Fill(gmom)
        print_progress(j + 1, l, prefix = 'Progress:', suffix = 'Complete')
    texts[n] = r.TLatex(-280,530,"entries: {}, #mu-: {}, #mu+: {}".format(int(hists[n].Integral()),muon,amuon))
    texts[n].SetTextSize(0.05)
    n=n+1

cweigh = r.TCanvas('cweigh', 'cweigh', 950,650)
cweigh.cd(1)
weight.GetXaxis().SetTitle('weights')
weight.Draw('BAR')
cweigh.Draw()
cweigh.SaveAs("hists/mc_weights.pdf")
cweigh.Close()

ctarg = r.TCanvas('ctarg', 'ctarg', 950,2*650)
ctarg.Divide(1,2)
ctarg.cd(1)
target_xy.GetXaxis().SetTitle('x')
target_xy.GetYaxis().SetTitle('y')
target_xy.Draw('COLZ')
ctarg.cd(2)
target_dist.GetXaxis().SetTitle('dist')
target_dist.Draw('HIST BAR')
ctarg.Draw()
ctarg.SaveAs("hists/target_dist.pdf")
ctarg.Close()

c01 = r.TCanvas('c01', 'c01', 2*950,3*650)
c01.Divide(2,3)

i = 1
for hist in hists:
    c01.cd(i)
    hist.GetXaxis().SetTitle('x')
    hist.GetYaxis().SetTitle('y')
    hist.Draw('COLZ')
    texts[i-1].Draw()
    i = i+1
c01.Draw()
c01.SaveAs("hists/fittracks.pdf")
c01.Close()

#----------------------------------------------------------------------------------

cproj = r.TCanvas('cproj', 'cproj', 2*950,6*650)
cproj.Divide(2,6)

bfields = [0.0,0.2,0.4,0.6,1.0,1.2]
infox = [0,0,0,0,0,0,0,0,0,0,0,0]
infoy = [0,0,0,0,0,0,0,0,0,0,0,0]

m = 1
k = 0
for hist in hists:
    X = hist.ProjectionX()
    X.SetTitle('x projection {}'.format(bfields[k]))
    X.SetStats(False)
    meanx = X.GetMean()
    infox[k] = r.TLatex(.2, .75,"mean: {:.4g}".format(meanx))
    infox[k].SetNDC(r.kTRUE)
    Y = hist.ProjectionY()
    Y.SetTitle('y projection {}'.format(bfields[k]))
    Y.SetStats(False)
    meany = Y.GetMean()
    infoy[k] = r.TLatex(.2, .75,"mean: {:.4g}".format(meany))
    infoy[k].SetNDC(r.kTRUE)
    cproj.cd(m)
    X.Draw('BAR')
    infox[k].Draw()
    m = m+1
    cproj.cd(m)
    Y.Draw('BAR')
    infoy[k].Draw()
    cproj.Update()
    m = m+1
    k = k+1
cproj.Draw()
cproj.SaveAs("hists/fittracks_proj.pdf")
cproj.Close()

#----------------------------------------------------------------------------------

cpz = r.TCanvas('cpz', 'cpz', 2*950,3*650)
cpz.Divide(2,3)

i = 1
for hist in pzhists:
    cpz.cd(i)
    hist.GetXaxis().SetTitle('p_{z} [GeV]')
    hist.Draw('HIST BAR')
    i = i+1
cpz.Draw()
cpz.SaveAs("hists/fittracks_pz.pdf")
cpz.Close()

cp = r.TCanvas('cp', 'cp', 2*950,650)
cp.Divide(2,1)

cp.cd(1)
loweight.GetXaxis().SetTitle('p [GeV]')
loweight.Draw('HIST BAR')
cp.cd(2)
hiweight.GetXaxis().SetTitle('p [GeV]')
hiweight.Draw('HIST BAR')
cp.Draw()
cp.SaveAs("hists/fittracks_pforweights.pdf")
cp.Close()


#----------------------------------------------------------------------------------

tz1 = r.TLatex(.15, .8,"Veto station")
tz2 = r.TLatex(.5, .85,"T1 entries: {:.7g}".format(z00.Integral()))
tz1.SetTextSize(0.05)
tz2.SetTextSize(0.05)
tz1.SetNDC(r.kTRUE)
tz2.SetNDC(r.kTRUE)

c01z = r.TCanvas('c01z', 'c01z', 950,650)

c01z.cd(1)
z00.GetXaxis().SetTitle('z')
z00.Draw('HIST BAR')
#tz1.Draw()
tz2.Draw()
c01z.Draw()
c01z.SaveAs("hists/nofield/z_500k.pdf")
c01z.Close()


'''
c01z = r.TCanvas('c01z', 'c01z', 2*950,3*650)
c01z.Divide(2,3)

n = 1
for zhist in zhists:
    c01z.cd(n)
    zhist.GetXaxis().SetTitle('z')
    zhist.Draw('BAR')
    zhist.SetFillColor(r.kBlue)
    tz1.Draw()
    tz2.Draw()
    c01z.Update()
    n=n+1
c01z.Draw()
c01z.SaveAs("hists/nofield/fittracks_z.pdf")
c01z.Close()'''
