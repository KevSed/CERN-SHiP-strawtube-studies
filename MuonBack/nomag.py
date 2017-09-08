# -*- coding: utf-8 -*-
from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys
import re

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

# Calculate impact parameter at target
def IPtoTarget(vtx, mom):
    p = np.sqrt(mom.Px()**2+mom.Py()**2+mom.Pz()**2)
    target = r.TVector3(0., 0., -7067.0)
    delta = 0.
    for i in range(3):
        delta += (target(i) - vtx(i)) * mom(i)/p
    ip = 0.
    for i in range(3):
        ip += (target(i) - vtx(i) - delta*mom(i)/p)**2.
    return r.TMath.Sqrt(ip)

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)

# read the rootfile
if '-old' in sys.argv:
    f00 = r.TFile("/eos/user/k/ksedlacz/files/00-100000-MuonBack_rec.root", "READ")
    directory='old'
    t00 = f00.cbmsim
if '-smallFieldOld' in sys.argv:
    f10 = r.TFile('/eos/user/k/ksedlacz/files/10mT-100000_rec.root', 'READ')
    directory='10mT'
    t00 = f10.cbmsim
if '-smallFieldNew' in sys.argv:
    f10 = r.TFile('/eos/user/k/ksedlacz/files/10mT-100000-MuonBack-new_rec.root', 'READ')
    directory='new/10mT'
    t00 = f10.cbmsim
if len(sys.argv)==1:
    f00 = r.TFile("/eos/user/k/ksedlacz/files/00-100000-MuonBack-new_rec.root", "READ")
    t00 = f00.cbmsim
    directory='new'
'''
geo = r.TFile("files/geofile_500k.root","read")
sGeo = geo.FAIRGeom
fGeo = r.gGeoManager'''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# initialising histograms

target_xy = r.TH2D('traget_xy','x:y at target (#mu^{-} and #mu^{+} reco) ',160,-600,600,160,-600,600)
target_xy.SetStats(False)

target_xy_mu = r.TH2D('traget_xy_mu','x:y at target (#mu^{-} reco) ',160,-600,600,160,-600,600)
target_xy_mu.SetStats(False)

target_xy_amu = r.TH2D('traget_xy_amu','x:y at target (#mu^{+} reco) ',160,-600,600,160,-600,600)
target_xy_amu.SetStats(False)

target_dist = r.TH1D('target_dist','#mu^{-} and #mu^{+} dist to target (reco) ',200,0,400)
target_dist.SetStats(False)

target_dist_mu = r.TH1D('target_dist_mu','#mu^{-} dist to target (reco) ',200,0,400)
target_dist_mu.SetStats(False)

target_dist_amu = r.TH1D('target_dist_amu','#mu^{+} dist to target (reco) ',200,0,400)
target_dist_amu.SetStats(False)

mc_trgt_xy_mu = r.TH2D('mc_trgt_xy_mu','x:y at target (#mu^{-} mc) ',160,-600,600,160,-600,600)
mc_trgt_xy_mu.SetStats(False)

mc_trgt_dist_mu = r.TH1D('mc_trgt_dist_mu','#mu^{-} dist to target (mc) ',100,0,100)
mc_trgt_dist_mu.SetStats(False)

mc_trgt_xy_amu = r.TH2D('mc_trgt_xy_amu','x:y at target (#mu^{+} mc) ',160,-600,600,160,-600,600)
mc_trgt_xy_amu.SetStats(False)

mc_trgt_dist_amu = r.TH1D('mc_trgt_dist_amu','#mu^{+} dist to target (mc) ',100,0,100)
mc_trgt_dist_amu.SetStats(False)

mc_pvtx_z = r.TH1D('mc_pvtx_z','z distr. of production vertices (mc)',200,-8000,5000)
mc_pvtx_z.SetStats(False)

rec_mom = r.TH1D('rec_mom','P distr. of tracks (reco)',200,0,200)
rec_mom.SetStats(False)

PDG = r.TH1D('PDG','Particle distribution (mc)',40,-20,20)
PDG.SetStats(False)

rPDG = r.TH1D('rPDG','Particle distribution (reco)',40,-20,20)
rPDG.SetStats(False)

targets = ['Target_1_1','Target_2_1','Target_3_1','Target_4_1','Target_5_1','Target_6_1','Target_7_1','Target_8_1','Target_9_1','Target_10_1','Target_11_1','Target_12_1','Target_13_1','Target_14_1','Target_15_1','Target_16_1','Target_17_1','Target_18_1','Target_19_1','Target_20_1','Target_1_1','Target_21_1','Target_22_1','Target_23_1']

#trees = [t00,t10]












n=0
strg = 'TargetArea'
items = list(range(0, t00.GetEntries()))
l = len(items)
print '\n' "Analysing tree ", n+1, "   "
muon=0
mc_muon=0
amuon=0
mc_amuon=0
j=0
target_z = -7067.0
for event in t00:
    j=j+1
    w = event.MCTrack[1].GetWeight()
    for candidate in event.MCTrack:
        w = candidate.GetWeight()
        pdg = candidate.GetPdgCode()
        x = candidate.GetStartX()
        y = candidate.GetStartY()
        z = candidate.GetStartZ()
        px = candidate.GetPx()
        py = candidate.GetPy()
        pz = candidate.GetPz()
        frac = (target_z-z)/pz
        if z<-7003.5 and z>-7130.5:
            xtarget = x#+frac*px
            ytarget = y#+frac*py
            dist = np.sqrt(xtarget**2+ytarget**2)
        PDG.Fill(pdg)
        if pdg==13:
            mc_muon=mc_muon+1
            mc_trgt_xy_mu.Fill(xtarget,ytarget,w)
            mc_trgt_dist_mu.Fill(dist,w)
        if pdg==-13:
            mc_amuon=mc_amuon+1
            mc_trgt_xy_amu.Fill(xtarget,ytarget,w)
            mc_trgt_dist_amu.Fill(dist,w)
        mc_pvtx_z.Fill(z,w)
        #part = fGeo.FindNode(mctrack.GetStartX(),mctrack.GetStartY(),mctrack.GetStartZ()).GetName()
        #if part in targets:
    for track in event.FitTracks:
        if track.getFitStatus().getTrackLen() >0:
            state = track.getFittedState()
            mom = state.getMom()
            pos = state.getPos()
            pdg = state.getPDG()
            Dir = state.getDir()
            #if pdg==13:
            #elif pdg==-13:
            x = pos[0]
            y = pos[1]
            z = pos[2]
            frac = (-z+target_z)/mom[2]
            xtarget = x+frac*mom[0]
            ytarget = y+frac*mom[1]
            dist = np.sqrt(xtarget**2+ytarget**2)
            MOM = np.sqrt(mom[0]**2+mom[1]**2+mom[2]**2)
            rec_mom.Fill(MOM)
            rPDG.Fill(pdg)
            if pdg==13:
                target_xy_mu.Fill(xtarget,ytarget)
                target_dist_mu.Fill(dist)
                muon=muon+1
            if pdg==-13:
                target_xy_amu.Fill(xtarget,ytarget)
                target_dist_amu.Fill(dist)
                amuon=amuon+1
            target_xy.Fill(xtarget,ytarget)
            target_dist.Fill(dist)
        print_progress(j + 1, l, prefix = 'Progress:', suffix = 'Complete')

comp = PDG.Integral()
rcomp = rPDG.Integral()






##############################################################################
#------------------ MONTE CARLO ----------------------------------------------
##############################################################################



flag = r.TLatex(.66, .85,"{}".format(directory))
flag.SetNDC(r.kTRUE)
flag.SetTextSize(0.04)

target1 = r.TLine(-7130.5,-1000,-7130.5,200e3)
target1.SetLineColor(r.kGreen)
target2 = r.TLine(-7003.5,-1000,-7003.5,200e3)
target2.SetLineColor(r.kGreen)
text_t = r.TLatex(.18, .85,"target")
text_t.SetNDC(r.kTRUE)
text_t.SetTextSize(0.04)
text_t.SetTextColor(r.kGreen)

shield1 = r.TLine(-1715.,-1000,-1715.,200e3)
shield1.SetLineColor(r.kBlue)
shield2 = r.TLine(1715.,-1000,1715.,200e3)
shield2.SetLineColor(r.kBlue)
text_s = r.TLatex(.6, .85,"shield")
text_s.SetNDC(r.kTRUE)
text_s.SetTextSize(0.04)
text_s.SetTextColor(r.kBlue)


decay1 = r.TLine(-5021.6,-1000,-5021.6,200e3)
decay1.SetLineColor(r.kRed)
decay2 = r.TLine(60.,-1000,60.,200e3)
decay2.SetLineColor(r.kRed)
text_d = r.TLatex(.33, .85,"decay volume")
text_d.SetNDC(r.kTRUE)
text_d.SetTextSize(0.03)
text_d.SetTextColor(r.kRed)

T1 = r.TLine(2580,-1000,2590,200e3)
T1.SetLineColor(r.kRed)
text_1 = r.TLatex(.8, .85,"T1")
text_1.SetNDC(r.kTRUE)
text_1.SetTextSize(0.04)
text_1.SetTextColor(r.kRed)

myon = 33/40
antimyon =7/40

t_mu = r.TLatex(13, 435e3,"#mu^{-}")

t_mu.SetTextSize(0.04)
t_mu.SetTextColor(r.kBlack)

t_amu = r.TLatex(-13, 140e3,"#mu^{+}")

t_amu.SetTextSize(0.04)
t_amu.SetTextColor(r.kBlack)

t_el = r.TLatex(11, 435e3,"e^{-}")

t_el.SetTextSize(0.04)
t_el.SetTextColor(r.kBlack)

t_ael = r.TLatex(-11, 140e3,"e^{+}")

t_ael.SetTextSize(0.04)
t_ael.SetTextColor(r.kBlack)

t_rel = r.TLatex(.13, .85,"e-: {:.3g}%  e+: {:.3g}%  mu-: {:.3g}%  mu+: {:.3g}% ".format((PDG.Integral(31,33)/comp)*100,(PDG.Integral(9,11)/comp)*100,(PDG.Integral(33,35)/comp)*100,(PDG.Integral(7,9)/comp)*100))
t_rel.SetNDC(r.kTRUE)
t_rel.SetTextSize(0.04)
t_rel.SetTextColor(r.kBlack)

t_rel2 = r.TLatex(.13, .75,"e-/e+: {:.3g}  mu-/mu+: {:.3g}".format(PDG.Integral(31,33)/PDG.Integral(9,11),PDG.Integral(33,35)/PDG.Integral(7,9)))
t_rel2.SetNDC(r.kTRUE)
t_rel2.SetTextSize(0.04)
t_rel2.SetTextColor(r.kBlack)


'''
shield1 = r.TLine(-1715.,-1000,-1715.,200e3)
shield1.SetLineColor(r.kBlue)
shield2 = r.TLine(1715.,-1000,1715.,200e3)
shield2.SetLineColor(r.kBlue)
text_s = r.TLatex(.55, .85,"shield")
text_s.SetNDC(r.kTRUE)
text_s.SetTextSize(0.035)'''

#----------------- PARTICLE DISTRIBUTION --------------------------------
pdg_canv = r.TCanvas('pdg_canv','pdg_canv',950,650)
pdg_canv.cd(1)
PDG.GetXaxis().SetTitle('pdg code')
PDG.Draw('HIST BAR')
t_mu.Draw()
t_amu.Draw()
t_el.Draw()
t_ael.Draw()
t_rel.Draw()
t_rel2.Draw()
pdg_canv.Draw()
pdg_canv.SaveAs('hists/nofield/{}/mc_pdg.pdf'.format(directory))
pdg_canv.Close()


#----------------- PRODUCTION VERTEX Z --------------------------------
z_canv = r.TCanvas('z_canv','z_canv',950,650)
z_canv.cd(1)
mc_pvtx_z.GetXaxis().SetTitle('z of production vertex')
mc_pvtx_z.Draw('HIST BAR')
target1.Draw()
target2.Draw()
text_t.Draw()
T1.Draw()
text_1.Draw()
z_canv.Draw()
z_canv.SaveAs('hists/nofield/{}/mc_pvtx_z.pdf'.format(directory))
z_canv.Close()


#----------------- MC PROJECTIONS MU ----------------------------------
mc_trgt_x_mu = mc_trgt_xy_mu.ProjectionX()
mc_trgt_x_mu.SetTitle('x projection of #mu^{-} positions at target')
mc_trgt_x_mu.SetFillColor(r.kGreen+2)
mc_trgt_y_mu = mc_trgt_xy_mu.ProjectionY()
mc_trgt_y_mu.SetTitle('y projection of #mu^{-} positions at target')
mc_trgt_y_mu.SetFillColor(r.kBlue+2)
mc_trgt_x_mu.SetStats(False)
mc_trgt_y_mu.SetStats(False)

#----------------- MC PROJECTIONS AMU ----------------------------------
mc_trgt_x_amu = mc_trgt_xy_amu.ProjectionX()
mc_trgt_x_amu.SetTitle('x projection of #mu^{+} positions at target')
mc_trgt_x_amu.SetFillColor(r.kGreen+2)
mc_trgt_y_amu = mc_trgt_xy_amu.ProjectionY()
mc_trgt_y_amu.SetTitle('y projection of #mu^{+} positions at target')
mc_trgt_y_amu.SetFillColor(r.kBlue+2)
mc_trgt_x_amu.SetStats(False)
mc_trgt_y_amu.SetStats(False)

#----------------- TEXTS MU -----------------------------------------
xtext_mu = r.TLatex(.13, .85,"RMS = {:.4g}cm                                       mean = {:.4g}cm #pm {:.4g}cm".format(mc_trgt_x_mu.GetRMS(),mc_trgt_x_mu.GetMean(),mc_trgt_x_mu.GetMeanError()))
xtext_mu.SetNDC(r.kTRUE)
xtext_mu.SetTextSize(0.035)
ytext_mu = r.TLatex(.13, .85,"RMS = {:.4g}cm                                       mean = {:.4g}cm #pm {:.4g}cm".format(mc_trgt_y_mu.GetRMS(),mc_trgt_y_mu.GetMean(),mc_trgt_y_mu.GetMeanError()))
ytext_mu.SetNDC(r.kTRUE)
ytext_mu.SetTextSize(0.035)
dtext2_mu = r.TLatex(.5, .65,"<1cm: {:.4g}%".format(mc_trgt_dist_mu.Integral(0,1)/mc_trgt_dist_mu.Integral()*100))
dtext2_mu.SetNDC(r.kTRUE)
mc_xytext_mu = r.TLatex(.13,.85,"Number of tracks: {}".format(mc_trgt_xy_mu.Integral()))
mc_xytext_mu.SetNDC(r.kTRUE)
mc_xytext_mu.SetTextSize(0.035)

#----------------- TEXTS AMU -----------------------------------------
xtext_amu = r.TLatex(.13, .85,"RMS = {:.4g}cm                                       mean = {:.4g}cm #pm {:.4g}cm".format(mc_trgt_x_amu.GetRMS(),mc_trgt_x_amu.GetMean(),mc_trgt_x_amu.GetMeanError()))
xtext_amu.SetNDC(r.kTRUE)
xtext_amu.SetTextSize(0.035)
ytext_amu = r.TLatex(.13, .85,"RMS = {:.4g}cm                                       mean = {:.4g}cm #pm {:.4g}cm".format(mc_trgt_y_amu.GetRMS(),mc_trgt_y_amu.GetMean(),mc_trgt_y_amu.GetMeanError()))
ytext_amu.SetNDC(r.kTRUE)
ytext_amu.SetTextSize(0.035)
dtext = r.TLatex(.5, .75,"distance to target centre")
dtext2_amu = r.TLatex(.5, .65,"<1cm: {:.4g}%".format(mc_trgt_dist_amu.Integral(0,1)/mc_trgt_dist_amu.Integral()*100))
dtext.SetNDC(r.kTRUE)
dtext2_amu.SetNDC(r.kTRUE)
mc_xytext_amu = r.TLatex(.13,.85,"Number of tracks: {}".format(mc_trgt_xy_amu.Integral()))
mc_xytext_amu.SetNDC(r.kTRUE)
mc_xytext_amu.SetTextSize(0.035)

#----------------------- DRAW MUON MONTE CARLO -----------------------------------
ctrgt_mu = r.TCanvas('ctrgt_mu', 'ctrgt_mu', 2*950,2*650)
ctrgt_mu.Divide(2,2)

el = r.TEllipse(0.,0.,15.,15.)
el.SetFillColor(0)

ctrgt_mu.cd(1)
mc_trgt_xy_mu.GetXaxis().SetTitle('x')
mc_trgt_xy_mu.GetYaxis().SetTitle('y')
mc_trgt_xy_mu.Draw('COLZ')
mc_xytext_mu.Draw()
flag.Draw()
#el.Draw()

ctrgt_mu.cd(2)
mc_trgt_dist_mu.GetXaxis().SetTitle('dist [cm]')
mc_trgt_dist_mu.Draw('HIST BAR')
dtext.Draw()
dtext2_mu.Draw()

ctrgt_mu.cd(3)
mc_trgt_x_mu.GetXaxis().SetTitle('x [cm]')
mc_trgt_x_mu.Draw('HIST BAR')
xtext_mu.Draw()

ctrgt_mu.cd(4)
mc_trgt_y_mu.GetXaxis().SetTitle('y [cm]')
mc_trgt_y_mu.Draw('HIST BAR')
ytext_mu.Draw()

ctrgt_mu.Draw()
ctrgt_mu.SaveAs("hists/nofield/{}/mc_target_dist_mu.pdf".format(directory))
ctrgt_mu.Close()

#----------------------- DRAW AMUON MONTE CARLO -----------------------------------
ctrgt_amu = r.TCanvas('ctrgt_amu', 'ctrgt_amu', 2*950,2*650)
ctrgt_amu.Divide(2,2)

ctrgt_amu.cd(1)
mc_trgt_xy_amu.GetXaxis().SetTitle('x')
mc_trgt_xy_amu.GetYaxis().SetTitle('y')
mc_trgt_xy_amu.Draw('COLZ')
mc_xytext_amu.Draw()
flag.Draw()
#el.Draw()

ctrgt_amu.cd(2)
mc_trgt_dist_amu.GetXaxis().SetTitle('dist [cm]')
mc_trgt_dist_amu.Draw('HIST BAR')
dtext.Draw()
dtext2_amu.Draw()

ctrgt_amu.cd(3)
mc_trgt_x_amu.GetXaxis().SetTitle('x [cm]')
mc_trgt_x_amu.Draw('HIST BAR')
xtext_amu.Draw()

ctrgt_amu.cd(4)
mc_trgt_y_amu.GetXaxis().SetTitle('y [cm]')
mc_trgt_y_amu.Draw('HIST BAR')
ytext_amu.Draw()

ctrgt_amu.Draw()
ctrgt_amu.SaveAs("hists/nofield/{}/mc_target_dist_amu.pdf".format(directory))
ctrgt_amu.Close()






##############################################################################
#------------------ RECONSTRUCTED --------------------------------------------
##############################################################################






rt_mu = r.TLatex(13, 1200,"#mu^{-}")

rt_mu.SetTextSize(0.04)
rt_mu.SetTextColor(r.kBlack)

rt_amu = r.TLatex(-13, 1550,"#mu^{+}")

rt_amu.SetTextSize(0.04)
rt_amu.SetTextColor(r.kBlack)

rt_rel = r.TLatex(.35, .85,"e-: {:.3g}%  e+: {:.3g}%  mu-: {:.3g}%  mu+: {:.3g}% ".format((rPDG.Integral(31,33)/rcomp)*100,(rPDG.Integral(9,11)/rcomp)*100,(rPDG.Integral(33,35)/rcomp)*100,(rPDG.Integral(7,9)/rcomp)*100))
rt_rel.SetNDC(r.kTRUE)
rt_rel.SetTextSize(0.04)
rt_rel.SetTextColor(r.kBlack)

rt_rel2 = r.TLatex(.35, .8,"e-/e+: {:.3g}  mu-/mu+: {:.3g}".format(rPDG.Integral(31,33)/(rPDG.Integral(9,11)+1),rPDG.Integral(33,35)/(rPDG.Integral(7,9)+1)))
rt_rel2.SetNDC(r.kTRUE)
rt_rel2.SetTextSize(0.04)

#----------------- PARTICLE DISTRIBUTION --------------------------------
rpdg_canv = r.TCanvas('pdg_canv','pdg_canv',950,650)
rpdg_canv.cd(1)
rPDG.GetXaxis().SetTitle('pdg code')
rPDG.SetFillColor(r.kBlue)
rPDG.Draw('HIST BAR')
rt_mu.Draw()
rt_amu.Draw()
rt_rel.Draw()
rt_rel2.Draw()
rpdg_canv.Draw()
rpdg_canv.SaveAs('hists/nofield/{}/reco_pdg.pdf'.format(directory))
rpdg_canv.Close()

#----------------- PROJECTIONS MU & AMU ----------------------------------
targ_x = target_xy.ProjectionX()
targ_x.SetTitle('x projection of #mu^{-} and #mu^{+} positions at target')
targ_x.SetFillColor(r.kGreen+2)
targ_y = target_xy.ProjectionY()
targ_y.SetTitle('y projection of #mu^{-} and #mu^{+} positions at target')
targ_y.SetFillColor(r.kBlue+2)
targ_x.SetStats(False)
targ_y.SetStats(False)

#----------------- PROJECTIONS MU ----------------------------------
targ_x_mu = target_xy_mu.ProjectionX()
targ_x_mu.SetTitle('x projection of #mu^{-} positions at target')
targ_x_mu.SetFillColor(r.kGreen+2)
targ_y_mu = target_xy_mu.ProjectionY()
targ_y_mu.SetTitle('y projection of #mu^{-} positions at target')
targ_y_mu.SetFillColor(r.kBlue+2)
targ_x_mu.SetStats(False)
targ_y_mu.SetStats(False)

#----------------- PROJECTIONS AMU ----------------------------------
targ_x_amu = target_xy_amu.ProjectionX()
targ_x_amu.SetTitle('x projection of #mu^{+} positions at target')
targ_x_amu.SetFillColor(r.kGreen+2)
targ_y_amu = target_xy_amu.ProjectionY()
targ_y_amu.SetTitle('y projection of #mu^{+} positions at target')
targ_y_amu.SetFillColor(r.kBlue+2)
targ_x_amu.SetStats(False)
targ_y_amu.SetStats(False)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S e t u p   m o d e l
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Declare variables x,mean,sigma with associated name, title, initial value and allowed range
top = 90
bot = -top
x_mu = ROOT.RooRealVar("x_mu","x_mu",bot,top)
x_mu_mean = ROOT.RooRealVar("x_mu_mean","mean in x",0,-30,30)
x_mu_sigma = ROOT.RooRealVar("x_mu_sigma","sigma x",50,1,130)

x_amu = ROOT.RooRealVar("x_amu","x_amu",bot,top)
x_amu_mean = ROOT.RooRealVar("x_amu_mean","mean in x",0,-30,30)
x_amu_sigma = ROOT.RooRealVar("x_amu_sigma","sigma x",50,1,130)

y_mu = ROOT.RooRealVar("y_mu","y_mu",bot,top)
y_mu_mean = ROOT.RooRealVar("y_mu_mean","mean in y",0,-30,30)
y_mu_sigma = ROOT.RooRealVar("y_mu_sigma","sigma y",50,1,130)

y_amu = ROOT.RooRealVar("y_amu","y_amu",bot,top)
y_amu_mean = ROOT.RooRealVar("y_amu_mean","mean in y",0,-30,30)
y_amu_sigma = ROOT.RooRealVar("y_amu_sigma","sigma y",50,1,130)

# Build gaussian p.d.f in terms of x,mean and sigma
x_mu_gauss = ROOT.RooGaussian("x_mu_gauss","gaussian PDF",x_mu,x_mu_mean,x_mu_sigma)
x_amu_gauss = ROOT.RooGaussian("x_amu_gauss","gaussian PDF",x_amu,x_amu_mean,x_amu_sigma)
y_mu_gauss = ROOT.RooGaussian("y_mu_gauss","gaussian PDF",y_mu,y_mu_mean,y_mu_sigma)
y_amu_gauss = ROOT.RooGaussian("y_amu_gauss","gaussian PDF",y_amu,y_amu_mean,y_amu_sigma)


# F i t   m o d e l   t o   d a t a
# -----------------------------

data_x_mu = ROOT.RooDataHist("data_x_mu", "data_x_mu", r.RooArgList(x_mu),r.RooFit.Import(targ_x_mu))
data_x_amu = ROOT.RooDataHist("data_x_amu", "data_x_amu", r.RooArgList(x_amu),r.RooFit.Import(targ_x_amu))
data_y_mu = ROOT.RooDataHist("data_y_mu", "data_y_mu", r.RooArgList(y_mu),r.RooFit.Import(targ_y_mu))
data_y_amu = ROOT.RooDataHist("data_y_amu", "data_y_amu", r.RooArgList(y_amu),r.RooFit.Import(targ_y_amu))


# Make plot of binned dataset showing Poisson error bars (ROOT.RooFit
# default)
frame1 = x_mu.frame(ROOT.RooFit.Title("gauss fit to x projection of reconstructed IP for #mu^{-} at 10mT"))
data_x_mu.plotOn(frame1)

frame2 = x_amu.frame(ROOT.RooFit.Title("gauss fit to x projection of reconstructed IP for #mu^{+} at 10mT"))
data_x_amu.plotOn(frame2)

frame3 = y_mu.frame(ROOT.RooFit.Title("gauss fit to y projection of reconstructed IP for #mu^{-} at 10mT"))
data_y_mu.plotOn(frame3)

frame4 = y_amu.frame(ROOT.RooFit.Title("gauss fit to y projection of reconstructed IP for #mu^{+} at 10mT"))
data_y_amu.plotOn(frame4)

# Fit a Gaussian p.d.f to the data

x_mu_gauss.fitTo(data_x_mu)
x_mu_gauss.plotOn(frame1)
x_amu_gauss.fitTo(data_x_amu)
x_amu_gauss.plotOn(frame2)
y_mu_gauss.fitTo(data_y_mu)
y_mu_gauss.plotOn(frame3)
y_amu_gauss.fitTo(data_y_amu)
y_amu_gauss.plotOn(frame4)

core1 = r.TLatex(.13, .85,"Inner {:.3g}%".format(targ_x_mu.Integral(bot,top)/targ_x_mu.Integral()*100))
core1.SetNDC(r.kTRUE)
core1.SetTextSize(0.044)

core2 = r.TLatex(.13, .85,"Inner {:.3g}%".format(targ_x_amu.Integral(bot,top)/targ_x_amu.Integral()*100))
core2.SetNDC(r.kTRUE)
core2.SetTextSize(0.044)

core3 = r.TLatex(.13, .85,"Inner {:.3g}%".format(targ_y_mu.Integral(bot,top)/targ_y_mu.Integral()*100))
core3.SetNDC(r.kTRUE)
core3.SetTextSize(0.044)

core4 = r.TLatex(.13, .85,"Inner {:.3g}%".format(targ_y_amu.Integral(bot,top)/targ_y_amu.Integral()*100))
core4.SetNDC(r.kTRUE)
core4.SetTextSize(0.044)

gauss_can = r.TCanvas('gauss_can','gauss_can',2*950,2*650)
gauss_can.Divide(2,2)

gauss_can.cd(1)
frame1.Draw()
core1.Draw()

gauss_can.cd(2)
frame2.Draw()
core2.Draw()

gauss_can.cd(3)
frame3.Draw()
core3.Draw()

gauss_can.cd(4)
frame4.Draw()
core4.Draw()

gauss_can.Draw()
gauss_can.SaveAs('hists/nofield/{}/gauss_fit.pdf'.format(directory))
gauss_can.Close()

#----------------- TEXTS MU & AMU -----------------------------------------

xytext = r.TLatex(.13, .85,"muons: {} antimuons: {}".format(muon,amuon))
xytext.SetNDC(r.kTRUE)
xytext.SetTextSize(0.044)
xtext = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_x.GetRMS(),targ_x.GetMean(),targ_x.GetMeanError()))
xtext.SetNDC(r.kTRUE)
xtext.SetTextSize(0.044)
ytext = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_y.GetRMS(),targ_y.GetMean(),targ_y.GetMeanError()))
ytext.SetNDC(r.kTRUE)
ytext.SetTextSize(0.044)
dtext2 = r.TLatex(.5, .65,"<15cm: {:.4g}%".format(target_dist.Integral(0,15)/target_dist.Integral()*100))
dtext2.SetNDC(r.kTRUE)

#----------------- TEXTS MU -----------------------------------------
xtext_mu = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_x_mu.GetRMS(),targ_x_mu.GetMean(),targ_x_mu.GetMeanError()))
xtext_mu.SetNDC(r.kTRUE)
xtext_mu.SetTextSize(0.045)
ytext_mu = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_y_mu.GetRMS(),targ_y_mu.GetMean(),targ_y_mu.GetMeanError()))
ytext_mu.SetNDC(r.kTRUE)
ytext_mu.SetTextSize(0.045)
dtext2_mu = r.TLatex(.5, .65,"<15cm: {:.4g}%".format(target_dist_mu.Integral(0,15)/target_dist_mu.Integral()*100))
dtext2_mu.SetNDC(r.kTRUE)
xytext_mu = r.TLatex(.13,.85,"Number of tracks: {}".format(target_xy_mu.Integral()))
xytext_mu.SetNDC(r.kTRUE)
xytext_mu.SetTextSize(0.045)

#----------------- TEXTS AMU -----------------------------------------
xtext_amu = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_x_amu.GetRMS(),targ_x_amu.GetMean(),targ_x_amu.GetMeanError()))
xtext_amu.SetNDC(r.kTRUE)
xtext_amu.SetTextSize(0.045)
ytext_amu = r.TLatex(.13, .85,"RMS = {:.4g}cm                         mean = {:.4g}cm #pm {:.4g}cm".format(targ_y_amu.GetRMS(),targ_y_amu.GetMean(),targ_y_amu.GetMeanError()))
ytext_amu.SetNDC(r.kTRUE)
ytext_amu.SetTextSize(0.045)
dtext = r.TLatex(.5, .75,"distance to target centre")
dtext2_amu = r.TLatex(.5, .65,"<15cm: {:.4g}%".format(target_dist_amu.Integral(0,15)/target_dist_amu.Integral()*100))
dtext.SetNDC(r.kTRUE)
dtext2_amu.SetNDC(r.kTRUE)
xytext_amu = r.TLatex(.13,.85,"Number of tracks: {}".format(target_xy_amu.Integral()))
xytext_amu.SetNDC(r.kTRUE)
xytext_amu.SetTextSize(0.045)


#----------------- DRAWING MU & AMU -----------------------------------------
ctarg = r.TCanvas('ctarg', 'ctarg', 2*950,2*650)
ctarg.Divide(2,2)
ctarg.cd(1)
target_xy.GetXaxis().SetTitle('x [cm]')
target_xy.GetYaxis().SetTitle('y [cm]')
target_xy.Draw('COLZ')
flag.Draw()
xytext.Draw()
ctarg.cd(2)
target_dist.GetXaxis().SetTitle('dist [cm]')
target_dist.Draw('HIST BAR')
dtext.Draw()
dtext2.Draw()
ctarg.cd(3)
targ_x.GetXaxis().SetTitle('x [cm]')
targ_x.Draw('HIST BAR')
xtext.Draw()
ctarg.cd(4)
targ_y.GetXaxis().SetTitle('y [cm]')
targ_y.Draw('HIST BAR')
ytext.Draw()
ctarg.Draw()
ctarg.SaveAs("hists/nofield/{}/target_dist.pdf".format(directory))
ctarg.Close()

#----------------- DRAWING MU -----------------------------------------
ctarg_mu = r.TCanvas('ctarg_mu', 'ctarg_mu', 2*950,2*650)
ctarg_mu.Divide(2,2)
ctarg_mu.cd(1)
target_xy_mu.GetXaxis().SetTitle('x [cm]')
target_xy_mu.GetYaxis().SetTitle('y [cm]')
target_xy_mu.Draw('COLZ')
flag.Draw()
xytext_mu.Draw()
ctarg_mu.cd(2)
target_dist_mu.GetXaxis().SetTitle('dist [cm]')
target_dist_mu.Draw('HIST BAR')
dtext.Draw()
dtext2_mu.Draw()
ctarg_mu.cd(3)
targ_x_mu.GetXaxis().SetTitle('x [cm]')
targ_x_mu.Draw('HIST BAR')
xtext_mu.Draw()
ctarg_mu.cd(4)
targ_y_mu.GetXaxis().SetTitle('y [cm]')
targ_y_mu.Draw('HIST BAR')
ytext_mu.Draw()
ctarg_mu.Draw()
ctarg_mu.SaveAs("hists/nofield/{}/target_dist_mu.pdf".format(directory))
ctarg_mu.Close()

#----------------- DRAWING AMU -----------------------------------------
ctarg_amu = r.TCanvas('ctarg_amu', 'ctarg_amu', 2*950,2*650)
ctarg_amu.Divide(2,2)
ctarg_amu.cd(1)
target_xy_amu.GetXaxis().SetTitle('x [cm]')
target_xy_amu.GetYaxis().SetTitle('y [cm]')
target_xy_amu.Draw('COLZ')
flag.Draw()
xytext_amu.Draw()
ctarg_amu.cd(2)
target_dist_amu.GetXaxis().SetTitle('dist [cm]')
target_dist_amu.Draw('HIST BAR')
dtext.Draw()
dtext2_amu.Draw()
ctarg_amu.cd(3)
targ_x_amu.GetXaxis().SetTitle('x [cm]')
targ_x_amu.Draw('HIST BAR')
xtext_amu.Draw()
ctarg_amu.cd(4)
targ_y_amu.GetXaxis().SetTitle('y [cm]')
targ_y_amu.Draw('HIST BAR')
ytext_amu.Draw()
ctarg_amu.Draw()
ctarg_amu.SaveAs("hists/nofield/{}/target_dist_amu.pdf".format(directory))
ctarg_amu.Close()

mom_canv = r.TCanvas('mom_canv', 'mom_canv', 950,650)

cut = r.TLine(10,-5,10,130)
cut.SetLineColor(r.kRed)

mom_canv.cd(1)
rec_mom.GetXaxis().SetTitle('P [GeV]')
rec_mom.Draw('HIST BAR')
flag.Draw()
cut.Draw()
mom_canv.Draw()
mom_canv.SaveAs("hists/nofield/{}/mom_reco_tracks.pdf".format(directory))

#----------------------------------------------------------------------------------
'''
cproj = r.TCanvas('cproj', 'cproj', 2*950,650)
cproj.Divide(2,1)

bfields = [0.0]
infox = [0,0]
infoy = [0,0]

m = 1
k = 0
X = h00.ProjectionX()
X.SetTitle('x projection {}'.format(bfields[k]))
X.SetStats(False)
meanx = X.GetMean()
infox[k] = r.TLatex(.2, .75,"mean: {:.4g}".format(meanx))
infox[k].SetNDC(r.kTRUE)
Y = h00.ProjectionY()
Y.SetTitle('y projection {}'.format(bfields[k]))
Y.SetStats(False)
meany = Y.GetMean()
infoy[k] = r.TLatex(.2, .75,"mean: {:.4g}".format(meany))
infoy[k].SetNDC(r.kTRUE)
cproj.cd(m)
X.Draw('HIST BAR')
infox[k].Draw()
m = m+1
cproj.cd(m)
Y.Draw('HIST BAR')
infoy[k].Draw()
cproj.Update()
m = m+1
k = k+1
cproj.Draw()
cproj.SaveAs("hists/nofield/fittracks_proj.pdf")
cproj.Close()

#----------------------------------------------------------------------------------

cpz = r.TCanvas('cpz', 'cpz', 950,650)

cpz.cd(1)
pz00.GetXaxis().SetTitle('p_{z} [GeV]')
pz00.Draw('HIST BAR')
cpz.Draw()
cpz.SaveAs("hists/nofield/fittracks_pz.pdf")
cpz.Close()

cz = r.TCanvas('cz', 'cz', 950,650)

cz.cd(1)
z00.GetXaxis().SetTitle('z')
z00.Draw('HIST BAR')
cz.Draw()
cz.SaveAs("hists/nofield/fittracks_z.pdf")
cz.Close()
'''
