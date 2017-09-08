import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)

# read the rootfile
f02 = r.TFile("files/ship.conical.MuonBack-TGeant4-02-100000.root", "READ")
f10 = r.TFile("files/ship.conical.MuonBack-TGeant4-10-100000.root", "READ")

t = f02.cbmsim
tree2 = f10.cbmsim

h1 = r.TH2D('T1','T1 B-field=0.2 (weights)',80,-300,300,160,-600,600)
h2 = r.TH2D('T2','T2 B-field=0.2 (weights)',80,-300,300,160,-600,600)
h3 = r.TH2D('T3','T3 B-field=0.2 (weights)',80,-300,300,160,-600,600)
h4 = r.TH2D('T4','T4 B-field=0.2 (weights)',80,-300,300,160,-600,600)


for i in range(t.GetEntries()):
    t.GetEvent(i)
    for p in t.strawtubesPoint:
        x = p.GetX()
        y = p.GetY()
        z = p.GetZ()
        for track in t.MCTrack:
            w = track.GetWeight()
            if z>2500 and z<2700:
                h1.Fill(x,y,w)
            if z>2700 and z<3000:
                h2.Fill(x,y,w)
            if z>3000 and z<3400:
                h3.Fill(x,y,w)
            if z>3400:
                h4.Fill(x,y,w)

h1.SetStats(False)
h2.SetStats(False)
h3.SetStats(False)
h4.SetStats(False)

t1 = r.TLatex(0,530,"entries: {}".format(int(h1.Integral())))
t2 = r.TLatex(0,530,"entries: {}".format(int(h2.Integral())))
t3 = r.TLatex(0,530,"entries: {}".format(int(h3.Integral())))
t4 = r.TLatex(0,530,"entries: {}".format(int(h4.Integral())))

lc1 = r.TCanvas('can', 'can', 1900,1300)
lc1.Divide(2,2)

lc1.cd(1)
h1.Draw('COLZ')
h1.GetXaxis().SetTitle('x-Point in strawtubes')
h1.GetYaxis().SetTitle('y-Point')
t1.Draw()

lc1.cd(2)
h2.Draw('COLZ')
h2.GetXaxis().SetTitle('x-Point in strawtubes')
h2.GetYaxis().SetTitle('y-Point')
t2.Draw()

lc1.cd(3)
h3.Draw('COLZ')
h3.GetXaxis().SetTitle('x-Point in strawtubes')
h3.GetYaxis().SetTitle('y-Point')
t3.Draw()

lc1.cd(4)
h4.Draw('COLZ')
h4.GetXaxis().SetTitle('x-Point in strawtubes')
h4.GetYaxis().SetTitle('y-Point')
t4.Draw()
lc1.Draw()
lc1.SaveAs('hists/B0.2/single_strawtubesPoint_xy_flat.pdf')
lc1.Close()

############################################################################

lc2 = r.TCanvas('can', 'can', 1900,1300)
lc2.Divide(2,2)

lc2.cd(1)
h1.Draw('Lego2')
h1.GetXaxis().SetTitle('x')
h1.GetYaxis().SetTitle('y')

lc2.cd(2)
h2.Draw('Lego2')
h2.GetXaxis().SetTitle('x')
h2.GetYaxis().SetTitle('y')

lc2.cd(3)
h3.Draw('Lego2')
h3.GetXaxis().SetTitle('x')
h3.GetYaxis().SetTitle('y')

lc2.cd(4)
h4.Draw('Lego2')
h4.GetXaxis().SetTitle('x')
h4.GetYaxis().SetTitle('y')
lc2.Draw()
lc2.SaveAs('hists/B0.2/single_strawtubesPoint_xy.pdf')
lc2.Close()

############################################################################

hist1 = r.TH2D('T1','T1 B-field=1.0 (weights)',80,-300,300,160,-600,600)
hist2 = r.TH2D('T2','T2 B-field=1.0 (weights)',80,-300,300,160,-600,600)
hist3 = r.TH2D('T3','T3 B-field=1.0 (weights)',80,-300,300,160,-600,600)
hist4 = r.TH2D('T4','T4 B-field=1.0 (weights)',80,-300,300,160,-600,600)


for i in range(tree2.GetEntries()):
    tree2.GetEvent(i)
    for p in tree2.strawtubesPoint:
        x = p.GetX()
        y = p.GetY()
        z = p.GetZ()
        for track in tree2.MCTrack:
            w = track.GetWeight()
            if z>2500 and z<2700:
                hist1.Fill(x,y,w)
            if z>2700 and z<3000:
                hist2.Fill(x,y,w)
            if z>3000 and z<3400:
                hist3.Fill(x,y,w)
            if z>3400:
                hist4.Fill(x,y,w)


hist1.SetStats(False)
hist2.SetStats(False)
hist3.SetStats(False)
hist4.SetStats(False)

tex1 = r.TLatex(0,530,"entries: {}".format(int(hist1.Integral())))
tex2 = r.TLatex(0,530,"entries: {}".format(int(hist2.Integral())))
tex3 = r.TLatex(0,530,"entries: {}".format(int(hist3.Integral())))
tex4 = r.TLatex(0,530,"entries: {}".format(int(hist4.Integral())))

lc3 = r.TCanvas('can', 'can', 1900,1300)
lc3.Divide(2,2)

lc3.cd(1)
hist1.Draw('COLZ')
hist1.GetXaxis().SetTitle('x-Point in strawtubes')
hist1.GetYaxis().SetTitle('y-Point')
tex1.Draw()

lc3.cd(2)
hist2.Draw('COLZ')
hist2.GetXaxis().SetTitle('x-Point in strawtubes')
hist2.GetYaxis().SetTitle('y-Point')
tex2.Draw()

lc3.cd(3)
hist3.Draw('COLZ')
hist3.GetXaxis().SetTitle('x-Point in strawtubes')
hist3.GetYaxis().SetTitle('y-Point')
tex3.Draw()

lc3.cd(4)
hist4.Draw('COLZ')
hist4.GetXaxis().SetTitle('x-Point in strawtubes')
hist4.GetYaxis().SetTitle('y-Point')
tex4.Draw()
lc3.Draw()
lc3.SaveAs('hists/B1.0/single_strawtubesPoint_xy_flat.pdf')
lc3.Close()

############################################################################

lc4 = r.TCanvas('can', 'can', 1900,1300)
lc4.Divide(2,2)

lc4.cd(1)
hist1.Draw('Lego2')
hist1.GetXaxis().SetTitle('x')
hist1.GetYaxis().SetTitle('y')

lc4.cd(2)
hist2.Draw('Lego2')
hist2.GetXaxis().SetTitle('x')
hist2.GetYaxis().SetTitle('y')

lc4.cd(3)
hist3.Draw('Lego2')
hist3.GetXaxis().SetTitle('x')
hist3.GetYaxis().SetTitle('y')

lc4.cd(4)
hist4.Draw('Lego2')
hist4.GetXaxis().SetTitle('x')
hist4.GetYaxis().SetTitle('y')
lc4.Draw()
lc4.SaveAs('hists/B1.0/single_strawtubesPoint_xy.pdf')
lc4.Close()
