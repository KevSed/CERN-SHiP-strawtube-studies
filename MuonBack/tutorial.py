# -*- coding: utf-8 -*-
from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)

f00 = r.TFile("files/00-100000-MuonBack.root", "READ")
t = f00.cbmsim

# Let's start to do things in loops. E.G.: find the reconstructed vertices.
for event in t:
    for candidate in event.Particles:
        vtx = r.TLorentzVector()
        candidate.ProductionVertex(vtx)
        print "vertex:", vtx.X(), vtx.Y(), vtx.Z(), vtx.T()

# The fourth component of the reconstructed vertex is actually the DOCA, that is made persistent this way.
h_doca = r.TH1F('h_doca', 'h_doca', 100, 0., 103.)

# Retrieve also the reconstructed momentum
for event in t:
    for candidate in event.Particles:
        vtx = r.TLorentzVector()
        candidate.ProductionVertex(vtx)
        print "vertex:", vtx.X(), vtx.Y(), vtx.Z(), vtx.T()
        h_doca.Fill(vtx.T())
        mom = r.TLorentzVector()
        candidate.Momentum(mom)
        print "momentum:", mom.X(), mom.Y(), mom.Z(), mom.T()

# Draw the DOCA
c = r.TCanvas('c', 'c', 950,650)
c.cd(1)
h_doca.GetXaxis().SetTitle('DOCA')
h_doca.Draw('HIST BAR')
c.Draw()
c.SaveAs('hists/geo_doca.pdf')
c.Close()

# Let's do something serious and compute the impact parameter to the target.
def IPtoTarget(vtx, mom):
    p = mom.P()
    target = r.TVector3(0., 0., ShipGeo.target.z0)
    delta = 0.
    for i in range(3):
        delta += (target(i) - vtx(i)) * mom(i)/p
    ip = 0.
    for i in range(3):
        ip += (target(i) - vtx(i) - delta*mom(i)/p)**2.
    return r.TMath.Sqrt(ip)

# Let's look at the IPs in order to define a selection.
for event in t:
        for candidate in event.Particles:
                vtx = r.TLorentzVector()
                candidate.ProductionVertex(vtx)
                print "vertex:", vtx.X(), vtx.Y(), vtx.Z(), vtx.T()
                mom = r.TLorentzVector()
                candidate.Momentum(mom)
                print "momentum:", mom.X(), mom.Y(), mom.Z(), mom.T()
                print 'IP:', IPtoTarget(vtx, mom)
