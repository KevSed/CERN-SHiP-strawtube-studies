from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)

# read the rootfile
#f01 = r.TFile("files/ship.conical.MuonBack-TGeant4-01-100000_rec.root", "READ")
#f02 = r.TFile("files/ship.conical.MuonBack-TGeant4-02-100000_rec.root", "READ")
f04 = r.TFile("files/04-100000-MuonBack_rec.root", "READ")
#f06 = r.TFile("files/ship.conical.MuonBack-TGeant4-06-100000_rec.root", "READ")
#f10 = r.TFile("files/ship.conical.MuonBack-TGeant4-10-100000_rec.root", "READ")
#f12 = r.TFile("files/ship.conical.MuonBack-TGeant4-12_rec.root", "READ")
#t01 = f01.cbmsim
#t02 = f02.cbmsim
t04 = f04.cbmsim
#t06 = f06.cbmsim
#t10 = f10.cbmsim
#t12 = f12.cbmsim
#geo = r.TFile("files/geofile_full.conical.MuonBack-TGeant4-04.root","read")
#geo.ls()
#sGeo = geo.FAIRGeom
#fGeo = r.gGeoManager

###############################################################################

def IPtoTarget(vtx, mom):
    p = mom.P()
    target = r.TVector3(0.,0.,ShipGeo.target.z0)
    delta = 0.
    for i in range(3):
        delta = delta+(target(i)-vtx(i))*mom(i)/p
    ip = 0.
    for i in range(3):
        ip=ip+(target(i)-vtx(i)-delta*mom(i)/p)**2
    return np.sqrt(ip)

'''
for i in range(t04.GetEntries()):
    t04.GetEvent(i)
    for candidate in t04.Particles:
        print candidate
        vtx = r.TLorentzVector()
        candidate.ProductionVertex(vtx)
        mom = r.TLorentzVector()
        candidate.Momentum(mom)
        print "IP: ", IPtoTarget(vtx,mom)'''
'''
n=0
m=0
for i in range(t02.GetEntries()):
    t02.GetEvent(i)
    for track in t02.FitTracks:
        w = track.hasFitStatus()
        if w==True:
            n=n+1
        if w==False:
            m=m+1
print "Fit Status for B=0.2: True ", n, " False ", m
print("------------------------------------------------------")
'''

'''
n=0
m=0
for i in range(t04.GetEntries()):
    t04.GetEvent(i)
    for track in t04.FitTracks:
        w = track.hasFitStatus()
        t = track.getPoint(1)
        if w==True:
            n=n+1
            print "point: "
        if w==False:
            m=m+1
print "Fit Status for B=0.4: True ", n, " False ", m
print("------------------------------------------------------")
'''

n=0
m=0
k=0
for event in t04:
    for track in event.FitTracks:
        w = track.hasFitStatus()
        p = track.getPoints()
        t = len(event.FitTracks)
        print "points: ", track
        if t==0:
            n=n+1
        if t==1:
            m=m+1
        if t>=2:
            k=k+1
print "Len of FitTracks for B=0.4: len=0: ", n, " len=1: ", m, " len>=2: ", k
print("------------------------------------------------------")


'''
n=0
m=0
for i in range(t06.GetEntries()):
    event = t06.GetEvent(i)
    for track in t06.FitTracks:
        w = track.hasFitStatus()
        if w==True:
            n=n+1
        if w==False:
            m=m+1
print "Fit Status for B=0.6: True ", n, " False ", m
print("------------------------------------------------------")


n=0
m=0
for i in range(t10.GetEntries()):
    t10.GetEvent(i)
    for track in t10.FitTracks:
        w = track.hasFitStatus()
        if w==True:
            n=n+1
        if w==False:
            m=m+1
print "Fit Status for B=1.0: True ", n, " False ", m
'''


#fGeo.FindNode(track.GetStartX(),track.GetStartY(),track.GetStartZ()).GetName()
