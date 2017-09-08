# -*- coding: utf-8 -*-
from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *
import sys
import re

# Supress display output
if not '-showCanvas' in sys.argv:
    r.gROOT.SetBatch(r.kTRUE)

f = r.TFile("files/no-mag_rec.root", "READ")
t = f.cbmsim

geo = r.TFile("files/geofile_no-mag.root","read")
sGeo = geo.FAIRGeom
fGeo = r.gGeoManager


#part = np.array(t.GetEntries(),dtype='string')
part = []

for event in t:
    for mctrack in event.MCTrack:
        part.append(fGeo.FindNode(mctrack.GetStartX(),mctrack.GetStartY(),mctrack.GetStartZ()).GetName())
np.savetxt('targets.txt',part,fmt="%s")
