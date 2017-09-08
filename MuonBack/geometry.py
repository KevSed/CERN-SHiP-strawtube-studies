from __future__ import division
import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
from decorators import *

f = r.TFile("files/no-mag_rec.root", "READ")
t = f.Get("cbmsim")

geo = r.TFile("files/geofile_no-mag.root","read")
sGeo = geo.FAIRGeom
fGeo = r.gGeoManager

# The fourth component of the reconstructed vertex is actually the DOCA, that is made persistent this way.
#h_doca = r.TH1F('h_doca', 'h_doca', 100, 0., 103.)

# Retrieve also the reconstructed momentum

#for event in t:
#    for track in event.MCTrack:
#        print fGeo.FindNode(track.GetStartX(),track.GetStartY(),track.GetStartZ()).GetName()

#for node in fGeo.GetTopVolume().GetNodes():
#   print node.GetName()

for node in fGeo.GetTopVolume().GetNodes():
    print node.GetName()

fGeo.GetTopVolume().GetNode('TargetArea_1').GetMatrix().GetTranslation()
xtarget_center = fGeo.GetTopVolume().GetNode('TargetArea_1').GetMatrix().GetTranslation()[0]
ytarget_center = fGeo.GetTopVolume().GetNode('TargetArea_1').GetMatrix().GetTranslation()[1]
ztarget_center = fGeo.GetTopVolume().GetNode('TargetArea_1').GetMatrix().GetTranslation()[2]
ztarget = fGeo.GetTopVolume().GetNode('TargetArea_1').GetMatrix().GetTranslation()
dztarget = fGeo.GetTopVolume().GetNode('TargetArea_1').GetVolume().GetShape().GetDZ()


print "centre of TargetArea_1: ", ztarget_center, "width: ", dztarget

'''
c = r.TCanvas('c', 'c', 950,650)

# Draw the DOCA
c.cd(1)
h_doca.Draw()
c.Draw()
c.SaveAs("hists/geo_DOCA.pdf")
c.Close()
'''
