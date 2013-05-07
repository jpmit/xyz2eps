# makeeps.py
# 24th January 2012
# just try this on an xyz file and see what you get!
# to do (21st March 2013) - clean this up massively

import readwrite
import numpy as np
import os
import inspect
import sys

# parameter for three column image
# XDIM 45.0
# YDIM 30.0
# XPTS 78
# YPTS 78
# CRAD 2.3

# settings for pos10_2
# XDIM=24.0,YDIM=18.0,XPTS=234,YPTS=156,XOFF=0,YOFF=-20,XSHIFT=0,
# YSHIFT=-35,CRAD=4.0,PLANE='YZ',DRAWBOX=False,inname='pos10_2.xyz'

XDIM = 24.0 # width and height of image 15.0 * 15.0 for three column image
YDIM = 18.0 
XPTS = 234 # 78 for three column image
YPTS = 156
XOFF = 0 # bottom left corner of BoundingBox is at XOFF,YOFF
YOFF = -20
XSHIFT = 0 # shift from centre by XSHIFT,YSHIFT
YSHIFT = -35
CRAD = 4.0 # circle radius (for both surface and overlayer pars) 2.3 for three
PLANE = 'YZ' # choose plane of 3d xyz image to map
DRAWBOX = False # draw white box in top right for text?
infol = ''#os.path.expanduser('~/awork/montecarlo/papers/flat/overlayerfiles/')
inname = 'pos10_2.xyz'
outname = 'mytest.eps'
outfile = open('%s%s' %(infol,outname),'w')

# To do - these assignments are pretty inefficient
positions,symbols = readwrite.rxyz(infol+inname,True)
rawsurfpos = np.array([positions[i] for i in range(len(positions)) if symbols[i] == 'O'])
rawlaypos = np.array([positions[i] for i in range(len(positions)) if symbols[i] == 'S'])
rawliquidpos = np.array([positions[i] for i in range(len(positions)) if symbols[i] == 'N'])
rawotherpos = np.array([positions[i] for i in range(len(positions)) if symbols[i] == 'C'])

# next, we shuffle around positions so that surfpos[0,;] and laypos[0,;] contain positions on first chosen axis
# and surfpos[1,;] and laypos[1,;] contain positions on second chosen axis
surfpos = np.empty(rawsurfpos.shape)
laypos = np.empty(rawlaypos.shape)
liquidpos = np.empty(rawliquidpos.shape)
otherpos = np.empty(rawotherpos.shape)
lookup = {'X': 0, 'Y': 1, 'Z': 2}

if len(rawsurfpos):
    surfpos[:,0] = rawsurfpos[:,lookup[PLANE[0]]]
    surfpos[:,1] = rawsurfpos[:,lookup[PLANE[1]]]
if len(rawlaypos):
    laypos[:,0] = rawlaypos[:,lookup[PLANE[0]]]
    laypos[:,1] = rawlaypos[:,lookup[PLANE[1]]]    
if len(rawliquidpos):
    liquidpos[:,0] = rawliquidpos[:,lookup[PLANE[0]]]
    liquidpos[:,1] = rawliquidpos[:,lookup[PLANE[1]]]    
if len(rawotherpos):
    otherpos[:,0] = rawotherpos[:,lookup[PLANE[0]]]
    otherpos[:,1] = rawotherpos[:,lookup[PLANE[1]]]
    
xcenter = sum(laypos[:,0])/len(laypos)
ycenter = sum(laypos[:,1])/len(laypos)
if len(surfpos):
    surfpos = surfpos - np.array([xcenter,ycenter,0]) + np.array([XDIM/2,YDIM/2,0])
if len(laypos):
    laypos = laypos - np.array([xcenter,ycenter,0]) + np.array([XDIM/2,YDIM/2,0])
if len(liquidpos):
    liquidpos = liquidpos - np.array([xcenter,ycenter,0]) + np.array([XDIM/2,YDIM/2,0])
if len(otherpos):
    otherpos = otherpos - np.array([xcenter,ycenter,0]) + np.array([XDIM/2,YDIM/2,0])
# center of circles
laycircs = laypos
surfcircs = surfpos
liquidcircs = liquidpos
othercircs = otherpos

if len(laycircs):
    laycircs[:,0] = (laypos[:,0]/XDIM)*XPTS + XOFF + XSHIFT
    laycircs[:,1] = (laypos[:,1]/YDIM)*YPTS + YOFF + YSHIFT
if len(surfcircs):
    surfcircs[:,0] = (surfpos[:,0]/XDIM)*XPTS + XOFF + XSHIFT
    surfcircs[:,1] = (surfpos[:,1]/YDIM)*YPTS + YOFF + YSHIFT
if len(liquidcircs):
    liquidcircs[:,0] = (liquidpos[:,0]/XDIM)*XPTS + XOFF + XSHIFT
    liquidcircs[:,1] = (liquidpos[:,1]/YDIM)*YPTS + YOFF + YSHIFT
if len(othercircs):
    othercircs[:,0] = (otherpos[:,0]/XDIM)*XPTS + XOFF + XSHIFT
    othercircs[:,1] = (otherpos[:,1]/YDIM)*YPTS + YOFF + YSHIFT
    
epshead = """%!PS-Adobe-3.0 EPSF-3.0
%%Title: none
%%Creator: JP Mithen
%%CreationDate: today
%%BoundingBox: {0} {1} {2} {3}
%%EndComments
%%BeginProlog
/filledcircle {{
6 dict begin
/b exch def
/g exch def
/r exch def
/radius exch def
/y exch def
/x exch def
newpath
x y radius 0 360 arc
closepath
gsave
r g b setrgbcolor fill
grestore
gsave
0.5 setlinewidth
0 setgray
stroke
grestore
end
}} bind def
%%EndProlog
""".format(XOFF,YOFF,XPTS+XOFF,YPTS+YOFF)
epsfoot = 'showpage\n'

# draw circles
outstr = ''
if len(liquidcircs): # liquid circles are blue
    for pos in liquidcircs[:,0:2]:
        if (pos[0] < XPTS+XOFF-CRAD) and (pos[1] < YPTS+YOFF-CRAD) and (pos[0] > XOFF+CRAD) and (pos[1] > YOFF+CRAD):    
            outstr = '%s%.2f %.2f %.2f 0 0 1 filledcircle\n' %(outstr,pos[0],pos[1],CRAD)
if len(surfcircs): # surface circles are red
    for pos in surfcircs[:,0:2]:
        if (pos[0] < XPTS+XOFF-CRAD) and (pos[1] < YPTS+YOFF-CRAD) and (pos[0] > XOFF+CRAD) and (pos[1] > YOFF+CRAD):
            outstr = '%s%.2f %.2f %.2f 1 0 0 filledcircle\n' %(outstr,pos[0],pos[1],CRAD)
if len(laycircs): # xtal circles are yellow
    for pos in laycircs[:,0:2]:
        if (pos[0] < XPTS+XOFF-CRAD) and (pos[1] < YPTS+YOFF-CRAD) and (pos[0] > XOFF+CRAD) and (pos[1] > YOFF+CRAD):    
            outstr = '%s%.2f %.2f %.2f 1 1 0 filledcircle\n' %(outstr,pos[0],pos[1],CRAD)
if len(othercircs): # other circles are grey
    for pos in othercircs[:,0:2]:
        if (pos[0] < XPTS+XOFF-CRAD) and (pos[1] < YPTS+YOFF-CRAD) and (pos[0] > XOFF+CRAD) and (pos[1] > YOFF+CRAD):    
            outstr = '%s%.2f %.2f %.2f 0.5 0.5 0.5 filledcircle\n' %(outstr,pos[0],pos[1],CRAD)

# draw box in top right corner
boxdims = [35,15]
boxpos = [XPTS+XOFF - boxdims[0], YPTS+YOFF - boxdims[1]]
boxstr = """newpath
%d %d moveto
0 %d rlineto
%d 0 rlineto
0 %d rlineto
%d 0 rlineto
closepath
gsave
1 setgray
fill
stroke
""" %(boxpos[0],boxpos[1],boxdims[1],boxdims[0],-boxdims[1],-boxdims[0])
## # put text in box
## textstr = """
## /Helvetica findfont
## 12 scalefont
## setfont
## newpath
## %d %d moveto
## (delta) show
## """ %(boxpos[0],boxpos[1])
if DRAWBOX:
    outstr = '%s%s' %(outstr,boxstr)

outfile.write(epshead+outstr+epsfoot)
outfile.close()
