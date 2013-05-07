#! /usr/bin/env python
# xyz2eps.py
# 25th April 2013
# James Mithen
# Note: we store all co-ords normalised i.e. between 0 and 1

import x2edata
import readwrite
import numpy as np
import ast
import sys

# rgb tuples for some colors
YELLOW = (1,1,0)
RED = (1,0,0)
BLUE = (0,0,1)
GREY = (0.5,0.5,0.5)

# mapping from symbols that appear in xyz files to default colors
# these are the same colors that jmol uses
DCOLORS = {'S': YELLOW, 'N': BLUE, 'O' : RED, 'OTHER': GREY}

# first column of xyz file is x-coord, 2nd is y, 3rd is z
AXMAP = {'X': 0, 'Y': 1, 'Z': 2}

class Eps(object):
    def __init__(self,params):
        """Create Eps object from info given in xyzconfig.py"""
        self.params = params
        self.positions, self.symbols = readwrite.rxyz(self.params['INFILE'],
                                                      True)
        # bounding box
        self.bbox = Box(self.params['BBOXX0'],self.params['BBOXY0'],
                        self.params['BBOXX'],self.params['BBOXY'])
        # centre coords in units of sigma
        self.center = self.getcenter()
        self.circles = self.createcircles()
        self.boxes = self.createboxes()

    def getcenter(self):
        """Return center (centerx,centery) as tuple"""
        # at the moment, the only option for centering the eps image
        # is as implemented below: The center point is the center of
        # mass of the overlayer atoms (the yellow atoms) in the relevant
        # plane.
        # Note that the parameters XSHIFT and YSHIFT
        # (see creatcircles) can be used to 'shift'
        # the center of the eps.
        spos = np.array([self.positions[i] for i in range(len(self.positions))
                         if self.symbols[i] == 'S'])
        xcenter = sum(spos[:,AXMAP[self.params['PLANE'][0]]])/len(spos)
        ycenter = sum(spos[:,AXMAP[self.params['PLANE'][1]]])/len(spos)
        return (xcenter,ycenter)
        
    def createcircles(self):
        """Return list of circle objects for eps"""
        circles = []
        wpts,hpts = self.bbox.getwidthheight()
        cx,cy = self.getcenter()
        for (p,s) in zip(self.positions,self.symbols):
            # coords of circle in pts
            # in x direction and y direction
            px,py = (p[AXMAP[self.params['PLANE'][0]]],
                     p[AXMAP[self.params['PLANE'][1]]])
            px = (px - cx) + self.params['DIM1']/2.0
            py = (py - cy) + self.params['DIM2']/2.0
            px = (px/self.params['DIM1']) * wpts + self.params['XSHIFT']
            py = (py/self.params['DIM2']) * hpts + self.params['YSHIFT']
            # colour
            ccol = DCOLORS.get(s,DCOLORS['OTHER'])
            # create circle and add to list if it 'fits' into the eps
            # print px,py
            circ = Circle(px,py,self.params['CRAD'],ccol)            
            if self.bbox.circlefits(circ):
                circles.append(circ)
        return circles

    def createboxes(self):
        """White boxes, generally for top right of image so can write legend"""
        # currently only a single white box is possible in top right
        if not self.params['BOX']:
            return []
        box = Box(self.params['BBOXX'] - self.params['BOXWIDTH'],
                  self.params['BBOXY'] - self.params['BOXHEIGHT'],
                  self.params['BBOXX'], self.params['BBOXY'])
        return [box]
    
    def getepsstring(self):
        """Return a string that is text of the eps file"""
        box = self.bbox.getbox()
        
        # header string
        hstr = x2edata.EPSHEAD.format(box[0],box[1],box[2],box[3])

        # circles string
        cstr = self._getcstr()

        # boxes string
        boxstr = self._getboxstr()

        # font string
        fontstr = x2edata.FONTSTR

        # text string
        textstr = gettextstr(self.params)

        # tail string
        tstr = x2edata.EPSFOOT

        return '%s%s%s%s%s%s' %(hstr,fontstr,cstr,boxstr,textstr,tstr)

    def _getcstr(self):
        """Return a string that is text for drawing the circles"""
        cstr = ''
        for circ in self.circles:
            col = ' '.join([str(i) for i in circ.color])
            cstr = '%s%.2f %.2f %.2f %s filledcircle\n' %(cstr,circ.posx,
                                                            circ.posy,circ.rad,
                                                            col)
        return cstr

    def _getboxstr(self):
        """Return a string that is text for drawing the box"""
        box = self.boxes[0]
        boxstr = x2edata.BOXSTR.format(box.lleftx,
                                       box.llefty,
                                       box.urighty-box.llefty,
                                       box.urightx-box.lleftx,
                                       box.llefty-box.urighty,
                                       box.lleftx-box.urightx)
        return boxstr

    def writeeps(self,outfile):
        outstr = self.getepsstring()
        fout = open(outfile,'w')
        fout.write(outstr)
        fout.close()
        return

class Box(object):
    """Box used as bounding box for the EPS, box for legend etc."""
    def __init__(self,lleftx,llefty,urightx,urighty):
        # note the syntax in eps file for bounding box is
        # BoundingBox: lowerleftx, lowerlefty, upperrightx,
        # upperrighty
        self.lleftx = lleftx
        self.llefty = llefty
        self.urightx = urightx
        self.urighty = urighty

    def getbox(self):
        """Return tuple with info for writing to eps file"""
        return (self.lleftx,self.llefty,self.urightx, self.urighty)

    def getwidthheight(self):
        return (self.urightx - self.lleftx,
                self.urighty - self.llefty)

    def circlefits(self,circ):
        """
        Work out if the circle fits in bounding box based on its
        coordinates and radius. posx and posy should be absolute
        positions given in pts."""
        if ((circ.posx < self.urightx - circ.rad) and
            (circ.posx > self.lleftx + circ.rad) and
            (circ.posy < self.urighty - circ.rad) and
            (circ.posy > self.llefty + circ.rad)):
            return True
        return False

class Circle(object):
    def __init__(self,posx,posy,rad,color):
        self.posx = posx # in pts
        self.posy = posy # in pts
        self.rad = rad # in pts
        self.color = color # normalised rgb tuple

def getfontstr(params):
    if not params['TEXT']:
        return ''
    else:
        return FONTSTR
    
def gettextstr(params):
    """get text (if any) to write on image"""
    if not params['TEXT']:
        return ''
    # mapping between characters and font
    charmap = {'delta': 'Cmmi10',
               'minus': 'Cmsy10',
               }
    default = 'Cmr10'
    textstr = ''
    # go though each character in textdata
    # dictionary
    for text in params['TEXTDATA']:
        char = text[0]
        # choose the font for the character
        font = charmap.get(char, default)
        textstr = ('%s/%s findfont\n'
                   '10 scalefont\n'
                   'setfont\n'
                   '%d %d moveto\n'
                   '/%s glyphshow\n\n' %(textstr,font,
                                         text[1],text[2], char)
                   )
        
    return textstr

def getparams(fname):
    """Get dictionary of parameters for file fname"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    params = {}
    # conversion functions for all non-string args
    INT = int
    FLOAT = float
    BOOL = lambda x : True if x == 'True' else False
    IDENTITY = lambda x: x # default conversion used for str
    convert = {'BBOXX0': INT,
               'BBOXY0': INT,
               'BBOXX': INT,
               'BBOXY':INT,
               'DIM1':FLOAT,
               'DIM2':FLOAT,
               'XSHIFT':INT,
               'YSHIFT':INT,
               'BOX':BOOL,
               'BOXWIDTH':INT,
               'BOXHEIGHT':INT,
               'CRAD':FLOAT,
               'TEXT':BOOL,
               'TEXTDATA': ast.literal_eval # trick
               }
    # go through each line of file,
    # adding parameters to dictionary
    for lin in lines:
        if '#' in lin:
            continue
        try:
            name,val = lin.split('=')
        except ValueError:
            continue
        val = val.strip() # remove \n
        if name:
            cfunc = convert.get(name, IDENTITY)
            nval = cfunc(val)
            params[name] = nval
    return params

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    params = getparams(infile)
    eps = Eps(params)
    eps.writeeps(outfile)
