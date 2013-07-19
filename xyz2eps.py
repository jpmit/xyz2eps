#! /usr/bin/env python
# xyz2eps.py
#
# James Mithen
# j.mithen@surrey.ac.uk
#
# Create an eps image from an 'XYZ' file that stores positions of a
# set of particles in 3d space.  See README and the examples/
# directory for how this works.

import x2edata
import numpy as np
import ast
import sys

class Xyz2EpsError(Exception): pass

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
        """Create Eps object from info given in dictionary params"""
        self.params = params
        self.positions, self.symbols = readxyz(self.params['INFILE'],
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
        # mass of the overlayer atoms (the yellow atoms) in the
        # relevant plane.  Note that the parameters XSHIFT and YSHIFT
        # (see creatcircles) can be used to 'shift' the center of the
        # eps.
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
            # coords of circle in pts in x direction and y direction
            px,py = (p[AXMAP[self.params['PLANE'][0]]],
                     p[AXMAP[self.params['PLANE'][1]]])
            px = (px - cx) + self.params['DIM1']/2.0
            py = (py - cy) + self.params['DIM2']/2.0
            px = (px/self.params['DIM1']) * wpts + self.params['XSHIFT']
            py = (py/self.params['DIM2']) * hpts + self.params['YSHIFT']
            # colour
            ccol = DCOLORS.get(s,DCOLORS['OTHER'])
            # create circle and add to list if it 'fits' into the eps
            # bounding box
            circ = Circle(px,py,self.params['CRAD'],ccol)            
            if self.bbox.circlefits(circ):
                circles.append(circ)
        return circles

    def createboxes(self):
        """White boxes, generally for top right of image so can write legend"""
        boxes = []
        # can create a box in top right and/or top left
        if self.params['RBOX']:
            box = Box(self.params['BBOXX'] - self.params['RBOXWIDTH'],
                      self.params['BBOXY'] - self.params['RBOXHEIGHT'],
                      self.params['BBOXX'], self.params['BBOXY'])
            boxes.append(box)
        if self.params['LBOX']:
            box = Box(0, self.params['BBOXY'] - self.params['LBOXHEIGHT'],
                      0 + self.params['LBOXWIDTH'], self.params['BBOXY'])
            boxes.append(box)
        return boxes
    
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

    def ordercircles(self):
        """Order circles so that we draw any yellow circles last"""
        for i in range(len(self.circles)):
            if self.circles[i].color == YELLOW:
                self.circles.append(self.circles.pop(i))

    def _getcstr(self):
        """Return a string that is text for drawing the circles"""
        cstr = ''
        self.ordercircles()
        for circ in self.circles:
            col = ' '.join([str(i) for i in circ.color])
            cstr = '%s%.2f %.2f %.2f %s filledcircle\n' %(cstr,circ.posx,
                                                            circ.posy,circ.rad,
                                                            col)
        return cstr

    def _getboxstr(self):
        """Return a string that is text for drawing the box"""
        boxstr = ''
        for box in self.boxes:
            boxstr += x2edata.BOXSTR.format(box.lleftx,
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
        # note the syntax in eps file for bounding box is BoundingBox:
        # lowerleftx, lowerlefty, upperrightx, upperrighty
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
    """get font string which contains glyphs for eps file."""
    if not params['TEXT']:
        # we are not writing any text on the eps image, hence we
        # dont want to write FONTSTR, which contains the glyphs
        return ''
    else:
        return FONTSTR
    
def gettextstr(params):
    """get text (if any) to write on image"""
    if not params['TEXT']:
        return ''
    # mapping between characters and font
    fontmap = {'delta': 'Cmmi10',
               'minus': 'Cmsy10',
               'equal': 'Cmr10',
               'zero': 'Cmr10',
               'one': 'Cmr10',
               'two': 'Cmr10',
               'three': 'Cmr10',
               'four': 'Cmr10',
               'five': 'Cmr10',
               'six': 'Cmr10',
               'seven': 'Cmr10',
               'eight': 'Cmr10',
               'nine': 'Cmr10',               
               }
    defaultfont = 'DejaVuSans'
    textstr = ''
    # go though each character in textdata
    # dictionary
    for text in params['TEXTDATA']:
        char = text[0]
        font = fontmap.get(char, defaultfont)
        textstr = ('%s/%s findfont\n'
                   '10 scalefont\n'
                   'setfont\n'
                   '%d %d moveto\n'
                   '/%s glyphshow\n\n' %(textstr,font,
                                         text[1],text[2], char)
                   )
        
    return textstr

def getparams(fname):
    """Get dictionary of parameters from file fname"""
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
               'RBOX':BOOL,
               'RBOXWIDTH':INT,
               'RBOXHEIGHT':INT,
               'LBOX':BOOL,
               'LBOXWIDTH':INT,
               'LBOXHEIGHT':INT,
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

def readxyz(fname, retsymbols=False, splines=1):
    """Read a .xyz coordinate file, return particle positions"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    npar = int(lines[0])
    positions = np.empty([npar,3])
    symbols = [None]*npar
    coordline = 1 + splines
    i = 0
    for line in lines[coordline:]:
        li = line.split()
        symbols[i] = li[0]
        positions[i,0] = float(li[1])
        positions[i,1] = float(li[2])
        positions[i,2] = float(li[3])
        i = i + 1
    if retsymbols:
        return positions,symbols
    else:
        return positions

def main():
    # name of the infile
    try:
        infile = sys.argv[1]
    except IndexError:
        raise Xyz2EpsError, 'no input file name supplied'
    # name of eps file to write output to
    try:
        outfile = sys.argv[2]
    except IndexError:
        raise Xyz2EpsError, 'no output file name supplied'
    
    # get parameters from input file
    params = getparams(infile)
    # create the Eps object
    eps = Eps(params)
    # write the Eps object to an eps file, which is the image we want
    eps.writeeps(outfile)

if __name__ == '__main__':
    main()
