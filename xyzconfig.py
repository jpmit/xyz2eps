# name of input file
INFILE = 'left.xyz'
OUTFILE = 'test.eps'

# bounding box params
BBOXX0 = 0
BBOXY0 = 0
BBOXX = 76 
BBOXY = 76

# dimension of layer in X and Y coords
# of xyz file
DIM1 = 15.0#24
DIM2 = 15.0#18
PLANE = 'XY'
XSHIFT = 0
YSHIFT = -10

# box in top right
BOX = True
BOXWIDTH = 35
BOXHEIGHT = 15

# properties of the circles (particles)
CRAD = 2.3

# text for writing delta = -13 etc
TEXT = True
textdata = {'delta' : (43,64),
            'equal' : (48,64),
            'minus' : (56,64),
            'one' : (62,64),
            'three' : (67,64)
            }
