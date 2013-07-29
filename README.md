xyz2eps.py
==========

James Mithen  
j.mithen@surrey.ac.uk

xyz2eps.py is a python script for producing 'publication quality'
encapsulated postscript (eps) images of systems of particles from data
in XYZ file format.

Usage
------

    $ xyz2eps.py infile outfile.eps
Where infile is a configuration file specifying details of the eps to
be drawn.

The best thing to do is probably to look at the examples, which are
included in the examples/ directory.

To produce eps images for the examples, navigate to the examples
directory and type
    $ ../xyz2eps.py example1.in example1.eps
    $ ../xyz2eps.py example2.in example2.eps
This will produce eps images named example1.eps and example2.eps

License
-------

GPLv3 - see http://www.gnu.org/licenses/gpl.html

Todo 
---- 

* Support a larger set of glyphs, current only characters
'(',')','delta','-','=',and digits 0-9 can be drawn (!)
* Allow for different 'orderings' of circles on the eps image
