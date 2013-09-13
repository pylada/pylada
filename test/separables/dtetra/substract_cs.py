###############################
#  This file is part of PyLaDa.
#
#  Copyright (C) 2013 National Renewable Energy Lab
# 
#  PyLaDa is a high throughput computational platform for Physics. It aims to make it easier to submit
#  large numbers of jobs on supercomputers. It provides a python interface to physical input, such as
#  crystal structures, as well as to a number of DFT (VASP, CRYSTAL) and atomic potential programs. It
#  is able to organise and launch computational jobs on PBS and SLURM.
# 
#  PyLaDa is free software: you can redistribute it and/or modify it under the terms of the GNU General
#  Public License as published by the Free Software Foundation, either version 3 of the License, or (at
#  your option) any later version.
# 
#  PyLaDa is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
#  Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along with PyLaDa.  If not, see
#  <http://www.gnu.org/licenses/>.
###############################

#! /usr/bin/python
from readpi import importPiStructure, SFtoCE_Structure
from Pylada import rMatrix3d, rVector3d, make_rMatrix3d, Lattice, \
                 tetragonalCS
from sys import exit
from numpy import array
from EquivStructure import makeRotationGroup, importLDA,\
                           importStructures


lattice = Lattice()
lattice.fromXML( "input.xml" )

# h = lattice.syms()
# # print "nb syms: ", len( lattice.syms() )
# # for x in lattice.syms():
# #   if x[1] == rVector3d( [0,0,0] ):
# #     print "sym:\n", x[0]

# G = [ x[0] for x in lattice.syms() if x[0].det() > 0 ]

# string = "array( [ [%2i,%2i,%i], [%2i,%2i,%2i], [%2i,%2i,%2i] ] ),"
# for x in G:
#   print  string % ( x[0,0], x[0,1], x[0,2],
#                     x[1,0], x[1,1], x[1,2],
#                     x[2,0], x[2,1], x[2,2] )

cs = tetragonalCS()
cs.fromXML( "input.xml" )
lda = importLDA( filename="originalLDAs.dat" )
for key in lda.keys():
  esStruct = importStructures( key, lda[key], "."  )
  structure = SFtoCE_Structure( esStruct )
  cs.define( structure )
  cs.vars()[:] = [ 2*x-1 for x in esStruct.entries ]
  print "%40s  %6.2f + %6.2f" % \
        ( esStruct.name,
          lda[key] - cs.evaluate(), cs.evaluate() )
  
# filename = "PIs2thru16"
# file = open(filename, 'r')
# oldcell= array( [ [0,0,0], [0,0,0], [0,0,0] ] ) 
# while( True ):

#   try:
#     esStruct = importPiStructure( file )
#   except IOError: 
#     print "End of %s\n" % ( file )
#     exit();
#   else:

#     structure = Structure()
#     if abs( array( oldcell - esStruct.period ).sum() ) > 0.001:
#       structure = SFtoCE_Structure( esStruct )
#       cs.define( structure )
#       oldcell = esStruct.period
#     cs.vars()[:] = [ 2*x-1 for x in esStruct.entries ]
#     print "%8.2f" % ( cs.evaluate() )
#     cs.vars()[:] = [ -x for x in cs.vars() ]
#     print "%8.2f" % ( cs.evaluate() )


#file.close()

