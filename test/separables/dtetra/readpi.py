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

def importPiStructure(file):
    """ load a Structure from Pi file

    We map the 1,2 indicator for atom type to 0,1
    We map the integer and halfinteger designation for BCC sites
       to even and odd integers

    filename and path are separated so that filename can be the
    structure name

    """
    import re
    from EquivStructure import Structure
    from numpy import array
    from scipy.linalg import inv

    # finds first instance of "NO."
    reNdec = re.compile('^\sNO\.')
    line = ""
    result = None
    for line in file:
      result = reNdec.match( line )
      if result is None: continue
      break
    if result is None: raise IOError 
    splitted = line.strip().split()
    # creates atomic species entries in structure
    N = eval(splitted[3])

    decoration = int( round( eval( splitted[4] ) ) )
    entries = [ int( (decoration >> n) %2 ) for n in range(N) ]
    cell = array( [ [ int( round( eval(x) ) ) for x in splitted[ 6: 9] ],
                    [ int( round( eval(x) ) ) for x in splitted[ 9:12] ],
                    [ int( round( eval(x) ) ) for x in splitted[12:15] ] ] )
    reciprocal=inv(cell)
    index = splitted[1]
    # now gets atomic positions
    reNdec = re.compile('^\sBASIS')
    current_pos = file.tell()
    for line in file:
      if reNdec.match( line ) is not None: break

    vectors = []
    oldline = line;

    splitted = line.strip().split()
    splitted = splitted[1:]
    for i in range( len(splitted) / 3 ):
      vectors.append(( int( round( eval( splitted[  i*3] ) ) ),
                       int( round( eval( splitted[i*3+1] ) ) ),
                       int( round( eval( splitted[i*3+2] ) ) ) ))

    if len( vectors ) < N:
      for line in file:
        splitted = line.strip().split()
        for i in range( len(splitted) / 3 ):
          vectors.append(( int( round( eval( splitted[  i*3] ) ) ),
                           int( round( eval( splitted[i*3+1] ) ) ),
                           int( round( eval( splitted[i*3+2] ) ) ) ))
        if len( vectors ) >= N: break
#   

    return Structure(vectors, entries,cell, reciprocal,index, 0.0)

def SFtoCE_Structure(structure):
    """ Converts EquivStructure.Structure to an atom.Structure
    """
    import Pylada
    from EquivStructure import Structure
    result = Pylada.Structure()
    dummy = Pylada.rMatrix3d(); dummy.diag( Pylada.rVector3d( [0.5, 0.5, 0.5 ] ) )
    result.cell = Pylada.make_rMatrix3d( structure.period ) * dummy
    
    for i in range( len( structure.entries ) ):
      result.atoms.append( Pylada.Atom( [ 0.5 * structure.vectors[i][0],
                                        0.5 * structure.vectors[i][1],
                                        0.5 * structure.vectors[i][2],
                                        structure.entries[i]*2-1 ] ) )
    result.scale = 1.0 
    return result
  
