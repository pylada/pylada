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


string = \
"""
date Tue Oct 23 01:47:26 BST 2007
hostname cn141.scarf.rl.ac.uk
system Linux cn141.scarf.rl.ac.uk 2.6.9-55.0.9.ELsmp #1 SMP Tue Sep 25 02:16:15 EDT 2007 x86_64 x86_64 x86_64 GNU/Linux
user scarf044
executable in /work/scratch/scarf044/pub/CRYSTAL06/20070510_OPTbugfix.ANHARM/Linux-mpichpgf_Dynamic/pgi_6.2_apps_mpi_gm_pgi_1.2.7..15
temporary directory /work/scratch/scarf044/3L_optatom2_29478
output in /home/scarf044/giuseppe/TiO2/slab_110/rutile/B3LYP/Ti_86411d41_O_8411d1_TI777714_IS888/3L_optatom2.out
input data /home/scarf044/giuseppe/TiO2/slab_110/rutile/B3LYP/Ti_86411d41_O_8411d1_TI777714_IS888/3L_optatom2.d12
rutile
CRYSTAL
0 0 0
136
4.63909875  2.97938395
2
22 0.0 0.0 0.0
8  3.061526467783E-01 3.061526467783E-01 0.0
SLAB
1 1 0
2 3
BREAKSYM
ATOMDISP
6
 1     0.0000000000     0.0000000000    0.0000000000 
 2     0.0000000000     0.0000000000    0.0000000000 
 3     0.0000000000     0.0000000000    0.0000000000 
 4     0.0000000000     0.0000000000    0.0000000000 
 5     0.0000000000     0.0000000000    0.0000000000 
 6     0.0000000000     0.0000000000    0.0000000000 
BREAKSYM
ATOMDISP
6
 1     0.0000000000     0.0000000000    0.0000000000 
 2     0.0000000000     0.0000000000    0.0000000000 
 3     0.0000000000     0.0000000000    0.0000000000 
 4     0.0000000000     0.0000000000    0.0000000000 
 5     0.0000000000     0.0000000000    0.0000000000 
 6     0.0000000000     0.0000000000    0.0000000000 
OPTGEOM
CVOLOPT
MAXCYCLE
10
ENDOPT
END
22 7
0 0 8 2. 1.
225338.0 0.000228
32315.0 0.001929
6883.61 0.011100
1802.14 0.05
543.063 0.17010
187.549 0.369
73.2133 0.4033
30.3718 0.1445
0 1 6 8. 1.
554.042 -0.0059 0.0085
132.525 -0.0683 0.0603
43.6801 -0.1245 0.2124
17.2243 0.2532 0.3902
7.2248 0.6261 0.4097
2.4117 0.282 0.2181
0 1 4 8. 1.
24.4975 0.0175 -0.0207
11.4772 -0.2277 -0.0653
4.4653 -0.7946 0.1919
1.8904 1.0107 1.3778
0 1 1 0. 1.
0.8126 1.0 1.0
0 1 1 0. 1.
0.3297 1.0 1.0
0 3 4 2. 1.
16.2685 0.0675
4.3719 0.2934
1.4640 0.5658
0.5485 0.5450
0 3 1 0. 1.
0.26 1.0
8 5
0 0 8 2. 1.
8020.0 0.00108
1338.0 0.00804
255.4 0.05324
69.22 0.1681
23.90 0.3581
9.264 0.3855
3.851 0.1468
1.212 0.0728
0 1 4 7. 1.
49.43 -0.00883 0.00958
10.47 -0.0915 0.0696
3.235 -0.0402 0.2065
1.217 0.379 0.347
0 1 1 0. 1.
0.4567 1.0 1.0
0 1 1 0. 1.
0.1843 1.0 1.0
0 3 1 0. 1.
 0.6 1.0
99 0
END
DFT
B3LYP
LGRID
END
TOLINTEG
7 7 7 7 14
SHRINK
8 8
LEVSHIFT
5 1
FMIXING
30
PPAN
TOLDEE
7
EXCHSIZE
6937578
BIPOSIZE
9701800
SCFDIR
END
  PROCESS             5  OF             8  WORKING
 PROCESS             1  OF             8  WORKING
 PROCESS             3  OF             8  WORKING
   PROCESS             2  OF             8  WORKING
 PROCESS             0  OF             8  WORKING

         
 *******************************************************************************
 *                                                                             *
 *                                CRYSTAL06                                    *
 *                              Release : 1.0                                  *
 *                      V1.0.2 - sequential executable                         *
 *                                                                             *
 *                                                                             *
 *                                                                             *
 *                              MAIN AUTHORS                                   *
 *                                                                             *
 *    R. DOVESI(1,10), V.R. SAUNDERS(2), C. ROETTI(1,10), R. ORLANDO (1,3),    *
 *  C.M. ZICOVICH-WILSON(1,4), F. PASCALE(5), B. CIVALLERI(1,10), K. DOLL(6),  *
 *       N.M. HARRISON(2,7), I. J. BUSH(2), Ph. D'ARCO(8), M. LLUNELL(9)       *
 *                                                                             *
 * (1) THEORETICAL CHEMISTRY GROUP - UNIVERSITA' DI TORINO - TORINO (ITALY)    *
 *      http://www.crystal.unito.it                                            *
 * (2) COMPUTATIONAL SCIENCE & ENGINEERING DEPARTMENT - CCLRC DARESBURY (UK)   *
 *     http://www.cse.clrc.ac.uk/cmg/CRYSTAL/                                  *
 * (3) UNIVERSITA' DEL PIEMONTE ORIENTALE - ALESSANDRIA (ITALY)                *
 * (4) UNIVERSIDAD AUTONOMA DEL ESTADO DE MORELOS - CUERNAVACA (MEXICO)        *
 * (5) UNIVERSITE' HENRI POINCARE' - NANCY (FRANCE)                            *
 * (6) MPI FUER FESTKOERPERFORSCHUNG - STUTTGART (GERMANY)                     *
 * (7) IMPERIAL COLLEGE - LONDON (UK)                                          *
 * (8) UNIVERSITE' PIERRE ET MARIE CURIE - PARIS (FRANCE)                      *
 * (9) UNIVERSIDAD DE BARCELONA - BARCELONA (SPAIN)                            *
 *(10) NIS - NANOSTRUCTURED INTERFACES AND SURFACES - TORINO (ITALY)           *
 *     http://www.crystal.unito.it                                             *
 *******************************************************************************
 EEEEEEEEEE STARTING  DATE 23 10 2007 TIME 01:47:27.0
 rutile                                                                          

 CRYSTAL CALCULATION
 (INPUT ACCORDING TO THE INTERNATIONAL TABLES FOR X-RAY CRYSTALLOGRAPHY)
 CRYSTAL FAMILY                       :  TETRAGONAL  
 CRYSTAL CLASS  (GROTH - 1921)        :  DITETRAGONAL DIPYRAMIDAL             

 SPACE GROUP (CENTROSYMMETRIC)        :  P 42/M N M      

 LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - CONVENTIONAL CELL
        A           B           C        ALPHA        BETA       GAMMA
     4.63910     4.63910     2.97938    90.00000    90.00000    90.00000


 NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL:    2

 INPUT COORDINATES

 ATOM AT. N.              COORDINATES
   1  22     0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
   2   8     3.061526467783E-01  3.061526467783E-01  0.000000000000E+00

 *******************************************************************************

 << INFORMATION >>: FROM NOW ON, ALL COORDINATES REFER TO THE PRIMITIVE CELL

 *******************************************************************************

 LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL
       A          B          C         ALPHA     BETA     GAMMA        VOLUME
    4.63910    4.63910    2.97938    90.00000  90.00000  90.00000     64.120029

 COORDINATES OF THE EQUIVALENT ATOMS (FRACTIONARY UNITS)

 N. ATOM EQUIV AT. N.          X                  Y                  Z

   1   1   1   22 TI    0.00000000000E+00  0.00000000000E+00  0.00000000000E+00
   2   1   2   22 TI   -5.00000000000E-01 -5.00000000000E-01 -5.00000000000E-01

   3   2   1    8 O     3.06152646778E-01  3.06152646778E-01  0.00000000000E+00
   4   2   2    8 O    -3.06152646778E-01 -3.06152646778E-01  0.00000000000E+00
   5   2   3    8 O    -1.93847353222E-01  1.93847353222E-01 -5.00000000000E-01
   6   2   4    8 O     1.93847353222E-01 -1.93847353222E-01 -5.00000000000E-01

 NUMBER OF SYMMETRY OPERATORS         :   16
 *******************************************************************************
 * GEOMETRY EDITING - INPUT COORDINATES ARE GIVEN IN ANGSTROM          
 *******************************************************************************

 GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP


 *******************************************************************************
  * CELL ROTATION 
 *******************************************************************************

 DIRECT-LATTICE FUNDAMENTAL VECTORS 
             X         Y         Z
 A1       4.6391    0.0000    0.0000
 A2       0.0000    4.6391    0.0000
 A3       0.0000    0.0000    2.9794
 **********************************************************************
 *            DEFINITION OF THE NEW LATTICE VECTORS                   *
 *            TWO OF WHICH BELONG TO THE SELECTED PLANE               *
 **********************************************************************

                       PLANE INDICES =   1   1   0

 NEW FUNDAMENTAL DIRECT-LATTICE VECTORS B1,B2,B3
 B1=      0 A1    0 A2    1 A3 
 B2=     -1 A1    1 A2    0 A3 
 B3=      0 A1   -1 A2    0 A3 

 LATTICE PARAMETERS  (ANGSTROM  AND DEGREES)
        A           B           C        ALPHA        BETA       GAMMA
     2.97938     6.56068     4.63910   135.00000    90.00000    90.00000

 NUMBER OF PRIMITIVE CELLS CONTAINED IN THE NEW FUNDAMENTAL CELL :  1

 *** THE CARTESIAN FRAME IS ROTATED SO THAT
     THE SELECTED PLANE IS IN THE X-Y PLANE

 NEW FUNDAMENTAL DIRECT-LATTICE VECTORS B1,B2,B3
             X         Y         Z
 B1       2.9794    0.0000    0.0000
 B2       0.0000    6.5607    0.0000
 B3       0.0000   -3.2803    3.2803

 VOLUME OF THE 3D CELL      64.120
 AREA OF THE 2D CELL        19.547
 B3 MODULUS                  4.639
 B3 Z PROJECTION             3.280
 B3 X-Y PROJECTION           3.280

 UNIT CELL ATOM COORDINATES :
 LAB AT.NO.       CARTESIAN (ANG)               CRYSTALLOGRAPHIC
   1  22       0.000     0.000     0.000     0.000000   0.000000   0.000000
   2  22       1.490    -3.280     0.000     0.500000  -0.500000   0.000000
   3   8       0.000    -3.280     1.272     0.000000  -0.306153   0.387695
   4   8       0.000     0.000     2.009     0.000000   0.306153   0.612305
   5   8       1.490     1.272     0.000     0.500000   0.193847   0.000000
   6   8       1.490    -1.272     0.000     0.500000  -0.193847   0.000000
   1  22       0.000     0.000     0.000     0.000000   0.000000   0.000000
   2  22       1.490     3.280     0.000     0.500000   0.500000   0.000000
   3   8       0.000    -3.280     1.272     0.000000  -0.306153   0.387695
   4   8       0.000     0.000     2.009     0.000000   0.306153   0.612305
   5   8       1.490     1.272     0.000     0.500000   0.193847   0.000000
   6   8       1.490    -1.272     0.000     0.500000  -0.193847   0.000000

 ATOMS CLASSIFIED ACCORDING TO THE Z COORDINATE :
  LAYER   1 Z=   2.0086; LABEL AT.NO.              X,Y,Z (ANG.)
                           4     8        0.00000     0.00000     2.00857
  LAYER   2 Z=   1.2718; LABEL AT.NO.              X,Y,Z (ANG.)
                           3     8        0.00000    -3.28034     1.27177
  LAYER   3 Z=   0.0000; LABEL AT.NO.              X,Y,Z (ANG.)
                           1    22        0.00000     0.00000     0.00000
                           2    22        1.48969     3.28034     0.00000
                           5     8        1.48969     1.27177     0.00000
                           6     8        1.48969    -1.27177     0.00000

 *************************** CELL ROTATION COMPLETE ***************************


 *******************************************************************************
 * TWO DIMENSIONAL SLAB PARALLEL TO THE SELECTED PLANE
 *******************************************************************************
         SURFACE LAYER NO.  2
         TOTAL NUMBER OF LAYERS:   3

 ORIGIN SHIFT TO MAXIMIZE USE OF SYMMETRY=     0.000

 NUMBER OF SYMMOPS=  8

 COORDINATES OF THE ATOMS BELONGING TO THE SLAB
 LAB AT.N.     CARTESIAN (ANG)        FRACTIONARY (3D)      FRAC (2D)     Z(ANG)
   1   8  0.00000 -3.28034  1.27177  0.000 -0.306  0.388  0.000 -0.500   1.27177
   2  22  0.00000  0.00000  0.00000  0.000  0.000  0.000  0.000  0.000   0.00000
   3  22 -1.48969 -3.28034  0.00000  0.500  0.500  0.000 -0.500 -0.500   0.00000
   4   8 -1.48969  1.27177  0.00000  0.500  0.194  0.000 -0.500  0.194   0.00000
   5   8 -1.48969 -1.27177  0.00000  0.500 -0.194  0.000 -0.500 -0.194   0.00000
   6   8  0.00000  3.28034 -1.27177  0.000  0.306 -0.388  0.000  0.500  -1.27177


 ******************************* SLAB GENERATED *******************************
 *******************************************************************************
 *  GEOMETRY EDITING
 *  THE SYMMETRY OF THE CRYSTAL IS BROKEN BY THE PERTURBATION
 *******************************************************************************

 *******************************************************************************
 * DISPLACEMENT OF    6 ATOMS
 *******************************************************************************

 ATOM N.    1 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    2 AT. N.  22 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    3 AT. N.  22 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    4 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    5 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    6 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000
 *******************************************************************************PROCESS             6  OF             8  WORKING

 *  GEOMETRY EDITING
 *  THE SYMMETRY OF THE CRYSTAL IS BROKEN BY THE PERTURBATION
 *******************************************************************************

 *******************************************************************************
 * DISPLACEMENT OF    6 ATOMS
 *******************************************************************************

 ATOM N.    1 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    2 AT. N.  22 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    3 AT. N.  22 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    4 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000PROCESS             4  OF             8  WORKING
PROCESS             7  OF             8  WORKING


 ATOM N.    5 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000

 ATOM N.    6 AT. N.   8 DISPLACED BY (A)   0.00000   0.00000   0.00000
 INFORMATION **** CVOLOPT **** CONSTANT VOLUME OPTIMIZATION
 INFORMATION **** MAXCYCLE ****   MAXIMUM NUMBER OF OPTIMIZATION CYCLES        10
 *******************************************************************************
 ATOMIC POSITIONS OPTIMIZATION CONTROL

 MAXIMUM GRADIENT COMPONENT   0.00045 MAXIMUM DISPLACEMENT COMPONENT     0.00030
 R.M.S. OF GRADIENT COMPONENT 0.00180 R.M.S. OF DISPLACEMENT COMPONENTS  0.00120
 THRESHOLD ON ENERGY CHANGE 0.100E-06 EXTRAPOLATING POLYNOMIAL ORDER           2
 MAXIMUM ALLOWED NUMBER OF STEPS   10 SORTING OF ENERGY POINTS:               NO
 ANALYTICAL  GRADIENT                 HESSIAN UPDATING                      BFGS
 STEP SIZE NUMERICAL GRADIENT 0.00100
 INITIAL HESSIAN MATRIX: SCHLEGEL MODEL 2
 *******************************************************************************


 GCALCO - MAX INDICES DIRECT LATTICE VECTOR 123  56   0
 NO.OF VECTORS CREATED 6999 STARS 1801 RMAX   394.73687 BOHR

 GEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM    2
 (NON PERIODIC DIRECTION: LATTICE PARAMETER FORMALLY SET TO 500)
 *******************************************************************************
 LATTICE PARAMETERS (ANGSTROMS AND DEGREES) - BOHR = 0.5291772083 ANGSTROM
 PRIMITIVE CELL
         A              B              C           ALPHA      BETA       GAMMA 
     2.97938395     6.56067637   500.00000000    90.000000  90.000000  90.000000
 *******************************************************************************
 ATOMS IN THE ASYMMETRIC UNIT    4 - ATOMS IN THE UNIT CELL:    6
     ATOM              X/A                 Y/B             Z(ANGSTROM)
 *******************************************************************************
   1 T   8 O     0.000000000000E+00 -5.000000000000E-01  1.271769749560E+00
   2 T  22 TI    0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
   3 T  22 TI    5.000000000000E-01 -5.000000000000E-01  0.000000000000E+00
   4 T   8 O     5.000000000000E-01  1.938473532217E-01  0.000000000000E+00
   5 F   8 O     5.000000000000E-01 -1.938473532217E-01  0.000000000000E+00
   6 F   8 O     0.000000000000E+00 -5.000000000000E-01 -1.271769749560E+00

 T = ATOM BELONGING TO THE ASYMMETRIC UNIT

 ****   8 SYMMOPS - TRANSLATORS IN FRACTIONARY UNITS
  V INV                   ROTATION MATRICES                       TRANSLATOR
  1  1  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00  1.00    0.00  0.00  0.00
  2  2  1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00 -1.00    0.00  0.00  0.00
  3  3 -1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00 -1.00    0.00  0.00  0.00
  4  4 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00  1.00    0.00  0.00  0.00
  5  5 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00 -1.00    0.00  0.00  0.00
  6  6 -1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00  1.00    0.00  0.00  0.00
  7  7  1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00  1.00    0.00  0.00  0.00
  8  8  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00 -1.00    0.00  0.00  0.00

 DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)
          X                    Y                    Z
   0.297938395000E+01   0.000000000000E+00   0.000000000000E+00
   0.000000000000E+00   0.656067636944E+01   0.000000000000E+00
   0.000000000000E+00   0.000000000000E+00   0.500000000000E+03


 CARTESIAN COORDINATES - PRIMITIVE CELL
 *******************************************************************************
 *      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)
 *******************************************************************************
   1     8 O     0.000000000000E+00  3.280338184719E+00  1.271769749560E+00
   2    22 TI    0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
   3    22 TI    1.489691975000E+00  3.280338184719E+00  0.000000000000E+00
   4     8 O     1.489691975000E+00  1.271769749560E+00  0.000000000000E+00
   5     8 O     1.489691975000E+00 -1.271769749560E+00  0.000000000000E+00
   6     8 O     0.000000000000E+00  3.280338184719E+00 -1.271769749560E+00

 *******************************************************************************
 LOCAL ATOMIC FUNCTIONS BASIS SET
 *******************************************************************************
   ATOM  X(AU)  Y(AU)  Z(AU)    NO. TYPE  EXPONENT  S COEF   P COEF   D/F/G COEF
 *******************************************************************************
   1 O   0.000  6.199  2.403
                                  1 S  
                                         8.020E+03 1.080E-03 0.000E+00 0.000E+00
                                         1.338E+03 8.040E-03 0.000E+00 0.000E+00
                                         2.554E+02 5.324E-02 0.000E+00 0.000E+00
                                         6.922E+01 1.681E-01 0.000E+00 0.000E+00
                                         2.390E+01 3.581E-01 0.000E+00 0.000E+00
                                         9.264E+00 3.855E-01 0.000E+00 0.000E+00
                                         3.851E+00 1.468E-01 0.000E+00 0.000E+00
                                         1.212E+00 7.280E-02 0.000E+00 0.000E+00
                             2-   5 SP 
                                         4.943E+01-8.830E-03 9.580E-03 0.000E+00
                                         1.047E+01-9.150E-02 6.960E-02 0.000E+00
                                         3.235E+00-4.020E-02 2.065E-01 0.000E+00
                                         1.217E+00 3.790E-01 3.470E-01 0.000E+00
                             6-   9 SP 
                                         4.567E-01 1.000E+00 1.000E+00 0.000E+00
                            10-  13 SP 
                                         1.843E-01 1.000E+00 1.000E+00 0.000E+00
                            14-  18 D  
                                         6.000E-01 0.000E+00 0.000E+00 1.000E+00
   2 TI  0.000  0.000  0.000
                                 19 S  
                                         2.253E+05 2.280E-04 0.000E+00 0.000E+00
                                         3.232E+04 1.929E-03 0.000E+00 0.000E+00
                                         6.884E+03 1.110E-02 0.000E+00 0.000E+00
                                         1.802E+03 5.000E-02 0.000E+00 0.000E+00
                                         5.431E+02 1.701E-01 0.000E+00 0.000E+00
                                         1.875E+02 3.690E-01 0.000E+00 0.000E+00
                                         7.321E+01 4.033E-01 0.000E+00 0.000E+00
                                         3.037E+01 1.445E-01 0.000E+00 0.000E+00
                            20-  23 SP 
                                         5.540E+02-5.900E-03 8.500E-03 0.000E+00
                                         1.325E+02-6.830E-02 6.030E-02 0.000E+00
                                         4.368E+01-1.245E-01 2.124E-01 0.000E+00
                                         1.722E+01 2.532E-01 3.902E-01 0.000E+00
                                         7.225E+00 6.261E-01 4.097E-01 0.000E+00
                                         2.412E+00 2.820E-01 2.181E-01 0.000E+00
                            24-  27 SP 
                                         2.450E+01 1.750E-02-2.070E-02 0.000E+00
                                         1.148E+01-2.277E-01-6.530E-02 0.000E+00
                                         4.465E+00-7.946E-01 1.919E-01 0.000E+00
                                         1.890E+00 1.011E+00 1.378E+00 0.000E+00
                            28-  31 SP 
                                         8.126E-01 1.000E+00 1.000E+00 0.000E+00
                            32-  35 SP 
                                         3.297E-01 1.000E+00 1.000E+00 0.000E+00
                            36-  40 D  
                                         1.627E+01 0.000E+00 0.000E+00 6.750E-02
                                         4.372E+00 0.000E+00 0.000E+00 2.934E-01
                                         1.464E+00 0.000E+00 0.000E+00 5.658E-01
                                         5.485E-01 0.000E+00 0.000E+00 5.450E-01
                            41-  45 D  
                                         2.600E-01 0.000E+00 0.000E+00 1.000E+00
   3 TI  2.815  6.199  0.000
   4 O   2.815  2.403  0.000
   5 O   2.815 -2.403  0.000
   6 O   0.000  6.199 -2.403
 INFORMATION **** TOLINTEG **** COULOMB AND EXCHANGE SERIES TOLERANCES MODIFIED
 INFORMATION **** PPAN ****  MULLIKEN POPULATION ANALYSIS AT THE END OF SCF
 INFORMATION **** TOLDEE **** SCF TOL ON TOTAL ENERGY SET TO         7
 INFORMATION **** EXCHSIZE **** EXCHANGE BIPOLAR BUFFER SIZE SET TO   6937578
 INFORMATION **** BIPOSIZE **** COULOMB BIPOLAR BUFFER SET TO   9701800
 INFORMATION **** READM2 **** FULL DIRECT SCF (MONO AND BIEL INT) SELECTED
 *******************************************************************************
 N. OF ATOMS PER CELL        6  COULOMB OVERLAP TOL         (T1) 10** -7
 NUMBER OF SHELLS           34  COULOMB PENETRATION TOL     (T2) 10** -7
 NUMBER OF AO              126  EXCHANGE OVERLAP TOL        (T3) 10** -7
 N. OF ELECTRONS PER CELL   76  EXCHANGE PSEUDO OVP (F(G))  (T4) 10** -7
 CORE ELECTRONS PER CELL    44  EXCHANGE PSEUDO OVP (P(G))  (T5) 10**-14
 N. OF SYMMETRY OPERATORS    8  POLE ORDER IN MONO ZONE                4
 *******************************************************************************
 TYPE OF CALCULATION :  RESTRICTED CLOSED SHELL
 KOHN-SHAM HAMILTONIAN

 (EXCHANGE)[CORRELATION] FUNCTIONAL:(BECKE)[LEE-YANG-PARR]

 NON-LOCAL WEIGHTING FACTOR (EXCHANGE) =           0.9000
 NON-LOCAL WEIGHTING FACTOR [CORRELATION] =        0.8100

 HYBRID EXCHANGE - PERCENTAGE OF FOCK EXCHANGE    20.0000
 EIGENVALUE LEVEL SHIFTING OF  0.500 HARTREE
 LOCKING - FERMI ENERGY ALTERED BY LEVEL SHIFTER

 CAPPA: IS1=  8 IS2=  8 IS3=  1K POINTS MONKHORST NET   25 SYMMOPS K SPACE  8 SYMMOPS G SPACE  8


 CAPPA1: ISJ1=  8 ISJ2=  8 ISJ3=  1K POINTS GILAT NET   25 SYMMOPS K SPACE  8 SYMMOPS G SPACE  8

 *******************************************************************************
 MAX NUMBER OF SCF CYCLES      50  CONVERGENCE ON DELTAP        10**-17
 WEIGHT OF F(I) IN F(I+1)      30% CONVERGENCE ON ENERGY        10**- 7
 EIGENVALUE LEVEL SHIFTING OF  0.500 HARTREE
 LOCKING - FERMI ENERGY ALTERED BY LEVEL SHIFTER
 SHRINK. FACT.(MONKH.)    8  8  1  NUMBER OF K POINTS IN THE IBZ     25
 SHRINKING FACTOR(GILAT NET)    8  NUMBER OF K POINTS(GILAT NET)     25
 *******************************************************************************
 *** K POINTS COORDINATES (OBLIQUE COORDINATES IN UNITS OF IS =  8)
   1-R(  0  0  0)   2-C(  1  0  0)   3-C(  2  0  0)   4-C(  3  0  0)
   5-R(  4  0  0)   6-C(  0  1  0)   7-C(  1  1  0)   8-C(  2  1  0)
   9-C(  3  1  0)  10-C(  4  1  0)  11-C(  0  2  0)  12-C(  1  2  0)
  13-C(  2  2  0)  14-C(  3  2  0)  15-C(  4  2  0)  16-C(  0  3  0)
  17-C(  1  3  0)  18-C(  2  3  0)  19-C(  3  3  0)  20-C(  4  3  0)
  21-R(  0  4  0)  22-C(  1  4  0)  23-C(  2  4  0)  24-C(  3  4  0)
  25-R(  4  4  0)

 DIRECT LATTICE VECTORS COMPON. (A.U.)     RECIP. LATTICE VECTORS COMPON. (A.U.) 
        X            Y            Z              X            Y            Z
    5.6302197    0.0000000    0.0000000      1.1159752    0.0000000    0.0000000
    0.0000000   12.3978816    0.0000000      0.0000000    0.5067951    0.0000000
    0.0000000    0.0000000  944.8630670      0.0000000    0.0000000    0.0066498

 DISK SPACE FOR EIGENVECTORS (FTN 10)     730296 REALS

 SYMMETRY ADAPTION OF THE BLOCH FUNCTIONS ENABLED
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT gordsh1     TELAPSE        0.40 TCPU        0.07

 DIMENSIONS  P(G)=   96172 F(G)=   23448 P(G),F(G) (IRR)   13009
 MAX G-VECTOR INDEX FOR 1- AND 2-ELECTRON INTEGRALS  37

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INPUT       TELAPSE        0.40 TCPU        0.07

 NEIGHBORS OF THE NON-EQUIVALENT ATOMS

 N = NUMBER OF NEIGHBORS AT DISTANCE R
    ATOM  N     R/ANG      R/AU   NEIGHBORS (ATOM LABELS AND CELL INDICES)
   1 O    2     1.9587     3.7014   3 TI   0 0 0   3 TI  -1 0 0
   1 O    1     2.5435     4.8066   6 O    0 0 0
   1 O    4     2.8055     5.3017   4 O    0 0 0   4 O   -1 0 0   5 O   -1 1 0
                                    5 O    0 1 0
   1 O    2     2.9794     5.6302   1 O   -1 0 0   1 O    1 0 0
   1 O    2     3.5182     6.6485   2 TI   0 0 0   2 TI   0 1 0
   1 O    2     3.9174     7.4029   6 O    1 0 0   6 O   -1 0 0

   2 TI   4     1.9587     3.7014   4 O    0 0 0   4 O   -1 0 0   5 O    0 0 0
                                    5 O   -1 0 0
   2 TI   2     2.9794     5.6302   2 TI  -1 0 0   2 TI   1 0 0
   2 TI   4     3.5182     6.6485   1 O    0 0 0   1 O    0-1 0   6 O    0-1 0
                                    6 O    0 0 0
   2 TI   4     3.6027     6.8082   3 TI   0 0 0   3 TI  -1 0 0   3 TI   0-1 0
                                    3 TI  -1-1 0
   2 TI   8     4.6103     8.7122   1 O   -1 0 0   1 O    1 0 0   1 O   -1-1 0
                                    1 O    1-1 0   6 O   -1-1 0   6 O    1-1 0
                                    6 O   -1 0 0   6 O    1 0 0
   2 TI   4     4.6465     8.7806   4 O    1 0 0   4 O   -2 0 0   5 O    1 0 0
                                    5 O   -2 0 0

   3 TI   4     1.9587     3.7014   1 O    0 0 0   1 O    1 0 0   6 O    0 0 0
                                    6 O    1 0 0
   3 TI   2     2.0086     3.7956   4 O    0 0 0   5 O    0 1 0
   3 TI   2     2.9794     5.6302   3 TI  -1 0 0   3 TI   1 0 0
   3 TI   4     3.5932     6.7902   4 O    1 0 0   4 O   -1 0 0   5 O    1 1 0
                                    5 O   -1 1 0
   3 TI   4     3.6027     6.8082   2 TI   0 0 0   2 TI   1 0 0   2 TI   0 1 0
                                    2 TI   1 1 0
   3 TI   2     4.5521     8.6022   4 O    0 1 0   5 O    0 0 0

   4 O    2     1.9587     3.7014   2 TI   0 0 0   2 TI   1 0 0
   4 O    1     2.0086     3.7956   3 TI   0 0 0
   4 O    1     2.5435     4.8066   5 O    0 0 0
   4 O    4     2.8055     5.3017   1 O    0 0 0   1 O    1 0 0   6 O    1 0 0
                                    6 O    0 0 0
   4 O    2     2.9794     5.6302   4 O   -1 0 0   4 O    1 0 0
   4 O    2     3.5932     6.7902   3 TI  -1 0 0   3 TI   1 0 0


 SYMMETRY ALLOWED INTERNAL DEGREE(S) OF FREEDOM:   2
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SYMM        TELAPSE        0.40 TCPU        0.07
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_SCREEN  TELAPSE        0.40 TCPU        0.07
 *******************************************************************************
 COORDINATE OPTIMIZATION - POINT    1
 INFORMATION **** EXCBUF **** EXCH. BIPO BUFFER: WORDS USED =    468954

 DFT PARAMETERS

     ATOM       ELECTRONS   NET CHARGE   R(ANGSTROM)
   1   8  O       9.0000     -1.0000     1.07000000
   2  22  TI     20.0000      2.0000     0.86000000

 SIZE OF GRID=      13081
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MAKEGRID    TELAPSE        0.72 TCPU        0.37
 BECKE WEIGHT FUNCTION
 RADSAFE =     2.00
 TOLERANCES - DENSITY:10**- 6; POTENTIAL:10**- 9; GRID WGT:10**-14

 RADIAL INTEGRATION  - INTERVALS (POINTS,UPPER LIMIT):            1( 75,  4.0*R)

 ANGULAR INTEGRATION - INTERVALS (ACCURACY LEVEL [N. POINTS] UPPER LIMIT):
  1(  2[  50]   0.2)  2(  6[ 146]   0.5)  3(  8[ 194]   0.9)  4( 13[ 434]   3.5)
  5(  8[ 194]9999.0)
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_CALC    TELAPSE        0.72 TCPU        0.37

 *******************************************************************************
 rutile                                                                          
 CRYSTAL - SCF - TYPE OF CALCULATION :  RESTRICTED CLOSED SHELL
 *******************************************************************************

 CAPPA: IS1=  8 IS2=  8 IS3=  1K POINTS MONKHORST NET   25 SYMMOPS K SPACE  8 SYMMOPS G SPACE  8


 CAPPA1: ISJ1=  8 ISJ2=  8 ISJ3=  1K POINTS GILAT NET   25 SYMMOPS K SPACE  8 SYMMOPS G SPACE  8

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SDIK        TELAPSE        0.74 TCPU        0.38

 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         ATOMIC WAVEFUNCTION(S)


 NUCLEAR CHARGE  8.0  SYMMETRY SPECIES            S    P
 N. ELECTRONS    9.0  NUMBER OF PRIMITIVE GTOS   14    6
                      NUMBER OF CONTRACTED GTOS   4    3
                      NUMBER OF CLOSED SHELLS     2    0
                      OPEN SHELL OCCUPATION       0    5

  ZNUC SCFIT  TOTAL HF ENERGY   KINETIC ENERGY   VIRIAL THEOREM ACCURACY
   8.0  24   -7.476280579E+01  7.472109493E+01 -2.000558221E+00  3.9E-06

 NUCLEAR CHARGE 22.0  SYMMETRY SPECIES            S    P    D
 N. ELECTRONS   20.0  NUMBER OF PRIMITIVE GTOS   20   12    5
                      NUMBER OF CONTRACTED GTOS   5    4    2
                      NUMBER OF CLOSED SHELLS     3    2    0
                      OPEN SHELL OCCUPATION       0    0    2

  ZNUC SCFIT  TOTAL HF ENERGY   KINETIC ENERGY   VIRIAL THEOREM ACCURACY
  22.0  37   -8.476906516E+02  8.474944897E+02 -2.000231461E+00  4.8E-06

 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0000000  20.0000000  20.0000000   9.0000000   9.0000000   9.0000000
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE        0.77 TCPU        0.40
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE        3.74 TCPU        3.37
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE        5.60 TCPU        5.23
 NUMERICALLY INTEGRATED DENSITY     76.0000046416
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE        5.88 TCPU        5.49
 CYC   0 ETOT(AU) -2.000159471162E+03 DETOT -2.00E+03 tst  0.00E+00 PX  1.00E+00
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE        5.93 TCPU        5.54
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -2.8777787E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE        5.95 TCPU        5.57
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.5606732  18.8056073  18.9509933   9.5610265   9.5610265   9.5606732
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE        5.96 TCPU        5.58
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE        9.07 TCPU        8.65
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       11.02 TCPU       10.61
 NUMERICALLY INTEGRATED DENSITY     76.0000068932
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       11.45 TCPU       11.02
 CYC   1 ETOT(AU) -1.999065827080E+03 DETOT  1.09E+00 tst  0.00E+00 PX  1.00E+00
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       12.01 TCPU       11.05
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -6.0254838E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       12.05 TCPU       11.07
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   8.9986684  20.1205811  19.9925833   8.9447495   8.9447495   8.9986684
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       12.06 TCPU       11.08
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       15.20 TCPU       14.16
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       17.16 TCPU       16.11
 NUMERICALLY INTEGRATED DENSITY     76.0000072965
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       17.59 TCPU       16.53
 CYC   2 ETOT(AU) -1.999740347977E+03 DETOT -6.75E-01 tst  4.16E+00 PX  3.03E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       17.74 TCPU       16.57
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.6476804E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       17.76 TCPU       16.59
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0649325  19.8300222  19.8717512   9.0841808   9.0841808   9.0649325
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       17.77 TCPU       16.60
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       20.90 TCPU       19.71
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       22.79 TCPU       21.59
 NUMERICALLY INTEGRATED DENSITY     76.0000066483
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       23.22 TCPU       22.01
 CYC   3 ETOT(AU) -1.999922940386E+03 DETOT -1.83E-01 tst  9.25E-03 PX  1.03E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       23.55 TCPU       22.30
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4048302E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       23.57 TCPU       22.33
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0700892  19.8449142  19.8549593   9.0799741   9.0799741   9.0700892
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       23.59 TCPU       22.34
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       26.72 TCPU       25.45
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       28.60 TCPU       27.33
 NUMERICALLY INTEGRATED DENSITY     76.0000063589
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       29.02 TCPU       27.73
 CYC   4 ETOT(AU) -1.999924970743E+03 DETOT -2.03E-03 tst  5.28E-03 PX  1.46E-02
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       29.08 TCPU       27.76
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3046629E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       29.11 TCPU       27.78
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0426452  19.9007875  19.8956040   9.0591590   9.0591590   9.0426452
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       29.12 TCPU       27.80
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       32.21 TCPU       30.87
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       34.13 TCPU       32.79
 NUMERICALLY INTEGRATED DENSITY     76.0000062677
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       34.55 TCPU       33.21
 CYC   5 ETOT(AU) -1.999926448486E+03 DETOT -1.48E-03 tst  2.39E-02 PX  1.24E-02
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       34.60 TCPU       33.24
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3819954E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       34.62 TCPU       33.26
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0433052  19.8876492  19.8899010   9.0679197   9.0679197   9.0433052
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       34.63 TCPU       33.27
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       37.73 TCPU       36.35
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       39.65 TCPU       38.27
 NUMERICALLY INTEGRATED DENSITY     76.0000061497
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       40.07 TCPU       38.67
 CYC   6 ETOT(AU) -1.999927238395E+03 DETOT -7.90E-04 tst  2.08E-03 PX  8.56E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       40.27 TCPU       38.73
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3782032E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       40.30 TCPU       38.76
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0409098  19.8914844  19.8898266   9.0684347   9.0684347   9.0409098
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       40.31 TCPU       38.77
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       43.44 TCPU       41.88
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       45.31 TCPU       43.76
 NUMERICALLY INTEGRATED DENSITY     76.0000060579
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       45.75 TCPU       44.18
 CYC   7 ETOT(AU) -1.999927510606E+03 DETOT -2.72E-04 tst  2.22E-03 PX  3.66E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       45.85 TCPU       44.21
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3806742E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       45.87 TCPU       44.23
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0376926  19.8950272  19.8917937   9.0688969   9.0688969   9.0376926
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       45.88 TCPU       44.24
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       48.98 TCPU       47.32
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       50.90 TCPU       49.25
 NUMERICALLY INTEGRATED DENSITY     76.0000059915
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       51.33 TCPU       49.66
 CYC   8 ETOT(AU) -1.999927687161E+03 DETOT -1.77E-04 tst  1.78E-03 PX  3.15E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       51.66 TCPU       49.71
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3885174E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       51.69 TCPU       49.73
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0360330  19.8955548  19.8917929   9.0702931   9.0702931   9.0360330
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       51.70 TCPU       49.74
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       54.80 TCPU       52.82
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       56.71 TCPU       54.74
 NUMERICALLY INTEGRATED DENSITY     76.0000059357
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       57.14 TCPU       55.15
 CYC   9 ETOT(AU) -1.999927810815E+03 DETOT -1.24E-04 tst  1.24E-03 PX  2.27E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       57.30 TCPU       55.28
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3921899E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       57.33 TCPU       55.31
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0345482  19.8966362  19.8919677   9.0711498   9.0711498   9.0345482
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       57.34 TCPU       55.32
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       60.68 TCPU       58.40
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       62.60 TCPU       60.31
 NUMERICALLY INTEGRATED DENSITY     76.0000058895
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       63.02 TCPU       60.72
 CYC  10 ETOT(AU) -1.999927897023E+03 DETOT -8.62E-05 tst  9.89E-04 PX  2.03E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       63.09 TCPU       60.74
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3954903E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       63.11 TCPU       60.77
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0332900  19.8974449  19.8921803   9.0718975   9.0718975   9.0332900
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       63.13 TCPU       60.78
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       66.47 TCPU       63.87
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       68.42 TCPU       65.81
 NUMERICALLY INTEGRATED DENSITY     76.0000058514
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       68.85 TCPU       66.22
 CYC  11 ETOT(AU) -1.999927959672E+03 DETOT -6.26E-05 tst  7.86E-04 PX  1.72E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       68.90 TCPU       66.25
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.3983756E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       68.93 TCPU       66.27
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0323009  19.8979894  19.8922810   9.0725639   9.0725639   9.0323009
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       68.94 TCPU       66.28
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       72.05 TCPU       69.37
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       74.00 TCPU       71.31
 NUMERICALLY INTEGRATED DENSITY     76.0000058196
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       74.43 TCPU       71.72
 CYC  12 ETOT(AU) -1.999928006358E+03 DETOT -4.67E-05 tst  6.30E-04 PX  1.47E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       74.47 TCPU       71.75
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4006139E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       74.49 TCPU       71.77
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0314804  19.8984553  19.8923602   9.0731119   9.0731119   9.0314804
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       74.50 TCPU       71.78
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       77.60 TCPU       74.87
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       79.52 TCPU       76.80
 NUMERICALLY INTEGRATED DENSITY     76.0000057928
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       79.95 TCPU       77.21
 CYC  13 ETOT(AU) -1.999928041778E+03 DETOT -3.54E-05 tst  5.14E-04 PX  1.28E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       80.05 TCPU       77.27
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4024718E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       80.07 TCPU       77.30
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0307982  19.8988240  19.8924299   9.0735749   9.0735749   9.0307982
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       80.09 TCPU       77.31
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       83.20 TCPU       80.40
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       85.11 TCPU       82.30
 NUMERICALLY INTEGRATED DENSITY     76.0000057700
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       85.53 TCPU       82.71
 CYC  14 ETOT(AU) -1.999928069269E+03 DETOT -2.75E-05 tst  4.24E-04 PX  1.12E-03
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       85.62 TCPU       82.76
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4040006E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       85.65 TCPU       82.79
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0302249  19.8991277  19.8924896   9.0739664   9.0739664   9.0302249
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       85.66 TCPU       82.80
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       88.80 TCPU       85.90
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       90.69 TCPU       87.80
 NUMERICALLY INTEGRATED DENSITY     76.0000057506
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       91.12 TCPU       88.22
 CYC  15 ETOT(AU) -1.999928090945E+03 DETOT -2.17E-05 tst  3.70E-04 PX  9.93E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       91.20 TCPU       88.27
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4052661E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       91.23 TCPU       88.29
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0297480  19.8993810  19.8925375   9.0742927   9.0742927   9.0297480
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       91.24 TCPU       88.30
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE       94.58 TCPU       91.39
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE       96.55 TCPU       93.35
 NUMERICALLY INTEGRATED DENSITY     76.0000057338
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE       96.98 TCPU       93.77
 CYC  16 ETOT(AU) -1.999928108301E+03 DETOT -1.74E-05 tst  3.25E-04 PX  8.86E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE       97.33 TCPU       93.82
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4063133E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE       97.36 TCPU       93.84
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0293397  19.8995955  19.8925819   9.0745716   9.0745716   9.0293397
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE       97.37 TCPU       93.85
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      100.47 TCPU       96.93
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      102.44 TCPU       98.90
 NUMERICALLY INTEGRATED DENSITY     76.0000057193
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      102.86 TCPU       99.31
 CYC  17 ETOT(AU) -1.999928122410E+03 DETOT -1.41E-05 tst  2.88E-04 PX  7.97E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      102.93 TCPU       99.34
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4072031E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      102.95 TCPU       99.37
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0289953  19.8997748  19.8926218   9.0748064   9.0748064   9.0289953
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      102.96 TCPU       99.38
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      106.07 TCPU      102.47
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      108.03 TCPU      104.43
 NUMERICALLY INTEGRATED DENSITY     76.0000057066
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      108.46 TCPU      104.84
 CYC  18 ETOT(AU) -1.999928134012E+03 DETOT -1.16E-05 tst  2.56E-04 PX  7.19E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      109.17 TCPU      105.18
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4079394E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      109.25 TCPU      105.24
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0286957  19.8999352  19.8926621   9.0750056   9.0750056   9.0286957
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      109.26 TCPU      105.25
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      112.39 TCPU      108.34
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      114.30 TCPU      110.24
 NUMERICALLY INTEGRATED DENSITY     76.0000056955
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      114.72 TCPU      110.65
 CYC  19 ETOT(AU) -1.999928143673E+03 DETOT -9.66E-06 tst  2.30E-04 PX  6.55E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      114.75 TCPU      110.67
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4085644E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      114.78 TCPU      110.70
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0284400  19.9000731  19.8926953   9.0751758   9.0751758   9.0284400
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      114.79 TCPU      110.71
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      118.13 TCPU      113.80
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      120.04 TCPU      115.70
 NUMERICALLY INTEGRATED DENSITY     76.0000056857
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      120.47 TCPU      116.12
 CYC  20 ETOT(AU) -1.999928151779E+03 DETOT -8.11E-06 tst  2.07E-04 PX  5.98E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      120.51 TCPU      116.15
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4090962E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      120.53 TCPU      116.18
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0282169  19.9001965  19.8927251   9.0753223   9.0753223   9.0282169
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      120.55 TCPU      116.19
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      123.68 TCPU      119.26
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      125.57 TCPU      121.16
 NUMERICALLY INTEGRATED DENSITY     76.0000056770
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      126.00 TCPU      121.57
 CYC  21 ETOT(AU) -1.999928158649E+03 DETOT -6.87E-06 tst  1.87E-04 PX  5.49E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      126.07 TCPU      121.60
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4095589E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      126.09 TCPU      121.62
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0280219  19.9003023  19.8927555   9.0754492   9.0754492   9.0280219
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      126.10 TCPU      121.63
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      129.27 TCPU      124.71
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      131.24 TCPU      126.68
 NUMERICALLY INTEGRATED DENSITY     76.0000056692
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      131.67 TCPU      127.09
 CYC  22 ETOT(AU) -1.999928164503E+03 DETOT -5.85E-06 tst  1.70E-04 PX  5.06E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      131.73 TCPU      127.12
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4099496E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      131.76 TCPU      127.14
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0278508  19.9003991  19.8927815   9.0755589   9.0755589   9.0278508
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      131.77 TCPU      127.16
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      134.86 TCPU      130.24
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      136.83 TCPU      132.20
 NUMERICALLY INTEGRATED DENSITY     76.0000056623
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      137.26 TCPU      132.62
 CYC  23 ETOT(AU) -1.999928169524E+03 DETOT -5.02E-06 tst  1.54E-04 PX  4.68E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      137.53 TCPU      132.70
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4102923E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      137.59 TCPU      132.73
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0276984  19.9004845  19.8928097   9.0756545   9.0756545   9.0276984
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      137.60 TCPU      132.74
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      140.73 TCPU      135.85
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      142.62 TCPU      137.74
 NUMERICALLY INTEGRATED DENSITY     76.0000056561
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      143.05 TCPU      138.15
 CYC  24 ETOT(AU) -1.999928173854E+03 DETOT -4.33E-06 tst  1.41E-04 PX  4.33E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      143.15 TCPU      138.19
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4105882E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      143.18 TCPU      138.22
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0275638  19.9005626  19.8928313   9.0757392   9.0757392   9.0275638
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      143.19 TCPU      138.23
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      146.32 TCPU      141.34
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      148.23 TCPU      143.25
 NUMERICALLY INTEGRATED DENSITY     76.0000056505
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      148.65 TCPU      143.66
 CYC  25 ETOT(AU) -1.999928177611E+03 DETOT -3.76E-06 tst  1.29E-04 PX  4.03E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      148.71 TCPU      143.69
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4108491E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      148.73 TCPU      143.71
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0274430  19.9006345  19.8928529   9.0758134   9.0758134   9.0274430
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      148.74 TCPU      143.72
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      151.84 TCPU      146.80
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      153.81 TCPU      148.77
 NUMERICALLY INTEGRATED DENSITY     76.0000056454
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      154.23 TCPU      149.18
 CYC  26 ETOT(AU) -1.999928180883E+03 DETOT -3.27E-06 tst  1.18E-04 PX  3.75E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      154.32 TCPU      149.22
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4110804E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      154.34 TCPU      149.24
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0273354  19.9006988  19.8928726   9.0758789   9.0758789   9.0273354
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      154.35 TCPU      149.26
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      157.48 TCPU      152.33
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      159.42 TCPU      154.27
 NUMERICALLY INTEGRATED DENSITY     76.0000056408
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      159.84 TCPU      154.68
 CYC  27 ETOT(AU) -1.999928183750E+03 DETOT -2.87E-06 tst  1.09E-04 PX  3.50E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      159.92 TCPU      154.72
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4112814E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      159.94 TCPU      154.75
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0272381  19.9007582  19.8928905   9.0759376   9.0759376   9.0272381
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      159.96 TCPU      154.76
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      163.31 TCPU      157.85
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      165.22 TCPU      159.75
 NUMERICALLY INTEGRATED DENSITY     76.0000056367
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      165.65 TCPU      160.16
 CYC  28 ETOT(AU) -1.999928186267E+03 DETOT -2.52E-06 tst  9.98E-05 PX  3.27E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      165.72 TCPU      160.19
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4114607E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      165.75 TCPU      160.22
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0271496  19.9008138  19.8929090   9.0759890   9.0759890   9.0271496
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      165.76 TCPU      160.23
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      168.87 TCPU      163.31
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      170.77 TCPU      165.21
 NUMERICALLY INTEGRATED DENSITY     76.0000056328
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      171.19 TCPU      165.61
 CYC  29 ETOT(AU) -1.999928188488E+03 DETOT -2.22E-06 tst  9.19E-05 PX  3.06E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      171.28 TCPU      165.67
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4116206E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      171.31 TCPU      165.69
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0270701  19.9008642  19.8929244   9.0760356   9.0760356   9.0270701
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      171.32 TCPU      165.70
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      174.43 TCPU      168.78
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      176.37 TCPU      170.72
 NUMERICALLY INTEGRATED DENSITY     76.0000056294
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      176.80 TCPU      171.13
 CYC  30 ETOT(AU) -1.999928190448E+03 DETOT -1.96E-06 tst  8.48E-05 PX  2.88E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      176.89 TCPU      171.17
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4117648E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      176.92 TCPU      171.20
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0269961  19.9009103  19.8929399   9.0760788   9.0760788   9.0269961
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      176.93 TCPU      171.21
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      180.02 TCPU      174.29
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      181.93 TCPU      176.20
 NUMERICALLY INTEGRATED DENSITY     76.0000056262
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      182.36 TCPU      176.61
 CYC  31 ETOT(AU) -1.999928192187E+03 DETOT -1.74E-06 tst  7.83E-05 PX  2.70E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      182.42 TCPU      176.65
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4118949E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      182.44 TCPU      176.67
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0269289  19.9009546  19.8929534   9.0761171   9.0761171   9.0269289
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      182.46 TCPU      176.68
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      185.56 TCPU      179.77
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      187.48 TCPU      181.68
 NUMERICALLY INTEGRATED DENSITY     76.0000056233
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      187.91 TCPU      182.09
 CYC  32 ETOT(AU) -1.999928193731E+03 DETOT -1.54E-06 tst  7.23E-05 PX  2.54E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      187.99 TCPU      182.14
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4120108E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      188.02 TCPU      182.16
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0268678  19.9009962  19.8929642   9.0761520   9.0761520   9.0268678
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      188.03 TCPU      182.18
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      191.16 TCPU      185.25
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      193.12 TCPU      187.22
 NUMERICALLY INTEGRATED DENSITY     76.0000056206
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      193.55 TCPU      187.63
 CYC  33 ETOT(AU) -1.999928195110E+03 DETOT -1.38E-06 tst  6.69E-05 PX  2.39E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      193.67 TCPU      187.70
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4121154E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      193.70 TCPU      187.72
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0268111  19.9010338  19.8929749   9.0761845   9.0761845   9.0268111
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      193.71 TCPU      187.73
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      197.01 TCPU      190.81
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      198.97 TCPU      192.76
 NUMERICALLY INTEGRATED DENSITY     76.0000056181
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      199.39 TCPU      193.17
 CYC  34 ETOT(AU) -1.999928196344E+03 DETOT -1.23E-06 tst  6.19E-05 PX  2.25E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      199.48 TCPU      193.21
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4122129E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      199.50 TCPU      193.23
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0267589  19.9010697  19.8929852   9.0762137   9.0762137   9.0267589
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      199.51 TCPU      193.24
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      202.65 TCPU      196.36
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      204.53 TCPU      198.23
 NUMERICALLY INTEGRATED DENSITY     76.0000056159
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      204.96 TCPU      198.65
 CYC  35 ETOT(AU) -1.999928197442E+03 DETOT -1.10E-06 tst  5.73E-05 PX  2.12E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      205.03 TCPU      198.68
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4123026E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      205.06 TCPU      198.71
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0267099  19.9011031  19.8929953   9.0762409   9.0762409   9.0267099
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      205.07 TCPU      198.72
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      208.43 TCPU      201.82
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      210.32 TCPU      203.71
 NUMERICALLY INTEGRATED DENSITY     76.0000056138
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      210.76 TCPU      204.13
 CYC  36 ETOT(AU) -1.999928198435E+03 DETOT -9.93E-07 tst  5.31E-05 PX  2.01E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      211.05 TCPU      204.15
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4123816E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      211.08 TCPU      204.18
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0266641  19.9011352  19.8930054   9.0762656   9.0762656   9.0266641
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      211.09 TCPU      204.19
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      214.19 TCPU      207.27
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      216.09 TCPU      209.17
 NUMERICALLY INTEGRATED DENSITY     76.0000056118
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      216.52 TCPU      209.58
 CYC  37 ETOT(AU) -1.999928199326E+03 DETOT -8.92E-07 tst  4.92E-05 PX  1.90E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      216.61 TCPU      209.64
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4124569E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      216.63 TCPU      209.67
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0266224  19.9011645  19.8930134   9.0762887   9.0762887   9.0266224
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      216.64 TCPU      209.68
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      219.74 TCPU      212.76
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      221.65 TCPU      214.67
 NUMERICALLY INTEGRATED DENSITY     76.0000056100
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      222.08 TCPU      215.08
 CYC  38 ETOT(AU) -1.999928200126E+03 DETOT -8.00E-07 tst  4.56E-05 PX  1.79E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      222.12 TCPU      215.11
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4125244E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      222.14 TCPU      215.13
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0265829  19.9011929  19.8930214   9.0763099   9.0763099   9.0265829
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      222.15 TCPU      215.14
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      225.27 TCPU      218.22
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      227.18 TCPU      220.13
 NUMERICALLY INTEGRATED DENSITY     76.0000056084
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      227.61 TCPU      220.54
 CYC  39 ETOT(AU) -1.999928200848E+03 DETOT -7.22E-07 tst  4.23E-05 PX  1.70E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      228.00 TCPU      220.82
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4125861E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      228.03 TCPU      220.85
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0265458  19.9012195  19.8930296   9.0763296   9.0763296   9.0265458
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      228.04 TCPU      220.86
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      231.16 TCPU      223.96
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      233.05 TCPU      225.85
 NUMERICALLY INTEGRATED DENSITY     76.0000056068
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      233.48 TCPU      226.26
 CYC  40 ETOT(AU) -1.999928201507E+03 DETOT -6.59E-07 tst  3.93E-05 PX  1.61E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      233.63 TCPU      226.30
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4126466E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      233.66 TCPU      226.32
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0265120  19.9012438  19.8930351   9.0763486   9.0763486   9.0265120
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      233.67 TCPU      226.33
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      236.77 TCPU      229.42
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      238.69 TCPU      231.34
 NUMERICALLY INTEGRATED DENSITY     76.0000056054
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      239.12 TCPU      231.75
 CYC  41 ETOT(AU) -1.999928202101E+03 DETOT -5.94E-07 tst  3.65E-05 PX  1.53E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      239.46 TCPU      231.78
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4126988E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      239.48 TCPU      231.80
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0264792  19.9012675  19.8930423   9.0763659   9.0763659   9.0264792
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      239.50 TCPU      231.81
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      242.63 TCPU      234.92
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      244.51 TCPU      236.79
 NUMERICALLY INTEGRATED DENSITY     76.0000056040
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      244.94 TCPU      237.21
 CYC  42 ETOT(AU) -1.999928202638E+03 DETOT -5.37E-07 tst  3.39E-05 PX  1.45E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      245.01 TCPU      237.24
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4127497E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      245.03 TCPU      237.26
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0264492  19.9012901  19.8930480   9.0763817   9.0763817   9.0264492
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      245.05 TCPU      237.28
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      248.16 TCPU      240.38
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      250.09 TCPU      242.30
 NUMERICALLY INTEGRATED DENSITY     76.0000056028
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      250.51 TCPU      242.71
 CYC  43 ETOT(AU) -1.999928203127E+03 DETOT -4.89E-07 tst  3.15E-05 PX  1.37E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      250.59 TCPU      242.73
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4127964E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      250.61 TCPU      242.75
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0264210  19.9013113  19.8930530   9.0763969   9.0763969   9.0264210
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      250.62 TCPU      242.76
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      253.71 TCPU      245.84
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      255.63 TCPU      247.75
 NUMERICALLY INTEGRATED DENSITY     76.0000056016
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      256.05 TCPU      248.16
 CYC  44 ETOT(AU) -1.999928203571E+03 DETOT -4.43E-07 tst  2.92E-05 PX  1.31E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      256.09 TCPU      248.18
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4128403E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      256.11 TCPU      248.20
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0263947  19.9013315  19.8930571   9.0764110   9.0764110   9.0263947
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      256.12 TCPU      248.22
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      259.24 TCPU      251.32
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      261.14 TCPU      253.22
 NUMERICALLY INTEGRATED DENSITY     76.0000056006
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      261.56 TCPU      253.63
 CYC  45 ETOT(AU) -1.999928203976E+03 DETOT -4.05E-07 tst  2.72E-05 PX  1.24E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      261.68 TCPU      253.70
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4128802E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      261.71 TCPU      253.72
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0263694  19.9013507  19.8930623   9.0764241   9.0764241   9.0263694
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      261.72 TCPU      253.73
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      264.87 TCPU      256.84
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      266.74 TCPU      258.71
 NUMERICALLY INTEGRATED DENSITY     76.0000055995
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      267.17 TCPU      259.12
 CYC  46 ETOT(AU) -1.999928204347E+03 DETOT -3.71E-07 tst  2.53E-05 PX  1.18E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      267.26 TCPU      259.18
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4129192E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      267.28 TCPU      259.20
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0263460  19.9013688  19.8930658   9.0764367   9.0764367   9.0263460
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      267.29 TCPU      259.21
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      270.39 TCPU      262.30
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      272.30 TCPU      264.20
 NUMERICALLY INTEGRATED DENSITY     76.0000055986
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      272.73 TCPU      264.62
 CYC  47 ETOT(AU) -1.999928204685E+03 DETOT -3.38E-07 tst  2.35E-05 PX  1.12E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      273.07 TCPU      264.64
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4129540E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      273.10 TCPU      264.67
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0263235  19.9013860  19.8930698   9.0764486   9.0764486   9.0263235
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      273.11 TCPU      264.68
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      276.45 TCPU      267.75
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      278.38 TCPU      269.68
 NUMERICALLY INTEGRATED DENSITY     76.0000055977
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      278.80 TCPU      270.09
 CYC  48 ETOT(AU) -1.999928204993E+03 DETOT -3.08E-07 tst  2.18E-05 PX  1.07E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      278.87 TCPU      270.12
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4129874E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      278.89 TCPU      270.14
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0263028  19.9014027  19.8930728   9.0764595   9.0764595   9.0263028
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      278.91 TCPU      270.16
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      281.99 TCPU      273.23
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      283.97 TCPU      275.21
 NUMERICALLY INTEGRATED DENSITY     76.0000055969
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      284.39 TCPU      275.61
 CYC  49 ETOT(AU) -1.999928205276E+03 DETOT -2.83E-07 tst  2.03E-05 PX  1.02E-04
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK        TELAPSE      284.44 TCPU      275.63
 INSULATING STATE - TOP OF VALENCE BANDS (A.U.) -7.4130187E-01
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG        TELAPSE      284.46 TCPU      275.65
 CHARGE NORMALIZATION FACTOR   1.00000000
 TOTAL ATOMIC CHARGES:
   9.0262828  19.9014183  19.8930762   9.0764700   9.0764700   9.0262828
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MOQGAD      TELAPSE      284.47 TCPU      275.67
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELLX      TELAPSE      287.56 TCPU      278.74
 +++ ENERGIES IN A.U. +++
 ::: EXT EL-POLE                                  -1.6287525950891E+03
 ::: EXT EL-SPHEROPOLE                             0.0000000000000E+00
 ::: BIELET ZONE E-E                               1.7391762888654E+03
 ::: TOTAL E-E                                     1.1042369377639E+02
 ::: TOTAL E-N + N-E                              -3.2891090217371E+03
 ::: TOTAL N-N                                    -7.2125837783154E+02
 ::: KINETIC ENERGY                                1.9953134357279E+03
 ::: PSEUDO TOTAL   ENERGY                        -1.9046302700643E+03
 ::: VIRIAL COEFFICIENT                            1.0232524293951E+00
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MONMO3      TELAPSE      289.49 TCPU      280.68
 NUMERICALLY INTEGRATED DENSITY     76.0000055961
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT      TELAPSE      289.92 TCPU      281.09
 CYC  50 ETOT(AU) -1.999928205535E+03 DETOT -2.59E-07 tst  1.89E-05 PX  1.02E-04

 == SCF ENDED - TOO MANY CYCLES            E(AU) -1.9999282055351E+03 CYCLES  50


 ENERGY EXPRESSION=HARTREE+FOCK EXCH*0.20000+(BECKE EXCH)*0.80000+LYP   CORR


 TOTAL ENERGY(DFT)(AU)( 50) -1.9999282055351E+03 DE-2.6E-07 tester 1.9E-05
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT EDFT        TELAPSE      289.92 TCPU      281.09
 INFORMATION **** OPTGEOM **** WF DATA SAVED IN fort.9
 ERROR **** OPTGEOM **** SCF FAILED - OPTGEOM STOP
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ERR         TELAPSE      289.94 TCPU      281.09
[3] MPI Abort by user Aborting program !
[3] Aborting program!
[4] MPI Abort by user Aborting program !
[4] Aborting program!
[7] MPI Abort by user Aborting program !
[7] Aborting program!
[1] MPI Abort by user Aborting program !
[1] Aborting program!
[2] MPI Abort by user Aborting program !
[2] Aborting program!
[6] MPI Abort by user Aborting program !
[6] Aborting program!
[5] MPI Abort by user Aborting program !
[5] Aborting program!
[0] MPI Abort by user Aborting program !
[0] Aborting program!
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Killed by signal 15.
Job  /apps/lsf/6.2/linux2.6-glibc2.3-x86_64/bin/gmmpirun_wrapper /work/scratch/scarf044/pub/CRYSTAL06/20070510_OPTbugfix.ANHARM/Linux-mpichpgf_Dynamic/pgi_6.2_apps_mpi_gm_pgi_1.2.7..15/Pcrystal

TID   HOST_NAME   COMMAND_LINE            STATUS            TERMINATION_TIME
===== ========== ================  =======================  ===================
00002 cn134.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00000 cn134.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00001 cn134.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00003 cn134.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00006 cn141.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00005 cn141.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00007 cn141.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
00004 cn141.scar /work/scratch/sc  Signaled (SIGIOT)        10/23/2007 01:52:16
Tue Oct 23 01:52:20 BST 2007
total 14160
-rw-------  1 scarf044 cseg 2109440 Oct 23 01:52 core.29489
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 dffit3.dat
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe1
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe2
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe3
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe4
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe5
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe6
-rw-r--r--  1 scarf044 cseg     153 Oct 23 01:52 ERROR.pe7
-rw-r--r--  1 scarf044 cseg  290288 Oct 23 01:52 fort.10.pe0
-rw-r--r--  1 scarf044 cseg  278288 Oct 23 01:52 fort.10.pe1
-rw-r--r--  1 scarf044 cseg  290288 Oct 23 01:52 fort.10.pe2
-rw-r--r--  1 scarf044 cseg  178144 Oct 23 01:52 fort.10.pe3
-rw-r--r--  1 scarf044 cseg  283184 Oct 23 01:52 fort.10.pe4
-rw-r--r--  1 scarf044 cseg  295184 Oct 23 01:52 fort.10.pe5
-rw-r--r--  1 scarf044 cseg  283184 Oct 23 01:52 fort.10.pe6
-rw-r--r--  1 scarf044 cseg  316864 Oct 23 01:52 fort.10.pe7
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe0
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe1
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe2
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe3
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe4
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe5
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe6
-rw-r--r--  1 scarf044 cseg  104080 Oct 23 01:52 fort.11.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.12
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.17.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe0
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.18.pe7
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe0
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe1
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe2
-rw-r--r--  1 scarf044 cseg    9216 Oct 23 01:48 fort.19.pe3
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe4
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe5
-rw-r--r--  1 scarf044 cseg    9888 Oct 23 01:48 fort.19.pe6
-rw-r--r--  1 scarf044 cseg   12288 Oct 23 01:48 fort.19.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe0
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.1.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.20
-rw-r--r--  1 scarf044 cseg     430 Oct 23 01:47 fort.33
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe0
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.38.pe7
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe0
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe1
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe2
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe3
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe4
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe5
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe6
-rw-r--r--  1 scarf044 cseg  769384 Oct 23 01:48 fort.3.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe0
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.40.pe7
-rw-r--r--  1 scarf044 cseg  290288 Oct 23 01:52 fort.8.pe0
-rw-r--r--  1 scarf044 cseg  278288 Oct 23 01:52 fort.8.pe1
-rw-r--r--  1 scarf044 cseg  290288 Oct 23 01:52 fort.8.pe2
-rw-r--r--  1 scarf044 cseg  178144 Oct 23 01:52 fort.8.pe3
-rw-r--r--  1 scarf044 cseg  283184 Oct 23 01:52 fort.8.pe4
-rw-r--r--  1 scarf044 cseg  295184 Oct 23 01:52 fort.8.pe5
-rw-r--r--  1 scarf044 cseg  283184 Oct 23 01:52 fort.8.pe6
-rw-r--r--  1 scarf044 cseg  316864 Oct 23 01:52 fort.8.pe7
-rw-r--r--  1 scarf044 cseg  523460 Oct 23 01:52 fort.9
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe0
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe1
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe2
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe3
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe4
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe5
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe6
-rw-r--r--  1 scarf044 cseg  112096 Oct 23 01:52 fort.95.pe7
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe0
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe1
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe2
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe3
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe4
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe5
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe6
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 fort.96.pe7
-rw-r--r--  1 scarf044 cseg    1822 Oct 23 01:47 INPUT
-rw-r--r--  1 scarf044 cseg       0 Oct 23 01:47 SCFOUT.LOG
wave function binary file /home/scarf044/giuseppe/TiO2/slab_110/rutile/B3LYP/Ti_86411d41_O_8411d1_TI777714_IS888/3L_optatom2.f9
Tue Oct 23 01:52:20 BST 2007
file fort.33 saved as /home/scarf044/giuseppe/TiO2/slab_110/rutile/B3LYP/Ti_86411d41_O_8411d1_TI777714_IS888/3L_optatom2.xyz
/work/scratch/scarf044/3L_optatom2_29478 removed
"""

def test():
  """ Makes sure that structure and input structure are read the same way. 

      Fractional and cartesian coordinates do not point necessarily to the same
      place in real space. More specifically, crystal will randomly assign -0.5
      or 0.5 to fractional coordinates. There does not seem to be a convention.
  """
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import join
  from pylada.dftcrystal import Extract

  directory = mkdtemp()
  try:
    with open(join(directory, 'OUTCAR'), 'w') as file: file.write(string)

    result = Extract(join(directory, 'OUTCAR'))
    for i in xrange(6):
      assert all(abs(result.structure[i].pos - result.input_structure[i].pos) < 1e-8)
  finally:
    try: rmtree(directory)
    except: pass

if __name__ == '__main__':
  test()
