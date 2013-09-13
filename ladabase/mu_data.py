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

""" Elemental mus (FERE). """
__docstring__ = "restructuredtext end"
from quantities import eV

elemental_mus = { 'Ag': -0.82700958541595615*eV,
                  'Au': -2.2303410086960551*eV,
                  'Ba': -1.3944992462870172*eV,
                  'Be': -3.3972092621754264*eV,
                  'Bi': -4.3853003286558812*eV,
                  'Cl': -1.6262437135301639*eV,
                  'Co': -4.7543486260270402*eV,
                  'Cr': -7.2224146752384204*eV,
                  'Cu': -1.9725806522979044*eV,
                  'F' : -1.7037867766570287*eV,
                  'Fe': -6.1521343161090325*eV,
                  'Ge': -4.137439286830797*eV,
                  'Hf': -7.397695761161847*eV,
                  'Hg': -0.12361566177444684*eV,
                  'Ir': -5.964577394407752*eV,
                  'K' : -0.80499202755075006*eV,
                  'La': -3.6642174822805287*eV,
                  'Li': -1.6529591953887741*eV,
                  'Mn': -6.9965778258511993*eV,
                  'Na': -1.0640326227725869*eV,
                  'Nb': -6.6867516375690608*eV,
                  'Ni': -3.5687859474688026*eV,
                  'Pd': -3.1174044624888873*eV,
                  'Pt': -3.9527597082085424*eV,
                  'Rb': -0.6750560483522855*eV,
                  'Rh': -4.7622899695820369*eV,
                  'Sb': -4.2862260747305099*eV,
                  'Sc': -4.6302422200922519*eV,
                  'Si': -4.9927748122726356*eV,
                  'Sn': -3.7894939351245469*eV,
                  'Sr': -1.1674559193419329*eV,
                  'Ta': -8.8184831379805324*eV,
                  'Te': -3.2503408197224912*eV,
                  'Ti': -5.5167842601434147*eV,
                  'V' : -6.4219725884764864*eV,
                  'Y' : -4.812621315561298*eV,
                  'Zr': -5.8747056261113126*eV,
                  'O' :-4.76*eV,
                  'S' :-4.00*eV,
                  'Se':-3.55*eV,
                  'N' :-8.51*eV,
                  'P' :-5.64*eV,
                  'As':-5.06*eV,
                  'Mg':-0.99*eV,
                  'Ca':-1.64*eV,
                  'Zn':-0.84*eV,
                  'Cd':-0.56*eV,
                  'Ga':-2.37*eV,
                  'Al':-3.02*eV,
                  'In':-2.31*eV}
""" mus determined from original and FERE paper. """


original_paper = ['O', 'S', 'Se', 'N', 'P', 'As', 'Mg', 'Ca', 'Zn', 'Cd', 'Ga', 'Al', 'In']
""" Atoms determined from original Stephan Lany paper. """

def enthalpy(extract):
  """ Returns enthalpy if compuational parameters are ok, else None. 
  
      A FERE calculation is performed in DFT with PBE functional and PAW
      pseudo-potentials.  It is static (no ionic or cell-shape relaxation). A
      given set of transition metals must have U = 3eV, and another two U = 5eV
      (see code below).  

      :returns: None if this is not a FERE calculation, and the enthalpy in eV
                per atom if it is.
  """
  if not extract.is_dft: return None
  if extract.relaxation != "static": return None 
  if extract.pseudopotential != 'PAW_PBE': return None
  if len(extract.species) == 1: return None
  for specie in extract.species:
    if specie not in elemental_mus: return None
    U = None
    if specie in ["Cu", "Ag"]: U = 5 # in eV
    elif specie in [ "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",\
                     "Y", "Zr", "Nb", "Rh", "Pd", \
                     "La", "Hf", "Ta", "Ir", "Pt", "Au" ]: U = 3 # eV
    if U == None and specie in extract.HubbardU_NLEP: return None
    elif U != None:
      if specie not in extract.HubbardU_NLEP: return None
      if len(extract.HubbardU_NLEP[specie]) != 1: return None
      params = extract.HubbardU_NLEP[specie][0]
      if params["type"] != 2: return None # "dudarev"
      if params["l"] != 2: return None
      if params["func"] != "U": return None
      if abs(params["J"]) > 1e-3: return None
      if abs(params["U"]-U) > 1e-3: return None
  units = extract.energy.units
  result = -sum([ n * elemental_mus[s].rescale(units).magnitude.tolist() \
                  for n, s in zip(extract.stoichiometry, extract.species) ])
  result += extract.energy.magnitude.tolist()
  return result / float(sum(extract.stoichiometry))

def delta_muO(P,T):
  """ Evolution in pressure and temperature of mu oxygen.
  
      Scavenged from Stephan Lany.

      :Parameters:
        P
          Pressure in atomospheres.
        T 
          Temperature in Kelvin.
  """
  from numpy import log

  t0=298.15

  h0o=8.68e0*1.0364e-2
  s0o=205.04e0*1.0364e-5
  
  k=8.6174e-5
  cp2=3.5e0*k
  
  dp2=0.5e0*k*T*log(P)
  dmo=0.5e0*(h0o+cp2*(T-t0)-T*(s0o+cp2*log(T/t0)))

  return dmo+dp2

def partial_pressureO(delta_mu, T):
  """ Partial pressure of oxygen with respect to chemical potential and temperature. """
  from numpy import exp, log

  t0=298.15

  h0o=8.68e0*1.0364e-2
  s0o=205.04e0*1.0364e-5
  
  k=8.6174e-5
  cp2=3.5e0*k
  
  dmo=0.5e0*(h0o+cp2*(T-t0)-T*(s0o+cp2*log(T/t0)))
  return exp( -2e0*(dmo - delta_mu) / (k*T) )
