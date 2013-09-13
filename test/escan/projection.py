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

import pickle
# import matplotlib.pyplot as plt 

class Projections(object):
  def __init__(self): object.__init__(self)
  def __call__(self, outdir, comm=None, **kwargs): 
    from os.path import join
    from pylada.escan import KExtract
    extract = KExtract(outdir, comm=comm)
    for value in extract:
      filename = join(value.directory, "PROJECT_BS")
      print "where: ", value.directory, filename
      a = self._project(value, filename, **kwargs)

  def _project(self, extract, filename = "projections", alpha = 1e0, **kwargs):
    """ Computes projected densities around each atom for a calculation. """
    from os.path import exists
    from sys import exit, getrefcount
    if exists(filename): 
      from pickle import load
      with open(filename, "r") as file: return load(file)

    from pickle import dump
    from numpy import zeros, dot, array, exp, sum
    from numpy.linalg import norm
    from quantities import angstrom
    from pylada import periodic_table as table
    from pylada.crystal import gaussian_projector
    species = set([u.type for u in extract.structure.atoms])
    species = sorted(list(species))
    result = {}
    for key in species:
      if (key not in kwargs) and (key not in table.__dict__): continue
      radius = kwargs.get(key, getattr(table, key)).atomic_radius
      sigma = -alpha / radius / radius

      # create projection operator.
      proj = zeros(extract.rvectors.shape[:-1])
      cell = extract.structure.cell * extract.structure.scale 
      for atom in extract.structure.atoms:
        if atom.type != key: continue
        pos = atom.pos * extract.structure.scale 
        proj += gaussian_projector(extract.rvectors, cell * angstrom, pos * angstrom, sigma )
      result[key] = [(w.eigenvalue, w.expectation_value(proj)) for w in extract.rwfns]
    
    n = norm(array(result.values()))
    for key in result.keys(): result[key] /= n
    with open(filename, "w") as file: dump(result, file)
    return result

  def iter(self, outdir, **kwargs):
    from os.path import join
    from pylada.escan import KExtract
    extract = KExtract(outdir)
    for value in extract:
      filename = join(value.directory, "PROJECT_BS")
      yield self._project(value, filename, **kwargs)

def compute_bs():
  from numpy import array
  from pylada.escan import read_input, exec_input, ReducedBPoints
  from pylada.vff import Vff

  # reads input file.
  input = read_input("input.py")

  # creating unrelaxed structure.
  structure = input.vff.lattice.to_structure()
  structure.atoms[0].type = "Si"
  structure.atoms[1].type = "Ge"
  structure.scale = 5.65

  # some kpoints + associated name
  X = array( [1,0,0], dtype="float64" )
  G = array( [0,0,0], dtype="float64" )
  L = array( [0.5,0.5,0.5], dtype="float64" )
  W = array( [0, 0.5,1], dtype="float64" )

  # Each job is performed for a given kpoint (first argument), at a given
  # reference energy (third argument). Results are stored in a specific directory
  # (second arguement). The expected eigenvalues are given in the fourth argument.
  input = read_input('input.py')
  kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

  kescan.fft_mesh = 14, 14, 14
  kescan.kpoints = ReducedBPoints(density=20) + (X, G) + (G, L)
  result = kescan( structure, outdir='results/projections', 
                   nbstates = len(structure.atoms) * 4 + 4,
                   eref = None )

def plot_bands(extractor, offset = 0.05, labels = None, **kwargs):
  """ Tries and plots band-structure. """
  from numpy import dot, array, min, max, sqrt
  from numpy.linalg import norm
  from matplotlib import pyplot as plt, rcParams

  old = extractor.unreduce
  extractor.unreduce = True

  bandcolor = kwargs.pop('bandcolor', 'blue')
  edgecolor = kwargs.pop('edgecolor', 'red')
  edgestyle = kwargs.pop('edgestyle', '--')
  awidth    = kwargs.pop('awidth', 5.)
  type      = kwargs.pop('type', 'Si')
  alpha     = kwargs.pop('alpha', 1.)
  facecolor = kwargs.pop('facecolor', bandcolor)
  linewidth = kwargs.pop('linewdith', rcParams['lines.linewidth'])
  bandwidth = kwargs.pop('bandwidth', linewidth)
  c = float(len([0 for atom in extractor.structure.atoms if atom.type==type]))
  c /= float(len(extractor.structure.atoms))

  # first finds breaking point.
  kpoints = extractor.kpoints
  delta = kpoints[1:] - kpoints[:-1]
  norms = [norm(delta[i,:]) for i in range(delta.shape[0])]
  bk, offsets, x, lims = [0], [False], [], [0]
  for i, d in enumerate(norms[1:]):
    if abs(norms[i]-d) > 1e-6: 
      if i == bk[-1] and i != 0:
        offsets[-1] = True
        bk[-1] += 1
        lims.append(bk[-1])
      else:
        bk.append(i+1)
        offsets.append(False)
      x.append(0 if offsets[-1] else norms[i])
    else: x.append(d)
  bk.append(None)
  offsets.append(False)
  lims.append(None)
  x.append(x[-1])

  # then plot bands.
  x = array([sum(x[:i]) for i in range(len(x)+1)])
  y = array(extractor.eigenvalues)
  projs = Projections()
  si_projs = array([ u["Si"] for u in projs.iter(extractor.directory) ])
  si_projs = si_projs[:,:,1]
  ge_projs = array([ u["Ge"] for u in projs.iter(extractor.directory) ])
  ge_projs = ge_projs[:,:,1]

  # then line markers.
  offset *= x[-1]
  lines = []
  for first, last, add_offset in zip(bk[:-1], bk[1:], offsets):
    if first != 0: lines.append(x[first])
    if add_offset:
      x[first:] += offset
      lines.append(x[first])
  fermi = 0.5 * (extractor.vbm + extractor.cbm)
  pldict = dict(kwargs)
  pldict['color'] = bandcolor
  pldict['facecolor'] = facecolor
  pldict['linewidth'] = bandwidth
  for first, last in  zip(lims[:-1], lims[1:]):
    for band, si, ge in zip(y[first:last].T, si_projs[first:last].T.real, ge_projs[first:last].T.real):
      width = (si if type == "Si" else ge)  / sqrt(si * si + ge * ge) / c
      if float(band[0]) < float(fermi):
        plt.fill_between(x[first:last], band, band - awidth * width, **pldict)
  for first, last in  zip(lims[:-1], lims[1:]):
    for band, si, ge in zip(y[first:last].T[::-1], si_projs[first:last].T.real[::-1], ge_projs[first:last].T.real[::-1]):
      width = (si if type == "Si" else ge)  / sqrt(si * si + ge * ge) / c
      if float(band[0]) > float(fermi):
        plt.fill_between(x[first:last], band + awidth * width, band, **pldict)
#   plt.plot(x[first:last], y[first:last], color=bandcolor, **kwargs)
  for l in lines: plt.axvline(l, color='black', linewidth=rcParams['axes.linewidth'], **kwargs)

  # then plot vbm and cbm.
  kwargs['linestyle'] = linestyle
  kwargs['color'] = edgecolor
  kwargs['linewidth'] = rcParams['axes.linewidth']
  plt.axhline(extractor.vbm, **kwargs)
  plt.axhline(extractor.cbm, **kwargs)



  plt.xlim((x[0], x[-1]))
  ylims = min(y) - (max(y) - min(y))*0.05, max(y) + (max(y) - min(y))*0.05
  plt.ylim(ylims)
  axes = plt.gca()
  axes.yaxis.set_ticks_position('both')
  if labels is None: axes.xaxis.set_ticks([])
  else:
    lines.insert(0, 0)
    lines.append(x[-1])
    assert len(labels) <= len(lines),\
           ValueError("Could not find as many breaking points as were provided labels.")
    axes.xaxis.set_ticks(lines[:len(labels)])
    axes.xaxis.set_ticklabels(labels)

  extractor.unreduce = old

# compute_bs()
# p = Projections()
# p(".")
# plot_bands(extract_bs, awidth=1)
