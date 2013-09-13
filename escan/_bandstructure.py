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

""" Bandstructure plotting tools """
__docformat__  = 'restructuredtext en'
__all__ = ['BPoints', 'ReducedBPoints', 'InnerBPoints', 'ReducedInnerBPoints', 'plot_bands', 'plot_alloybands']

from .kpoints import KPoints, _reduced_grids_factory
from pylada import try_import_matplotlib

class BPoints(KPoints):
  """ *Grid* of kpoints for bands-structure calculation. """
  def __init__(self, lines=None, density=10, relax=True, mink=3):
    KPoints.__init__(self)
    self.lines = lines
    """ Lines of k-points for which to perform calculations. """
    if self.lines is None: self.lines = []
    self.density = density
    """ Requested density of kpoints. """
    self.mink = mink
    """ Minimum number of kpoints per line. """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  def _mnk(self, input, output):
    """ Yields lines of k-points with appropriate density. """
    from numpy import any, abs, dot
    from numpy.linalg import inv

    if len(self.lines) == 0: return
    if self.relax: deformation = dot(inv(output.cell.T), input.cell.T)
    
    last_end = self.lines[0][0] + 1e0
    for start, end in self.lines:
      if any(abs(last_end-start) > 1e-12):
        yield 1e0, (dot(deformation, start) if self.relax else start)
      last_end = end
      for kpoint in self._line(start, end):
        yield 1e0, (dot(deformation, kpoint) if self.relax else kpoint)

  def _line(self, start, end):
    """ Yields kpoints on a line. """
    from numpy.linalg import norm
    distance = norm(end-start)
    nbkpt = max(self.mink, int(float(self.density) * distance - 1))
    stepsize = 1e0 / float(nbkpt)
    for i in range(1, nbkpt+1):
      yield start + float(i) * stepsize * (end-start)

  def __iadd__(self, values):
    """ Adds segment to band-structure line. """
    from numpy import array
    assert hasattr(values, '__iter__'), ValueError('Value should be a sequence.')
    values = [v for v in values]
    if len(values) == 2 and len(values[0]) == 3 and len(values[1]) == 3: values = [values]
    for start, end in values:
      assert hasattr(start, '__iter__'), ValueError('Values should be two tuples of 3d-vectors.')
      start = array(start, dtype='float64')
      assert hasattr(end, '__iter__'), ValueError('Values should be two tuples of 3d-vectors.')
      end = array(end, dtype='float64')
      self.lines.append((start, end))
    return self

  def __add__(self, other):
    """ Returns new instance with summed segment. """
    result = self.__class__(self.lines, self.density, self.relax, self.mink)
    result += other.lines if hasattr(other, 'lines') else other
    return result

  def __repr__(self):
    """ Returns representation of this object. """
    string = '{0.__class__.__name__}(None, {0.density}, {0.relax}, {0.mink})'.format(self)
    for start, end in self.lines:
      string += '+([{0[0]},{0[1]},{0[2]}], [{1[0]},{1[1]},{1[2]}])'.format(start, end)
    return string

class InnerBPoints(BPoints):
  """ Cuts off band-structure directions to fit inside first Brillouin Zone. """
  def __init__(self, *args, **kwargs):
    """ See BPoints for parameters. """
    self.tolerance = kwargs.pop("tolerance", 1e-12)
    """ Convergence criteria for determining inside/outside BZ criteria. """
    super(InnerBPoints, self).__init__(*args, **kwargs)

  def _mnk(self, input, output):
    """ Yields lines of k-points with appropriate density. """
    from numpy import any, abs, dot, array
    from numpy.linalg import inv

    if len(self.lines) == 0: return
    if self.relax: deformation = dot(inv(output.cell.T), input.cell.T)
    kcell = inv(input.cell).T
    
    last_end = self.lines[0][0] + 10e0 
    for start, end in self.lines:
      start, end = array(start, dtype="float64"), array(end, dtype="float64")
      start, end = self._to_first_bz(start, end, kcell, self.tolerance)
      if any(abs(last_end-start) > 1e-12):
        yield 1e0, (dot(deformation, start) if self.relax else start)
      last_end = end
      for kpoint in self._line(start, end):
        yield 1e0, (dot(deformation, kpoint) if self.relax else kpoint)

  def _halfline(self, start, direction, kcell, tolerance):
    """ Finds intersection with BZ if start is in BZ. """
    from numpy.linalg import norm
    from ..crystal import to_voronoi
    is_inside = lambda x: norm(x - to_voronoi(x, kcell)) < tolerance
    if not is_inside(start):
      raise ValueError("Starting point is not inside first BZ. {0}, {1}.".format(start, to_voronoi(start, kcell)))

    direction = direction / norm(direction)
    xin, xout = 0, sum([norm(k) for k in kcell.T])
    result = start + direction * xout
    last_inside = is_inside(result)
    if last_inside: raise RuntimeError("Could not find outside point.")

    while abs(xout - xin) > tolerance:
      if is_inside(start + direction * 0.5 * (xout + xin)): xin  = 0.5 * (xout + xin)
      else:                                                 xout = 0.5 * (xout + xin)
    return start + direction * xin

  def _to_first_bz(self, start, last, kcell, tolerance=1e-12):
    """ Computes endpoints of directions inside first bz. """
    from numpy.linalg import norm
    from ..crystal import to_voronoi
    is_inside = lambda x: norm(x - to_voronoi(x, kcell)) < tolerance

    startin, lastin = is_inside(start), is_inside(start)
    if not (startin or lastin): # neither point is inside the first BZ
      # translate start into the cell, then figure out other point, keeping direction constant.
      direction = last - start
      start = to_voronoi(start, kcell)
      last = self._halfline(start, direction, kcell, tolerance)
      # then figure out first point as well.
      start = self._halfline(last, -direction, kcell, tolerance)
    elif startin and not lastin: # one point is in first BZ, figure out the other.
      last = self._halfline(start, last-start, kcell, tolerance)
    elif (not startin) and lastin: # one point is in first BZ, figure out the other.
      start = self._halfline(last, start-last, kcell, tolerance)
    return start, last

ReducedBPoints = _reduced_grids_factory('ReducedBPoints', BPoints)
ReducedInnerBPoints = _reduced_grids_factory('ReducedInnerBPoints', InnerBPoints)

if try_import_matplotlib:
  try: import matplotlib.pyplot as plt 
  except: 
    def plot_bands(extractor, **kwargs):
      """ Plots band-structure. """
      raise ImportError("Cannot use plot_bands without matplotlib. """)
    def plot_alloybands(extractor, multicell, tolerance=1e-6, **kwargs):
      """ Plots alloy band-structure using the majority representation. """
      raise ImportError("Cannot use plot_alloybands without matplotlib. """)
  else:
    def plot_bands(extractor, offset = 0.05, labels = None, **kwargs):
      """ Tries and plots band-structure. """
      from numpy import array, min, max
      from numpy.linalg import norm
      from matplotlib import rcParams

      old = extractor.unreduce
      extractor.unreduce = True

      bandcolor = kwargs.pop('bandcolor', 'blue')
      edgecolor = kwargs.pop('edgecolor', 'red')
      edgestyle = kwargs.pop('edgestyle', '-')
      linewidth = kwargs.pop('linewidth', rcParams['lines.linewidth'])
      bandwidth = kwargs.pop('bandwidth', linewidth)

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

      # then line markers.
      offset *= x[-1]
      lines = []
      for first, last, add_offset in zip(bk[:-1], bk[1:], offsets):
        if first != 0: lines.append(x[first])
        if add_offset:
          x[first:] += offset
          lines.append(x[first])
      for first, last in  zip(lims[:-1], lims[1:]):
        plt.plot(x[first:last], y[first:last], color=bandcolor, linewidth=bandwidth,**kwargs)

      # then plot vbm and cbm.
      kwargs['linewidth'] = rcParams['axes.linewidth']
      for l in lines: plt.axvline(l, color='black', **kwargs)
      kwargs['linestyle'] = edgestyle
      plt.axhline(extractor.vbm, color=edgecolor, **kwargs)
      plt.axhline(extractor.cbm, color=edgecolor, **kwargs)



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
        for i in axes.xaxis.get_major_ticks():
          i.tick1On=False
          i.tick2On=False

      extractor.unreduce = old

    def plot_alloybands(extractor, multicell, tolerance=1e-6, **kwargs):
      """ Plots alloy band-structure using the majority representation. """
      from numpy import array, min, max
      from numpy.linalg import norm
      from . import majority_representation

      edgecolor = kwargs.pop('edgecolor', 'red')
      edgestyle = kwargs.pop('edgestyle', '-')

      # first finds breaking point.
      istr = extractor.input_structure
      ostr = extractor.structure
      kpoints = array([u[1] for u in extractor.functional.kpoints.unreduced(istr, ostr)])
      delta = kpoints[1:] - kpoints[:-1]
      norms = [norm(delta[i,:]) for i in range(delta.shape[0])]
      bk = []
      for i, d in enumerate(norms[1:]):
        if abs(norms[i]-d) > 1e-6: bk.append(i+1)

      # then plot bands.
      xvalues = array([sum(norms[:i]) for i in range(len(norms)+1)])
      args = [[],[],[]]
      maj = majority_representation(extractor, multicell)
      for values, x in zip(maj, xvalues):
        for v in values: 
          args[0].append(x)
          args[1].append(v[0])
          args[2].append(v[1])

      # then line markers.
      plt.hexbin(*args, **kwargs)
      x = array([sum(norms[:i]) for i in range(len(norms)+1)])
      for i in bk: plt.axvline(x[i], color='black', **kwargs)

      kwargs.pop('linestyle', None) 
      plt.axhline(extractor.vbm, color=edgecolor, linestyle=edgestyle, **kwargs)
      plt.axhline(extractor.cbm, color=edgecolor, linestyle=edgestyle, **kwargs)

      y = array(extractor.eigenvalues)
      plt.xlim((x[0], x[-1]))
      ylims = min(y) - (max(y) - min(y))*0.05, max(y) + (max(y) - min(y))*0.05
      plt.ylim(ylims)
else:
  def plot_bands(extractor, **kwargs):
    """ Plots band-structure. """
    raise ImportError("Cannot use plot_bands without matplotlib. """)
  def plot_alloybands(extractor, multicell, tolerance=1e-6, **kwargs):
    """ Plots alloy band-structure using the majority representation. """
    raise ImportError("Cannot use plot_alloybands without matplotlib. """)
