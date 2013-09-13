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

"""" Creates diagram of formation enthalpies with respect to Fermi energy. """
from pylada.jobs import AbstractMassExtract
from pylada.opt.decorators import make_cached

class PointDefectExtactor(object):
  """ Extracts output across all charged states of a defect. """
  def __init__(self, extract, epsilon = 1e0, host = None, pa_maxdiff=-8):
    """ Initializes an enthalpy function. """
    # extraction object.
    self.extract = extract.copy(unix_re=False)
    if self.extract.excludes is None: self.extract.excludes = [".*relax_*"]
    else: self.extract.excludes.append(".*relax_*$")

    self.epsilon = epsilon
    """ Dimensionless dielectric constant. """
    self._host = host
    """ Host extraction object. 

        If None, will determine it from the directory structure.
    """
    self.pa_maxdiff = pa_maxdiff
    """ Potential alignment parameter. """

  @property
  def rootdir(self):
    """ Root directory of defects. """
    return self.extract.rootdir

  def _all_jobs(self):
    """ Loops over all jobs in special way. """
    for child in self.extract.children: # Each child is a different charge state.
      for job in child.itervalues(): yield job

  @property
  @make_cached
  def _charge_correction(self):
    """ Returns the charge corrections.
    
        Tries and minimizes the number of calculations by checking if performed
        in smae cell.
    """
    from numpy.linalg import inv, det
    from numpy import array
    from quantities import eV
    result, cells  = [], []
    # loops of all jobs.
    for job in self._all_jobs():
      cell = job.structure.cell * job.structure.scale
      invcell = inv(cell)
      found = None
      # looks if already exists.
      for i, other in enumerate(cells):
        rotmat = other * cell
        d = abs(det(rotmat))
        if abs(d - 1e0) < 1e-8: continue
        invrotmat = inv(rotmat)
        if all( abs(rotmat.T - invrotmat) < 1e-8 ): found = i; break
      if found is None: 
        cells.append(inv(cell))
        result.append( job.charge_corrections / self.epsilon )
      else: result.append(result[found])
    return array(result) * eV

  @property
  @make_cached
  def _potential_alignment(self):
    """ Potential alignments for all jobs. """
    from numpy import array
    from quantities import eV
    from pylada.crystal.point_defects import potential_alignment
    return array([ potential_alignment(state, self.host, self.pa_maxdiff) \
                   for state in self._all_jobs() ]) * eV

  @property
  @make_cached
  def _band_filling(self):
    """ Band-filling for all jobs. """
    from numpy import array
    from quantities import eV
    from pylada.crystal.point_defects import band_filling
    return array([ band_filling(state, self.host, maxdiff=self.pa_maxdiff) \
                   for state in self._all_jobs() ]) * eV

  @property
  @make_cached
  def _uncorrected(self):
    """ Uncorrected formation enthalpy. """
    from numpy.linalg import det
    from numpy import array
    from quantities import eV
    energies = []
    for state in self._all_jobs():
      n = int(det(state.structure.cell)/det(self.host.structure.cell) + 1.e-3) + 0.
      energies.append(state.total_energy - self.host.total_energy * n)
    return array(energies) * eV 

  @property
  def _corrected(self):
    """ Corrected formation enthalpy. """
    return   self._uncorrected \
           + self._charge_correction\
           + self._potential_alignment\
           + self._band_filling


  @property
  @make_cached
  def _charged_states(self):
    """ Yields extraction routine toward each charge states. 

        If a charge state has children, then only the lowest energy calculation
        is returned.
    """
    from os.path import basename
    from operator import itemgetter
    alles = {}
    names = [child.directory for child in self._all_jobs()]
    for n, u, c, p, b, corr in zip( names, self._uncorrected, \
                                    self._charge_correction, \
                                    self._potential_alignment, \
                                    self._band_filling,\
                                    self._corrected ):
      alles[n] = u, c, p, b, corr
    corrected = self._corrected
    charges  = [u.charge for u in self._all_jobs()]
    children  = [u for u in self._all_jobs()]
    result = []
    for charge in sorted(list(set(charges))):
      sequence = [(child, u) for child, u, c in zip(children, corrected, charges) if c == charge]
      child = sorted(sequence, key=itemgetter(1))[0][0].copy()
      child.__dict__['raw_DeltaH']          = alles[child.directory][0]
      child.__dict__['charge_corrections']  = alles[child.directory][1]
      child.__dict__['potential_alignment'] = alles[child.directory][2]
      child.__dict__['band_filling']        = alles[child.directory][3]
      child.__dict__['DeltaH']              = alles[child.directory][4]
      result.append(child)
    return result

  @property
  def _all_energies(self):
    """ Dictionary with all energies. """
    return result


  @property
  def host(self):
    """ Returns extraction object towards the host. """
    if self._host is None: 
      host = self.extract['../..' if self._is_site is None else '../../..']
      host = self.copy(excludes=[".*PointDefects"], naked_end=False)
      host.excludes.extend(host.excludes)
      lowest = sorted(child.total_energies.iteritems(), key=itemgetter(1))[0][0]
      self._host = [u for u in self.extract[lowest].itervalues()]
      assert len(self._host) == 1
      self._host = self._host[0]
    return self._host

  @property 
  def _site(self):
    """ Returns site number or None. """
    from re import match
    regex = match(r"site_(\d+)", self.extract.view.split('/')[-1])
    return int(regex.group(1)) if regex is not None else None

  @property 
  def name(self):
    """ Name of the defect. """
    return self.extract.view.split('/')[-2 if self._site is not None else -1]

  @property
  def is_vacancy(self):
    """ True if this is a vacancy. """
    from re import match
    return match("vacancy_[A-Z][a-z]?", self.name) is not None

  @property
  def is_interstitial(self):
    """ True if this is an interstitial. """
    from re import match
    return match("[A-Z][a-z]?_interstitial_\S+", self.name) is not None

  @property
  def is_substitution(self):
    """ True if this is a substitution. """
    from re import match
    return match("[A-Z][a-z]?_on_[A-Z][a-z]?", self.name) is not None

  @property
  def n(self):
    """ Number of atoms added/removed from system.
    
        This is a dictionary.
    """
    from re import match
    if self.is_vacancy:
      return {match("vacancy_([A-Z][a-z])?", self.name).group(1): -1}
    elif self.is_interstitial:
      return {match("[A-Z][a-z]?_interstitial_(\S+)", self.name).group(1): -1}
    else: 
      found = match("([A-Z][a-z])?_on_([A-Z][a-z])?", self.name)
      return {match.group(1): 1, match.group(2): -1}

  def uncache(self):
    """ Uncaches result. """
    from opt import uncache as opt_uncache
    opt_uncache(self)
    self.extract.uncache()
    self.host.unchache()


class PointDefectExtractor(PointDefectExtractorImpl):
  """ Properties of a single defect. """
  def __init__(self, extract, epsilon = 1e0, host = None, pa_maxdiff=-8):
    """ Initializes an enthalpy function. """
    super(PointDefectExtractor, self).__init__(extract, epsilon, host, pa_maxdiff)

  def chempot(self, mu):
    """ Computes sum of chemical potential from dictionary ``mu``. 
    
        :Param mu: Dictionary of chemical potentials. If no units, assumes eV.
        :return: Chemical potential of this defect. Value is always in eV.
    """
    from quantities import eV
    if mu is None: return 0 * eV 
    result = 0e0 * eV
    n = self.n
    for specie, value in self.n:
      assert specie in mu,\
             ValueError("Specie {0} not in input chemical potentials {1}.".format(specie, mu))
      chem = mu[specie]
      if not hasattr(chem, 'units'): chem = chem * eV
      result += value * chem
    return result.rescale(eV)

  def _lines(self):
    """ Returns lines composed by the different charge states. """
    from numpy import array
    from quantities import elementary_charge as e, eV
    lines = []
    states = set()
    for state in self._charged_states:
      assert state.charge not in states,\
             RuntimeError("Found more than one calculation for the same charge state.")
      states.add(state.charge)
      lines.append((state.DeltaH.rescale(eV), state.charge))
    return lines

  def _all_intersections(self, _lines):
    """ Returns all intersection points between vbm and cbm, ordered. """
    from numpy import array
    from quantities import eV
    vbm = 0.*eV
    cbm = (self.host.cbm - self.host.vbm).rescale(eV)
    result = []
    for i, (b0, a0) in enumerate(_lines[:-1]):
      for b1, a1 in _lines[i+1:]: result.append( (b0 - b1) / (a1 - a0) )
    result = [u for u in sorted(result) if u - 1e-6 * eV > vbm]
    result = [u for u in result if u + 1e-6 * eV < cbm]
    result.append(cbm)
    result.insert(0, vbm)
    return array([array(u.rescale(eV)) for u in result]) * eV

  def lines(self):
    """ Lines forming the formation enthalpy diagram. 
    
        :return: A list of 2-tuples with the first item b and the second a (a*x+b).
    """
    from numpy import array
    from quantities import eV
    _lines = self._lines()
    intersections = self._all_intersections(_lines)
    # adds line before vbm
    func  = lambda x: x[0] + (intersections[0]-eV)*x[1] 
    lines = [ min(_lines, key=func)  ]

    # now look for lines up to cbm
    for i, intersection in enumerate(intersections[1:]):
      func  = lambda x: x[0] + (intersection-intersections[i])*x[1] 
      min_line = min(lines, key=func)
      if    abs(min_line[0] - lines[-1][0]) > 1e-12*eV \
         or abs(min_line[1] - lines[-1][1]) > 1e-12:
        lines.append([min_line[0].rescale(eV), min_line[1]])

    # adds line after cbm
    func  = lambda x: x[0] + (intersections[-1]+eV)*x[1] 
    min_line = min(_lines, key=func)
    if    abs(min_line[0] - lines[-1][0]) > 1e-12*eV \
       or abs(min_line[1] - lines[-1][1]) > 1e-12:
      lines.append([min_line[0].rescale(eV), min_line[1]])
    return lines

  def enthalpy(self, fermi, mu = None):
    """ Point-defect formation enthalpy. 
    
        :Parameters:
          fermi  
            Fermi energy with respect to the host's VBM. If without
            units, assumes eV. 
          mu : dictionary or None
            Dictionary of chemical potentials. If without units, assumes eV.
            If None, chemical potential part of the formation enthalpy is
            assumed zero.

        :return: Lowest formation enthalpy for all charged states.
    """
    from quantities import eV
    if hasattr(fermi, "rescale"): fermi = fermi.rescale(eV)
    else: fermi = fermi * eV
    return (min(x[0]+fermi*x[1] for x in self.lines()) + self.chempot(mu)).rescale(eV)

  @property
  def latex_label(self):
    """ A label in LaTex format. """
    from re import match
    if self.is_interstitial:
      site = self._site
      if site is None:
        found = match("([A-Z][a-z]?)_interstitial_(.+)$", self.name) 
        return r"{0}$^{{(i)}}_{{ \mathrm{{ {1} }} }}$"\
               .format(found.group(1), found.group(2).replace('_', r"\_"))
      else:
        found = match("([A-Z][a-z]?)_interstitial_(.+)$", self.name) 
        return r"{0}$^{{(i,{2})}}_{{ \mathrm{{ {1} }} }}$"\
               .format(found.group(1), found.group(2).replace('_', r"\_"), site)
    if self.is_substitution:
      found = match("([A-Z][a-z]?)_on_([A-Z][a-z]?)", self.name) 
      site = self._site
      if site is None:
        return r"{0}$_{{ \mathrm{{ {1} }} }}$".format(found.group(1), found.group(2))
      else:
        return r"{0}$_{{ \mathrm{{ {1} }}_{{ {2} }} }}$"\
               .format(found.group(1), found.group(2), site)
    if self.is_vacancy:
      found = match("vacancy_([A-Z][a-z]?)", self.name) 
      site = self._site
      if site is None:
        return r"$\square_{{ \mathrm{{ {0} }} }}$".format(found.group(1))
      else:
        return r"$\square_{{ \mathrm{{ {0} }}_{{{1}}} }}$".format(found.group(1), site)
  
  def __str__(self):
    """ Energy and corrections for each charge defect. """
    from operator import itemgetter
    from os.path import relpath
    from numpy import array
    from numpy.linalg import det
    from quantities import eV
    result = "{0}: \n".format(self.name)
    states = sorted(((c, c.charge) for c in self._charged_states),  key = itemgetter(1))
    for extract, charge in states:
      n = int(det(extract.structure.cell)/det(self.host.structure.cell) + 1.e-3) + 0.
      a = float(extract.raw_DeltaH.rescale(eV))
      b = float(extract.charge_corrections.rescale(eV))
      c = float(extract.potential_alignment.rescale(eV))
      d = float(extract.band_filling.rescale(eV))
      e = relpath(extract.directory, extract.directory + "/../../")
      result += "  - charge {0:>3}: DeltaH = {1:8.4f} + {2:8.4f} + {3:8.4f}"\
                "+ {4:8.4f} = {5:8.4f} eV # {6}.\n"\
                .format(int(charge), a, b, c, d, a+b+c+d, e)
    return result


class PointDefectMassExtractorImpl(AbstractMassExtract):
  """ Enthalpy for a series of defects for a given material and lattice. """
  def __init__(self, path=None, epsilon = 1e0, pa_maxdiff=0.5, Extractor=None, **kwargs):
    """ Initializes an enthalpy function. """
    from pylada.vasp import MassExtract
    super(PointDefectMassExtractor, self).__init__(**kwargs)

    self.Extractor = Extractor
    """ Class for extracting data from a single defect. """
    if self.Extractor is None: self.Extractor = PointDefectExtractor
    self.massextract = MassExtract(path, unix_re=False, excludes=[".*relax_*"])
    """ Mass extraction object from which all results are pulled. """
    self.host = self._get_host()
    """ Result of host calculations. """

    # must be last. Can't use porperty setter.
    self._epsilon = epsilon
    self._pa_maxdiff = pa_maxdiff


  @property 
  def epsilon(self): 
    """ Dimensionless dielectric constant. """
    return self._epsilon
  @epsilon.setter
  def epsilon(self, value):
    self._epsilon = value 
    for v in self.itervalues(): v.epsilon = self._epsilon

  @property 
  def pa_maxdiff(self): 
    """ Dimensionless dielectric constant. """
    return self._pa_maxdiff
  @pa_maxdiff.setter
  def pa_maxdiff(self, value):
    self._pa_maxdiff = value 
    for v in self.itervalues(): v.pa_maxdiff = self._pa_maxdiff


  @property
  def rootdir(self): 
    """ Path to the root-directory containing the poin-defects. """
    return self.massextract.rootdir
  @rootdir.setter
  def rootdir(self, value): self.massextract.rootdir = value
    
  def _get_host(self):
    """ Returns extraction object towards the host. """
    from operator import itemgetter
    host = self.massextract.copy(excludes=[".*PointDefects"])
    host.excludes.extend(self.massextract.excludes)
    lowest = sorted(host.total_energies.iteritems(), key=itemgetter(1))[0][0]
    host = [u for u in host[lowest].itervalues()]
    assert len(host) == 1
    return host[0]


  def __iter_alljobs__(self):
    """ Walks through point-defects only. """
    for child in self.massextract["PointDefects"].children:
      # looks for site_n
      if len(child["site_\d+"].jobs) != 0:
        assert len(child["site_\d+"].jobs) == len(child.jobs),\
               RuntimeError("Don't understand directory structure of {0}.".format(child.view))
        for site in child.children: # should site specific defects.
          result = self.Extractor(site, self.epsilon, self.host, self.pa_maxdiff)
          # checks this is a point-defect.
          if result.is_interstitial or result.is_vacancy or result.is_substitution:
            yield site.view, result
      else:
        result = self.Extractor(child, host=self.host, pa_maxdiff=self.pa_maxdiff,\
                                epsilon = self.epsilon)
        # checks if this is a valid point-defect.
        if result.is_interstitial or result.is_vacancy or result.is_substitution:
          yield child.view, result

  def ordered_items(self):
    """ Returns items ordered by substitution, vacancy, and interstitial. """
    from operator import itemgetter
    interstitials = (u for u in self.iteritems() if u[1].is_interstitial)
    substitution  = (u for u in self.iteritems() if u[1].is_substitution)
    vacancy       = (u for u in self.iteritems() if u[1].is_vacancy)
    result = sorted(substitution, key = itemgetter(0)) 
    result.extend(sorted(vacancy, key = itemgetter(0)))
    result.extend(sorted(interstitials, key = itemgetter(0)))
    return result
  def ordered_keys(self):
    """ Returns keys ordered by substitution, vacancy, and interstitial. """
    return [u[0] for u in self.ordered_items()]
  def ordered_values(self):
    """ Returns values ordered by substitution, vacancy, and interstitial. """
    return [u[1] for u in self.ordered_items()]

  def __str__(self): 
    """ Prints out all energies and corrections. """
    return "".join( str(value) for value in self.ordered_values() )
      
        


class PointDefectMassExtractor(PointDefectMassExtractImpl):
  """ Enthalpy for a series of defects for a given material and lattice. """
  def __init__(self, **kwargs):
    """ Initializes an enthalpy function. """
    super(PointDefectMassExtractor, self).__init__(**kwargs)

  def enthalpies(self, fermi, mu=None):
    """ Dictionary of point-defect formation enthalpies. 
    
        :Parameters:
          fermi  
            Fermi energy with respect to the host's VBM. If without
            units, assumes eV. 
          mu : dictionary or None
            Dictionary of chemical potentials. If without units, assumes eV.
            If None, chemical potential part of the formation enthalpy is
            assumed zero.

        :return: Dictionary where keys are the name of the defects, and the
          values the formation enthalpy.
    """
    from quantities import eV
    results = {}
    for name, defect in self.iteritems():
      results[name] = defect.enthalpy(fermi, mu).rescale(eV)
    return results
  
  def plot_enthalpies(self, mu=None, **kwargs):
    """ Plots diagrams using matplotlib. """
    from quantities import eV
    try: import matplotlib.pyplot as plt
    except ImportError: 
      print "No matplotlib module."
      return
    from operator import itemgetter

    # sets up some stuff for legends.
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble="\usepackage{amssymb}")
    # finds limits of figure
    xlim = 0., float( (self.host.cbm-self.host.vbm).rescale(eV) ) 
    all_ys = [float(val.rescale(eV)) for x in xlim for val in self.enthalpies(x, mu).itervalues()]
    ylim = min(all_ys), max(all_ys)
    # creates figures and axes.
    figure = plt.figure()
    axes = figure.add_subplot(111, xlim=(self.host.vbm, self.host.cbm), ylim=ylim)

    # loop over defects.
    for defect in self.ordered_values():
      # finds intersection points. 
      x = [-5e0*eV]
      lines = defect.lines()
      for i in range(len(lines)-1):
        (b0, a0), (b1, a1) = lines[i], lines[i+1]
        x.append( ((b0-b1)/(a1-a0)).rescale(eV) - self.host.vbm)
      x.append(5e0*eV)

      # Now draws lines. 
      lines.append(lines[-1])
      y = [u[0] + u[1] * xx for u, xx in zip(lines, x)]
      axes.plot(x, y, label=defect.latex_label, **kwargs)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    plt.legend()
    plt.draw()

      
        
















