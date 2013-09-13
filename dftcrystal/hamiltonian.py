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

__docformat__ = "restructuredtext en"
__all__ = ['Dft']
from .input import AttrBlock, BoolKeyword
from ..tools.input import BaseKeyword, ChoiceKeyword, TypedKeyword

class Exchange(ChoiceKeyword):
  keyword = 'exchange'
  """ Corresponding CRYSTAL keyword. """
  values = ['becke', 'lda', 'pbe', 'pbesol', 'pwgga', 'sogga', 'vbh', 'wcgga']
  """" Set of values from which to choose. """
  def __init__(self, value=None): 
    super(Exchange, self).__init__(value)
    self.excludes = ['b3lyp', 'b3pw', 'pbe0', 'soggaxc']
    """ keywords incompatible with this one. """
  def __set__(self, instance, value): 
    """ Makes sure that excluded keywords are set to False. """
    from ..error import ValueError
    if value is None:
      raise ValueError('Exchange and Correlations cannot be None.')
    super(Exchange, self).__set__(instance, value)
    for u in self.excludes: 
      if getattr(instance, u, False): setattr(instance, u, False)
  def output_map(self, **kwargs):
    """ Print EXCHANGE keyword. 

        Also checks whether one of global keywords is set, in which case does
        not print.
    """
    if self._value is None: return None
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    return {self.keyword.upper(): self.value.upper()}

class Correlation(Exchange):
  keyword = 'correlat'
  """ Corresponding CRYSTAL keyword. """
  values = [ 'lyp', 'p86', 'pbe', 'pbesol', 'pwgga', 'pwlsd', 'pz',
             'vbh', 'wl', 'vwn']
  """ Set of values from which to choose. """
  def __init__(self, value=None):
    super(Correlation, self).__init__(value=value)

class GlobalExc(BaseKeyword):
  """ Takes care of global keywords for exchange-correlation. """
  def __init__(self, keyword, values):
    super(GlobalExc, self).__init__(keyword=keyword)
    self.values = values
    """ Values against which to test. """

  def __get__(self, instance, owner=None):
    """ True if global keyword is set. """
    def test(a, b):
      if a is None:
        if b is not None: return False
      else: 
        if b is None: return False
        if abs(a-b) > 1e-8: return False
      return True
    return instance.exchange == self.values[0]                                 \
           and instance.correlat == self.values[1]                             \
           and test(self.values[2], instance.hybrid)                           \
           and test(self.values[3], instance.nonlocal.exchange)                \
           and test(self.values[4], instance.nonlocal.correlation)
  def __set__(self, inst, value):
    """ Sets to global if value is True. 
    
        If False, ignores.
    """
    if value == True:
      inst.exchange, inst.correlat, inst.hybrid,                               \
        inst.nonlocal.exchange, inst.nonlocal.correlation = self.values
  def output_map(self, **kwargs):
    if self.__get__(kwargs['crystal'].dft): 
      return {self.keyword.upper(): None}
  def read_input(self, tree, owner=None, **kwargs):
    """ True if ever read. """
    self.__set__(owner, True)

class NonLocal(BaseKeyword):
  """ Non-local weight parameters. """
  keyword = 'NONLOCAL'
  def __init__(self, xc=None, corr=None):
    super(NonLocal, self).__init__()
    self.exchange = None
    """ Exchange weight. """
    self.correlation = None
    """ Correlation weight. """
  @property
  def raw(self):
    x = 0 if self.exchange is None else self.exchange
    c = 0 if self.correlation is None else self.correlation
    return '{0} {1}\n'.format(x, c)
  @raw.setter
  def raw(self, value):
    value = value.split()
    self.exchange, self.correlation = float(value[0]), float(value[1])
  def disable(self): 
    """ Sets exchange and correlation to None. """
    self.exchange, self.correlation = None, None
  def output_map(self, **kwargs):
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    x = 0 if self.exchange is None else self.exchange
    c = 0 if self.correlation is None else self.correlation
    if abs(x) < 1e-8 and abs(c) < 1e-8: return None
    return super(NonLocal, self).output_map(**kwargs)
  def __getitem__(self, index):
    """ list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: return self.exchange
    elif index == 1 or index == -1: return self.correlation
    raise IndexError('nonlocal can be indexed with 0, 1, or -1 only.')
  def __setitem__(self, index, value):
    """ sets as list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: self.exchange = value
    elif index == 1 or index == -1: self.correlation = value
    else: raise IndexError('nonlocal can be indexed with 0, 1, or -1 only.')
  def __len__(self): return 2
  def __set__(self, value):
    """ Sets from (x, c) tuple. """
    from ..error import ValueError
    if len(value) != 2: raise ValueError('Expected (xchange, corr) tuple')
    self.exchange = value[0]
    self.correlat = value[1]
  def __repr__(self):
    """ Dumps representation of self. """
    args = []
    if self.exchange is None: 
      if self.correlation is not None:
        args.append('corr={0.correlation!r}'.formas(self))
    else:
      args.append(repr(self.exchange))
      if self.correlation is not None: args.append(repr(self.correlation))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))
  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    from ..tools.uirepr import add_to_imports
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    if defaults is None or type(defaults) != type(self):
      add_to_imports(self, imports)
      return {name: self.__repr__()}
    if self.correlation is None: issame = defaults.correlation is None
    elif defaults.correlation is None: issame = False
    else: issame = abs(defaults.correlation - self.correlation) < 1e-8 
    if issame:
      if self.exchange is None: issame = defaults.exchange is None
      elif defaults.exchange is None: issame = False
      else: issame = abs(defaults.exchange - self.exchange) < 1e-8 
    if issame: return {}
    return {name: '{0.exchange!r}, {0.correlation!r}'.format(self)}
      
            


class Hybrid(TypedKeyword):
  keyword = 'hybrid'
  """ Crystal input keyword. """
  type = float
  """ Keyword must conform to this type, if not None. """
  def __init__(self, value=None):
    super(Hybrid, self).__init__(value=value)
  def output_map(self, **kwargs):
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    return super(Hybrid, self).output_map(**kwargs)

class Radial(BaseKeyword):
  """ Defines radial CRYSTAL keyword. """
  keyword = 'radial'
  def __init__(self, intervals=None, nbpoints=None):
    """ Creates RADIAL keyword """
    super(Radial, self).__init__()
    self.intervals = intervals
    self.nbpoints  = nbpoints
  @property
  def intervals(self): 
    """ Upper limits of the intervals in ascending order """
    return self._intervals
  @intervals.setter
  def intervals(self, value): 
    from numpy import array, all
    from ..error import input
    if value is None: self._intervals = None; return
    self._intervals = [float(u) for u in value]
    # sanity check
    if not all(array(self.intervals[1:]) > array(self.intervals[:-1])):
      raise input('Interval limits are not given in ascending order')
  @property
  def nbpoints(self): 
    """ Number of integration points per interval """
    return self._nbpoints
  @nbpoints.setter
  def nbpoints(self, value): 
    if value is None: self._nbpoints = None; return
    self._nbpoints = [int(u) for u in value]

  @property
  def raw(self):
    """ Raw CRYSTAL input """
    if self.intervals is None or self.nbpoints is None: return ''
    n = min(len(self.intervals), len(self.nbpoints))
    a = (str(u) for u in self.intervals[:n])
    b = (str(u) for u in self.nbpoints[:n])
    return '{0}\n{1}\n{2}'.format(n, ' '.join(a), ' '.join(b))
  @raw.setter
  def raw(self, value):
    """ Sets value from CRYSTAL input """
    value = value.split()
    n = int(value[0])
    self.intervals = value[1:n+1]
    self.nbpoints = value[n+1:]
  def output_map(self, **kwargs):
    """ Prints CRYSTAL input """
    from ..error import input
    if self.intervals is None or self.nbpoints is None: return None
    crystal = kwargs['crystal'].dft
    if crystal.lgrid or crystal.xlgrid or crystal.xxlgrid: return None
    # sanity check
    if len(self.intervals) != len(self.nbpoints):
      raise input('Different number of interval limits and points per interval')
    return super(Radial, self).output_map()
  def __repr__(self):
    """ Dumps representation of self. """
    args = []
    if self.intervals is None: 
      if self.nbpoints is not None:
        args.append('nbpoints={0.nbpoints!r}'.formas(self))
    else:
      args.append(repr(self.intervals))
      if self.nbpoints is not None: args.append(repr(self.nbpoints))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))
  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    from numpy import array, all, abs
    from ..tools.uirepr import add_to_imports
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    if defaults is None or type(defaults) != type(self):
      add_to_imports(self, imports)
      return {name: self.__repr__()}
    if self.nbpoints is None: issame = defaults.nbpoints is None
    elif defaults.nbpoints is None: issame = False
    else: issame = all(array(defaults.nbpoints) == self.nbpoints)
    if issame:
      if self.intervals is None: issame = defaults.intervals is None
      elif defaults.intervals is None: issame = False
      else:
        issame = all(abs(array(defaults.intervals) - self.intervals) < 1e-8)
    if issame: return {}
    return { '{0}.intervals'.format(name): repr(self.intervals), 
             '{0}.nbpoints'.format(name): repr(self.nbpoints) }

class Angular(BaseKeyword):
  """ Defines angular CRYSTAL keyword. """
  keyword = 'angular'
  def __init__(self, intervals=None, levels=None):
    """ Creates RADIAL keyword """
    super(Angular, self).__init__()
    self.intervals = intervals
    self.levels  = levels
  @property
  def intervals(self): 
    """ Upper limits of the intervals in ascending order """
    return self._intervals
  @intervals.setter
  def intervals(self, value): 
    from numpy import array, all
    from ..error import input
    if value is None: self._intervals = None; return
    self._intervals = [float(u) for u in value]
    # sanity check
    if not all(array(self.intervals[1:]) > array(self.intervals[:-1])):
      raise input('Interval limits are not given in ascending order')
  @property
  def levels(self): 
    """ Accuracy level of each interval """
    return self._levels
  @levels.setter
  def levels(self, value): 
    if value is None: self._levels = None; return
    self._levels = [int(u) for u in value]

  @property
  def raw(self):
    """ Raw CRYSTAL input """
    if self.intervals is None or self.levels is None: return ''
    n = min(len(self.intervals), len(self.levels))
    a = (str(u) for u in self.intervals[:n])
    b = (str(u) for u in self.levels[:n])
    return '{0}\n{1}\n{2}'.format(n, ' '.join(a), ' '.join(b))
  @raw.setter
  def raw(self, value):
    """ Sets value from CRYSTAL input """
    value = value.split()
    n = int(value[0])
    self.intervals = value[1:n+1]
    self.levels = value[n+1:]
  def output_map(self, **kwargs):
    """ Prints CRYSTAL input """
    from ..error import input
    if self.intervals is None or self.levels is None: return None
    crystal = kwargs['crystal'].dft
    if crystal.lgrid or crystal.xlgrid or crystal.xxlgrid: return None
    # sanity check
    if len(self.intervals) != len(self.levels):
      raise input('Different number of interval limits and points per interval')
    return super(Angular, self).output_map()
  def __repr__(self):
    """ Dumps representation of self. """
    args = []
    if self.intervals is None: 
      if self.levels is not None:
        args.append('levels={0.levels!r}'.formas(self))
    else:
      args.append(repr(self.intervals))
      if self.levels is not None: args.append(repr(self.levels))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))
  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    from numpy import array, all, abs
    from ..tools.uirepr import add_to_imports
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    if defaults is None or type(defaults) != type(self):
      add_to_imports(self, imports)
      return {name: self.__repr__()}
    if self.levels is None: issame = defaults.levels is None
    elif defaults.levels is None: issame = False
    else: issame = all(array(defaults.levels) == self.levels)
    if issame:
      if self.intervals is None: issame = defaults.intervals is None
      elif defaults.intervals is None: issame = False
      else:
        issame = all(abs(array(defaults.intervals) - self.intervals) < 1e-8)
    if issame: return {}
    return { '{0}.intervals'.format(name): repr(self.intervals), 
             '{0}.levels'.format(name): repr(self.levels) }

class GlobalGridKeyword(BaseKeyword):
  """ Defines global grid keywords """
  def __init__(self, values=None, keyword=None):
    from copy import deepcopy
    super(GlobalGridKeyword, self).__init__(keyword=keyword)
    self.values = deepcopy(values)
    """ Values against which to test. """

  def __get__(self, instance, owner=None):
    """ True if global keyword is set. """
    def test(a, b):
      from numpy import array, abs, any
      if a is None:
        if b is not None: return False
      else: 
        if b is None: return False
        if len(a) != len(b): return False
        if any(abs(array(a)-array(b)) > 1e-8): return False
      return True
    return test(self.values[0], instance.radial.intervals)                     \
           and test(self.values[1], instance.radial.nbpoints)                  \
           and test(self.values[2], instance.angular.intervals)                \
           and test(self.values[3], instance.angular.levels)
  def __set__(self, inst, value):
    """ Sets to global if value is True. 
    
        If False, ignores.
    """
    if value == True:
      inst.radial.intervals, inst.radial.nbpoints,                             \
        inst.angular.intervals, inst.angular.levels =  self.values
  def output_map(self, **kwargs):
    if self.__get__(kwargs['crystal'].dft): 
      return {self.keyword: True}
  def read_input(self, tree, owner=None, **kwargs):
    """ True if ever read. """
    self.__set__(owner, True)

class Dft(AttrBlock):
  """ DFT attribute block. 
  
      Holds parameters relative to the Hamiltonian used in the calculation.
      These parameters can be accessed indirectly in tha functional:

      >>> functional.dft.b3lyp = True
      >>> functional.dft.xxlgrid = True
     
      The above would set the Hamiltonian to the B3LYP exchange-correlation
      functional and the integration grid to 64 fl.oz.
  """ 
  keyword = 'DFT'
  """ CRYSTAL keyword. """
  def __init__(self):
    """ Creates the DFT attribute block. """
    super(Dft, self).__init__()
    self.exchange = Exchange() 
    """ Exchange functional. 
    
        Name of the exchange functional to use.  It should be one of the
        following: 'becke', 'lda', 'pbe', 'pbesol', 'pwgga', 'sogga', 'vbh',
        'wcgga', or None(default).
    """
    self.correlat = Correlation() 
    """ Correlation functional. 
    
        Name of the correlation functional to use. It can be one of the
        following: 'lyp', 'p86', 'pbe', 'pbesol', 'pwgga', 'pwlsd', 'pz',
        'vbh', 'wl', 'vwn', or None. 
    """
    self.hybrid   = Hybrid()
    """ Amount of exchange to add to functional. 
    
        It should be a floating point or None.
    """
    self.nonlocal = NonLocal()
    """ Non-local weights on exchange-correlation. 
    
        A tuple of two floating point which sets the weights of the non-local
        part of the exchange (first) and of the correlations (second). It can
        also be None.
    """
    self.spin     = BoolKeyword(keyword='spin')
    """ If True, then perform spin-polarized calculation. 
    
        It should be None(default), True, or False.
    """
    self.b3lyp    = GlobalExc('b3lyp', ['becke', 'lyp', 20, 0.9, 0.81])
    """ B3LYP global keyword.
    
        Sets :py:attr:`exchange`, :py:attr:`correlat`, :py:attr:`nonlocal`
        correctly for this functional. Since it acts on/checks upon other
        attributes, it can only be set to True.
    """
    self.b3pw     = GlobalExc('b3pw', ['becke', 'pwgga', 20, 0.9, 0.81])
    """ B3PW global keyword.
    
        Sets :py:attr:`exchange`, :py:attr:`correlat`, :py:attr:`nonlocal`
        correctly for this functional. Since it acts on/checks upon other
        attributes, it can only be set to True.
    """
    self.pbe0     = GlobalExc('pbe0', ['pbe', 'pbe', None, None, None])
    """ PBE0 global keyword. 
    
        Sets :py:attr:`exchange`, :py:attr:`correlat`, :py:attr:`nonlocal`
        correctly for this functional. Since it acts on/checks upon other
        attributes, it can only be set to True.
    """
    self.soggaxc  = GlobalExc('soggaxc', ['sogga', 'pbe', None, None, None])
    """ SOGGAXC global keyword. 
    
        Sets :py:attr:`exchange`, :py:attr:`correlat`, :py:attr:`nonlocal`
        correctly for this functional. Since it acts on/checks upon other
        attributes, it can only be set to True.
    """
    self.angular  = Angular()
    """ Angular integration grid.
    
        Contains two attributes ``intervals`` and ``levels`` which can be used
        to set the angular grid.
    """
    self.radial   = Radial()
    """ Radial integration grid.
    
        Contains two attributes ``intervals`` and ``nbpoints`` which can be used
        to set the radial integration grid.
    """
    self.lgrid    = GlobalGridKeyword(( [4], [75],
                                        [0.1667, 0.5, 0.9, 3.05, 9999.0], 
                                        [2, 6, 8, 13, 8] ))
    """ Preset large integration grid. """
    self.xlgrid   = GlobalGridKeyword(( [4], [75],
                                        [0.1667, 0.5, 0.9, 3.5, 9999.0], 
                                        [2, 8, 12, 16, 12] ))
    """ Preset extra large integration grid. """
    self.xxlgrid  = GlobalGridKeyword(( [4], [99],
                                        [0.1667, 0.5, 0.9, 3.5, 9999.0], 
                                        [6, 8, 14, 18, 14] ))
    """ Preset extra extra large integration grid. """
    self.tollgrid = TypedKeyword(type=int)
    """ DFT grid weight tolerance.
    
        Should be None(default) or an integer.
    """
    self.tolldens = TypedKeyword(type=int)
    """ DFT density tolerance 

        Should be None(default) or an integer.
    """

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    exclude = [] if exclude is None else list(exclude)
    if not self.b3lyp: exclude.append('b3lyp')
    if not self.b3pw: exclude.append('b3pw')
    if not self.pbe0: exclude.append('pbe0')
    if not self.soggaxc: exclude.append('soggaxc')
    if self.b3lyp or self.b3pw or self.pbe0 or self.soggaxc: 
      exclude.extend(['exchange', 'correlat', 'hybrid', 'nonlocal'])
    if not self.lgrid: exclude.append('lgrid')
    if not self.xlgrid: exclude.append('xlgrid')
    if not self.xxlgrid: exclude.append('xxlgrid')
    if self.lgrid or self.lgrid or self.xxlgrid:
      exclude.extend(['angular', 'radial'])
    return super(Dft, self).__ui_repr__(imports, name, defaults, exclude)

