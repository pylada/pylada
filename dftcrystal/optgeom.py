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

from .input import AttrBlock, BoolKeyword
from ..tools.input import BoolKeyword as BaseBoolKeyword, TypedKeyword

class GeometryOpt(BaseBoolKeyword):
  """ Geometry optimization keyword. 

      When set to True, makes sure that:

        1. MAXCYCLE is non-zero. If it is, sets it to 100.
        2. only one of FULLOPTG, ITATOCELL, ITREDUN, CELLONLY is True.
        3. CVOLOPT is False unless this is a FULLOPTG or CELLONLY calculation.

      \(2\) and (3) are decided by the keyword argument in input.
  """
  def __init__(self, keyword, value=False):
    """ Initializes a geometry optimization keyword. """
    self.keyword = keyword
    """ CRYSTAL keyword. """
    super(GeometryOpt, self).__init__(value=value)

  def __set__(self, instance, value):
    """ Sets keyword to appear or not. 

        Also clears attendant conditions described in :py:class:`GeometryOpt`.
    """
    super(GeometryOpt, self).__set__(instance, value)
    if self.value:
      if instance.maxcycle is not None and instance.maxcycle < 1:
        instance.maxcycle = 100
      if self.keyword != 'fulloptg':  instance.fulloptg  = False
      if self.keyword != 'cellonly':  instance.cellonly  = False
      if self.keyword != 'itatocell': instance.itatocell = False
      if self.keyword != 'intredun':  instance.intredun  = False
      if self.keyword != 'fulloptg' and self.keyword != 'cellonly':
        instance.cvolopt = False

class CVolOpt(BoolKeyword):
  """ Constant volume optimization keyword. 
  
      Only appears if FULLOPTG or CELLONLY exist.
  """
  keyword = 'cvolopt'
  """ Crystal input keyword """
  def output_map(self, **kwargs):
    """ Only prints if FULLOPTG or CELLONLY exist. """
    optgeom = kwargs['crystal'].optgeom
    if optgeom.fulloptg == False and optgeom.cellonly == False:
      return None
    return super(CVolOpt, self).output_map(**kwargs)

class FixCell(BoolKeyword):
  """ Constant volume optimization keyword. 
  
      Only appears for INTREDUN relaxations.
  """
  keyword = 'fixcell'
  """ Crystal input keyword """
  def output_map(self, **kwargs):
    optgeom = kwargs['crystal'].optgeom
    if optgeom.intredun == False: return None
    return super(FixCell, self).output_map(**kwargs)

class MaxCycle(TypedKeyword):
  """ Maxcycle input.

      It is expected to be an integer or None.
  """
  keyword = 'maxcycle'
  """ CRYSTAL keyword. """ 
  type = int
  """ Type of the keyword. """

class ExclAttrBlock(AttrBlock):
  """ An attribute block set up to exclude others. 
  
      Expects both "keyword" and "excludegroup" attributes (or class
      attributes) to exist in derived instances. This class makes it convenient
      to describe mutually exclusive blocks such as optgeom and freqcalc.
  """
  excludegroup = 'optgeom', 'freqcalc', 'anharm', 'confcnt', 'cphf',           \
                 'elastcon', 'eos'
  """ Groups of blocks of which only one should be enabled. """
  def __init__(self):
    """ Initializes the exclusive attribute block. """
    super(ExclAttrBlock, self).__init__()
    self.enabled = False
    self._parent = None
    """ Weak-reference to parent instance. """

  @property
  def enabled(self):
    """ True  if this block is enabled. """
    return self._enabled
  @enabled.setter
  def enabled(self, value):
    self._enabled = value == True
    # disable other instances. 
    if self._enabled and self._parent is not None                               \
       and len(getattr(self, 'excludegroup', ())) > 0:
      parent = self._parent()
      for u in self.excludegroup:
        if u == self.keyword: continue
        if not hasattr(parent, u): continue
        inst = getattr(parent, u)
        if not hasattr(inst, 'enabled'): continue
        inst.enabled = False
  def __get__(self, instance, owner=None):
    """ Sets up weak ref to parent. """
    from weakref import ref
    self._parent = ref(instance)
    return self

  def read_input(self, tree, owner=None, **kwargs):
    """ Reads from input. """
    if owner is not None:
      from weakref import ref
      self._parent = ref(owner)
    self.enabled = True
    return super(ExclAttrBlock, self).read_input(tree, owner, **kwargs)

  def output_map(self, **kwargs):
    """ Does not print if disabled. """
    from ..tools.input import Tree
    if not self.enabled: return None
    result = super(ExclAttrBlock, self).output_map(**kwargs)
    if result is None:
      result = Tree()
      result[self.keyword] = Tree()
    return result

  def __getstate__(self):
    d = self.__dict__.copy()
    d['_parent'] = None
    return d
  def __setstate__(self, value):
    self.__dict__.update(value)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly output. """
    # make sure enabled is printed.
    if defaults is not None: defaults.enabled = not self.enabled
    return super(ExclAttrBlock, self).__ui_repr__( imports, name,              \
                                                   defaults, exclude )

class OptGeom(ExclAttrBlock):
  """ Geometry Optimization block. 
  
      Defines the input block for geometry optimization. Geometry optimization
      must be explicitely enabled:

      >>> functional.optgeom.enabled = True

      It is disabled by default. When enabled, other sub-blocks within
      CRYSTAL_'s first input block are automatically disabled (e.g. ``freqcalc`` [*]_). 
      The full set of geometry optimization parameters can be accessed by the
      user within the :py:attr:`~pylada.dftcrystal.functional.Functional.optgeom`
      attribute of an :py:class:`~pylada.dftcrystal.functional.Functional`
      instance:

      >>> functional.optgeom.fulloptg = True
      >>> functional.optgeom.maxcycle = 300
      >>> functional.optgeom.enabled  = True

      The above would set the optimization method to ``FULLOPTG`` and the
      number of geometry optimization to 300. The third line would enable the
      block so that the geomtry optimization can happen.

      .. note::

        Inner keywords can be modified when the block is disabled. However, the
        block will not appear in the input until is is explicitely enabled.

      .. [*]
      
         Currently, :py:class:`~pylada.dftcrystal.functional.Functional` only
         contains the :py:attr:`~pylada.dftcrystal.functional.Functional.optgeom`
         attribute. If implemented, other sub-blocks should derive from
         :py:class:`~pylada.dftcrystal.optgeom.ExclAttrBlock`.
  """
  keyword = 'optgeom'
  """ CRYSTAL input keyword (class-attribute) """
  def __init__(self, breaksym=None, keepsymm=None): 
    """ Creates an optimization block. """
    from quantities import UnitQuantity, hartree, angstrom
    from ..tools.input import QuantityKeyword
    from ..error import ValueError

    super(OptGeom, self).__init__()

    self.breaksym = True
    """ Whether to keep or break symmetries.
    
        This value can be True, False, or None. If True, then symmetries will
        be broken. If False, symmetries will be kept during the minimization.
        By default, symmetries are broken.
    """
    if breaksym is None and keepsymm is not None:
      self.breaksym = keepsymm == False
    elif breaksym is not None and keepsymm is None:
      self.breaksym = breaksym == True
    elif breaksym is not None and keepsymm is not None                         \
         and breaksym == keepsymm:
      raise ValueError( 'OptGeom: breaksym and keepsymm '                      \
                        'give opposite instructions.' )
 
    self.maxcycle   = MaxCycle()
    """ Maxium number of iterations in geometry optimization loop. """
    self.fulloptg   = GeometryOpt('fulloptg')
    """ Full optimization. 

        Unless :py:attr:`cellonly` is True, this will optimize both cell
        internal and cell external degrees of freedom. This option is not
        compatible with other optimization methods. Setting it will disable
        other optimization methods:
    
        >>> functional.optgeom.fulloptg = True
        >>> functional.intredun  # disables intredun automatically.
        False
    """
    self.cellonly   = GeometryOpt('cellonly')
    """ Cell-shape optimization only. 
    
        Constrains optimization to cell-external degrees of freedom only. The
        fractional coordinates of the atomic sites are kept constant.
    """
    self.itatocell  = GeometryOpt('itatocell')
    """ Iterative cell-shape/atom optimization.
    
        Alternates between relaxing cell-internal and cell-external degrees of
        freedom. This option is not compatible with other optimization methods.
        Setting it to ``True`` auttomatically disables other optimization
        methods.

        >>> functional.optgeom.itatocell = True
        >>> functional.fulloptg, functional.intredun 
        False, False
    """
    self.intredun   = GeometryOpt('intredun')
    """ Constrained optimization. 
    
        Relaxes the strain energy by optimizing a complete and redundant set of
        chemically intuitive parameters, such as bond-length and bond-angles. 
        This option is not compatible with other optimization methods.

        >>> functional.optgeom.intredun = True
        >>> functional.fulloptg, functional.itatocell 
        False, False
    """
    self.cvolopt    = CVolOpt()
    """ Constant volume optimization keyword. 
    
        Only appears if FULLOPTG or CELLONLY exist.
    """
    self.fixcell    = FixCell()
    """ Constant volume optimization keyword. 
    
        When used in conjunction with :py:attr:`intredun` optimization method,
        keeps the cell-external degrees of freedom constant. Only appears for
        :py:attr:`intredun` relaxations.
    """
    self.toldee     = TypedKeyword('toldee', int)
    """ Electronic structure minimization convergence criteria. 

        The criteria is with respect to the total energy. It is logarithmic and
        should be an integer: :math:`|\\Delta E| < 10^{-\\mathrm{toldee}}`. 
    """
    self.toldeg     = TypedKeyword('toldeg', float)
    """ Structural relaxation convergence criteria.

        Convergence criteria for the root mean square gradient of the total energy
        with respect to the displacements. Should be a floating point.
    """ 
    self.toldex     = TypedKeyword('toldex', float)
    """ Structural relaxation convergence criteria.

        Convergence criteria as the root mean square of the displacements
        during structural relaxation. Should be a floating point.
    """
    bohr = UnitQuantity('crystal_bohr', 0.5291772083*angstrom, symbol='bohr')
    """ Bohr radius as defined by CRYSTAL """
    self.extpress   = QuantityKeyword(units=hartree/bohr**3)
    """ Hydrostatic pressure in :math:`\\frac{\\text{hartree}}{\\text{bohr}^{3}}`.
    
        Sets the hydrostatic pressure for which an structural relaxation is
        performed.
    """
    self.extstress  = QuantityKeyword(units=hartree/bohr**3, shape=(3,3))
    """ External stress in :math:`\\frac{\\text{hartree}}{\\text{bohr}^{3}}`. 
    
        Sets the stress for which a structural relaxation is performed.
    """
    self.onelog     = BoolKeyword()
    """ Whether to print electronic minimization at each geometry step. 

        The electronic minimization steps can be printed to the same file as
        the geometry optimization step (e.g. the standard output). By default,
        Pylada prefers to print to the same file, in contrast to CRYSTAL_'s
        default.
    """
    self.noxyz       = BoolKeyword()
    """ Whether to print cartesian coordinates at the end of optimization. """
    self.nosymmops   = BoolKeyword()
    """ Whether to print symmetry operations at the end of optimization. """
    self.printforces = BoolKeyword()
    """ Whether to print forces at the end of optimization. """
    self.printhess   = BoolKeyword()
    """ Whether to print Hessian at the end of optimization. """
    self.printopt    = BoolKeyword()
    """ Whether to print extended optimization information. """
    self.verbose     = BoolKeyword(keyword='print')
    """ Verbose printint. """

  @property
  def keepsymm(self): 
    """ Alias to the opposite of breaksym. """
    return self.breaksym != True
  @keepsymm.setter
  def keepsymm(self, value): self.breaksym = value != True

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly output. """
    exclude = ['keepsym', 'keepsymm']
    return super(OptGeom, self).__ui_repr__(imports, name, defaults, exclude)

  def output_map(self, **kwargs):
    """ Reads input tape. """
    from ..tools.input import Tree
    if self.breaksym != kwargs.get('breaksym', True): 
      result = Tree()
      result['breaksym' if self.breaksym else 'keepsymm'] = True
      result.update(super(OptGeom, self).output_map(**kwargs))
      if len(result) == 0: return None
      return result
    return super(OptGeom, self).output_map(**kwargs)
  
  def read_input(self, value, **kwargs):
    """ Reads input tape. """
    self.breaksym = kwargs.get('breaksym', True)
    return super(OptGeom, self).read_input(value, **kwargs)
