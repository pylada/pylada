""" Subpackage grouping single-electron properties. """
__docformat__ = "restructuredtext en"
__all__ = ['Band', 'NewK', 'Rdfmwf']
from ..input import BaseKeyword, BoolKeyword
from ..electronic import Shrink

class Band(BaseKeyword):
  """ Creates a band-structure object. """
  keyword = 'band'
  """ CRYSTAL keyword. """
  def __init__( self, title=None, nsub=None, iss=None, minband=None,
                maxband=None, lines=None ):
    super(Band, self).__init__()

    self.title = title
    """ Title to give the band-structure. 

        This is only of interest to the BAND.DAT file.
    """
    self.nsub = nsub
    """ Number of points along full path. """
    self.iss = iss
    """ Parameter to avoid giving floating points. 

        If an integer, the line should be given in fractional coordinates in
        units of this point, and integers. If None, then the lines should be
        given in fractional coordinates. 
    """
    self.minband = minband
    """ Minimum band to output. """
    self.maxband = maxband
    """ Maximum band to ouput. """
    self.lines = lines

  def _check_values(self, value):
    """ Check input values are correct. """
    from numpy import array
    type = 'float64' if self.iss is None else 'int64'
    value = array(value, dtype=type)
    if value.shape[-1] is None:
      raise ValueError('Last dimension of array should be 1, 2, or 3.')
    if value.shape[-1] < 1 or value.shape[-1] > 3:
      raise ValueError('Last dimension of array should be 1, 2, or 3.')
    try: value.reshape(-1, 2, value.shape[-1])
    except: raise ValueError("Could not reshape array to n by 2 by 3.")
    return value.reshape(-1, value.shape[-1])

  def __iadd__(self, value):
    """ Adds points to the band-structure. """
    from numpy import concatenate
    value = self._check_values(value)
    return value if self._lines is None else concatenate((self._lines, value))
  def __set__(self, instance, value):
    """ Sets k-point path. """
    value = self._check_values(value)
    self._lines = value

  @property
  def lines(self): 
    """ Lines for which to perform calculation. """
    return self._lines
  @lines.setter
  def lines(self, value):
    if value is None: self._lines = None; return
    value = self._check_values(value)
    self._lines = value.copy()

  def output_map(self, structure=None, properties=None, input=None, **kwargs):
    """ Returns an output map. """
    if self._lines is None or len(self.lines) == 0: return None
    formats = {}
    formats["title"] = self.title if self.title is not None                    \
                       else getattr(structure, 'name', "")
    formats["iss"] = int(1e8) if self.iss is None else int(self.iss)
    formats["nsub"] = len(self.lines) * 10 if self.nsub is None else self.nsub
    formats["inzb"] = 0 if self.minband is None else int(self.minband)
    nelectrons = properties.nbelectrons(structure, input)
    formats["ifnb"] = nelectrons + 6  if self.maxband is None                  \
                      else int(self.maxband)
    formats["iplo"] = 1
    formats["lpr66"] = 0
    formats["nlines"] = len(self.lines) // 2
    string = "{title}\n{nlines} {iss} {nsub} {inzb} {ifnb} {iplo} {lpr66}\n"   \
             .format(**formats)
    for a, b in zip(self.lines[::2], self.lines[1::2]):
      if self.iss is None:
        a = (formats['iss'] * a).astype('int64')
        b = (formats['iss'] * b).astype('int64')
      string += "{0[0]: <10} {0[1]: <10} {0[2]: <10}   "                       \
                "{1[0]: <10} {1[1]: <10} {1[2]: <10}\n".format(a, b)
    return {self.keyword: string}
      
  def read_input(self, value, owner=None, **kwargs):
    """ Reads from an input string. """
    from numpy import array
    self.title = value[:value.find('\n')].rstrip().lstrip()
    value = value[:value.find('\n')], value[value.find('\n')+1:].split()
    nlines = int(value[0])
    self.iss = int(self.value[1])
    if self.iss == int(1e6): self.iss = None
    self.nsub = int(value[2])
    self.minband = int(value[3])
    if self.minband == 0: self.minband = None
    self.maxband = int(value[4])
    self._lines = array(value[5:3*2*nlines], dtype='int64')
    if self.iss is not None:
      self._lines = self._lines.astype('float64') / float(self.iss)

  def save_files(self, **kwargs):
    """ Copies files appropriately. 

        Should be called with same arguments as output_map.
    """
    from os.path import exists, join
    from ...misc import copyfile
    from ... import CRYSTAL_propnames as propnames
    if self.output_map(**kwargs) is None: return 
    workdir = kwargs['workdir']
    outdir = kwargs['outdir']
    if exists(join(workdir, 'BANDS.DAT')):
      copyfile( join(workdir, 'BANDS.DAT'), 
                join(outdir, propnames['BANDS.DAT'].format('crystal')),
                nothrow='same' )

  def __repr__(self):
    """ Dumps representation to string. """
    args = []
    if self.title is not None and len(self.title) > 0: 
      args.append(repr(self.title))
    if self.nsub is not None: 
      if len(args) != 1: args.append("nsub={0.nsub!r}".format(self))
      else: args.append("{0.nsub!r}".format(self))
    if self.iss is not None: 
      if len(args) != 2: args.append("iss={0.iss!r}".format(self))
      else: args.append("{0.iss!r}".format(self))
    if self.minband is not None: 
      if len(args) != 3: args.append("minband={0.minband!r}".format(self))
      else: args.append("{0.minband!r}".format(self))
    if self.maxband is not None: 
      if len(args) != 4: args.append("maxband={0.maxband!r}".format(self))
      else: args.append("{0.maxband!r}".format(self))
    if self._lines is not None and len(self._lines): 
      if len(args) != 5: args.append("lines={0!r}".format(self._lines.tolist()))
      else: args.append("{0!r}".format(self._lines.tolist()))
    return "{0.__class__.__name__}({1})".format(self, ", ".join(args)) 

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    from ...tools.uirepr import add_to_imports
    if defaults is not None: 
      if type(defaults) is not type(self): 
        add_to_imports(self, imports)
        return {name: repr(self)}
      if self.title == defaults.title and self.nsub == defaults.nsub           \
         and self.iss == defaults.iss and self.minband == defaults.minband     \
         and self.maxband == defaults.maxband:
        if self._lines is None or len(self._lines) == 0: return {}
        else: return {name: repr(self._lines.tolist())}
      return {name: repr(self)}
    elif name is None:
      add_to_imports(self, imports)
      return {None: 'band = {0!r}'.format(self)}
    add_to_imports(self, imports)
    return {name: self.__repr__()}

class NewK(Shrink):
  """ Recomputes on new k-point grid. """
  keyword = 'newk'
  """ CRYSTAL keyword. """
  def __init__(self, mp=None, gallat=None, recompute_fermi=False, printing=None):
    from ..input import SetPrint
    super(NewK, self).__init__(mp=mp, gallat=gallat)
    self.recompute_fermi = recompute_fermi
    """ Whether or not to recompute fermi energy. """
    self._printing = SetPrint(printing)
    """ Map of printing options, if any. """
  
  @property
  def printing(self): 
    """ Map of printing options, if any. """
    return self._printing.__get__(self)
  @printing.setter
  def printing(self, value):
    return self._printing.__set__(self)
  @printing.deleter
  def printing(self):
    return self._printing.__delete__(self)

  def output_map(self, **kwargs):
    result = super(NewK, self).output_map(**kwargs)
    if result is None: return None
    string = result[self.keyword].rstrip()
    string +=  '\n1' if self.recompute_fermi else '\n0'
    if len(self.printing) != 0:
      other = self.printing.output_map(**kwargs)
      string +=  ' ' + other[self.printing.keyword].rstrip().lstrip() + '\n'
    else: string += ' 0\n'
    result[self.keyword] = string
    return result
  def read_input(self, value, **kwargs):
    self.raw = value
    value = " ".join(value.split()[2 if isinstance(self.mp, int) else 5:])
    self.printing.read_input(value, **kwargs)

  def __repr__(self):
    """ Dumps representation of this object. """
    args = []
    if self.mp is not None: args.append(repr(self.mp))
    if self.gallat is not None:
      if len(args) == 0: args.append('gallat={0.gallat!r}'.format(self))
      else: args.append('{0.gallat!r}'.format(self))
    if self.recompute_fermi: 
      args.append('recompute_fermi=True' if len(args) != 2 else 'True')
    if len(self.printing) != 0:
      if len(args) != 3:
        args.append('printing={0.printing.options!r}'.format(self))
      else: args.append('{0.printing.options!r}'.format(self))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creaters user-friendly representation. """
    from ...tools.uirepr import add_to_imports
    if name is None: name = 'newk'
    result = super(NewK, self).__ui_repr__( imports, name=name, 
                                            defaults=defaults, exclude=exclude )
    if len(self.printing) == 0 and not self.recompute_fermi: return result
    add_to_imports(self, imports)
    return {result.keys()[0]: repr(self)}


class Rdfmwf(BoolKeyword):
  """ Whether to read formated wavefunctions. """
  keyword = 'rdfmwf'
  """ CRYSTAL keyword. """
  def __init__(self, value=None):
    super(BoolKeyword, self).__init__(value=value)
  def output_map(self, **kwargs):
    from os.path import join
    from ...misc import copyfile, latest_file
    from ...error import IOError

    value = self.value
    if value is None: 
      which = latest_file( join(kwargs['input'].directory, 'crystal.f9'), 
                           join(kwargs['input'].directory, 'crystal.f98'), 
                           join(kwargs['outdir'], 'crystal.f9'), 
                           join(kwargs['outdir'], 'crystal.f98') )
    elif value:
      which = latest_file( join(kwargs['input'].directory, 'crystal.f98'), 
                           join(kwargs['outdir'], 'crystal.f98') )
    else:
      which = latest_file( join(kwargs['input'].directory, 'crystal.f9'), 
                           join(kwargs['outdir'], 'crystal.f9') )
    if which is not None:
      value = which[which.rfind('.'):] == '.f98'
      if kwargs.get('filework', False):
        copyfile(which, kwargs['workdir'], nothrow='same')
    elif kwargs.get('filework', False): 
      raise IOError('Could not find input density.')
    return None if value is False or value is None else {self.keyword: True}
