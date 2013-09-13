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

""" Extraction wrapper around database OUTCAR. """
__all__ = ['VaspExtract']

from ..opt import AbstractExtractBase, OutcarSearchMixin
from ..vasp.extract._common import Extract as ExtractCommonBase
from ..vasp.extract._dft import Extract as ExtractDFTBase
from ..vasp.extract._gw import Extract as ExtractGWBase


class StringStream(object):
  """ Converts string to file-like object. """
  def __init__(self, *args, **kwargs):
    """ Initializes the object using the stream. """
    from cStringIO import StringIO
    super(StringStream, self).__init__()
    self._stringio = StringIO(*args, **kwargs)
    """ String-as-a-filestream. """

  def __enter__(self):
    """ Makes this a context returning a cStringIO. """
    return self._stringio
  def __exit__(self, *args, **kwargs):
    """ Makes this a context returning a cStringIO. """
    self._stringio.close()

  def __getattr__(self, name):
    """ Forwards to stringio. """
    return getattr(self._stringio, name)

class ExtractCommon(AbstractExtractBase, ExtractCommonBase, OutcarSearchMixin):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    from bson import ObjectId
    from .misc import get_ladabase
    AbstractExtractBase.__init__(self, comm=comm)
    ExtractCommonBase.__init__(self)
    OutcarSearchMixin.__init__(self)
    # figure out how the filter has been specified.
    filter = filter if isinstance(filter, dict) else {'_id': ObjectId(filter)}
    # retrives first element with that filter.
    iterator = get_ladabase().files.find(filter)
    try: element = iterator.next()
    except StopIteration:
      raise ValueError("Could not find database document answering "\
                       "to the following filter:\n{0}.\n".format(filter))
    try: iterator.next()
    except StopIteration: pass
    else:
      raise ValueError("Found more than one document answering "\
                       "to the following filter:\n{0}.".format(filter))
    self._element = element
    """ Elements in database associated with current OUTCAR. """

  def __outcar__(self):
    """ File-like object for database outcars. """
    from .misc import get_ladabase
    if self.compression == "bz2":  
      from bz2 import decompress
      return StringStream(decompress(get_ladabase().outcars.get(self._element['_id']).read()))
    elif self.compression is None: 
      return StringStream(get_ladabase().outcars.get(self._element['_id']).read())
    elif self.compression.lower() == "none":
      return StringStream(get_ladabase().outcars.get(self._element['_id']).read())
  def __funccar__(self):
    raise IOError('FUNCCARs are not stored in the database.')
  def __contcar__(self):
    raise IOError('CONTCARs are not stored in the database.')

  def __getattr__(self, name):
    """ Adds read-only properties corresponding to database metadata. """
    if name == "comment": return self._element["comment"]["text"]
    if name in self._element: return self._element[name] 
    raise AttributeError()

  def __dir__(self):
    """ Attributes in this object. """
    result = set([u for u in self._element.keys() if u[0] != '_'] + ['_id', 'success']) \
             | set([u for u in dir(self.__class__) if u[0] != '_'])
    return list(result)

  @property
  def success(self):
    """ True if calculation was successfull. """
    return ExtractCommonBase.success.__get__(self)

  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}(\"{0._id}\")".format(self)

class ExtractDFT(ExtractCommon, ExtractDFTBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, filter, comm, **kwargs)
    ExtractDFTBase.__init__(self)

class ExtractGW(ExtractCommon, ExtractGWBase):
  """ Extracts GW data from an OUTCAR. """
  def __init__(self, filter, comm=None, **kwargs):
    """ Initializes extraction object. """
    ExtractCommon.__init__(self, filter, comm, **kwargs)
    ExtractGWBase.__init__(self)

def VaspExtract(filter, *args, **kwargs): 
  """ Chooses between DFT or GW extraction object, depending on OUTCAR. """
  a = ExtractCommon(filter, *args, **kwargs)
  try: which = ExtractDFT if a.is_dft else ExtractGW
  except: which = ExtractCommon
  return which(filter, *args, **kwargs)
