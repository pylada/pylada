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

""" PyMongo interface. """
__all__ = ['Manager', 'MassExtract', 'VaspExtract']

from .vasp import VaspExtract
from .massextract import MassExtract

class Manager(object): 
  """ Holds data regarding database management. """
  def __init__(self, host=None, port=None, database=None, username=None, prefix=None): 
    """ Initializes a connection and database. """
    from pymongo import Connection
    from gridfs import GridFS
    from .. import pymongo_host, pymongo_port, vasp_database_name,  OUTCARS_prefix, pymongo_username
    super(Manager, self).__init__()

    self._host = host if host is not None else pymongo_host
    """ Host where the database is hosted. """
    self._port = port if port is not None else pymongo_port
    """ Port of the host where the database is hosted. """
    self._vaspbase_name = database if database is not None else vasp_database_name
    """ Name of the vasp database. """
    self._outcars_prefix = prefix if prefix is not None else OUTCARS_prefix
    """ Username for database. """
    self._username = username if username is not None else pymongo_username
    """ Name of the OUTCAR database. """
    self.connection = Connection(host=self._host, port=self._port)
    """ Holds connection to pymongo. """
    self.database = getattr(self.connection, self._vaspbase_name)
    """ Database within pymongo. """
    self.outcars = GridFS(self.database, self.outcars_prefix)
    """ GridFS object for OUTCARs. """
    self.comments = self.database["{0}.comments".format(self.outcars_prefix)]
    """ Collection of comments attached when adding OUTCAR's to a file. """
    self.files = self.database["{0}.files".format(self.outcars_prefix)]
    """ OUTCAR files collection. """
    self.extracted = self.database["extracted"]
    """ Collection with pre-extracted values from the outcar. """
    self.fere = self.database["fere_summary"]
    """ Collection with FERE ground-state analysis. """

  @property
  def host(self):
    """ Host where the database is hosted. """
    return self._host
  @property
  def port(self):
    """ Port of the host where the database is hosted. """
    return self._port
  @property
  def vasp_database_name(self):
    """ Name of the vasp database. """
    return self._vaspbase_name
  @property
  def outcars_prefix(self):
    """ Name of the OUTCAR GridFS collection. """
    return self._outcars_prefix

  def push(self, filename, outcar, comment, **kwargs):
    """ Pushes OUTCAR to database. 

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.
    """
    from hashlib import sha512
    from os import uname
    from .misc import get_username

    assert len(comment.replace(' ', '').replace('\n', '')) != 0,\
           ValueError("Cannot push file if comment is empty.")
    
    try: kwargs["comment"] = self.comments.find({'text': comment}).next()
    except StopIteration: # add comment to database
      kwargs["comment"] = self.comments.insert({'text': comment})
      print kwargs['comment']

    hash = sha512(outcar).hexdigest()
    if self.outcars.exists(sha512=hash): 
      print "{0} already in database. Please use 'ladabase.update'.".format(filename)
      return 
    if 'filename' not in kwargs: kwargs['filename'] = filename
    if 'uploader' not in kwargs: kwargs['uploader'] = get_username()
    if 'host'     not in kwargs: kwargs['host']     = uname()[1]
    compression = kwargs.get('compression', None)
    kwargs['compression'] = compression
    if compression == "bz2": 
      from bz2 import compress
      return self.outcars.put(compress(outcar), sha512=hash, **kwargs)
    elif compression is None: return self.outcars.put(outcar, sha512=hash, **kwargs)
    else: raise ValueError("Invalid compression format {0}.".format(compression))

  def find_fromfile(self, path):
    """ Returns the database object corresponding to this file.

        :raise ValueError:  if the corresponding object is not found.
        :raise IOError:  if the path does not exist or is not a file.

        Finds the corresponding file using sha512 hash. 
    """
    from os.path import exists, isfile
    from hashlib import sha512
    from ..opt import RelativeDirectory

    ipath = RelativeDirectory(path).path
    assert exists(ipath), IOError('{0} does not exist.'.format(path))
    assert isfile(ipath), IOError('{0} is not a file.'.format(path))

    with open(ipath, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()
   
    assert self.outcars.exists(sha512=hash),\
           ValueError('{0} could not be found in database.'.format(path))

    return self.files.find_one({'sha512': hash})

  def __contains__(self, path):
    """ True if path already in database. """
    from os.path import exists
    from hashlib import sha512
    from ..opt import RelativeDirectory
    path = RelativeDirectory(path).path
    if not exists(path): ValueError("File {0} does not exist.".format(path))
    with open(path, 'r') as file: string = file.read()
    hash = sha512(string).hexdigest()
    return self.outcars.exists(sha512=hash)
