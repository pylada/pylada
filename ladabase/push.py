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

""" Push/Pulls from a remote computer. """
from contextlib import contextmanager

def _getcomment(self=None, cmdl=None):
  """ Gets comment from user. """
  from sys import stderr
  from os import remove
  import re
  from tempfile import NamedTemporaryFile
  from IPython.ipapi import TryNext, get as get_ipy
  try: from .. import fullname
  except ImportError:
    print >>stderr, "Could not import fullname with which to tag files in database.\n"\
                    "Please add `fullname = 'my full name'` in ~/.pylada.\n"
    return
  if len(fullname) == 0:
    print >>stderr, "Username with which to tag files in database is empty.\n"\
                    "Please add `fullname = 'my full name'` in ~/.pylada.\n"
    return

  if self == None: self = get_ipy()

  with NamedTemporaryFile(mode='w', delete=False) as file: 
    file.write("\n# operator: {0}\n".format(fullname))
    if cmdl is not None: file.write("# command-line: {0}\n".format(cmdl))
    filename = file.name
  try: self.shell.hooks.editor(filename, 0)
  except TryNext:
    print "Could not open editor."
    return
  with open(filename, 'r') as file: comment = file.read()
  remove(filename)
  stripped = re.sub('#.*(?:\n|$)', '', comment, re.M)
  stripped = stripped.replace('\n','').replace(' ', '')
  if len(stripped) == 0: 
    print "Empty comment. Aborting."
    comment = None
  return comment

def _get_push_parser():
  """ Creates and returns CLI parser for push. """
  import argparse
  import textwrap
  
  description ="Push OUTCAR and directories of OUTCAR to database."
  epilog = textwrap.wrap(textwrap.dedent(\
           """\
                 This function will pull remote files and push them to the database. It will
              also push local files to the database.  To do so 'fullname' and
              'pymongo_username' should be declared in your ~/.pylada file. The first shoudl
              be your full name and will be attached within the database to files you have
              pushed. The second is you NREL username. Your ~/.pylada file should read
              something like:
          """)) + [ "", "fullname = \"John Doe\"", "pymongo_username = \"jdoe\"", "" ] \
                + textwrap.wrap(textwrap.dedent(\
          """
                Before pushing files to the database, you should know what kind of computation
             you have been doing. Currently, only FERE calculations are accepted. These
             should conform to a given standard such that FERE enthalpies can be computed.
             To push local files one can simply run the following line: 
         """)) + ["", "%push fere directoryA directoryB directoryC/OUTCAR-fere_*", "" ] \
         + textwrap.wrap(textwrap.dedent(\
         """\
            This will look into directoryA and  directoryB for any file containing OUTCAR
            in its name, as well as any file in directoryC starting with OUTCAR-fere_. To
            push files from a remote computer, add the following option:
         """)) + ["", "%push --hostname jcdoe@greenmesa.sandia.gov fere directoryA\n\n", ""] \
         + textwrap.wrap(textwrap.dedent(\
        """\
           In this case, OUTCAR files from directoryA on the greenmesa computer will be
           first transfered to the local computer, and then pushed to the database.
        """))
           
  parser = argparse.ArgumentParser( prog='%push', description=description, epilog='\n'.join(epilog),\
                                    formatter_class=argparse.RawDescriptionHelpFormatter )
  parser.add_argument( 'algo', choices=['fere', 'gw'],
                       help='Type of data to push.' )
  parser.add_argument( '--hostname', type=str, default=None, 
                       help='Remote computer from which to pull. '\
                            'Should be given as username@hostname. '\
                            'Ignore if pushing local files. ')
  parser.add_argument( '--norecurrence', action="store_true", help='Do not walk into subdirectory.')
  parser.add_argument( '--pattern', type=str, nargs=1, default="*OUTCAR*",
                       help='When searching within directories, examine only files '\
                            'fitting this pattern. Defaults to "*OUTCAR*". Can include bash wildcards.')
  parser.add_argument( '--exclude', type=str, nargs='*', help="Files and directories matching "\
                       "these patterns are excluded from the search. 'relax_cellshape' and "\
                       "'relax_ions' are always excluded. Can include bash wildcards. "\
                       "More than one --exclude argument can be given. ")
  parser.add_argument( 'directories', metavar='DIR/OUTCAR', type=str, nargs='*',
                       help='OUTCAR file, or root directory.' )
  return parser

@contextmanager
def repatriate_file(path, sftp):
  """ Repatriates a file from a remote computer.
  
      This context returns None if this is not a FERE calculation. If it is a
      FERE calculation, it returns a valid extraction object pointing to a
      local file. If sftp is not None, a local file is first repatriated, and
      then destroyed on leaving the context.
  """
  from tempfile import NamedTemporaryFile
  from os import remove

  localfile = NamedTemporaryFile(delete=False)
  localfile.close()
  try: 
    sftp.get(path, localfile.name)
    yield localfile.name
  except: yield None
  remove(localfile.name)

def local_iglob(path):
  """ Globs local files. """
  from glob import iglob
  from os import stat
  from stat import S_ISDIR
  from ..opt import RelativeDirectory
  for result in iglob(RelativeDirectory(path).path):
    yield result, S_ISDIR(stat(result).st_mode)

def remote_iglob(path, ssh, sftp):
  """ Globs remote files. """
  from os.path import dirname, basename, join
  from stat import S_ISDIR
  from fnmatch import filter as fnfilter
  dir, suffix = dirname(path), basename(path)
  for dir in ssh.exec_command("ls -d " + dir)[1].read()[:-1].split('\n'):
    paths = [join(dir, u) for u in sftp.listdir(dir)]
    if len(suffix) != 0:
      if '~' not in path: paths = fnfilter(paths, path)
      else: paths = fnfilter(paths, path.replace('~/', '*/'))
    for path in paths:
      try: attr = sftp.stat(path)
      except: continue
      yield path, S_ISDIR(attr.st_mode)

def remote_walk(path, sftp):
  """ Walk over remote files and directories. """
  from os.path import join
  from stat import S_ISDIR
  all = sftp.listdir(path), sftp.listdir_attr(path)
  dirs, files = [], []
  for u, a in zip(*all):
    if S_ISDIR(a.st_mode): dirs.append(u)
    else: files.append(u)
  yield path, dirs, files
  for dir in dirs:
    for dummy in remote_walk(join(path, dir), sftp): yield dummy

def walk_calc_files(args, context, iglob, walk):
  """ Walks over acceptable result files on a remote computer. """
  from os.path import basename, join, dirname, relpath
  from fnmatch import translate as fntranslate
  from re import compile, search

  excludes = [compile(fntranslate('relax_cellshape')), compile(fntranslate('relax_ions'))]
  if args.exclude is not None:
    for i in args.exclude: excludes.append(compile(fntranslate(i)))
  pattern = compile(fntranslate(args.pattern))

  for input in args.directories:
    for path, isdir in iglob(input):
      if any(search(i, path) is not None for i in excludes): continue
      if isdir: 
        for root, dirs, files in walk(path):
          if args.norecurrence: dirs[:] = []
          else: dirs[:] = [d for d in dirs if all(search(i, d) is None for i in excludes)]
          for file in files:
            if search(pattern, basename(file)) == None: continue
            with context(join(root, file)) as result:
              if result is not None: yield result, relpath(join(root, file), dirname(dirname(path)))
        continue
      with context(path) as result:
        if result is not None: yield result, relpath(path, dirname(path))

def push(self, cmdl):
  """ Pulls files from a remote directory. """
  from getpass import getpass
  from paramiko import SSHClient, AutoAddPolicy
  from pylada.ladabase import Manager
  from pylada.ladabase.extracted import generate_extracted
  from hashlib import sha512
  try: from .. import fullname
  except ImportError:
    print "Could not import fullname with which to tag files in database.\n"\
          "Please add `fullname = 'my full name'` in ~/.pylada.\n"
    return
  if len(fullname) == 0:
    print "Username with which to tag files in database is empty.\n"\
          "Please add `fullname = 'my full name'` in ~/.pylada.\n"
    return

  if __name__ == "__main__":
    from sys import argv
    cmdl = " ".join(argv)
    try: args = _get_push_parser().parse_args()
    except SystemExit: return None
  else: 
    try: args = _get_push_parser().parse_args(cmdl.split())
    except SystemExit: return None

  if args.algo == "gw": 
    print "Pushing GW data is no yet implemented."
    return
 
  # gets comment. 
  comment = _getcomment(self, cmdl)
  if comment is None: return

  # try connecting to host if requested.
  if getattr(args, 'hostname', None) is not None:
    client = SSHClient()
    client.set_missing_host_key_policy(AutoAddPolicy())
    try: username, hostname = args.hostname.split('@')
    except: 
      print "Remote hostname should be given an username@hostname."
      return 
    found = False
    for i in range(3):
      client.password = getpass("Pass for {0}@{1}:".format(username, hostname))
      try: client.connect(hostname, username=username, password=client.password, timeout=5)
      except KeyboardInterrupt:
        print "Aborting."
        return
      except Exception as e: print 'error', e; continue
      else: found = True; break
    if not found: 
      print "Could not connect to {0}".format(args.remote)
      return
    # sets up environment for iterating over files.
    client_sftp = client.open_sftp()
    iglob = lambda x: remote_iglob(x, client, client_sftp)
    walk = lambda x: remote_walk(x, client_sftp)
    def context(other): 
      @contextmanager
      def _context(path):
        with repatriate_file(path, client_sftp) as filepath:
          with other(filepath) as result: yield result
      return _context

  # otherwise, look for local files.
  else: 
    # sets up environment for iterating over files.
    from os import walk as local_walk
    iglob = local_iglob
    walk = local_walk
    context = lambda x: x



  # Now performs work.
  manager = Manager()
  outcardb = manager.files
  if args.algo == "fere": 
    from pylada.ladabase.fere import check_fere_context, generate_fere_summary
    found = False
    for extract, path in walk_calc_files(args, context(check_fere_context), iglob, walk):
      hash = sha512(extract.__outcar__().read()).hexdigest()
      if outcardb.find_one({'sha512': hash}) != None: 
        print path, "is already in the database."
        continue
      with extract.__outcar__() as file: outcar = file.read()
      item = manager.push( path, outcar, comment, compression="bz2",\
                           is_dft=extract.is_dft, is_gw=extract.is_gw, uploader=username )
      found = True
      print "Pushing", path, "."
      try: generate_extracted(filter={'_id': item})
      except:
        print "Could not extract values from ", path, "."
        pass
    if not found:
      print "No new OUTCAR found. "
      return
    generate_fere_summary(2)
  


if __name__ == "__main__":
  push(None, None)
