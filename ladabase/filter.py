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

""" IPython interface to filter through the database. """

def filters(self, cmdl):
  """ Provides a stack of filters through which to filter the database. """
  import argparse
  import pylada
  from .massextract import MassExtract

  if not hasattr(pylada, "filters_stack"):
    pylada.filters_stack = [("True", MassExtract(namespace=ip.user_ns, mongofilter={'groundstates': {'$exists:1'}}))]
  stack = pylada.filters_stack

  # empty commandline. Print filters.
  if len(cmdl.rstrip()) == 0: cmdl = "--view"

  parser = argparse.ArgumentParser(prog='%filter',
                     formatter_class=argparse.RawDescriptionHelpFormatter,
                     epilog = \
"""
Creates stack of filters to navigate database. Only those elements satisfying
all filters can be analyzed with "collect". Each call to this magic function
will add a further filter. Filters are specified as below:

  >>> %filters 'Ba' in extract.species

Where *extract* is an object used to access job-parameters from the outcar.
To view which parameters are accessible, type:

  >>> collect.[TAB]

including the \".\" (dot) and hit TAB.
""" )

  group = parser.add_mutually_exclusive_group()
  group0 = group.add_argument_group()
  group0.add_argument( 'exprs', type=str, nargs="*", default="",
                        help='Python expression returning True or False.')
  group0.add_argument( '--mongofilter', action="store_true",
                       dest="mongofilter", 
                       help='Mongo filter (fast) defining documents over which '  \
                            'further filters (slow) trawl. Changing this value '  \
                            'resets the stack of (slow) filters. It should be a ' \
                            'dictionary (or a string defining a dictionary) '     \
                            'honoring the pymongo syntax. Google "mongodb '       \
                            'query". Most values available with collect are '     \
                            'also available as direct attributes of the '         \
                            'database, eg: {"species": {"%all": ["Ni", "In", '    \
                            '"O"]}}. Note that "$" has a special significance '   \
                            'for both mongodb and ipython. As such, they must be '\
                            'replaced with "%" signs. The above is a valid query '\
                            'which will match all systems containing only the '   \
                            'relevant element (including binaries). ')
  group.add_argument( '--iterate', action="store_true", 
                       dest="iterate", help='Iterates over current job.' )
  group.add_argument( '--view', action="store_true", 
                       dest="view", help='Prints all filters.' )
  group.add_argument( '--pop', action="store_true", 
                       dest="pop", help='Remove last filter.' )
  group.add_argument( '--outcar', action="store_true", 
                       dest="outcar", help='View current OUTCAR, if there is only one.' )



  try: args = parser.parse_args(cmdl.replace("%", "$").split())
  except SystemExit as e: return None

  if args.view: 
    print "Filters:\n{0:>2}. {1}".format(0, repr(stack[0][1].mongofilter))
    for i, (filter, value) in enumerate(stack[1:]):
      if i == 0: continue
      else: print "{0:>2}. {1}".format(i, filter)
  elif args.iterate: 
    if not hasattr(pylada, "iterstack"): pylada.iterstack = 0
    else: stack.pop(-1)
    stopiter = True
    for i, key in enumerate(stack[-1][1]):
      if i == pylada.iterstack: stopiter = False; break
    if stopiter == True:
      print "Reached last filtered element."
      key = stack[-1][1].keys()[-1]
    else: pylada.iterstack += 1
    stack.append(("id == {0}".format(repr(key)), stack[-1][1]["id == {0}".format(repr(key))]))
  elif args.pop and len(stack) > 1:
    if hasattr(pylada, "iterstack"): del pylada.iterstack
    stack.pop(-1)
  elif args.outcar:
    if len(stack[-1][1]) != 1:
      print "outcar only valid if database is filtered to a single outcar."
    else: 
      from tempfile import NamedTemporaryFile
      try: 
        filename = None
        with NamedTemporaryFile("w", delete=False) as file:
          filename = file.name
          with stack[-1][1].values()[0].raw.__outcar__() as foutcar:
            file.write(foutcar.read())
        self.api.system("less {0}\n".format(filename))
      except: raise
      finally:
        from os import remove
        if filename is not None:
         try: remove(filename)
         except: pass
  elif args.mongofilter: 
    mongofilter = eval(" ".join(args.exprs).replace("%%", "$"), self.api.user_ns)
    dictionary = self.api.user_ns.copy()
    dictionary.pop("all", None)
    dictionary.pop("$all", None)
    stack = [["True", MassExtract(namespace=self.api.user_ns, mongofilter=mongofilter)]]
  else: 
    if hasattr(pylada, "iterstack"): del pylada.iterstack
    cmdl = " ".join(args.exprs)
    try:
      add = stack[-1][1][cmdl]
      print "Found {0} calculations statisfying the filters.".format(len(add))
    except: 
      print "Could not make sense of expression {0}.".format(cmdl)
      if len(stack[-1][1]): 
        self.api.user_ns["extract"] = stack[-1][1].values()[0]
        print "I have put an \"extract\" variable in your namespace for you to try it out."
      else: print "But in any case, there is no jobs to filter anymore."
    else:
      if len(add) == 0: 
        print "This filter has not been added to the list."
      else: stack.append((cmdl, add))

  pylada.filters_stack = stack
  self.api.user_ns["collect"] = stack[-1][1]

def completer(self, event):
  """ Completer for filter magic function. """
  return ["--iterate", "--view", "--pop", "--outcar", '--mongofilter']

def init(ip):
  """ Initializes filters stuff. """
  from .massextract import MassExtract
  import pylada
  ip.expose_magic("filters", filters)
  ip.set_hook('complete_command', completer, re_key = '\s*%?filter')
  if not hasattr(pylada, "filters_stack"):
    pylada.filters_stack = [("True", MassExtract(namespace=ip.user_ns, mongofilter={'groundstate': {'$exists': 1}}))]
  ip.user_ns["collect"] = pylada.filters_stack[-1][1]

