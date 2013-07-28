""" Runs one job from the jobfolder. """
def main():
  import re 
  from sys import path as python_path
  from os.path import exists
  from argparse import ArgumentParser
  from pylada import jobfolder
  from pylada.process.mpi import create_global_comm
  import pylada

  # below would go additional imports.

  parser = ArgumentParser( prog="runone", description = re.sub("\\s+", " ", __doc__[1:]))

  parser.add_argument('--bugLev', dest="bugLev", default=0, type=int,\
                      help="Debug level.")
  parser.add_argument( "--jobid", dest="names", nargs='+', type=str, \
                       help="Job name", metavar="N" )
  parser.add_argument( "--ppath", dest="ppath", default=None, \
                       help="Directory to add to python path",
                       metavar="Directory" )
  parser.add_argument('--nbprocs', dest="nbprocs", default=pylada.default_comm['n'], type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--ppn', dest="ppn", default=pylada.default_comm['ppn'], type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--timeout', dest="timeout", default=300, type=int,\
                      help="Time to wait for job-dictionary to becom available "
                           "before timing out (in seconds). A negative or null "
                           "value implies forever. Defaults to 5mn.")
  parser.add_argument('pickle', metavar='FILE', type=str, help='Path to a job-folder.')

  try: options = parser.parse_args()
  except SystemExit: return

  from pylada.misc import setBugLev
  setBugLev( options.bugLev)         # set global debug level
  from pylada.misc import bugLev     # must import after calling setBugLev

  # additional path to look into.
  if options.ppath is not None: python_path.append(options.ppath)

  if not exists(options.pickle): 
    print "Could not find file {0}.".format(options.pickle)
    return

  # Set up mpi processes.
  pylada.default_comm['ppn'] = options.ppn
  pylada.default_comm['n'] = options.nbprocs
  create_global_comm(options.nbprocs)

  timeout = None if options.timeout <= 0 else options.timeout
  
  jobfolder = jobfolder.load(options.pickle, timeout=timeout)
  print '  ipy/lau/scattered_script: jobfolder: ', jobfolder
  print '  ipy/lau/scattered_script: options: ', options
  for name in options.names:
    if bugLev >= 1:
      print '  ipy/lau/scattered_script: name: %s' % ( name,)
      print '  ipy/lau/scattered_script: jobfolder[name]: %s' \
        % ( jobfolder[name],)
      print '  ipy/lau/scattered_script: type(jobfolder[name]): %s' \
        % ( type(jobfolder[name]),)
      print '  ipy/lau/scattered_script: jobfolder[name].compute: %s' \
        % ( jobfolder[name].compute,)
      print '  ipy/lau/scattered_script: type(jobfolder[name].compute): %s' \
        % ( type(jobfolder[name].compute),)
      print '  ipy/lau/scattered_script: before compute for name: %s' \
        % ( name,)

    jobfolder[name].compute(comm=pylada.default_comm, outdir=name)
    if bugLev >= 1:
      print '  ipy/lau/scattered_script: after compute for name: %s' \
        % ( name,)

if __name__ == "__main__": main()
