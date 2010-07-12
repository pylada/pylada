#! /usr/bin/python
""" Module with different methods for vasp 

    Most of the generators in this module will not recompute the result if an
    appropriate output directory is found. In other words, this function can be
    called once to perform calculation, and then to perform output extractions.

    >>> # First create method.
    >>> some_method = SomeMethod(vasp)
    >>> # Then call the method.
    >>> some_method(structure)
    >>> # do something else
    >>>
    >>> # Then  analyze each step.
    >>> for extract in some_method.generator(structure):
    >>>   print "Where: ", extract.directory
    >>>   print extract.total_energy

    The output extraction object will be the output the vasp callable.
"""
__docformat__ = "restructuredtext en"

class RelaxCellShape(object):
  """ Functor for cell-shape relaxation.
  
      Since vasp is a plane-wave code, cell-relaxation are never quite accurate.
      This functional keeps working until convergence is achieved, and then
      creates a static calculation.
  """
  def __init__(self, vasp, relaxation="volume ionic cellshape", first_trial=None, maxiter=0):
    """ Initializes a cell-shape relaxation. 
        
        :Parameters:
          vasp : lada.vasp.Vasp
            The functional with which to perform cell relaxation.
          relaxation
            Degrees of freedom to relax. Can be:
            
            - volume
            - ionic
            - cellshape

            Note that using this method with only *ionic* relaxation is a waste
            of computer ressources.
          first_trial : dict or None
            Does nothing if *None*.  A dictionary of parameters to apply for
            the first run of the functional. It may be interesting to first get a
            qualitative relaxation with less computationally intensive
            parameters of the vasp functional. 
          maxiter : int
            Maximum number of iterations before bailing out. Zero or negative
            number means infinit iterations.
    """
    super(RelaxCellShape, self).__init__()
    self.vasp = vasp
    """ The functional with which to perform cell relaxation. see `__init__`. """
    self.relaxation = relaxation
    """ Degrees of freedom to relax. see `__init__`. """
    self.first_trial = first_trial
    """ Dictionary of parameters for the first run of the functional. see `__init__`. """
    self.maxiter = maxiter
    """ Maximum number of iterations before bailing out. """

  def generator(self, structure, outdir=None, comm=None, **kwargs ):
    """ Performs a vasp relaxation, yielding each result.
    
        The result from the `vasp` calculations are yielded at each step. This
        makes it easy to analyse results during or after the run.

        >>> relaxor = RelaxCellShape(vasp)
        >>> for output in relaxor.generator(structure): 
        >>>   print output.energy

        :Parameters:
          structure
            The structure to relax with `vasp`.
          outdir
            Output directory passed on to the `vasp` functional.
          comm : boost.mpi.communicator or None
            MPI communicator passed on to the `vasp` functional.
          kwargs 
            Other keywords will overide attributes of this instance of
            `RelaxCellShape` (though for this run only. This function is
            stateless) if they are named after attributes of `RelaxCellShape`.
            Otherwise, the keywords are passed on to the `vasp` functional.
    """
    from copy import deepcopy
    from math import fabs 
    from os import getcwd
    from os.path import join, exists

    # does not run code. Just creates directory.
    if kwargs.pop("norun", False): 
      yield self._norun(structure, outdir, comm, **kwargs)
      return

    # make this function stateless.
    vasp = kwargs.pop("vasp", self.vasp)
    structure = deepcopy(structure)
    set_relaxation = kwargs.pop("relaxation", self.relaxation)
    set_relaxation = kwargs.pop("set_relaxation", set_relaxation)
    first_trial = kwargs.pop("first_trial", self.first_trial)
    maxiter = kwargs.pop("maxiter", self.maxiter)
    outdir = kwargs.pop("outdir", getcwd())


    # number of restarts.
    nb_steps, olde = 0, None
   
    # sets parameter dictionary for first trial.
    if first_trial != None and nb_steps == 0:
      params = kwargs.copy()
      params.update(first_trial)
    else: param = kwargs
    
    # performs relaxation calculations.
    while maxiter <= 0 or nbsteps < maxiter:
      # performs initial calculation.
      output = vasp\
               (\
                 structure, \
                 outdir = join(outdir, "relax_cellshape_" + str(nb_steps)),\
                 comm=comm,\
                 set_relaxation = set_relaxation,\
                 **params
               )
      yield output
      
      nb_steps += 1
      params = kwargs
      olde = output.energy
      if nb_steps == 1: continue
      if fabs(output.energy - olde) < float(len(structure.atoms)) * self.vasp.ediff: break

    # Does not perform static calculation if convergence not reached.
    if fabs(output.energy - olde) > float(len(structure.atoms)) * self.vasp.ediff: 
      yield output 
    # performs final calculation outside relaxation directory. 
    output = vasp\
             (\
               structure, \
               outdir = outdir,\
               comm=comm,\
               set_relaxation = "static",\
               **kwargs\
             )
    yield output

  def __call__(*args, **kwargs):
    """ Performs a vasp relaxation. 

        Arguments are the same as `generator`.

        The plane wave basis depends upon cell-shape. Hence, the basis goes out
        of sync with the actual structure during cell-shape relaxation.
        Convergence can only be achieved by restarting the calculation. And
        eventually performing a static calculation.

        It is wastefull to use this functional unless a cell-shape relaxation
        is wanted.

        On top of the usual vasp parameter, this functional accepts the keyword
        agument *relaxation* as a shorthand for *set_relaxation* attribute of
        the vasp functional. If both are passed to this method,
        *set_relaxation* takes precedence.

        If you want to examine the result of each and every vasp calculation,
        use `generator` instead.
    """ 
    if kwargs.pop("norun", False): 
      return self._norun(structure, outdir, comm, **kwargs)

    for output in self.generator(*args, **kwargs): pass
    return output
  
  def _norun(structure, outdir, comm, **kwargs):
    """ Just creates directory for debugging. """
    from ..opt.changedir import Changedir

    is_root = True if comm == None else comm.rank == 0
    if is_root:
      # creates a file describing the relaxation parameters.
      with Changedir(outdir) as pwd:
        with open("relaxe_cell_shape_parameters") as file:
          file.write( "self.relaxation = %s\n" % (repr(self.relaxation)) )
          file.write( "self.maxiter = %s\n" % (repr(self.maxiter)) )
          file.write( "self.first_trial = %s\n" % (repr(self.first_trial)) )

    # Now have vasp do a fake run to create anything it does create.
    return vasp(structure, outdir, comm, norun=True, **kwargs)
