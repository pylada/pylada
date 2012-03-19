""" Classes to manipulate job-dictionaries. """
__docformat__ = "restructuredtext en"
__all__ = ['JobParams']

from .extract import AbstractMassExtract

class JobParams(AbstractMassExtract):
  """ Get and sets job parameters for a job-dictionary. """
  def __init__(self, jobdict=None, only_existing=True, **kwargs):
    """ Initializes job-parameters.

        :param jobdict :
            :py:class:`JobDict <lada.jobs.jobdict.JobDict>` instance for which
            to get/set parameters. If None, will look for ipython's
            current_jobdict.
        :param bool only_existing:
            If true (default), then only existing parameters can be modified:
            non-existing parameters will not be added.
        :param kwargs:
            Variable length keyword argument passed on to
            :py:meth:`AbstractMassExtract.__init__`.
    """
    from .. import jobparams_only_existing

    self.__dict__["_jobdict"] = jobdict
    self.__dict__['only_existing'] = jobparams_only_existing\
                                       if only_existing is None\
                                       else only_existing
    """ Only modifies parameter which already exist. """

    if 'dynamic' not in kwargs: kwargs['dynamic'] = True
    super(JobParams, self).__init__(**kwargs)


  @property
  def addattr(self):
    """ Returns manipulator with ability to *add new* attributes. """
    return self.copy(only_existing=False)
    
  @property
  def jobdict(self):
    """ Jobdictionary for which to get/set parameters. """
    if self._jobdict is None: 
      from lada import is_interactive
      try: from lada import current_jobdict
      except ImportError:
        if is_interactive: print "No current job-dictionary."
        else: raise RuntimeError("No current job-dictionary.")
      return current_jobdict.root
    return self._jobdict.root
  @jobdict.setter
  def jobdict(self, value): self._jobdict = value
  @jobdict.deleter
  def jobdict(self, value): self._jobdict = None

  @property
  def onoff(self):
    """ Dictionary with calculations which will run.

	Whereas other properties only report untagged jobs, this will report
        both. Effectively checks wether a job is tagged or not. Calculations which 
    """
    result = {}
    for name, job in self.iteritems():
      result[name] = "off" if job.is_tagged else "on"
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return result

  @onoff.setter
  def onoff(self, value):
    """ Dictionary with tagged and untagged jobs.

        Whereas other properties only report untagged jobs, this will report
        both.
    """
    if hasattr(value, 'iteritems'):
      for key, value in value.iteritems():
        try: job = self[key]
        except: continue
        else: job.onoff = value
    elif value == "on" or value == True:
      for name, job in self.iteritems(): job.untag()
    elif value == "off" or value == False:
      for name, job in self.iteritems(): job.tag()

  @property
  def extractors(self):
    """ Returns dictionary of extrators. """
    from .forwarding_dict import ForwardingDict
    result = self.dicttype()
    for k, j in self.iteritems(): result[k] = j
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return ForwardingDict( dictionary=result, naked_end=self.naked_end, \
                           only_existing=self.only_existing, readonly=False)
    

  @property
  def view(self):
    """ A regex pattern which the name of extracted jobs should match.

        If None, then no match required. Should be a string, not an re object.
    """
    return self.jobdict.name if self._view is None else self._view
  @view.setter
  def view(self, value): self._view = value

  def __iter_alljobs__(self):
    """ Loops through all correct jobs. """
    for name, job in self.jobdict.iteritems(): yield job.name, job

  def __getattr__(self, name): 
    """ Returns extracted values. """
    from .forwarding_dict import ForwardingDict
    result = self.dicttype()
    for key, value in self.iteritems():
      if value.is_tagged: continue
      try: result[key] = getattr(value, name)
      except: result.pop(key, None)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    if len(result) == 0: 
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )
    return ForwardingDict( dictionary=result, naked_end=self.naked_end, \
                           only_existing=self.only_existing, readonly=False)

  def __setattr__(self, name, value):
    """ Returns dictionary with job parameters for each job. """
    from re import match
    # initialization not done yet.
    if "only_existing" not in self.__dict__: super(JobParams, self).__setattr__(name, value)
    # some cached attributes.
    if match("_cached_attr\S+", name): super(JobParams, self).__setattr__(name, value)
    # Look for other attriubtes in current instance.
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: pass
    else:
      super(JobParams, self).__setattr__(name, value)
      return 

    found = False
    for jobname, job in self.iteritems():
      if job.is_tagged: continue
      if hasattr(job, name):
        setattr(job, name, value)
        found = True
      elif not self.only_existing: 
        job.jobparams[name] = value
        found = True
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  def __delattr__(self, name):
    try: super(JobParams, self).__getattribute__(name)
    except AttributeError: pass
    else:
      super(JobParams, self).__delattr__(name, value)
      return

    found = False
    for jobname, job in self.iteritems():
      if job.is_tagged: continue
      if hasattr(job, name):
        delattr(job, name)
        found = True
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  @property
  def _attributes(self):
    """ Attributes which already exist. """
    result = set()
    for name, job in self.iteritems():
      if not job.is_tagged: result |= set([u for u in dir(job) if u[0] != '_'])
    return result
 
