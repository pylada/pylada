def raw_input(*args): return 'y'

def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os import makedirs
  from os.path import exists, join
  from IPython.core.interactiveshell import InteractiveShell
  from pylada.jobfolder import JobFolder
  import pylada
  from dummy import functional
  import __builtin__ 
  __builtin__.raw_input = raw_input

  self = InteractiveShell.instance()

  root = JobFolder()
  for type, trial, size in [('this', 0, 10), ('this', 1, 15), ('that', 2, 20), ('that', 1, 20)]:
    jobfolder = root / type / str(trial)
    jobfolder.functional = functional
    jobfolder.params['indiv'] = size
    if type == 'that': jobfolder.params['value'] = True

  directory =  '/tmp/test' # mkdtemp() # 
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    self.user_ns['jobfolder'] = root
    self.magic("explore jobfolder")
    jobfolder = pylada.interactive.jobfolder
    assert 'this/0' in jobfolder and 'this/1' in jobfolder and 'that/2' in jobfolder and 'that/1'
    assert '0' in jobfolder['this'] and '1' in jobfolder['this']
    assert '1' in jobfolder['that'] and '2' in jobfolder['that']
    assert 'other' not in jobfolder
    for job in jobfolder.values():
      assert repr(job.functional) == repr(functional)
    assert getattr(jobfolder['this/0'], 'indiv', 0) == 10
    assert getattr(jobfolder['this/1'], 'indiv', 0) == 15
    assert getattr(jobfolder['that/1'], 'indiv', 0) == 20
    assert getattr(jobfolder['that/2'], 'indiv', 0) == 20
    assert not hasattr(jobfolder['this/0'], 'value')
    assert not hasattr(jobfolder['this/1'], 'value')
    assert getattr(jobfolder['that/1'], 'value', False)
    assert getattr(jobfolder['that/2'], 'value', False)
    assert pylada.interactive.jobfolder_path is None
    assert 'jobparams' in self.user_ns
    assert jobfolder is self.user_ns['jobparams'].jobfolder

    self.magic("savefolders {0}/dict".format(directory))
    pylada.interactive.jobfolder = None
    pylada.interactive.jobfolder_path = None
    self.magic("explore {0}/dict".format(directory))
    jobfolder = pylada.interactive.jobfolder
    assert 'this/0' in jobfolder and 'this/1' in jobfolder and 'that/2' in jobfolder and 'that/1'
    assert '0' in jobfolder['this'] and '1' in jobfolder['this']
    assert '1' in jobfolder['that'] and '2' in jobfolder['that']
    assert 'other' not in jobfolder
    for job in jobfolder.values():
      assert repr(job.functional) == repr(functional)
    assert getattr(jobfolder['this/0'], 'indiv', 0) == 10
    assert getattr(jobfolder['this/1'], 'indiv', 0) == 15
    assert getattr(jobfolder['that/1'], 'indiv', 0) == 20
    assert getattr(jobfolder['that/2'], 'indiv', 0) == 20
    assert not hasattr(jobfolder['this/0'], 'value')
    assert not hasattr(jobfolder['this/1'], 'value')
    assert getattr(jobfolder['that/1'], 'value', False)
    assert getattr(jobfolder['that/2'], 'value', False)
    assert pylada.interactive.jobfolder_path is not None
    assert 'jobparams' in self.user_ns
    assert jobfolder is self.user_ns['jobparams'].jobfolder
    assert jobfolder is self.user_ns['collect'].jobfolder

    for name, job in root.iteritems():
      if name == 'this/1': continue
      job.compute(outdir=join(directory, name))

    self.magic("explore results".format(directory))
    assert set(['/this/0/', '/that/1/', '/that/2/']) == set(self.user_ns['collect'].iterkeys())
    self.magic("explore errors".format(directory))
    assert len(self.user_ns['collect']) == 0
    self.magic("explore all".format(directory))
    self.magic("explore errors".format(directory))
    assert set(self.user_ns['collect'].keys()) == set(['/this/1/'])
    
  finally: 
    if directory != '/tmp/test': rmtree(directory)
    pass
  


if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1: path.extend(argv[1:])
  test()
