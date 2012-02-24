def test():
  from pickle import loads, dumps
  from os import chdir, getcwd, remove
  from os.path import join, exists
  from shutil import rmtree
  from tempfile import mkdtemp
  from lada.vasp.incar._params import Restart
  from lada.vasp.files import WAVECAR, CHGCAR
  from restart_class import Extract, Vasp

  cwd = getcwd()
  directory = mkdtemp()
  try: 
    # no prior run.
    for v, istart, icharg in [(Vasp(True), None, 12), (Vasp(False), None, None)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(None).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(None)))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)

    # no prior run.
    for v, istart, icharg in [(Vasp(True), None, 12), (Vasp(False), None, None)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(Extract(directory, False)).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(Extract(directory, False))))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert r.value.directory == directory 
        assert r.value.success == False
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)
    
    # empty prior run.
    with open(join(directory, CHGCAR), "w") as file: pass
    with open(join(directory, WAVECAR), "w") as file: pass
    for v, istart, icharg in [(Vasp(True), 0, 12), (Vasp(False), 0, 2)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(Extract(directory, True)).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(Extract(directory, True))))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert r.value.directory == directory 
        assert r.value.success == True
        assert not exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)
    
    # prior run with charge only.
    with open(join(directory, CHGCAR), "w") as file: file.write('hello')
    for v, istart, icharg in [(Vasp(True), 0, 11), (Vasp(False), 0, 1)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(Extract(directory, True)).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))
        remove(join(playdir, CHGCAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(Extract(directory, True))))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert r.value.directory == directory 
        assert r.value.success == True
        assert exists(join(playdir, CHGCAR))
        assert not exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)

    # prior run with charge and wavecar.
    with open(join(directory, WAVECAR), "w") as file: file.write('hello')
    for v, istart, icharg in [(Vasp(True), 1, 11), (Vasp(False), 1, 1)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(Extract(directory, True)).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert exists(join(playdir, CHGCAR))
        assert exists(join(playdir, WAVECAR))
        remove(join(playdir, CHGCAR))
        remove(join(playdir, WAVECAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(Extract(directory, True))))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert r.value.directory == directory 
        assert r.value.success == True
        assert exists(join(playdir, CHGCAR))
        assert exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)

    # prior run with wavecar only.
    remove(join(directory, CHGCAR))
    with open(join(directory, CHGCAR), "w") as file: pass
    for v, istart, icharg in [(Vasp(True), 1, 10), (Vasp(False), 1, 0)]:
      playdir = mkdtemp()
      try: 
        chdir(playdir)
        v.istart, v.icharge = None, None
        assert Restart(Extract(directory, True)).incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert not exists(join(playdir, CHGCAR))
        assert exists(join(playdir, WAVECAR))
        remove(join(playdir, WAVECAR))

        v.istart, v.icharge = None, None
        r = loads(dumps(Restart(Extract(directory, True))))
        assert r.incar_string(vasp=v) is None
        assert v.istart == istart and v.icharg == icharg
        assert r.value.directory == directory 
        assert r.value.success == True
        assert not exists(join(playdir, CHGCAR))
        assert exists(join(playdir, WAVECAR))
      finally:
        chdir(cwd)
        rmtree(playdir)

  finally: rmtree(directory)

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

