def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import join, exists
  from os import mkdir
  from lada.dftcrystal import Crystal, Functional, Shell
  from lada.dftcrystal.properties import Properties
  from lada import default_comm

  functional = Functional()
  functional.basis['Si'] = [
      Shell('s', a0=[16120.0, 0.001959],
                 a1=[2426.0, 0.01493], 
                 a2=[553.9, 0.07285],
                 a3=[156.3, 0.2461], 
                 a4=[50.07, 0.4859],
                 a5=[17.02, 0.325]),
      Shell('sp', a0=[292.7, -0.002781, 0.004438],
                  a1=[69.87, -0.03571, 0.03267],
                  a2=[22.34, -0.115, 0.1347],
                  a3=[8.15, 0.09356, 0.3287],
                  a4=[3.135, 0.603, 0.4496]), 
      Shell('sp', 4.0, a0=[1.22, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.55, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.27, 1.0, 1.0]) ]

  functional.dft.pbe0 = True
  functional.fmixing = 30
  functional.shrink = 8, 16
  functional.levshift = 5, True
  functional.maxcycle = 600
  functional.dft.spin = True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  directory =  mkdtemp() 
  if directory == '/tmp/test/':
    rmtree(directory)
    mkdir(directory)
  firstdir, seconddir = join(directory, '0'), join(directory, '1')
  try: 
     si = functional(crystal, outdir=firstdir, comm=default_comm)
     assert exists(join(firstdir, 'crystal.f9'))
     assert not exists(join(directory, 'crystal.f98'))
     properties = Properties(si)
     properties.rdfmwf = True
     result = properties(outdir=directory, workdir=join(directory, 'try0'))
     assert result.success
     assert exists(join(directory, 'crystal.f9'))
     assert not exists(join(directory, 'crystal.f98'))
     properties = Properties(directory)
     properties.fmwf = True
     result = properties(outdir=directory, workdir=join(directory, 'try1'))
     assert not exists(join(directory, 'crystal.f98'))
     result = properties(outdir=directory, workdir=join(directory, 'try1'), overwrite=True)
     assert result.success
     assert exists(join(directory, 'crystal.f98'))

  finally: 
    if directory != '/tmp/test/': rmtree(directory)
if __name__ == '__main__':
  test()