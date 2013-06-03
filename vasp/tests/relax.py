def test():
  from pickle import loads, dumps
  from pylada.vasp import Vasp
  from pylada.error import ValueError
 
  a = Vasp()

  for key in a._input.keys():
    if key not in ['isif', 'nsw', 'ibrion', 'relaxation']: 
      del a._input[key]
  assert a.relaxation == 'static'
  assert len(a.output_map(vasp=a)) == 1
  assert a.output_map(vasp=a)['ibrion'] == str(-1)

  a.relaxation = 'cellshape'
  assert a.relaxation == 'cellshape'
  assert a.isif == 5
  assert a.nsw == 50
  assert a.ibrion == 2
  assert len(a.output_map(vasp=a)) == 3
  assert a.output_map(vasp=a)['ibrion'] == str(2)
  assert a.output_map(vasp=a)['isif'] == str(5)
  assert a.output_map(vasp=a)['nsw'] == str(50)
  a = loads(dumps(a))
  assert len(a.output_map(vasp=a)) == 3
  assert a.output_map(vasp=a)['ibrion'] == str(2)
  assert a.output_map(vasp=a)['isif'] == str(5)
  assert a.output_map(vasp=a)['nsw'] == str(50)

  a.relaxation = 'cellshape volume'
  a.nsw = 25
  assert a.relaxation == 'cellshape volume'
  assert a.isif == 6
  assert a.nsw == 25
  assert a.ibrion == 2
  assert len(a.output_map(vasp=a)) == 3
  assert a.output_map(vasp=a)['ibrion'] == str(2)
  assert a.output_map(vasp=a)['isif'] == str(6)
  assert a.output_map(vasp=a)['nsw'] == str(25)

  a.relaxation = 'ions'
  assert a.relaxation == 'ionic'
  assert a.isif == 2
  a.relaxation = 'ionic'
  assert a.relaxation == 'ionic'
  assert a.isif == 2

  a.relaxation = 'cellshape, volume ions'
  assert a.relaxation == 'cellshape ionic volume'
  assert a.isif == 3

  a.relaxation = 'cellshape, ionic'
  assert a.relaxation == 'cellshape ionic'
  assert a.isif == 4

  a.relaxation = 'volume'
  assert a.relaxation == 'volume'
  assert a.isif == 7

  a.relaxation = 'static'
  assert a.ibrion == -1
  assert a.nsw == 0
  assert a.isif == 2

  try: a.relaxation = 'ions, volume'
  except ValueError: pass
  else: raise Exception

if __name__ == "__main__": test()

