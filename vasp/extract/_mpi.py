""" Decorator to make any extraction method mpi-aware. """
from decorator import decorator

@decorator
def _bound_mpi_extraction(method, *args, **kwargs): 
  """ Reads from root, and broadcasts to all 
  
      Decorator for bound methods.
  """
  assert len(args) >= 1, "Decorator for bound methods only.\n"
  return _cache_method_result(0, method, *args, **kwargs)
  
@decorator
def _attribute_mpi_extraction(method, *args, **kwargs): 
  """ Reads from root, and broadcasts to all 
  
      Decorator for properties.
  """
  assert len(args) >= 2, "Decorator for bound methods only.\n"
  return _cache_method_result(1, method, *args, **kwargs)

def _cache_method_result(which, method, *args, **kwargs):
  """ Saves grep result for later reuse """
  cache_name = "_cached_attr%s" % (method.__name__)
  if not hasattr(args[which], cache_name): 
    setattr(args[which], cache_name, _extract_which_is_self(which, method, *args, **kwargs))
  return getattr(args[which], cache_name)

def _uncache(ob):
  for key in ob.__dict__.keys():
    if key[:len("_cached_attr")] == "_cached_attr": del ob.__dict__[key]

def _extract_which_is_self(which, method, *args, **kwargs):
  from boost.mpi import broadcast
  
  mpicomm = args[which].mpicomm
  if mpicomm == None: return method(*args, **kwargs)
  if mpicomm.size == which: return method(*args, **kwargs)
  if mpicomm.rank != 0:
    result = broadcast(mpicomm, root = 0)
    assert broadcast(mpicomm, root = 0) == "Am in sync", "Processes not in sync"
    return result

  args[which].mpicomm = None
  result = method(*args, **kwargs)
  args[which].mpicomm = mpicomm
  assert mpicomm != None
  assert args[which].mpicomm != None
  broadcast(mpicomm, result, root = 0)
  assert broadcast(mpicomm, "Am in sync", root = 0) == "Am in sync", "Processes not in sync"
  return result
