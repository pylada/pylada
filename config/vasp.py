""" VASP parameters for pylada. """
if "vasp" in globals()["pyladamodules"]:

  #vasp_program = "/projects/nrel/cid/bin/vasp"
  vasp_program = "/nopt/nrel/apps/vasp/5.2.2-bench/bin/vasp"

  """ Path of vasp binary executable (if launching as external program). """
  vasp_has_nlep = True
  """ Should be set to True if one wants to use NLEP. """

