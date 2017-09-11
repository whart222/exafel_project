# LIBTBX_SET_DISPATCHER_NAME xfel.map_height_at_atoms

from __future__ import division
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import sys

master_phil_str = """
map_type = anom
  .type = str
exclude_free_r_reflections = False
  .type = bool
fill_missing_f_obs = False
  .type = bool
resolution_factor = 0.25
  .type = float
selection = element CA or element ZN
  .type = atom_selection
"""

def master_phil () :
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string=master_phil_str,
    enable_automatic_twin_detection=False)

def run (args, out=sys.stdout) :
  usage_string = """\
xfel.map_height_at_atoms model.pdb data.mtz [map_type=MAP_TYPE] \\
    [selection=ATOM_SELECTION]
"""
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=False,
    prefer_anomalous=True,
    usage_string=usage_string,
    #set_wavelength_from_model_header=True,
    #set_inelastic_form_factors="sasaki",
    out=out)
  params = cmdline.params
  fmodel = cmdline.fmodel
  xray_structure = fmodel.xray_structure
  pdb_hierarchy = cmdline.pdb_hierarchy
  sel_cache = pdb_hierarchy.atom_selection_cache()
  selection = sel_cache.selection(params.selection).iselection()
  if (len(selection) == 0) :
    raise Sorry("No atoms selected!")
  params = cmdline.params
  map_coeffs = fmodel.map_coefficients(
    map_type=params.map_type,
    exclude_free_r_reflections=params.exclude_free_r_reflections,
    fill_missing=params.fill_missing_f_obs)
  fft_map = map_coeffs.fft_map(
    resolution_factor=params.resolution_factor).apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  make_sub_header("Map analysis", out=out)
  from scitbx.array_family import flex
  print >> out, "Maximum grid point value: %6.2f sigma" % \
    flex.max(real_map.as_1d())
  print >> out, ""
  for i_seq in selection :
    sc = xray_structure.scatterers()[i_seq]
    map_value = real_map.tricubic_interpolation(sc.site)
    print >> out, "%s : %6.2f sigma" % (sc.label, map_value)

if (__name__ == "__main__") :
  run(sys.argv[1:])
