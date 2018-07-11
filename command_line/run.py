# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit2
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT


from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
  
from iotbx import file_reader
import iotbx.pdb
import libtbx
from libtbx import phil
from libtbx.utils import Sorry
from mmtbx.dynamics import simulated_annealing as sa
import mmtbx.model
import mmtbx.utils
import os
from scitbx.array_family import flex
import sys

#master_params_str are used for default values of options in GUI
master_params_str = """
cryo_fit2 {
include scope libtbx.phil.interface.tracking_params
Input{
  model_file_name = None
    .type = path
    .short_caption = Starting model file 
    .multiple = False
    .help = map_to_model derived model / homology model / model from different organism/experimental method (available format: .cif/.pdb)
    .style = file_type:pdb bold input_file
  map_file_name = None
    .type = path
    .short_caption = Target map file
    .help = Cryo-EM map file (available format: .ccp4/.map)
    .style = bold input_file
}
Output
{
  output_file_name_prefix = None
    .type = str
    .short_caption = Output prefix
    .help = Prefix for output filename
}
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
}
"""
master_params = master_params_str
master_phil = phil.parse(master_params_str, process_includes=True)
# This sentence works before main function

def calculate_cc(map_data, model, resolution):
  xrs = model.get_xray_structure()
  fc = xrs.structure_factors(d_min = resolution).f_calc()
  f_map = fc.structure_factors_from_map(
    map            = map_data,
    use_scale      = True,
    anomalous_flag = False,
    use_sg         = False)
  return fc.map_correlation(other = f_map)
# end of calculate_cc fn


def clean_pdb_for_phenix(input_pdb_file_name):
  if "_RNA_cleaned.pdb" in input_pdb_file_name:
    return input_pdb_file_name
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_RNA_cleaned.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    if line[18:20] == "RA":
      new_line = line[:18] + " A" + line[20:]
      f_out.write(new_line)
    elif line[18:20] == "RU":
      new_line = line[:18] + " U" + line[20:]
      f_out.write(new_line)
    elif line[18:20] == "RG":
      new_line = line[:18] + " G" + line[20:]
      f_out.write(new_line)
    elif line[18:20] == "RC":
      new_line = line[:18] + " C" + line[20:]
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  return output_pdb_file_name
#end of clean_pdb_for_phenix fn


def write_ccp4_map(crystal_symmetry, file_name, map_data):
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=crystal_symmetry.unit_cell(),
      space_group=crystal_symmetry.space_group(),
      map_data=map_data,
      labels=flex.std_string([""]))
# end of write_ccp4_map fn

def validate_params(params): # validation for GUI
  if (params.cryo_fit2.Input.model_file_name is None):
    raise Sorry("Model file should be given")
  if (params.cryo_fit2.Input.map_file_name is None):
    raise Sorry("Map file should be given")
  # check if file type is OK
  
  file_reader.any_file(file_name = params.cryo_fit2.Input.model_file_name).check_file_type(expected_type = 'pdb')
  
  print('validate_params pass')
  return True
# end of validate_params function

def run(args, prefix="tst_00", validated=False):
  user_input_pdb = ''
  user_input_map = ''
  
  # very simple parsing of model and map
  for i, arg in enumerate(args):
    if arg.endswith('.cif') or arg.endswith('.ent') or arg.endswith('.pdb'): # EMD-3981 has 6exv.ent instead of .pdb
      user_input_pdb = arg
      if arg.find('=')==-1:
        args[i]='model=%s' % arg
    elif arg.endswith('.ccp4') or arg.endswith('.map'):
      user_input_map = arg
      if arg.find('=')==-1:
        args[i]='map=%s' % arg
  
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="cryo_fit2",
  )
  
  user_input_pdb = clean_pdb_for_phenix(user_input_pdb)
  
  pdbs = []
  maps = []
  phils = []
  phil_args = []
  for arg in args:
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
      elif arg.endswith('.ccp4') or arg.endswith('.map'): # not the smartest
        maps.append(arg)
      else:
        try :
          file_phil = phil.parse(file_name=arg)
        except RuntimeError :
          pass
        else :
          phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  working_phil.show()
  working_params = working_phil.extract()
  
  if (not validated):
    validate_params(working_params)
    
  # Compute a target map
  from iotbx import ccp4_map
  ccp4_map = ccp4_map.map_reader(user_input_map)
  print('Map read from %s' %(user_input_map))
  target_map_data = ccp4_map.map_data()
    
  # initial atomic model that we want to fit to an EM-map
  pdb_inp = iotbx.pdb.input(file_name=user_input_pdb)
  model = mmtbx.model.manager(model_input = pdb_inp)
  
  # Initialize states accumulator
  states = mmtbx.utils.states(
    pdb_hierarchy  = model.get_hierarchy(),
    xray_structure = model.get_xray_structure())
  states.add(sites_cart = model.get_xray_structure().sites_cart())
  #
  params = sa.master_params().extract()
  params.start_temperature=2000
  params.final_temperature=0
  params.cool_rate = 100
  params.number_of_steps = 1000
  params.update_grads_shift = 0.
  params.interleave_minimization=False #Pavel will fix the error that occur when params.interleave_minimization=True
  
  print('CC: %s' %(calculate_cc(map_data=target_map_data, model=model, resolution=3.)))
  #STOP()
  result = sa.run(
    params = params,
    xray_structure     = model.get_xray_structure(),
    restraints_manager = model.get_restraints_manager(),
    target_map         = target_map_data,
    real_space         = True,
    wx                 = 100, # wx=5 broke helix conformation of tst_00_poor.pdb, wx=100 kept helix well
    wc                 = 1,
    states_collector   = states)
  states.write(file_name = "all_states.pdb")
  model.set_xray_structure(result.xray_structure)
  with open("refined.pdb", "w") as f:
    f.write(model.model_as_pdb())
# end of run function

if (__name__ == "__main__"):
  args = sys.argv[1:]
  if len(args) == 0:
    print ('Please provide user.pdb user.map')
    print ('Example: python cryo_fit2 user.pdb user.map')
    exit(1)
  run(args)
  print ('OK')
  