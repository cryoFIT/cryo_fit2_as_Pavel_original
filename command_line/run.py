from __future__ import division
import iotbx.pdb
from mmtbx.dynamics import simulated_annealing as sa
import mmtbx.model
import mmtbx.utils
from scitbx.array_family import flex
import sys

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

def run(args, prefix="tst_00"):
  user_provided_pdb = False
  user_provided_map = False
  
  user_input_pdb = ''
  user_input_map = ''
  
  # very simple parsing of model and map
  for i, arg in enumerate(args):
    if arg.endswith('.pdb') or arg.endswith('.cif'):
      user_provided_pdb = True
      user_input_pdb = arg
      
    elif arg.endswith('.ccp4') or arg.endswith('.map') or arg.endswith('.sit'):
      user_provided_map = True
      user_input_map = arg
     
  # Compute a target map
  target_map_data = ''
  if user_provided_map == False:
    xrs_answer = pdb_inp.xray_structure_simple()
    fc = xrs_answer.structure_factors(d_min=3.0).f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
    write_ccp4_map(crystal_symmetry=fc.crystal_symmetry(), file_name="map.ccp4", \
                   map_data=target_map_data)
  else: # user_provided_map == True
    from iotbx import ccp4_map
    print 'user_input_map:', user_input_map
    ccp4_map = ccp4_map.map_reader(user_input_map)
    print "Map read from %s" %(user_input_map)
    target_map_data = ccp4_map.map_data()
    
  # initial atomic model that we want to fit to an EM-map
  pdb_inp = '' 
  if user_provided_pdb == False:
    pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_poor)
    pdb_inp.write_pdb_file(file_name="%s_poor.pdb"%prefix)
  else: #user_provided_pdb = True
    print "use user input pdb:", user_input_pdb
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
  #
  print "CC:", calculate_cc(map_data=target_map_data, model=model, resolution=3.)
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
    print "Please provide user.pdb user.map"
    print "Example: python dynamics.py user.pdb user.map"
    exit(1)
  run(args)
  print "OK"
  