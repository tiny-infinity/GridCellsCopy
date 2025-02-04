"""
Setup for single simulations
This file,
1. Gets input params from specs file
2. Generates full parameter dict from specs file 
    and default parameters.
3. Builds connectivity matrix if required.
4. Saves params file and matrix in cache/
5. Calls s_run.py to run simulation
6. Clears cache files
7. Saves specs file and params dict in 
    data directory for future reference
"""

import subprocess
import os
os.environ["NEURON_MODULE_OPTIONS"] = "-nogui" #Stops no gui warnings in output file
import time
import importlib
import logging
import sim_utils as s_utils
from param import Param
import shutil
import sys

tstart = time.perf_counter()
args,unk = s_utils.sim_setup_arg_parser().parse_known_args()
#Process specs file
fname = args.specs_file
mod_name = s_utils.get_module_from_path(fname)

#load input parameters
param_file = importlib.import_module(mod_name)
input_params = param_file.generate_input_params()

#check if sim_id declared
if not input_params.get("sim_id",None):
    raise ValueError("sim_id not set in specs file")

#Initialize params dictionary
params = Param()
params.update_params(input_params)
params.save_curr_time()

sim_id = params["sim_id"]
n_cpus = params["n_cpus"] if params["n_cpus"] else os.cpu_count()//2

#start logger
logging.basicConfig(handlers=[logging.FileHandler(f"logs/setup_{sim_id}.log",mode="w"),
        logging.StreamHandler()], encoding='utf-8', level=args.verbose,
        format=f'%(asctime)s:%(levelname)s:{sim_id}:Sim 0: %(message)s')

# save parameters to cache for running the simulation
s_utils.json_save(params,f"cache/params_{sim_id}.json")

#create data directory in data_root/
data_root = s_utils.process_data_root(params["data_root"])
data_loc = data_root+f"{sim_id}/"


try:
    if args.overwrite_data:
        logging.info(f"Overwriting data at {data_loc}")
        if os.path.exists(data_loc):
            shutil.rmtree(data_loc)
        os.makedirs(data_loc)
    else:
        os.makedirs(data_loc)
except FileExistsError as err:
    logging.error(f"Simulation data already exists at {data_loc}")
    err.add_note(f"Simulation data already exists at {data_loc}. Use -o to overwrite data")
    raise err

logging.debug(f"Created data location at {data_loc}")

# copy specs file to data_loc for future reference
os.system(f"cp {fname} {data_loc}{sim_id}.py")
 


#build matrix



if params["build_conn_matrix"]:
    t2 = time.perf_counter()
    logging.info(f"Building connectivity matrix")
    cmd = f"{sys.executable} network_configs/connections/{params['conn_id']}_config.py -i {sim_id}"
    proc=subprocess.run(cmd.split(),check=True)
    logging.info(f"Connectivity matrix built in {round(time.perf_counter()-t2,2)}s")
else:
    logging.info(f"Skipping matrix build. Using matrix_{params['conn_id']}.hdf5")

#run simulation
logging.debug(f"Launching s_run")
verbosity = "-v" if args.verbose == logging.DEBUG else ""

if params["split_sim"][0]:
    proc=subprocess.run(["mpiexec","-n", f"{n_cpus}",f"{sys.executable}", "s_run_split.py", "--sim_id",f"{sim_id}",f"{verbosity}"],check=True)
else:
    proc=subprocess.run(
        ["mpiexec","-n", f"{n_cpus}",f"{sys.executable}", "s_run.py", "--sim_id",f"{sim_id}",f"{verbosity}"],check=True)


# save simulation time
t_simulation=round(time.perf_counter() - tstart, 2)
params["t_simulation"]= t_simulation


#remove cache file
try:
    os.remove(f"cache/params_{sim_id}.json")
except FileNotFoundError as err:
    logging.debug(f"Cannot remove cached params file: {err}")


try:
    os.remove(f"cache/matrix_{params['conn_id']}_{sim_id}_0.hdf5")
except FileNotFoundError as err:
    logging.debug(f"Cannot remove cached matrix file: {err}")


            
# copy params file to data for future reference
param_fname = data_loc + f"{params['sim_id']}.json"
s_utils.json_save(params, param_fname)

logging.info(f"Total time (s): {round(time.perf_counter() - tstart, 2)}")