
import subprocess
import os
os.environ["NEURON_MODULE_OPTIONS"] = "-nogui" #Stops no gui warnings in output file
import sim_utils as s_utils
import time
import importlib
from param import mParam
import logging

tstart = time.perf_counter()
args = s_utils.sim_setup_arg_parser().parse_args()
#Process specs file
fname = args.specs_file
mod_name = os.path.split(fname)[0] + "." + \
    (os.path.split(fname)[1]).split(".")[0]

#load input parameters
param_file = importlib.import_module(mod_name)
mult_input_params = param_file.generate_mult_input_params()
#Initialize params dictionary
mult_params = mParam()
mult_params.load_update_mult_params(mult_input_params)

sim_id = mult_params["0"]["sim_id"]

n_cpus = 32 #mult_params["0"]["n_cpus"] if mult_params["0"]["n_cpus"] else os.cpu_count()//2

#start logger
logging.basicConfig(handlers=[logging.FileHandler(f"logs/setup_{sim_id}.log",mode="w"),
        logging.StreamHandler()], encoding='utf-8', level=args.verbose,
        format=f'%(asctime)s:%(levelname)s:{sim_id}:%(message)s')

# save parameters to cache for running the simulation
s_utils.json_save(mult_params,f"cache/params_{sim_id}.json")

#create data directory in data_root/
data_root = s_utils.process_data_root(mult_params["0"]["data_root"])
data_loc = data_root+f"{sim_id}/"
over_write_data=True #set to False to prevent overwriting data

try:
    os.makedirs(data_loc)
except FileExistsError as err:
    logging.error(f"Simulation data already exists at {data_loc}")
    raise err

logging.debug(f"Created data location at {data_loc}")

# copy specs file to data_loc for future reference
os.system(f"cp {fname} {data_loc}{sim_id}.py".format())

n_cpus = 3 #params["0"]["n_cpus"] if params["0"]["n_cpus"] \
    #else os.cpu_count()//2

#run simulation
logging.debug(f"Launching m_run")
verbosity = "-v" if args.verbose == logging.DEBUG else ""
proc=subprocess.run(["mpiexec","-n", f"{n_cpus}","python", "m_run.py", 
                     "--sim_id",f"{sim_id}",f"{verbosity}"],check=True)

# save simulation time
t_simulation=round(time.perf_counter() - tstart, 2)
mult_params["0"]["t_simulation"]= t_simulation

#remove cache file
try:
    os.remove(f"cache/params_{sim_id}.json")
except FileNotFoundError as err:
    logging.debug(f"Cannot remove cached params file: {err}")

# copy params file to data for future reference
param_fname = data_loc + f"{sim_id}.json"
s_utils.json_save(mult_params, param_fname)
logging.info(f"Total time (s): {round(time.perf_counter() - tstart, 2)}")
