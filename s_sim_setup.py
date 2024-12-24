import subprocess
import os
import sys
import warnings
import time
import importlib
import logging
import sim_hf as s_hf
from param import Param


tstart = time.perf_counter()

#Process specs file
try:
    fname = sys.argv[1]
except IndexError as err:
    err.add_note("Please pass specs file")
    raise err
mod_name = os.path.split(fname)[0] + "." + \
    (os.path.split(fname)[1]).split(".")[0]

#load input parameters
param_file = importlib.import_module(mod_name)
input_params = param_file.generate_input_params()

#Initialize params dictionary
params = Param().update_input_params(input_params)
params.save_curr_time()
sim_id = params["sim_id"]
n_cpus = params["n_cpus"] if params["n_cpus"] else os.cpu_count()//2

#start logger
logging.basicConfig(handlers=[logging.FileHandler(f"logs/setup_{sim_id}.log"),
        logging.StreamHandler()], encoding='utf-8', level=logging.DEBUG,
        format=f'%(asctime)s:%(levelname)s:{sim_id}: %(message)s')

# save parameters to cache for running the simulation
s_hf.json_save(params,f"cache/s_params_merged_{sim_id}.json")

#create data directory in data_root/
data_root = params["data_root"]+"/" if  params["data_root"][-1]!="/" else  params["data_root"] #ensure trailing slash
data_loc = data_root+f"{sim_id}/"
over_write_data=True #set to False to prevent overwriting data
if over_write_data:
    os.makedirs(data_loc, exist_ok=True) 
else:
    try:
        os.makedirs(data_loc)
    except FileExistsError as err:
        logging.error(f"Simulation data already exists at {data_loc}")
        raise err

logging.debug(f"Created data location at {data_loc}")

# copy specs file to data_loc for future reference
os.system(f"cp {fname} {data_loc}{sim_id}.py".format())
os.environ["NEURON_MODULE_OPTIONS"] = "-nogui" #Stops no gui warnings in output file
 


#build matrix
logging.debug(f"Building connectivity matrix")

t2 = time.perf_counter()
if params["build_conn_matrix"]:
    cmd = f"python networks/connections/{params['conn_id']}_conf.py -i {sim_id} -s {params["split_matrx"]}"
    proc=subprocess.run(cmd.split(),check=True)
logging.debug(f"Connectivity matrix built in {round(time.perf_counter()-t2,2)}s. saved in network_configs/connections/saved_matrices/")

#run simulation
# logging.debug(f"Starting simulation")
# if params["split_sim"][0]:
#     proc=subprocess.run(["mpiexec","-n", f"{n_cpus}","python", "s_run_split.py",f"{sim_id}"],check=True)
# else:
#     proc=subprocess.run(
#         ["mpiexec","-n", f"{n_cpus}","python", "s_run.py"],check=True)


# # save simulation time
# t_simulation=round(time.perf_counter() - tstart, 2)
# params["t_simulation"]= t_simulation


# #remove cache file
# try:
#     os.remove(f"cache/s_params_merged_{sim_id}.json")
#     os.remove(f"cache/conn_matrix_{sim_id}.json")

# except:
#     warnings.warn("Could not remove cache file")
#     logging.warning(f"Could not remove cache file")
        
# # copy params file to data for future reference
# param_fname = data_loc + "{}.json".format(params['sim_id'])
# s_hf.json_save(params, param_fname)

# logging.INFO(f"Total time (s): {round(time.perf_counter() - tstart, 2)}")