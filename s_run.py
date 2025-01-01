"""
Main simulation run. 
This file,
1. Initializes network based on params by calling 
    network_initialize from network_init.py
2. Runs simulation.
3. Saves data in hdf5 format.
"""

import numpy as np
from neuron import h
from neuron.units import ms, mV
import time
from network_init import network_intialize
import sim_hf as s_hf
import h5py
h5py.get_config().track_order = True
import logging


h.nrnmpi_init()
pc = h.ParallelContext()

tstart = time.perf_counter()
args = s_hf.sim_run_arg_parser().parse_args()
sim_id = args.sim_id
if pc.id() == 0:
    logging.basicConfig(handlers=[logging.FileHandler(f"logs/sim_{sim_id}.log"),
            logging.StreamHandler()], encoding='utf-8', level=logging.DEBUG,
            format=f'%(asctime)s:%(levelname)s: {sim_id}: %(message)s')
    logger = logging.getLogger()
else:
    logger = None
#load params file
params = s_hf.json_read(f"cache/s_params_merged_{sim_id}.json")

#initialize network
t2 = time.perf_counter()
s_hf.log_from_rank_0(logger,pc.id(),"Initializing network")
network = network_intialize(params)
tinit = time.perf_counter()
s_hf.log_from_rank_0(logger,pc.id(),f"Network intialized in {round(tinit-t2,2)}s")
sim_dur = network.params["sim_dur"]
sim_num = str(network.params["sim_num"])
sim_id = str(network.params["sim_id"])
data_root = network.params["data_root"]+"/" if  network.params["data_root"][-1]!="/" else  network.params["data_root"]
data_loc = data_root+f"{sim_id}/"

#run simulation
t3 = time.perf_counter()
h.celsius = 37
t = h.Vector().record(h._ref_t)
pc.set_maxstep(10 * ms)
s_hf.log_from_rank_0(logger,pc.id(),f"Simulation {sim_num} started")
if pc.id() == 0:
    if params["progress"]:
        pbar=s_hf.ProgressBar(total=int(sim_dur))
if not params["progress"]:
    #main simulatiion run without progress bar
    h.finitialize(-65 * mV)
    pc.psolve(sim_dur * ms)
else:
    h.finitialize(-65 * mV)
    for i in range(int(sim_dur)):
        pc.psolve(i * ms)
        if pc.id()==0:
            pbar.increment(int(pc.t(0)))
    if pc.id()==0:
        pbar.finish()        
tsim = round(time.perf_counter()-t3, 2)
s_hf.log_from_rank_0(logger,pc.id(),f"Simulation {sim_num} completed")
#save data for which recording is enabled
for param_to_record,states in network.params['record_handle_stell'].items():
    if states['state']==True:
        local_data={}
        for cell in network.stellate_cells:
            try:
                local_data[cell._gid]= list(cell.recorder[f'{param_to_record}'])
            except:
                continue            
        all_data = pc.py_alltoall([local_data] + [None] * (pc.nhost() - 1)) 
        pc.barrier()
        if pc.id() == 0:
            data = {}
            for process_data in all_data:
                data.update(process_data)
            data = dict(sorted(data.items()))
            data_arr = s_hf.list_to_numpy(data.values())
            with h5py.File(data_loc + f"{param_to_record}_{sim_id}.hdf5", "w") as file:
                group = file.create_group(sim_num)
                group.create_dataset(f"{param_to_record}", data=data_arr, compression="gzip",dtype=np.float32)

#save data for internuerons
for param_to_record,states in network.params['record_handle_intrnrn'].items():
    if states['state']==True:
        local_data={}
        for cell in network.interneurons:
            try:
                local_data[cell._gid]= list(cell.recorder[f'{param_to_record}'])
            except:
                raise Exception
                continue
        all_data = pc.py_alltoall([local_data] + [None] * (pc.nhost() - 1))
        pc.barrier()
        if pc.id() == 0:
            data = {}
            for process_data in all_data:
                data.update(process_data)
            data = dict(sorted(data.items()))
            data_arr = s_hf.list_to_numpy(data.values())
            with h5py.File(data_loc + f"{param_to_record}_{sim_id}.hdf5", "w") as file:
                group = file.create_group(sim_num)
                group.create_dataset(f"{param_to_record}", data=data_arr, compression="gzip",dtype=np.float32)

s_hf.log_from_rank_0(logger,pc.id(),f"Data saved in {data_loc}")
pc.barrier()
pc.done()
h.quit()
