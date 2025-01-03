import sys
import numpy as np
from neuron import h
from neuron.units import ms, mV
import time
import sim_hf as s_hf
import h5py
h5py.get_config().track_order = True
import logging
import subprocess
import os
h.nrnmpi_init()
pc = h.ParallelContext()
tstart = time.perf_counter()
args, unknown= s_hf.sim_run_arg_parser().parse_known_args()
sim_id = args.sim_id
if pc.id() == 0:
    logging.basicConfig(handlers=[logging.FileHandler(f"logs/sim_{sim_id}.log",mode='w'),
            logging.StreamHandler()], encoding='utf-8', level=args.verbose,
            format=f'%(asctime)s:%(levelname)s:{sim_id}:%(message)s')
    logger = logging.getLogger()
else:
    logger = None

#load params file
mult_params = s_hf.json_read(f"cache/params_{sim_id}.json")

for sim_num, params in mult_params.items():
    tsim = time.perf_counter()
    s_hf.log_from_rank_0(logger,pc.id(),
                         f"Sim {sim_num}: Building connectivity matrix",
                         level=logging.DEBUG)
    if params["build_conn_matrix"] and pc.id()==0:
        cmd = f"python network_configs/connections/{params['conn_id']}_config.py -i {sim_id} -n {sim_num}"
        proc=subprocess.run(cmd.split(),check=True)    
    tmat = time.perf_counter()
    s_hf.log_from_rank_0(logger,pc.id(),
                         f"Sim {sim_num}: Connectivity matrix built in {round(tmat-tsim,2)}s",
                         level=logging.DEBUG)
    pc.barrier()
    #initialize network
    s_hf.log_from_rank_0(logger,pc.id(),
                         f"Sim {sim_num}: Initializing network",level=logging.DEBUG)
    
    network = s_hf.network_intialize(params)
    tinit = time.perf_counter()
    s_hf.log_from_rank_0(logger,pc.id(),
                         f"Sim {sim_num}: Network intialized in {round(tinit-tmat,2)}s",
                         level=logging.DEBUG)
    
    #load params
    sim_dur = network.params["sim_dur"]
    data_root = s_hf.process_data_root(network.params["data_root"])
    data_loc = data_root+f"{sim_id}/"
    
    #run simulation
    h.celsius = 37
    t = h.Vector().record(h._ref_t)
    pc.set_maxstep(10 * ms)
    h.finitialize(-65 * mV)
    pc.psolve(sim_dur * ms)
    
    #save data for stellate cells
    for param_to_record,states in network.params['record_handle_stell'].items():
        if states['state']==True:
            local_data={}
            
            for cell in network.stellate_cells:
                if cell.recorder.get(f'{param_to_record}',None):
                    local_data[cell._gid]= list(cell.recorder[f'{param_to_record}'])                        
            all_data = pc.py_alltoall([local_data] + [None] * (pc.nhost() - 1))
            pc.barrier()
            if pc.id() == 0:
                data = {}
                for process_data in all_data:
                    data.update(process_data)
                data = dict(sorted(data.items()))
                data_arr = s_hf.list_to_numpy(data.values())
                with h5py.File(data_loc + f"{param_to_record}_{sim_id}.hdf5", "a") as file:
                    group = file.create_group(sim_num)
                    group.create_dataset(f"{param_to_record}", data=data_arr, compression="gzip",dtype=np.float32)
    
    #save data for internuerons
    for param_to_record,states in network.params['record_handle_intrnrn'].items():
        if states['state']==True:
            local_data={}
            for cell in network.interneurons:
                if cell.recorder.get(f'{param_to_record}',None):
                    local_data[cell._gid]= list(cell.recorder[f'{param_to_record}'])
                    
            all_data = pc.py_alltoall([local_data] + [None] * (pc.nhost() - 1))
            pc.barrier()
            if pc.id() == 0:
                data = {}
                for process_data in all_data:
                    data.update(process_data)
                data = dict(sorted(data.items()))
                data_arr = s_hf.list_to_numpy(data.values())
                with h5py.File(data_loc + f"{param_to_record}_{sim_id}.hdf5", "a") as file:
                    group = file.create_group(sim_num)
                    group.create_dataset(f"{param_to_record}", data=data_arr, compression="gzip",dtype=np.float32)
    
    if pc.id() == 0:
        try:
            os.remove(f"cache/matrix_{params["conn_id"]}_{sim_id}_{sim_num}.hdf5")
        except FileNotFoundError as err:
            logging.warning(f"Cannot not remove cached matrix file")
            logging.debug(f"{err}")
    s_hf.log_from_rank_0(logger,pc.id(), f"Sim {sim_num}: Completed in {round(time.perf_counter()-tsim, 2)}s")
    pc.gid_clear(0)
    del network
    pc.barrier()
    pc.done()
h.quit()
