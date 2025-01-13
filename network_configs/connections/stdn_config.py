import numpy as np
import h5py
import conn_utils

args = conn_utils.file_args().parse_args()
conn_utils.add_project_to_sys_path() #to access parent directories
params = conn_utils.find_params(args)


conn_id = params['conn_id']
sim_id = params['sim_id']
sim_num = params['sim_num']
n_intrnrn = params['N_intrnrn']
n_stell = params['N_stell']
is_peak = params["is_peak"]
si_peak = params["si_peak"]
    
#--------------------------------------------------------------------------------#

cnnct_ss = np.zeros((n_stell, n_stell))

#--------------------------------------------------------------------------------#
cnnct_si = np.zeros((n_stell, n_intrnrn))
cnnct_si[0,0]=si_peak


#---------------------------------------------------------------------------------#
cnnct_is = np.zeros((n_intrnrn, n_stell))
cnnct_is[0,0]=is_peak

#----------------------------------------------------------------------------------#

cnnct_ii = np.zeros((n_intrnrn, n_intrnrn))

#----------------------------------------------------------------------------------#

adj_matrix = np.hstack(
    (np.vstack((cnnct_ss, cnnct_is)), np.vstack((cnnct_si, cnnct_ii))))

if params["save_conn_matrix"]:
    with h5py.File(f"network_configs/connections/saved_matrices/matrix_{conn_id}_{params['matrix_id']}.hdf5", "w") as f:
        z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip')

with h5py.File(f"cache/matrix_{conn_id}_{sim_id}_{sim_num}.hdf5", "w") as f:
    z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip') 