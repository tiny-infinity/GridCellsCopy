import sys
import time
from multiprocessing import Pool
import numpy as np
import conn_utils
import h5py

#build all at once,toroidal offsets
args = conn_utils.file_args().parse_args()
conn_utils.add_project_to_sys_path() #to access parent directories
params = conn_utils.find_params(args)

conn_id = params['conn_id']
sim_id = params['sim_id']
sim_num = params['sim_num']
N_per_sheet = params['N_per_sheet']
si_peak = params['si_peak']
ii_peak = params['ii_peak']
is_peak = params['is_peak']
N_per_sheet = params['N_per_sheet']
N_stell = params['N_per_sheet']*4
N_cells = N_stell + N_per_sheet
N_per_axis = params['N_per_axis']
si_mean = params['si_mean']
si_stdev = params['si_stdev']
is_stdev = params['is_stdev']
is_mean = params['is_mean']
ii_stdev = params['ii_stdev']
ii_mean = params['ii_mean']
threshold = 1e-7
N_per_axis = int(np.sqrt(N_per_sheet))
w_h = np.array((N_per_axis,N_per_axis))[:,np.newaxis]
coords_sheet = conn_utils.assign_positions_rect(N_per_axis)

#--------------------------------------------------------------------------------#

cnnct_ss = np.zeros((N_stell,N_stell))
cnnct_si = np.zeros((N_stell,N_per_sheet))
cnnct_is = np.zeros((N_per_sheet,N_stell))
cnnct_ii = np.zeros((N_per_sheet,N_per_sheet))

def task(i,j):
    vec=((coords_sheet[j]-coords_sheet[i]+N_per_axis/2)%N_per_axis)-N_per_axis/2
    dist=np.linalg.norm(vec)
    is_val= conn_utils.gaussian(dist,is_peak,is_mean,is_stdev)
    si_val_E= conn_utils.gaussian_2D(vec[0],vec[1],si_mean[0][0],si_mean[0][1],si_stdev,si_peak)
    si_val_W= conn_utils.gaussian_2D(vec[0],vec[1],si_mean[1][0],si_mean[1][1],si_stdev,si_peak)
    si_val_N= conn_utils.gaussian_2D(vec[0],vec[1],si_mean[2][0],si_mean[2][1],si_stdev,si_peak)
    si_val_S= conn_utils.gaussian_2D(vec[0],vec[1],si_mean[3][0],si_mean[3][1],si_stdev,si_peak)

    ii_val= conn_utils.gaussian(dist,ii_peak,ii_mean,ii_stdev)
    return is_val,ii_val,si_val_E,si_val_W,si_val_N,si_val_S
if __name__=='__main__':
    with Pool() as pool:
        # prepare arguments
        args = [(i,j) for i in range(N_per_sheet) for j in range(N_per_sheet)]
        result =  pool.starmap(task, args)

    for i,(j,k) in enumerate(args):
        cnnct_is[j,k]=cnnct_is[j,k+N_per_sheet]=cnnct_is[j,k+2*N_per_sheet]=cnnct_is[j,k+3*N_per_sheet] = result[i][0]
        cnnct_ii[j,k]= result[i][1]

        cnnct_si[j,k]= result[i][2]
        cnnct_si[j+N_per_sheet,k]= result[i][3]
        cnnct_si[j+2*N_per_sheet,k]= result[i][4]
        cnnct_si[j+3*N_per_sheet,k]= result[i][5]


        

    np.fill_diagonal(cnnct_ii,0)
    cnnct_si[cnnct_si<threshold] = 0
    cnnct_is[cnnct_is<threshold] = 0
    cnnct_ii[cnnct_ii<threshold] = 0
    adj_matrix=np.hstack((np.vstack((cnnct_ss,cnnct_is)),np.vstack((cnnct_si,cnnct_ii))))


    if params["save_conn_matrix"]:
        with h5py.File(f"network_configs/connections/saved_matrices/matrix_{conn_id}.hdf5", "w") as f:
            z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip')

    with h5py.File(f"cache/matrix_{conn_id}_{sim_id}_{sim_num}.hdf5", "w") as f:
        z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip') 
    

