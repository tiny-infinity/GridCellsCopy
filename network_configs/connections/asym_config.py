
import numpy as np
import h5py
import conn_utils

args = conn_utils.file_args().parse_args()
conn_utils.add_project_to_sys_path()
params = conn_utils.find_params(args)


def asymmetric_gaussian(x, N,mu, sigma_left, sigma_right):
    return N*np.where(
        x < mu,
        np.exp(-((x - mu) ** 2) / (2 * sigma_left ** 2)),
        np.exp(-((x - mu) ** 2) / (2 * sigma_right ** 2)))

conn_id = params['conn_id']
sim_id = params['sim_id']
sim_num = params['sim_num']
n_intrnrn = N_per_sheet = params['N_intrnrn']
n_stell = params['N_stell']
n_stell_per_ring = n_stell//2
si_peak = params['si_peak']
ii_peak = params['ii_peak']
is_peak = params['is_peak']

si_mean = params['si_mean']*((2*np.pi)/N_per_sheet)
si_stdev = params['si_stdev']
is_stdev = params['is_stdev']
is_mean = params['is_mean']*((2*np.pi)/N_per_sheet)
ii_stdev = params['ii_stdev']
ii_mean = params['ii_mean']*((2*np.pi)/N_per_sheet)
threshold = 1e-7
si_asym_factor = params["si_asym_factor"]
is_asym_factor = params["is_asym_factor"]
is_peak_asym_fact = params["is_peak_asym_fact"]
is_mean_asym_factor = params["is_mean_asym_factor"]*((2*np.pi)/N_per_sheet)

#--------------------------------------------------------------------------------#

cnnct_ss = np.zeros((n_stell, n_stell))


#--------------------------------------------------------------------------------#
cnnct_si = np.zeros((n_stell, n_intrnrn))

for i in range(N_per_sheet):
    for j in range(N_per_sheet):
        d_theta = (2*np.pi/N_per_sheet)*(j-i)
        i_rad = (2*np.pi/N_per_sheet)*i
        j_rad = (2*np.pi/N_per_sheet)*j
        dtheta=conn_utils.signed_arc_length(i_rad,j_rad)
        cnnct_si[i, j] = asymmetric_gaussian(dtheta, si_peak, si_mean, si_stdev*si_asym_factor[0],si_stdev*si_asym_factor[1])


temp_mat = np.zeros((N_per_sheet, N_per_sheet))
for i in range(N_per_sheet):
    for j in range(n_intrnrn):
        d_theta = (2*np.pi/N_per_sheet)*(j-i)
        i_rad = (2*np.pi/N_per_sheet)*i
        j_rad = (2*np.pi/N_per_sheet)*j
        dtheta=conn_utils.signed_arc_length(i_rad,j_rad)

        
        temp_mat[i, j] = asymmetric_gaussian((-2*si_mean-dtheta), si_peak, -si_mean, si_stdev*si_asym_factor[0],si_stdev*si_asym_factor[1])
cnnct_si[N_per_sheet:, :] = temp_mat
cnnct_si[cnnct_si < threshold] = 0


#---------------------------------------------------------------------------------#
cnnct_is = np.zeros((n_intrnrn, n_stell))

for i in range(n_intrnrn):
    for j in range(N_per_sheet):
        i_rad = (2*np.pi/N_per_sheet)*i
        j_rad = (2*np.pi/N_per_sheet)*j
        dtheta=conn_utils.signed_arc_length(i_rad,j_rad)
        
        if dtheta >= 0:
            cnnct_is[i, j] = asymmetric_gaussian(dtheta, is_peak*is_peak_asym_fact, is_mean+is_mean_asym_factor, is_stdev*is_asym_factor[0],is_stdev*is_asym_factor[1])
        else:
            cnnct_is[i, j] = asymmetric_gaussian(dtheta, is_peak*is_peak_asym_fact, -is_mean+is_mean_asym_factor, is_stdev*is_asym_factor[0],is_stdev*is_asym_factor[1])
temp_mat = np.zeros((N_per_sheet, N_per_sheet))

for i in range(n_intrnrn):
    for j in range(N_per_sheet):
        i_rad = (2*np.pi/N_per_sheet)*i
        j_rad = (2*np.pi/N_per_sheet)*j
        dtheta=conn_utils.signed_arc_length(i_rad,j_rad)
        
        if dtheta >= 0:
            temp_mat[i, j] = asymmetric_gaussian((2*is_mean-dtheta), is_peak*is_peak_asym_fact, is_mean+is_mean_asym_factor, is_stdev*is_asym_factor[0],is_stdev*is_asym_factor[1])
        else:
            temp_mat[i, j] = asymmetric_gaussian((-2*is_mean-dtheta), is_peak*is_peak_asym_fact, -is_mean+is_mean_asym_factor, is_stdev*is_asym_factor[0],is_stdev*is_asym_factor[1])
cnnct_is[:, N_per_sheet:] = temp_mat
            
cnnct_is[cnnct_is < threshold] = 0


#----------------------------------------------------------------------------------#

cnnct_ii = np.zeros((n_intrnrn, n_intrnrn))

for i in range(n_intrnrn):
    for j in range(n_intrnrn):
        i_rad = (2*np.pi/N_per_sheet)*i
        j_rad = (2*np.pi/N_per_sheet)*j
        dtheta=conn_utils.signed_arc_length(i_rad,j_rad)
        
        if dtheta >= 0:
            cnnct_ii[i, j] = conn_utils.gaussian(dtheta, ii_peak, ii_mean, ii_stdev)
        else:
            cnnct_ii[i, j] = conn_utils.gaussian(dtheta, ii_peak, -ii_mean, ii_stdev) 
            
cnnct_ii[cnnct_ii < threshold] = 0
np.fill_diagonal(cnnct_ii,0)



#----------------------------------------------------------------------------------#

adj_matrix = np.hstack(
    (np.vstack((cnnct_ss, cnnct_is)), np.vstack((cnnct_si, cnnct_ii))))

if params["save_conn_matrix"]:
    with h5py.File(f"network_configs/connections/saved_matrices/matrix_{conn_id}_{params['matrix_id']}.hdf5", "w") as f:
        z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip')

with h5py.File(f"cache/matrix_{conn_id}_{sim_id}_{sim_num}.hdf5", "w") as f:
    z=f.create_dataset(str('matrix'),data=adj_matrix,compression='gzip') 

