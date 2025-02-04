"""Utility functions for simulation setup, runs and handling data.
"""

import numpy as np
from itertools import zip_longest
import h5py
import json
import os
os.environ["NEURON_MODULE_OPTIONS"] = "-nogui" #Stops no gui warnings
import argparse
import importlib
import logging

def network_intialize(params):
    """Initialize network and setup instrumentation

    Args:
        params (dict or Param): Parameter dictionary

    Returns:
        network: Network object that includes cells and instrumentations

    """
    from network import Network
    #if asked to build matrix, load from cache. Otherwise load a previously saved matrix.
    if params['build_conn_matrix']:
        file = h5py.File(
            f"cache/matrix_{params['conn_id']}_{params['sim_id']}_{params['sim_num']}.hdf5", "r"
        )
        adj_matrix = file["matrix"]
    else:
        try:
            file = h5py.File(
                f"network_configs/connections/saved_matrices/matrix_"
                f"{params['conn_id']}_{params['matrix_id']}.hdf5", "r"
            )
        except FileNotFoundError as err:
            err.add_note(f"Cannot find saved matrix network_configs/connections/"
                         f"saved_matrices/matrix_{params['conn_id']}_{params['matrix_id']}.hdf5")
            raise err
        adj_matrix = file["matrix"]
    network = Network(0, adj_matrix, params)  # initialize grid cells
    file.close()
    
    #Add instrumentation
    setup_instrumentation = importlib.import_module(f"network_configs.instrument"
                                                    f"ations.{params['instr_id']}_instr").setup_instrumentation
    setup_instrumentation(network)
    return network

def load_sim_params(sim_id: str, file_path: str = None) -> dict:
    """Load simulation parameters.

    Args:
        sim_id (str): Simulation ID to load the parameters for.
        file_path (str, optional: Direct path to the JSON file containing 
            the simulation parameters. If not provided, the parameters are 
            located for the given `sim_id`.

    Returns:
        dict: The simulation parameters loaded from the JSON file.

    """

    if file_path == None:
        data_dir = locate_data_dir(sim_id)
        data_loc = f'{data_dir}{sim_id}/'
        file_path = data_loc+f'{sim_id}.json'
    return json_read(file_path)

def locate_data_dir(sim_id: str) -> str:
    """Find the location of the data directory for a given simulation ID.


    The function checks for the existence of the data directory in several 
    predefined locations.
    
    1. "data/" (local)
    2. "/data/{user}/data/" (global)

    The {user} placeholder is replaced with the current user's username obtained 
    from the environment variables "USERNAME" or "USER".
    
    Args:
        sim_id (str): The sim ID to locate the data directory for.

    Returns:
        str: The location of the data directory.
    """

    user= os.getenv("USERNAME") or os.getenv("USER")
    

    #predefined data locations
    data_locations ={"local":"data/","global":f"/data/{user}/data/","ada":f"/data/{user}/"}
    for key,val in data_locations.items():
        if os.path.exists(f"{val}{sim_id}"):
            return val
    raise FileNotFoundError(f"Cannot locate data for sim ID: {sim_id}")

def json_save(obj: dict, fname: str):
    """Wrapper function to save a dictionary object to a JSON file.
    
    Args:
        obj (dict): The dictionary object to save.
        fname (str): The file name to save the dictionary to.
    """

    with open(fname, "w") as file:
        json.dump(obj, file, indent=0)

def json_read(fname: str) -> dict:
    """Wrapper function to read a JSON file and return the dictionary object.
    
    Args:
        fname (str): The file name to read the dictionary from.
    Returns:
        dict: The dictionary object read from the JSON file.
    """
    with open(fname, "r") as file:
        obj = json.load(file)
    return obj

def json_modify(obj: dict, fname: str):
    """Modify an entry in JSON file using the given dictionary.
    
    If the file exists, it updates the file with the key-value pairs from the object.
    If the file does not exist, it creates a new file with the given object.

    Args:
        obj (dict): The dictionary object containing key-value pairs to be 
            added or updated in the JSON file.
        fname (str): The file name (including path) of the JSON file to be modified.
    """
    if os.path.isfile(fname):
        file_dict = json_read(fname)
        for key, val in obj.items():
            file_dict[key] = val
        json_save(file_dict, fname)
    else:
        json_save(obj, fname)


def list_to_numpy(LoL: list, fill:float = np.nan) -> np.ndarray:
    """Converts a list of lists (LoL) to a NumPy array, filling missing values with a specified fill value.
    
    Used for to convert list spike times with non-homogeneous lengths to a NumPy array to save in hdf5 format.

    Args:
        LoL (list of lists): The input list of lists to be converted to a NumPy array.
        fill (float, optional): The value to use for filling missing values. Defaults to np.nan.
    
    Returns
        numpy.ndarray: A NumPy array with the contents of the input list of lists, with missing values filled.
    """
    return np.array(list(zip_longest(*LoL, fillvalue=fill))).T


def check_sim_dup(sim_id: str, sim_num: int)->bool:
    """Checks if a simulation number exists in the HDF5 file for a given simulation ID.

    Args:
        sim_id (int): The ID of the simulation.
        sim_num (int): The number of the simulation to check for duplication.
    Returns:
        bool: True if the simulation number exists in the file, False otherwise.
    """

    fname = "data/sim_spikes_data_m_{}.hdf5".format(sim_id)
    if os.path.exists(fname):
        with h5py.File("data/sim_spikes_data_m_{}.hdf5".format(sim_id), "a") as file:
            if str(sim_num) in file.keys():
                return True
    return False



def find_sim_num(params: dict, param_check: dict) -> dict:
    """Find simulation numbers that was run with a given set of paramters.

    Useful for analysis. For e.g Finding the simulation number that was run with 
    ``si_peak=1.0``

    Args:
        params (dict or Param): Dictionary of all simulations generated after a 
            run of multi simulation.
        param_check (dict): A dictionary of parameters to check for. E.g: {si_peak: 1}

    Returns:
        matches (dict): A subset of global dictoinary containing only the simulations 
            that match the given parameters in param_check
    """
    
    matches = {}
    for sim_num, sim_param in params.items():
        correct_sim = True
        for key, value in sim_param.items():
            for check_cond_key, check_cond_val in param_check.items():
                if key == check_cond_key:
                    if str(value).startswith(str(check_cond_val)):
                        correct_sim = correct_sim and True
                    else:
                        correct_sim = correct_sim and False
        if correct_sim:
            matches[sim_num] = params[sim_num]
    return matches




def get_sim_num(iters: tuple, n_iters: tuple) -> int:
    """Calculate the simulation number for the running iterartor indices.

    Useful for analysis.
    
    Args:
        iters (tuple): A tuple of current iteration indices.
        n_iters (tuple): A tuple of the total size for each iterator.

    Returns:
        int: simulation number matching the given iterator indices.

    """

    assert len(iters)==len(n_iters) #The inputs should have equal dimensions
    sim_num=0
    n_iters = np.array(n_iters)
    for i,iter in enumerate(iters):
        sim_num+=iter*np.prod(n_iters[i+1:])
    return sim_num


def sim_setup_arg_parser()->argparse.ArgumentParser:
    """Set up the argument parser for the simulation.

    Argument parser is initialzed and processed here to avoid cluttering
    the main setup file.
    
    Returns:
        argparse.ArgumentParser: The argument parser with the arguments added.
    """
    
    parser = argparse.ArgumentParser(description="Run a single simulation")
    parser.add_argument("specs_file",
                    help="specificatons file",
                    type=str)
    parser.add_argument("-v","--verbose",
                    help="show verbose output",
                    action="store_const",const=logging.DEBUG,default=logging.INFO)
    
    parser.add_argument("-o","--overwrite_data",
                    help="overwrite data",
                    action="store_true")
    return parser

def sim_run_arg_parser()->argparse.ArgumentParser:
    """Set up the argument parser for the simulation run.

    Argument parser is initialzed and processed here to avoid cluttering
    the main run file.
    
    Returns:
        argparse.ArgumentParser: The argument parser with the arguments added.
    """
    
    parser = argparse.ArgumentParser(description="Run a single simulation")
    parser.add_argument("-i","--sim_id",
                    help="simulation ID",
                    type=str,
                    required=True)
    parser.add_argument("-v","--verbose",
                        help="show verbose output",
                        action="store_const",const=logging.DEBUG,default=logging.INFO)
    return parser

def log_from_rank_0(logger:logging.Logger,rank:int,msg:str,level:int=logging.INFO):
    """Logs a message if the rank is 0.
    
    Args:
        logger (logging.Logger): The logger object to use for logging.
        rank (int): The rank of the process.
        msg (str): The message to be logged.
    """
    if rank==0:
        logger.log(level,msg)

def process_data_root(data_root:str)->str:
    """Processes the given data root path to ensure it ends with a slash.
    
    Args:
        data_root (str): The root path to the data directory.
    Returns:
        str: The data root path ending with a slash.
    """
    return data_root+"/" if data_root[-1]!="/" else data_root


def load_spikes(sim_id:str,sim_num:int=0)->tuple:

    """Load spike data for a given simulation ID.
    
    Args:
        sim_id (str): Simulation ID to load spikes from.
        sim_num (int, optional): Simulation number to load spikes from, 
            by default 0 (single sim).
    
    Returns
        tuple
            A tuple containing two lists of lists:
            - stell_spikes_l: List of spike times for stellate cells.
            - intrnrn_spikes_l: List of spike times for interneurons.
    """
    
    data_dir = locate_data_dir(sim_id)
    # get spikes of a sim from hdf5 file
    sim_num = str(sim_num)
    data_loc = f"{data_dir}{sim_id}/"
    file_path_stell = data_loc + f"stell_spks_{sim_id}.hdf5"
    file_path_intrnrn = data_loc + f"intrnrn_spks_{sim_id}.hdf5"
    with h5py.File(file_path_stell, "r") as file:
        stellate_spks_arr = np.array(file[f"{sim_num}/stell_spks"][:])
        stell_spikes_l = [list(cell[~np.isnan(cell)]) for cell in stellate_spks_arr]
    with h5py.File(file_path_intrnrn, "r") as file:
        intrnrn_spks_arr = np.array(file[f"{sim_num}/intrnrn_spks"][:])
        intrnrn_spikes_l = [list(cell[~np.isnan(cell)]) for cell in intrnrn_spks_arr]
    
    return stell_spikes_l, intrnrn_spikes_l

class ProgressBar:
    """Progress bar for simulations.
    
    Not tested for multiple simulations.
    
    :meta private:
    """
    def __init__(self,total,pc=None):
        if os.name == 'nt':
            try:
                from colorama import just_fix_windows_console
                import sys
                sys.stdout.reconfigure(encoding='utf-8')
                sys.stderr.reconfigure(encoding='utf-8')
                just_fix_windows_console()
            except ModuleNotFoundError:
                print("colorama package not found. progress bar might not work")
        rank0=self._check_rank(pc)
        if rank0:
            self.marker='\x1b[31m█\x1b[39m'
            self.total= total
            self.length=50
            self.curr_progress=0
        

    def finish(self,pc=None):
        rank0=self._check_rank(pc)
        if rank0:
            self.iteration=self.total
            percent = 100 * (self.iteration / float(self.total))
            filled_length = int(self.length * self.iteration // self.total)
            bar = self.marker * filled_length + '-' * (self.length - filled_length)
            print(f"Progress ({self.iteration} of {self.total} ms): |{bar}| {percent:.2f}%",end="\n",flush=True)

    def increment(self,iteration,pc=None,flush=False):
        rank0=self._check_rank(pc)
        if rank0:
            self.iteration=iteration
            percent = 100 * (self.iteration / float(self.total))
            filled_length = int(self.length * self.iteration // self.total)
            bar = self.marker * filled_length + '-' * (self.length - filled_length)
            print(f"Progress ({self.iteration} of {self.total} ms): |{bar}| {percent:.2f}%",end="\r",flush=flush)

    def _check_rank(self,pc):
        if pc is not None:
            if pc.id()==0:
                return True
            return False
        else:
            return True



def get_module_from_path(file_path: str) -> str:
    """Convert a file path to a module name.

    Used to import specs file that is passed as as argument to 
    the simulation setup scripts.

    Args:
        file_path (str): The relative file path (e.g., "specs/s_template.py")
    
    Returns:
        str: The corresponding module name (e.g., "specs.s_template")
    """
    # Remove the file extension
    module_name, _ = os.path.splitext(file_path)
    # Replace path separators with dots
    module_name = module_name.replace(os.path.sep, ".")

    return module_name

def load_data(sim_id:str,data_id:str,cell_n:int=0,sim_num:str=0)->np.ndarray:
    """Load Non-spiking data for a given simulation ID.
    
    Args:
        sim_id (str): Simulation ID to load data from.
        data_id (str): Data ID. e.g. 'stell_v'
        cell_n (int, optional): Cell number to load data for, by default 0.
        sim_num (str, optional): Simulation number to load data from, by default 0.
    
    Returns:
        np.ndarray: The data array loaded from the HDF5 file
    """
    sim_num = str(sim_num)
    data_dir = locate_data_dir(sim_id)
    with h5py.File(f'{data_dir}{sim_id}/{data_id}_{sim_id}.hdf5', 'r') as f:
            data= np.array(f[str(sim_num)]['{}'.format(data_id)][cell_n])
    return data



def load_full_data(sim_id,data_id,sim_num=0):
    """Load Non-spiking data of all cells for a given simulation ID and number.
    
    Args:
        sim_id (str): Simulation ID to load data from.
        data_id (str): Data ID. e.g. 'stell_v'
        cell_n (int, optional): Cell number to load data for, by default 0.
        sim_num (str, optional): Simulation number to load data from, by default 0.
        
    Returns:
        np.ndarray: The data array loaded from the HDF5 file.
    """
    
    data_dir = locate_data_dir(sim_id)
    with h5py.File(f'{data_dir}{sim_id}/{data_id}_{sim_id}.hdf5', 'r') as f:
            data= np.array(f[str(sim_num)]['{}'.format(data_id)])
    return data