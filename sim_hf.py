
import numpy as np
from itertools import zip_longest
import h5py
import json
from itertools import islice
import os
import argparse

"""
Helper functions for the simulation setup, runs and data handling
"""

def load_sim_params(sim_id: str, file_path: str = None) -> dict:
    """Load simulation parameters.

    Parameters:
        sim_id : str
            Simulation ID to load the parameters for.
        file_path : str, optional 
            Direct path to the JSON file containing the simulation parameters. 
            If not provided, the parameters are located for the given `sim_id`.

    Returns:
        dict 
            The simulation parameters loaded from the JSON file.

    """

    if file_path == None:
        data_dir = locate_data_dir(sim_id)
        data_loc = f'{data_dir}{sim_id}/'
        file_path = data_loc+f'{sim_id}.json'
    return json_read(file_path)

def locate_data_dir(sim_id: str) -> str:
    """Find the location of the data directory for a given simulation ID.


    The function checks for the existence of the data directory in several predefined locations.
    
    1. "data/" (local)
    2. "/data/{user}/data/" (global)

    The {user} placeholder is replaced with the current user's username obtained from the environment variables "USERNAME" or "USER".
    
    Parameters:
        sim_id : str
            The sim ID to locate the data directory for.

    Returns:
        str
            The location of the data directory.
    """

    user= os.getenv("USERNAME") or os.getenv("USER")
    

    #predefined data locations
    data_locations ={"local":"data/","global":f"/data/{user}/data/","ada":f"/data/{user}/"}
    for key,val in data_locations.items():
        if os.path.exists(f"{val}{sim_id}"):
            return val
    raise FileNotFoundError("cannot locate data. verify sim_id")

def json_save(obj: dict, fname: str):
    """Wrapper function to save a dictionary object to a JSON file.
    
    Parameters:
        obj : dict 
            The dictionary object to save.
        fname: str
            The file name to save the dictionary to.
    """

    with open(fname, "w") as file:
        json.dump(obj, file, indent=0)

def json_read(fname: str) -> dict:
    """Wrapper function to read a JSON file and return the dictionary object.
    
    Parameters:
        fname : str
            The file name to read the dictionary from.
    Returns:
        dict 
            The dictionary object read from the JSON file.
    """
    with open(fname, "r") as file:
        obj = json.load(file)
    return obj

def json_modify(obj: dict, fname: str):
    """Modify an entry in JSON file using the given dictionary.
    
    If the file exists, it updates the file with the key-value pairs from the object.
    If the file does not exist, it creates a new file with the given object.
    
    Parameters
        obj : dict
            The dictionary object containing key-value pairs to be added or updated in the JSON file.
        fname : str
            The file name (including path) of the JSON file to be modified.
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

    Parameters
        LoL : list of lists
            The input list of lists to be converted to a NumPy array.
        fill : float, optional
            The value to use for filling missing values. Defaults to np.nan.
    
    Returns
        numpy.ndarray: A NumPy array with the contents of the input list of lists, with missing values filled.
    """
    return np.array(list(zip_longest(*LoL, fillvalue=fill))).T


def check_sim_dup(sim_id: str, sim_num: int)->bool:
    """Checks if a simulation number exists in the HDF5 file for a given simulation ID.

    Parameters:
        sim_id : int 
            The ID of the simulation.
        sim_num : int
            The number of the simulation to check for duplication.
    Returns:
        bool: 
            True if the simulation number exists in the file, False otherwise.
    """

    fname = "data/sim_spikes_data_m_{}.hdf5".format(sim_id)
    if os.path.exists(fname):
        with h5py.File("data/sim_spikes_data_m_{}.hdf5".format(sim_id), "a") as file:
            if str(sim_num) in file.keys():
                return True
    return False



def find_sim_num(params: dict, param_check: dict) -> dict:
    """Find simulation numbers that was run with a given set of paramters.

    Useful for analysis. For e.g Finding the simulation number that was run with ``si_peak=1.0``

    Parameters:
        params : dict or Param
            Dictionary of all simulations generated after a run of multi simulation.
        param_check : dict
            A dictionary of parameters to check for. E.g: {si_peak: 1}

    Returns:
        matches : dict
            A subset of global dictoinary containing only the simulations that match the given parameters in param_check
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
    
    Parameters:
        iters : tuple
            A tuple of current iteration indices.
        n_iters : tuple
            A tuple of the total size for each iterator.

    Returns:
        int
            simulation number matching the given iterator indices.

    """

    assert len(iters)==len(n_iters) #The inputs should have equal dimensions
    sim_num=0
    n_iters = np.array(n_iters)
    for i,iter in enumerate(iters):
        sim_num+=iter*np.prod(n_iters[i+1:])
    return sim_num


def sim_setup_arg_parser():
    """Set up the argument parser for the simulation.

    Argument parser is initialzed and processed here to avoid cluttering
    the main setup file.
    
    Returns
    -------
    argparse.ArgumentParser
        The argument parser with the arguments added.
    """
    
    parser = argparse.ArgumentParser(description="Run a single simulation")
    parser.add_argument("-p","--specs_file",
                    help="specificatons file",
                    type=str,
                    required=True)
    return parser


class ProgressBar:
    """Progress bar for simulations.
    
    Not tested for multiple simulations.
    
    :meta private:
    """
    def __init__(self,total):
        # self.progress_bar = tqdm(total=n_sim,position=0,leave=True,ncols=100)
        self.marker='\x1b[31m█\x1b[39m'
        self.total= total
        self.length=50
        self.curr_progress=0
        
    def update(self):
        self.curr_progress=self.curr_progress+1

        self.iteration=self.curr_progress
        percent = 100 * (self.iteration / float(self.total))
        filled_length = int(self.length * self.iteration // self.total)
        bar = self.marker * filled_length + '-' * (self.length - filled_length)
        # sys.stdout.write(f'\rProgress ({self.curr_progress} of {self.total}): |{bar}| {percent:.2f}%')
        # sys.stdout.flush()
        print(f"Progress ({self.iteration} of {self.total}): |{bar}| {percent:.2f}%",end="\r")
        

    def finish(self):
        self.iteration=self.total
        percent = 100 * (self.iteration / float(self.total))
        filled_length = int(self.length * self.iteration // self.total)
        bar = self.marker * filled_length + '-' * (self.length - filled_length)
        # sys.stdout.write(f'\rProgress ({iteration} of {self.total}): |{bar}| {percent:.2f}%\n')
        # sys.stdout.flush()
        print(f"Progress ({self.iteration} of {self.total}): |{bar}| {percent:.2f}",end="\n",flush=True)

    def increment(self,iteration):
        self.iteration=iteration
        percent = 100 * (self.iteration / float(self.total))
        filled_length = int(self.length * self.iteration // self.total)
        bar = self.marker * filled_length + '-' * (self.length - filled_length)
        # sys.stdout.write(f'\rProgress ({self.iteration} of {self.total}): |{bar}| {percent:.2f}%')
        # sys.stdout.flush()
        print(f"Progress ({self.iteration} of {self.total}): |{bar}| {percent:.2f}%",end="\r")




def invrtd_gaussian(x, N, mean, stdv):
    return -N * np.exp((-((x - mean) ** 2)) / (2 * (stdv**2))) + N 


def chunks(data, SIZE=10000):
    it = iter(data.items())
    res = []
    for i in range(0, len(data), SIZE):
        res.append(dict(islice(it, SIZE)))
    return res

