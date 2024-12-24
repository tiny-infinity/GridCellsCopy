import argparse
import sys
from pathlib import Path
import os
import importlib

def file_args():
    """Creates an argument parser for connectivity configuration.

    Returns:
        argparse.ArgumentParser: The argument parser with the following arguments:
            - -i, --sim_id (str, optional): Simulation ID.
            - -c, --conn_id (str, optional): Connection ID.
            - -t, --sim_type (str, optional): Simulation type. Default is "single".
            - -n, --sim_num (int, optional): Simulation number. Default is 0.
            - -p, --specs_file (str, optional): Path to the specifications file.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--sim_id",
                        metavar='N',
                        help="Sim id",
                        type=str,
                        required=False)

    parser.add_argument("-c","--conn_id",
                        metavar='N',
                        help="Connection id",
                        type=str,
                        required=False)

    parser.add_argument("-t","--sim_type",
                        metavar='N',
                        help="simulation type",
                        type=str,
                        required=False,
                        default="single")

    parser.add_argument("-n","--sim_num",
                        metavar='N',
                        help="simulation number",
                        type=int,
                        required=False,
                        default=0)

    parser.add_argument("-p","--specs_file",
                        metavar='N',
                        help="specs file",
                        type=str,
                        required=False)
    return parser

def add_project_to_sys_path():
    """Adds the project's root directory to the system path.

    Python interpreter does not allow importing modules from parent directories.
    This function resolves the path of the current file, navigates two levels up 
    to reach the project's root directory, and inserts this directory at the 
    beginning of the system path. This allowes importing helper_functions, 
    param files, specs_file etc.

    """
    
    directory=Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(directory))


def find_params(args: argparse.Namespace) -> dict:
    """Finds and returns simulation parameters based on the provided arguments.
    
    Parameters:
        args (argparse.Namespace): The arguments namespace containing either a 
        simulation ID or a specifications file.
    Returns:
        dict: A dictionary containing the simulation parameters.
    Raises:
        FileNotFoundError: If neither sim_id nor specs_file is provided, or if 
        the specified files are not found.
    Notes:
        Parameter dict is generated based on following order of precedence:
        * If `args.sim_id` is provided, it attempts to read the parameters from 
            a JSON file in the cache directory.
            This is what occurs in a typical simulation run.
        * If the JSON file is not found, it loads the simulation parameters from 
            data directories using the provided simulation ID.
        * If `args.specs_file` is provided, it manually generates params dictionary.
    """

    import sim_hf as s_hf
    if args.sim_id:
        try:
            return s_hf.json_read(f"cache/s_params_merged_{args.sim_id}.json")
        except FileNotFoundError:
            
            return s_hf.load_sim_params(args.sim_id)
    if args.specs_file:
        from param import Param
        mod_name = os.path.split(args.specs_file)[0] + "." + \
            (os.path.split(args.specs_file)[1]).split(".")[0]
        param_file = importlib.import_module(mod_name)
        input_params = param_file.generate_input_params()

        #Initialize params dictionary
        return Param().update_input_params(input_params)
    raise  FileNotFoundError("No sim_id or specs_file provided")



        

