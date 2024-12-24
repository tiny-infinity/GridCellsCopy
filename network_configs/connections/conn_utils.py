import argparse
import sys
from pathlib import Path
import os
import importlib

def file_args():
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
    directory=Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(directory))

def find_params(args):
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



        

