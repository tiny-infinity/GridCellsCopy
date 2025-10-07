"""
Global Parameter object
- Load default parameters from default_model_params.json and 
   default_sim_params.json
- Update and return the default parameters with the input parameters

mParam:
   - Contains nested dictionaries of parameters for multiple simulations with keys 
   as simulation numbers
"""

import json
import time
import copy

class Param(dict):
   def __init__(self):
      super(Param, self).__init__()
      self.load_default_params()
        
   def load_default_params(self):
      # load default parameters from: default_model_params.json and default_sim_params.json
      with open('default_model_params.json', 'r') as file:
         self.update(json.load(file))

      with open('default_sim_params.json', 'r') as file:
         self.update(json.load(file))

   def update_params(self, input_params):
      """Update the default parameters with the input parameters."""
      for key, val in input_params.items():
         if key not in self.keys():
               raise NameError(f"Invalid input param: {key}")
         
         if key in {"record_handle_stell", "record_handle_intrnrn"}:
               for rec_param, rec_val in val.items():
                  if rec_param not in self[key]:
                     raise NameError(f"Invalid recorder param: {rec_param}")
                  
                  self[key][rec_param]["state"] = rec_val["state"]
                  self[key][rec_param]["cells_to_record"] = rec_val.get(
                      "cells_to_record", self[key][rec_param].get("cells_to_record"))
         else:
               self[key] = val
              

   def save_curr_time(self):
      self.update({"sim_date_time":time.strftime("%Y-%m-%d, %X",time.localtime())})

   def __repr__(self):
      return "Parameter dictionary"
   

class mParam(dict):
   def __init__(self):
      super(mParam, self).__init__()
        
   def get_default_params(self):
      # load default parameters from: default_model_params.json and default_sim_params.json
      _default_params={}
      with open('default_model_params.json', 'r') as file:
         _default_params.update(json.load(file))

      with open('default_sim_params.json', 'r') as file:
         _default_params.update(json.load(file))   
   
      return _default_params
   def update_mult_params(self, sim_num,input_params):
      """Update the default parameters with the input parameters."""
      for key, val in input_params.items():
         if key not in self[sim_num].keys():
               raise NameError(f"Invalid input param: {key}")
         
         if key in {"record_handle_stell", "record_handle_intrnrn"}:
               for rec_param, rec_val in val.items():
                  if rec_param not in self[sim_num][key]:
                     raise NameError(f"Invalid recorder param: {rec_param}")
                  
                  self[sim_num][key][rec_param]["state"] = rec_val["state"]
                  self[sim_num][key][rec_param]["cells_to_record"] = rec_val.get(
                      "cells_to_record", self[sim_num][key][rec_param].get("cells_to_record"))
         else:
               self[sim_num][key] = val

   def load_update_mult_params(self,mult_input_params):
      """Load default parameters and update the with input parameters
      
      """
      _default_params=self.get_default_params()
      
      for sim_num,input_params in mult_input_params.items():
         sim_num=str(sim_num)
         self[sim_num]=copy.deepcopy(_default_params)
         self.update_mult_params(sim_num,input_params)
         self[sim_num]["sim_date_time"]=time.strftime("%Y-%m-%d, %X",time.localtime())
   def __repr__(self):
      return "Multiple Parameter dictionary"