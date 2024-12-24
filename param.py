
import json
import time

"""
Global Parameter object to access all parameters in the simulation.
"""

class Param(dict):
   def __init__(self):
      super(Param, self).__init__()
      self.load_default_params_from_file()
        

   def load_default_params_from_file(self):
      # load default parameters from: default_model_params.json and default_sim_params.json
      with open('default_model_params.json', 'r') as file:
         self.update(json.load(file))

      with open('default_sim_params.json', 'r') as file:
         self.update(json.load(file))

   def update_input_params(self,input_params):
      # update the default parameters with the input parameters
      for key,val in input_params.items():
         if key in self.keys():
               if key=="record_handle_stell" or key=="record_handle_intrnrn":
                  for rec_param,rec_val in input_params[key].items():
                     if rec_param in self.params[key].keys():
                           self[key][rec_param]=rec_val
                     else:
                           raise NameError(f"Invalid recorder param: {rec_param}")
               else:
                  self[key]=val
         else:
               raise NameError(f"Invalid input param: {key}")
      return self
   
   def save_curr_time(self):
      self.update({"sim_date_time":time.strftime("%Y-%m-%d, %X",time.localtime())})

   def __repr__(self):
      return "Parameter dictionary"
   
   