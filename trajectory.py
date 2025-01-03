import numpy as np
import h5py as hdf
import pickle
from scipy.interpolate import BSpline
import sim_hf as s_hf
from neuron import h

class Trajectory:
    def __init__(self, params, save_mem=True):
        # self.max_vel = 35 * 1e-3  # cm/ms
        self.params = params
        self.sim_dur = params["sim_dur"]
        self.dt = params["dt"]
        self.t = np.arange(0, self.sim_dur + self.dt, self.dt)
        self.N_per_sheet = self.params["N_per_sheet"]
        self.save_mem = save_mem
        self.additional_inputs = {"allothetic_input":params["allothetic_input"]}
        if self.params["vel_type"] == "const":
            self.const_vel()
        elif self.params["vel_type"] == "input":
            self.input_velocity()
        elif self.params["vel_type"] == "ACVT-1DAC":
            self.figure_1_pulse_input()


    def vel_to_dc(self, vel_input, min_dc, max_dc):
        return ((max_dc - min_dc) / self.max_vel) * vel_input + min_dc

    def vel_to_dc_fit(self, input_vel):
        #spline
        spline_params=s_hf.json_read("input_data/vi_transform/spline_params.json")
        dc_out=BSpline(*list(spline_params.values()))(input_vel)
        zero_mask = np.abs(input_vel) < 0.0001
        dc_out[zero_mask] = self.params["vel_integ_zero"]
        other_ring_mask = np.isnan(input_vel)
        dc_out[other_ring_mask] = self.params["vel_integ_or"]

        return dc_out



    def decompose_vel(self):
        self.right_vel = self.vel_input.copy()
        self.left_vel = self.vel_input.copy()
        self.right_vel[self.right_vel < 0] = np.nan
        self.left_vel[self.left_vel > 0] = np.nan

    def const_vel(self):
        self.intrnrn_dc=np.full_like(self.t, self.params["intrnrn_dc_amp"])
        self.right_const_dc = self.params["stell_const_dc"][0]
        self.left_const_dc = self.params["stell_const_dc"][1]
        self.left_dc = np.full_like(self.t, self.left_const_dc)
        self.right_dc = np.full_like(self.t, self.right_const_dc)


    def input_velocity(self):

        with hdf.File("input_data/trajectories/traj_{}.hdf5".format(self.params["traj_id"]), "r") as file:
            self.vel_input = np.array(file["vel_rinb"][:])
            self.vel_input=self.vel_input*self.params["vel_integ_multiple"]   #*1.37
            self.pos_input = np.array(file["pos_rinb"][:])
            self.allothetic_dur = float(file.attrs["allothetic_dur"])
        self.decompose_vel()
        self.right_dc = self.vel_to_dc_fit(self.right_vel)
        self.left_dc = self.vel_to_dc_fit(-1 * self.left_vel)
        self.intrnrn_dc=np.full_like(self.t, self.params["intrnrn_dc_amp"])
        if self.params["allothetic_input"]:
            self.allothetic_input()
        if self.save_mem:
            del self.vel_input, self.right_vel, self.left_vel, self.pos_input

    def allothetic_input(self):
        # duration of allothetic input
        self.init_idx = int((self.allothetic_dur) / 0.025)
        self.init_position = self.pos_input[
            self.init_idx
        ]
        self.params["allothetic_nrn_n"]
        self.allothetic_left_dc = self.left_dc.copy()
        self.allothetic_right_dc = self.right_dc.copy()
        self.allothetic_left_dc[: self.init_idx] = np.linspace(1e-2,self.left_dc[self.init_idx],self.allothetic_left_dc[: self.init_idx].shape[0])
        self.allothetic_right_dc[: self.init_idx] = np.linspace(1e-2,self.right_dc[self.init_idx],self.allothetic_right_dc[: self.init_idx].shape[0])
        self.left_dc[: self.init_idx+int(0/0.025)] =self.params["allothetic_stell_dc"]
        self.right_dc[: self.init_idx+int(0/0.025)]=self.params["allothetic_stell_dc"]
        phi = np.round(
            (self.init_position*self.params["n_phases"]) / (self.params["lambda0"])
        ).astype("int")
        t_active_grids = ((np.arange(int(self.N_per_sheet /self.params["n_phases"])) * \
        self.params["n_phases"])[:,None] +np.arange(-self.params["allothetic_nrn_n"],self.params["allothetic_nrn_n"]+1)+phi)%self.N_per_sheet
        t_active_grids=t_active_grids.ravel()
        t_active_grids= t_active_grids + self.params["N_stell"]
        self.active_cells =t_active_grids

        self.allothetic_intrnrn_dc = self.intrnrn_dc.copy()
        self.allothetic_intrnrn_dc[: self.init_idx] =1.5e-3
        self.intrnrn_dc[: self.init_idx] =-1e-3
        self.ext_amp_right_allo = h.Vector(self.allothetic_right_dc)
        self.ext_amp_left_allo = h.Vector(self.allothetic_left_dc)
        self.ext_amp_intnrn_allo = h.Vector(self.allothetic_intrnrn_dc)
    def figure_1_pulse_input(self):
        "Used to generate a pulse input that is used for raster plot in figure one of the paper"
        # define step inputs
        x = [0, 5000, 10000, 15000]
        left_ring = [-1.5e-3, -3e-3, 0]
        right_ring = [-1.5e-3, 0, -3e-3]
        self.left_dc = self.create_piecewise(x, left_ring)
        self.right_dc = self.create_piecewise(x, right_ring)
        self.intrnrn_dc=np.full_like(self.t, self.params["intrnrn_dc_amp"])

    def create_piecewise(self, x, y):
        t_ = np.arange(0, self.sim_dur + self.dt, self.dt)
        l = []
        for i in range(len(x) - 1):
            l.append(np.logical_and(t_ >= x[i], t_ <= x[i + 1]))

        return np.piecewise(t_, l, y)

    def add_additional_inputs(self,cell):
        if self.additional_inputs["allothetic_input"]:
            if cell.name == "StellateCell" and cell._gid in self.active_cells:
                self.ext_amp_left_allo.play(cell.ext_dc._ref_amp, self.t, True)
            if cell.name == "Interneuron" and cell._gid in self.active_cells:
                self.ext_amp_intnrn_allo.play(cell.ext_dc._ref_amp, self.t, True)

