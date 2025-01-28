
"""
Initialize cells and connect them based on the connectivity matrix. 
Called by network_init()
"""
from neuron import h
import numpy as np
from neuron.units import ms, mV
from stellate import Stellate
from interneuron import Interneuron

h.nrnmpi_init()

pc = h.ParallelContext()


class Network():
    def __init__(self, id, adj_matrix, params):
        self.adj_matrix = adj_matrix
        self.params = params
        self.sim_dur = params['sim_dur']
        self.n_stell = params['N_stell']
        self.n_intrnrn = params['N_intrnrn']
        self._N = self.n_stell+self.n_intrnrn
        self.stell_gids = list(range(self.n_stell))
        self.intrnrn_gids = list(range(self.n_stell, self._N))

    
        self.exc_syn_ss_gmax = params['ss_syn_gmax'] 
        self.exc_syn_si_gmax = params['si_syn_gmax']
        self.inhb_syn_ii_gmax = params['ii_syn_gmax']
        self.inhb_syn_is_gmax = params['is_syn_gmax']
    
        
        self._set_gids()
        self._create_cells()
        self._connect_cells()

        self.id = id  # network id

    def __repr__(self):
        return 'network_{}'.format(self.id)

    def _create_cells(self):

        self.stellate_cells = []
        self.interneurons = []

        for i in self.gidlist:  # only create the cells that exist on this host
            if i < self.n_stell:
                self.stellate_cells.append(Stellate(i))

            else:
                self.interneurons.append(Interneuron(i))

        # associate the cell with this host and gid
        for cell in self.stellate_cells + self.interneurons:
            pc.cell(cell._gid, cell._spike_detector)

    def _connect_cells(self):
        for target in self.stellate_cells + self.interneurons:
            for i, weight in enumerate(self.adj_matrix[:, target._gid]):
                if weight != 0.:  # ignore zero weight (no synapse)
                    source_gid = i
                    if self.i_or_s(source_gid) == 'StellateCell' and target.name=="StellateCell":  #SS
                        nc = pc.gid_connect(source_gid, target.exc_syn)
                        nc.weight[0] = weight*self.exc_syn_ss_gmax
                    elif self.i_or_s(source_gid) == 'Interneuron' and target.name=="StellateCell":  #IS
                        nc = pc.gid_connect(source_gid, target.inhb_syn)
                        nc.weight[0] = weight*self.inhb_syn_is_gmax
                    elif self.i_or_s(source_gid) == 'StellateCell' and target.name=="Interneuron": #SI
                        nc = pc.gid_connect(source_gid, target.exc_syn)
                        nc.weight[0] = weight*self.exc_syn_si_gmax
                    elif self.i_or_s(source_gid) == 'Interneuron' and target.name=="Interneuron": #II
                        nc = pc.gid_connect(source_gid, target.inhb_syn)
                        nc.weight[0] = weight*self.inhb_syn_ii_gmax
                    nc.delay = self.params["netcon_delay"]
                    target._ncs.append(nc)


    def _set_gids(self):
        """Set the gidlist on this host."""
        self.gidlist = list(range(pc.id(), self._N, pc.nhost()))
        for gid in self.gidlist:
            pc.set_gid2node(gid, pc.id())

    def i_or_s(self, x):
        """Return the type of cell based on the gid."""
        if x < self.n_stell:
            return "StellateCell"
        elif x >= self.n_stell:
            return "Interneuron"


