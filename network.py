
"""
Initialize cells and connect them based on the connectivity matrix. Called by network_init.py
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
        self.N_per_sheet = params['N_per_sheet'] #or ring
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

        # synapse where all targets are stellate cells
        for target in self.stellate_cells:

            for i, weight in enumerate(self.adj_matrix[:, target._gid]):
                if weight != 0.:  # ignore zero weight (no synapse)
                    source_gid = i
                    if self.i_or_s(source_gid) == 's':  # excitatory

                        nc = pc.gid_connect(source_gid, target.exc_syn)
                        nc.weight[0] = weight*self.exc_syn_ss_gmax
                        nc.delay = 1  # important for parallel

                        target._ncs.append(nc)
                    if self.i_or_s(source_gid) == 'i':  # inhibitory
                        nc = pc.gid_connect(source_gid, target.inhb_syn)
                        nc.weight[0] = weight*self.inhb_syn_is_gmax

                        nc.delay = 1
                        target._ncs.append(nc)

        # synapse where all targets are stellate cells
        for target in self.interneurons:

            for i, weight in enumerate(self.adj_matrix[:, target._gid]):

                if weight != 0.:

                    source_gid = i
                    if self.i_or_s(source_gid) == 's':  # exc

                        nc = pc.gid_connect(source_gid, target.exc_syn)
                        nc.weight[0] = weight*self.exc_syn_si_gmax
                        nc.delay = 1

                        target._ncs.append(nc)
                    if self.i_or_s(source_gid) == 'i':  # inhib

                        nc = pc.gid_connect(source_gid, target.inhb_syn)
                        nc.weight[0] = weight*self.inhb_syn_ii_gmax
                        nc.delay = 1

                        target._ncs.append(nc)

    def _set_gids(self):
        """Set the gidlist on this host."""
        self.gidlist = list(range(pc.id(), self._N, pc.nhost()))
        for gid in self.gidlist:
            pc.set_gid2node(gid, pc.id())

    def i_or_s(self, x):
        """Return type of neuron based in gid"""
        if x < self.n_stell:
            return 's'
        elif x >= self.n_stell:
            return 'i'