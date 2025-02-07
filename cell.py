#Parent Class for all cells

from neuron import h
class Cell:
    def __init__(self,gid):
        self._gid = gid               #global id of the cell
        self._set_morphology()
        self.all = self.soma.wholetree()
        self._set_biophysics()
        self._default_instrumentation()

        #vectors to record spike times
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self._ncs = []

    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
