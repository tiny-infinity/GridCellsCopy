Installation
===================
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Grid Cells Conductance-based Model
A conductance-based model of [Grid cells](https://en.wikipedia.org/wiki/Grid_cell), implemented using the [NEURON](https://www.neuron.yale.edu/neuron/) simulator with parallelized functionality for efficient simulations. Grid cells are crucial for spatial navigation, and this model replicates their dynamics with biologically realistic conductances.

Key Features:
- Conductance-based implementation of Stellate cells and Inhibitory Interneurons of the medial entorhinal cortex.
- Parallelized Simulation: leveraging NEURON for multi-core simulations to handle large networks efficiently.
- User-Friendly Code: Clear structure with modular components for easy setup, modification, and analysis.

# Prerequisites

| To run simulations  | For analysis and plotting: |
| ------------- | ------------- |
| Python >= 3.12  | [SciPy>=1.13](https://nrn.readthedocs.io/en/latest/index.html)   |
| [NEURON>=8.2](https://nrn.readthedocs.io/en/latest/index.html) (with MPI support)  | [matplotlib>=3.9](https://nrn.readthedocs.io/en/latest/index.html) |
|[numpy >= 1.26.4](https://nrn.readthedocs.io/en/latest/index.html)|[seaborn>=0.13.2](https://nrn.readthedocs.io/en/latest/index.html) |
|[h5py >= 3.12](https://nrn.readthedocs.io/en/latest/index.html)||



# Running a simulation
To run a simulation, pass a `specs` file to `s_sim_setup.py`.
```
python s_sim_setup.py specs/s_template.py
```
A `specs` file contains a subset of parameters that override the default parameters to run a simulation. The default parameters are stored in `default_model_params.json` and `default_sim_params.json`. Data from the simulation is saved in `data/{sim_id}`, with sim_id specified in the `specs` file.

Check out `analysis\examples\base_model.ipynb` for some basic plots generated from the simulated data.

**Refer the [docs](https://inayath-sh.github.io/GridCellsCond/) for details on the project's structure and parameters.**

# Dependency Issues?
We use [`uv`](https://github.com/astral-sh/uv) package manager to manage dependencies in this project. You can use `uv.lock` to sync your local environment to the exact state in which these simulations were run. However, before running `uv` ensure, 

#### For macOS and Linux:
You have OpenMPI/MPICH
```bash
which mpiexec
```

#### For Windows:
Install [NEURON](https://nrn.readthedocs.io/en/latest/index.html) using the GUI
installer. 

Remove "NEURON>=8.2" from `pyproject.toml`.
<pre>
dependencies = [
    <strike>"NEURON>=8.2",</strike>
 "scipy>=1.13",
 "numpy>=1.26.4",
 ...
]
</pre>

### Install uv and sync environment:
In the project directory,

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync
```
### Activate env
```bash
source .venv/bin/activate
```
This must be executed in every instance of the terminal. You can configure [VS Code](https://code.visualstudio.com/docs/python/environments) to handle python environments.

# Cite
# Funding
