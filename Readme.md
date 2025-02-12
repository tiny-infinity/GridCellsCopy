[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Grid Cells Conductance-based Model
A conductance-based model of [Grid cells](https://en.wikipedia.org/wiki/Grid_cell), implemented using the [NEURON](https://www.neuron.yale.edu/) simulator with parallelized functionality for efficient simulations. Grid cells are crucial for spatial navigation, and this model replicates their dynamics with biologically realistic conductances.

Key Features:
- Conductance-based implementation of Stellate cells and Inhibitory Interneurons of the medial entorhinal cortex.
- Parallelized Simulation: leveraging NEURON for multi-core simulations to handle large networks efficiently.
- User-Friendly Code: Clear structure with modular components for easy setup, modification, and analysis.

## Prerequisites

| To run simulations  | For analysis and plots: |
| ------------- | ------------- |
| Python >= 3.12  | [SciPy>=1.13](https://scipy.org/install/)   |
| [NEURON>=8.2](https://nrn.readthedocs.io/en/latest/index.html) (with MPI support)  | [matplotlib>=3.9](https://matplotlib.org/stable/) |
|[numpy >= 1.26.4](https://numpy.org/install/)|[seaborn>=0.13.2](https://seaborn.pydata.org/installing.html) |
|[h5py >= 3.12](https://docs.h5py.org/en/latest/build.html)||



## Installing using uv
This repository uses [`uv`](https://github.com/astral-sh/uv) for dependency management. You can use `uv.lock` to 
sync your local environment to match the simulation requirements.

### Clone the repository:

```bash
$ git clone https://github.com/assisilab/GridCellsCond.git
$ cd GridCellsCond
```

### Install uv and sync environment:

In the project directory,

```bash
wget -qO- https://astral.sh/uv/install.sh | sh
uv sync
```
> [!IMPORTANT]
> For Windows, install NEURON through the GUI installer - [NEURON](https://nrn.readthedocs.io/en/latest/index.html)

### Activate environment:

For linux/macOS:

```bash
source .venv/bin/activate
```
For Windows:

```powershell
.venv\Scripts\activate
```
> [!WARNING]
> This must be executed in every instance of the terminal. You can configure [VS Code](https://code.visualstudio.com/docs/python/environments) to handle python environments.

## Running a simulation
To run a simulation, compile the mod files and pass a `specs` file to `s_sim_setup.py`.

```bash
nrnivmodl mod
```

```bash
python s_sim_setup.py specs/s_template.py # specs\s_template.py for windows
```

A `specs` file contains a subset of parameters that override the default parameters to run a simulation. The default parameters are stored in `default_model_params.json` and `default_sim_params.json`. Data from the simulation is saved in `data/{sim_id}`, with `sim_id` specified in the `specs` file.

`analysis\examples\BaseModel.ipynb` provides some basic plots generated from the simulated data.

**Refer the [docs](https://assisilab.github.io/GridCellsCond/) for details on the project's structure and parameters.**

# Cite
# Funding
