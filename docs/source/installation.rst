Installation
============

Prerequisites
-------------

.. csv-table::
   :file: prereqs.csv
   :widths: 50,50
   :header-rows: 1



Running a simulation
--------------------
To run a simulation, compile the mod files and pass a ``specs`` file to ``s_sim_setup.py``.

.. code-block:: console

    $ nrnivmodl mod

.. code-block:: console

    $ python s_sim_setup.py specs/s_template.py

A ``specs`` file contains a subset of parameters that override the default parameters to run a simulation. 
The default parameters are stored in ``default_model_params.json`` and ``default_sim_params.json``. 
Data from the simulation is saved in ``data/{sim_id}``, with sim_id specified in the ``specs`` file.

Check out :doc:`project structure <simulations>` for a detailed overview of the project. ``analysis\examples\BaseModel.ipynb`` 
contains some basic plots generated from the simulated data.

Handling Dependencies
---------------------
This repository is `uv <https://github.com/astral-sh/uv>`_ enabled to manage dependencies. You can use ``uv.lock`` to 
sync your local environment to the exact state in which these simulations were run in.


**Install uv and sync environment:**

In the project directory,

.. code-block:: console

    $ curl -LsSf https://astral.sh/uv/install.sh | sh
    $ uv sync

.. important::

    For Windows: Install NEURON independently through the GUI installer - `NEURON <https://nrn.readthedocs.io/en/latest/index.html>`_ 


**Activate env:**

.. TAB:: Linux/macOS
    
    .. code-block:: console
        
        $ source .venv/bin/activate

.. TAB:: Windows

    .. code-block:: powershell

        .venv\Scripts\activate

.. caution::

    This must be executed in every instance of the terminal. You can configure `VS Code <https://code.visualstudio.com/docs/python/environments>`_ to handle python environments.

Compile the mod files:

.. code-block:: console

    $ nrnivmodl mod

Cite
----

Funding
-------
