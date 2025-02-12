Installation
============

Prerequisites
-------------

.. csv-table::
   :file: prereqs.csv
   :widths: 50,50
   :header-rows: 1


Installation with uv
--------------------

This repository uses `uv <https://github.com/astral-sh/uv>`_ for dependency management. You can use ``uv.lock`` to 
sync your local environment to match the simulation requirements.

Clone the repository:
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    $ git clone https://github.com/assisilab/GridCellsCond.git
    $ cd GridCellsCond


Install uv
^^^^^^^^^^

.. TAB:: Linux
    
    .. code-block:: console
        
        $ wget -qO- https://astral.sh/uv/install.sh | sh

.. TAB:: macOS/Linux
    
    .. code-block:: console
        
        $ curl -LsSf https://astral.sh/uv/install.sh | sh

.. TAB:: Windows

    .. code-block:: powershell

        powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"


See the `uv installation guide <https://docs.astral.sh/uv/getting-started/installation/>`_ for alternative installation methods.        

Sync environment
^^^^^^^^^^^^^^^^

Running ``uv sync`` sets up a virtual environment and installs all the dependencies.

.. important::

    On Windows systems, NEURON must be installed separately through its GUI installer - `NEURON <https://nrn.readthedocs.io/en/latest/index.html>`_ 

In the project directory, run:

.. code-block:: console
    
    $ uv sync

Activate environment:
^^^^^^^^^^^^^^^^^^^^^

.. TAB:: Linux/macOS
    
    .. code-block:: console
        
        $ source .venv/bin/activate

.. TAB:: Windows

    .. code-block:: powershell

        .venv\Scripts\activate

.. caution::

    This must be executed in every instance of the terminal. You can configure `VS Code <https://code.visualstudio.com/docs/python/environments>`_ to handle python environments.


Running a simulation
--------------------
To run a simulation, compile the mod files and pass a ``specs`` file to ``s_sim_setup.py``.

.. code-block:: console

    $ nrnivmodl mod

.. TAB:: Linux/macOS
    
    .. code-block:: console

        $ python s_sim_setup.py specs/s_template.py

.. TAB:: Windows

    .. code-block:: powershell

        $ python s_sim_setup.py specs\s_template.py

A ``specs`` file contains a subset of parameters that override the default parameters to run a simulation. 
The default parameters are stored in ``default_model_params.json`` and ``default_sim_params.json``. 
Data from the simulation is saved in ``data/{sim_id}``, with ``sim_id`` specified in the ``specs`` file.

For a high-level overview of the repository, checkout the :doc:`project structure <simulations>`.  
``analysis/examples/BaseModel.ipynb`` provides some basic plots generated from the simulation data.
