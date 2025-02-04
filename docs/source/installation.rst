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

Handling Dependencies
---------------------
This repository uses `uv <https://github.com/astral-sh/uv>`_ for dependency management. You can use ``uv.lock`` to 
sync your local environment to match the simulation requirements.


**Install uv and sync environment:**

In the project directory,

.. TAB:: Linux/macOS
    
    .. code-block:: console
        
        $ curl -LsSf https://astral.sh/uv/install.sh | sh
        $ uv sync

.. TAB:: Windows

    .. code-block:: powershell

        powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
        uv sync

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

**Compile the mod files and run a simulation:**

.. code-block:: console

    $ nrnivmodl mod

.. TAB:: Linux/macOS
    
    .. code-block:: console

        $ python s_sim_setup.py specs/s_template.py

.. TAB:: Windows

    .. code-block:: powershell

        $ python s_sim_setup.py specs\s_template.py

Cite
----

Funding
-------
