===============
Getting Started
===============


Installation
------------
Prerequisites
=============

1. `python >= 3.10 <https://www.python.org/>`_
2. `neuron >=8.2 <https://neuron.yale.edu/neuron/>`_
3. `numpy <https://numpy.org/>`_
4. `h5py <https://www.h5py.org/>`_

Clone the repository:
=====================

.. code-block:: bash

    git clone https:://github.com/GridCellsCond.git

First simulation
----------------
1. Compile the mod files:
    .. code-block:: bash

        cd GridCellsCond
        nrnivmodl mod_files
2. Run the simulation:
    .. code-block:: bash

        python s_sim_setup.py specs/s_template.py

3. Analysis:
    The simulation results will be saved in the ``data/`` folder.
    Refer to the notebooks in the ``analysis`` folder for few examples of how to visualise and analyse the data.

4. Modify the simulation:
    Modify the parameters in the ``specs/s_template.py`` file and run the simulation again.
    Refer `Model Parameters <model_parameters>`_ for the list of availaible parameters.

5. Next steps:
    Checkout the `Introduction <introduction>`_ section for more details on the model and the simulation setup.
    


.. _model_parameters: model_parameters
.. _introduction: introduction