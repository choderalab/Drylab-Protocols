How to install Openbabel
************************

Requires a valid licence.

Installation
------------

For python toolkit, download using ``pip`` as stated in the here_.

.. _here: https://docs.eyesopen.com/toolkits/python/quickstart-python/index.html

.. code-block:: bash
    pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits

Once installed, create the directory ``.openeye``, preferably in your home directory. Copy the licence file to that directory.

In your ``.bashrc`` (Linux) or ``.bash_profile`` (Mac), include

.. code-block::
    export OE_LICENSE=$HOME/.openeye/oe_license.txt
    export OE_DIR=$HOME/.openeye
