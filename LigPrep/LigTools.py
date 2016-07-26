"""
Small molecule preparation toolkit.

Description
-----------

A set of wrapper functions to protonate, add partial charges, and generate conformers of small molecules using a variety of
computational chemistry software.

The aim of these functions is to simplify comp. chem. tools, help standardize set-up protocols, and reduce barriers when
hopping between different comp. chem. tools.

Notes
-----

Functions are still being added.

Libraries that are supported are currently: OpenEye, Openbabel

References
----------

[1] Insert reference here!

Examples
--------

ProtonateLig_openeye('lig.sdf','lig.pdb','sdf','pdb')
ProtonateLig_openeye('lig.sdf',lig.pdb,protonate=False)
ProtonateLig_obabel('lig.mol2','lig_protonated.mol2')

TODO
----
    * Replace print statements with logger.
    * Consider aggregating functions into ligand prep class.
    
Copyright and license
---------------------

@author Gregory A. Ross <gregoryross.uk@gmail.com>

"""


def ProtonateLig_obabel(file_in,file_out,protonate=True,type_in=None,type_out=None):
    """
    Function to protonate small molecule using Openbabel. Requires either Openbabel python bindings, or openbabel
    accessible from the command line with the alias 'obabel'.

    Also doubles as a function to perform chemical file conversions. Files types can be specified or detected via the
    file extension. File types supported by OpenEye listed here: https://docs.eyesopen.com/toolkits/python/oechemtk/molreadwrite.html

    Parameters
    ----------
    file_in: str
        Name of the input file containing the coordinates of the molecule. One molecule expected.
    file_out: str
        Name of the output file.
    type_in: str
        File type of input file. Typically the file extension.
    type_out: str
        File type of the output file. Typically the file extension.

    Returns
    -------
    The output is written to a file
    """
    if type_in == None:
        type_in = file_in.split('.')[-1]
    if type_out == None:
        type_out = file_out.split('.')[-1]

    try:
        import openbabel

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(type_in, type_out)
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, file_out)

        if protonate==True: mol.AddHydrogens()

        obConversion.WriteFile(mol, file_out)
        print "OpenBabel: Added hydrogens to {0}. Output = {1}".format(file_in,file_out)

    except:
        print 'Warning: openbabel python bindings not found. Using command line...'
        import subprocess

        cmd = 'obabel -i{0} {1} -o{2} -O {3} -h'.format(type_in,file_in,type_out,file_out)
        subprocess.call(cmd,shell=True)
        print "OpenBabel: Added hydrogens to {0}. Output = {1}".format(file_in,file_out)

def ProtonateLig_openeye(file_in,file_out,protonate=True,type_in=None,type_out=None):
    """
    Function to protonate small molecule using OpenEye's python libraries. Requires OpenEye licence.
    More information found https://docs.eyesopen.com/toolkits/python/quickstart-python/index.html

    Also doubles as a file conversion script. Files types can be specified or detected via the file extension. File types
    supported by OpenEye listed here: https://docs.eyesopen.com/toolkits/python/oechemtk/molreadwrite.html

    Parameters
    ----------
    file_in: str
        Name of the input file containing the coordinates of the molecule. One molecule expected.
    file_out: str
        Name of the output file.
    type_in: str
        File type of input file. Typically the file extension.
    type_out: str
        File type of the output file. Typically the file extension.

    Returns
    -------
    The output is written to a file
    """
    try:
        import openeye.oechem as oechem

        if type_in == None:
            type_in = file_in.split('.')[-1]
        if type_out == None:
            type_out = file_out.split('.')[-1]

        mol = oechem.OEGraphMol()           # initialising molecule object
        ifs = oechem.oemolistream()         # initialising input stream for reading in data
        ofs = oechem.oemolostream()         # initialising the OUTPUT stream for writing data

        ifs.SetFormat(eval('oechem.OEFormat_' + type_in.upper()))
        ofs.SetFormat(eval('oechem.OEFormat_' + type_out.upper()))

        ifs.open(file_in)
        if oechem.OEReadMolecule(ifs,mol):  # this function automatically returns True or False, to help spot for errors.
            pass
        else:
            print "Problem loading molecule!"

        mol = oechem.OEMol(mol)
        oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)

        if protonate==True: oechem.OEAddExplicitHydrogens(mol)

        ofs.open(file_out)
        oechem.OEWriteMolecule(ofs,mol)

    except ImportError:
        print 'OpenEye python bindings not found'