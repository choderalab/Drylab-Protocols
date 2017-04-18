""""

A script to run Schrodinger's Protein Prep, relies on Schrodinger 2016-4, on a specified PDB file

Written by Steven Albanese, with gratuitious borrowing from Openmoltools' Schrodinger package, courtesy of Andrea Rizzi

OpenMM Refinement code written by Gregory Ross

Modified to separate input and output directories by Mehtap Isik.
"""


#################
#     Import    #
#################

from openmoltools import utils
from openmoltools import schrodinger
import os
import argparse
import shutil
import csv
import subprocess
from openmoltools.schrodinger import need_schrodinger
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmoltools.forcefield_generators import gaffTemplateGenerator

#################
#     Parser    #
#################

parser = argparse.ArgumentParser(description="Automated script to search PDB by chemical ID")
parser.add_argument('--input', required=True, dest='pdb_file',
                    help='THe pdb file that you would like to run through this program')
parser.add_argument('--ph', required=False, default=7.4, type=float, dest='ph',
                    help='Use to set pH to something other than 7.0')
parser.add_argument('--output', required=False, default='output.pdb', dest='output',
                    help='the name of the output file')
parser.add_argument('--cap', required=False, action='store_true', dest='capterminal',
                    help='set flag to cap terminal amino acids')

args = parser.parse_args()

file_name = args.pdb_file
ph = args.ph
output = args.output
cap = args.capterminal

def write_file(filename, contents):
    """
    Little helper function to write the pdb files

    Args:
        filename: String, 4-letter PDB ID
        contents: string that will be written to the file

    Returns: Nothing, just writes the file

    """

    with open(filename, 'w') as outfile:
        outfile.write(contents)

@need_schrodinger
def protein_prep(input_file, output_file, cap, pH=7.4, fillsidechains=True, fillloops=True,
                 noepik=False, rehtreat=True, max_states=32, tolerance=0):

    # Locate PrepWizard executable
    prepwiz_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'prepwizard')

    # Normalize paths
    input_file_path = os.path.abspath(input_file)
    input_file_dir, input_name = os.path.split(input_file_path)
    output_file_path = os.path.abspath(output_file)
    output_file_dir, output_name = os.path.split(output_file_path)

    output_file_name = os.path.join(output_file_dir, output_name + '-prepped.pdb')
    working_dir = os.path.join(output_file_dir, '%s-prepped' % input_name)

    # Check for output file pathway
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    # Format arguments for PrepWizard command

    wiz_args = dict(ms=max_states, ph=pH)
    wiz_args['fillsidechains'] = '-fillsidechains' if fillsidechains else ''
    wiz_args['fillloops'] = '-fillloops' if fillloops else ''
    wiz_args['pht'] = tolerance
    wiz_args['rehtreat'] = '-rehtreat' if rehtreat else ''
    wiz_args['water_hbond_cutoff'] = 0
    wiz_args['noepik'] = '-noepik' if noepik else ''
    wiz_args['captermini'] = '-captermini' if cap else ''


    cmd = [prepwiz_path]
    cmd += '{captermini} -mse -propka_pH {ph} {fillsidechains} {fillloops} {rehtreat} {noepik} -delwater_hbond_cutoff {water_hbond_cutoff} ' \
           '-keepfarwat -disulfides -ms {ms} -minimize_adj_h -epik_pH {ph} -epik_pHt {pht} -fix -NOJOBID'.format(**wiz_args).split()

    cmd.append(input_file_path)
    cmd.append(output_file_name)

    with utils.temporary_cd(working_dir):
        log = schrodinger.run_and_log_error(cmd)
        write_file('protein_prep.log', log)


def discard_organic(model, verbose=True):
    """
    Return an OpenMM modeller object that doesn't contain small organic molecules from the mother liquor

    Parameter
    ---------
    model: simtk.openmm.app.modeller.Modeller
        The modeller object with the PDB file of interest

    Returns
    -------
    model: simtk.openmm.app.modeller.Modeller
        The same object as the input with certain residues discarded
    """
    unwanted = ['DTT', 'EDO', 'GOL', 'SO4', 'PO4', 'DMS']
    for junk in unwanted:
        atms = [atm for atm in model.topology.atoms() if atm.residue.name == junk]
        if len(atms) > 0:
            if verbose == True: print('Deleting {0} from topology'.format(junk))
            model.delete(atms)
    return model

def openmm_clean(input_file, output_file, solvate=False):
    """
    Removes unwanted molecules and residues, and performs a short minimization.

    If returns IOError about missing gaff.xml file, manually copy
    gaff.xml to OpenMM path indicated in the error message.
    """

    # Normalize paths
    input_file_path = os.path.abspath(input_file)
    input_file_dir, input_name = os.path.split(input_file_path)

    output_file_path = os.path.abspath(output_file)
    output_file_dir, output_name = os.path.split(output_file_path)

    input_name = os.path.join(input_file_dir, input_name + '-prepped.pdb')
    output_name = os.path.join(output_file_dir, output_name + '-prepped.pdb')

    output_file_name = os.path.join(output_file_dir, output_file + '-minimized.pdb')
    working_dir = os.path.join(output_file_dir, '%s-prepped' % output_name)

    # Initialize forcefield with small molecule capabilities
    forcefield = ForceField('gaff.xml', 'tip3p.xml', 'amber99sbildn.xml')
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)

    # Use modeller to remove unwanted residues
    pdb = PDBFile(output_name)
    model = Modeller(pdb.topology, pdb.positions)

    # Remove unwanted molecules
    model = discard_organic(model, verbose=False)

    # Add waters in a cubic box
    if solvate == True:
        model.addSolvent(forcefield, padding=1.0 * nanometers)

    # Create the system with a cheap electrostatic cutoff
    system = forcefield.createSystem(model.topology, nonbondedMethod=CutoffNonPeriodic)

    # Minimize system with a placeholder integrator
    integrator = VerletIntegrator(0.001 * picoseconds)
    simulation = Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)
    simulation.minimizeEnergy()

    # Print PDB
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions,
                      open(os.path.join(output_file_dir, output_file_name + '-minimized.pdb'), 'w'))

if __name__ == '__main__':
    protein_prep(file_name, output, cap, pH=ph)
    openmm_clean(file_name, output, solvate=True)
