"""
General setup pipeline for F@H.

@author Steven K. Albanese, John D. Chodera
@
@date 15 Nov 2016
"""

#
# IMPORTS
#

import os, os.path
import copy
import numpy
import shutil
import tempfile
from openmoltools.forcefield_generators import gaffTemplateGenerator
from pdbfixer import PDBFixer
from simtk import openmm, unit
from simtk.openmm import app
import urllib
import urllib.request
import os
from natsort import natsorted
import pypdb
import protprep

#
# PARAMETERS FROM COMMAND LINE
#

import argparse

parser = argparse.ArgumentParser(description='Solvating and equillibrating PDB')
parser.add_argument('--source_pdbid', dest='source_pdbfile', action='store', default=None,
                    help='the four letter pdbid of interest')
parser.add_argument('--output_directory', dest='output_directory', action='store', default=None,
                    help='output directory')
parser.add_argument('--run', dest='run_number', action='store', required=False, default=0,
                    help='number of run being set up. Default is 0. For use for multiple run projects')
parser.add_argument('--fix', required=False, action='store_true', dest='fix',
                    help='Setting flag fixes problems with the PDB')
parser.add_argument('--sch', required=False, action='store_true', dest='sch',
                    help='Uses Schrodinger PrepWizard to fix pdb problems. Must also use --fix')
parser.add_argument('--ph', required=False, default=7.4, type=float, dest='ph',
                    help='Use to set pH to something other than 7.4')
parser.add_argument('--biological_unit', required=False, action='store_true', dest='biological_unit',
                    help='Set flag to retrieve biological unit for all structures')
args = parser.parse_args()

if args.source_pdbfile==None or args.output_directory==None:
    parser.print_help()
    raise Exception("All arguments must be specified.")

run = args.run_number
pH = args.ph
fix = args.fix
use_schrodinger = args.sch
bunit = args.biological_unit

print("Source PDB filename: %s" % args.source_pdbfile)
print("Output directory: %s" % args.output_directory)

#
# PARAMETERS
#

# Path to put all output in
output_path = args.output_directory

# Source PDB
pdbfilename = args.source_pdbfile # This is just the 4 character ID
keep_crystallographic_water = False # keep crystallographic waters?
# Forcefield
ff_name = 'amber99sbildn'
water_name = 'tip3p'
solvate = True # if True, will add water molecules using simtk.openm.app.modeller
padding = 10.0 * unit.angstroms
nonbonded_cutoff = 10.0 * unit.angstroms
nonbonded_method = app.CutoffPeriodic
max_minimization_iterations = 5000
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
collision_rate = 5.0 / unit.picoseconds
barostat_frequency = 50
timestep = 2.0 * unit.femtoseconds
nsteps = 5000 # number of steps to take for testing
ionicStrength = 150 * unit.millimolar

# Verbosity level
verbose = True
trim = True

########################
# Function Definitions #
########################

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)


def get_pdb_biological_unit(pdb_id):
    """

    Args:
        pdb_id: 4-letter PDB code

    Returns: a string with the full

    """


    fullurl = 'http://www.rcsb.org/pdb/files/'+pdb_id+'.pdb1'
    req = urllib.request.Request(fullurl)
    f = urllib.request.urlopen(req)
    result = f.read()
    result = result.decode('unicode_escape')
    result = result.replace('XXXX', pdb_id)

    return result


def pdb_fix_schrodinger(pdbid, file_pathway, ph):

    print(pdbid)

    input_file = os.path.join(file_pathway, '%s.pdb' % pdbid)
    output_file_pathway = file_pathway
    protprep.protein_prep(input_file, output_file_pathway, pdbid, pH=ph)


def download_pdb(pdbid, file_pathway):

    """

    Args:
        pdbid: 4 letter string specifying the PDB ID of the file yoou want to fix
        file_pathway: a string containing the pathway specifying how you want to organize the PDB files once written

    Returns: nothing, but it does write the PDB file

    ***Note: this function does NOT fix any mistakes with the PDB file

    """

    if not os.path.exists(file_pathway):
        os.makedirs(file_pathway)

    if bunit == True:
        pdb = get_pdb_biological_unit(pdbid)

    else:
        pdb = pypdb.get_pdb_file(pdbid)

    write_file(os.path.join(file_pathway, '%s.pdb' % pdbid), pdb)


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
    unwanted = ['DTT', 'EDO', 'GOL', 'SO4', 'PO4', 'DMS', 'MAL']
    for junk in unwanted:
        atms = [atm for atm in model.topology.atoms() if atm.residue. == junk]
        if len(atms) > 0:
            if verbose == True: print('Deleting {0} from topology'.format(junk))
            model.delete(atms)
    return model

def pdb_fix_pdbfixer(pdbid, file_pathway, ph):
    """
    Args:
        pdbid: 4 letter string specifying the PDB ID of the file yoou want to fix
        file_pathway: a string containing the pathway specifying how you want to organize the PDB files once written
        ph: the pH at which hydrogens will be determined and added
        chains_to_remove: dictionary containing pdbs with chains to remove
    Returns: nothing, but it does right PDB files
    """
    print(pdbid)

    # Download the topology from rcsb based on pdbod
    fixer = PDBFixer(pdbid=pdbid)

    # Determine the first and last residue resolved in chain 0
    chains = [chain for chain in fixer.topology.chains()]
    resindices = [residue.index for residue in chains[0].residues()]
    resindices = natsorted(resindices)
    first_resindex = resindices[0]
    last_resindex = resindices[-1]

    # Find Missing residues and determine if they are C or N terminal fragments (which will be removed)

    fixer.findMissingResidues()
    if len(fixer.missingResidues) >0:
        if sorted(fixer.missingResidues.keys())[0][-1] <= first_resindex:
            fixer.missingResidues.pop((sorted(fixer.missingResidues.keys())[0]))

        if sorted(fixer.missingResidues.keys())[-1][-1] >= last_resindex:
            fixer.missingResidues.pop((sorted(fixer.missingResidues.keys())[-1]))

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    # Write fixed PDB file, with all of the waters and ligands
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                            '%s_fixed.pdb' % (pdbid, ph)), 'w'), keepIds=keepNumbers)

    # Remove the ligand and write a pdb file
    fixer.removeHeterogens(True)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                            '%s_fixed_ph%s_apo.pdb' % (pdbid, ph)), 'w'), keepIds=keepNumbers)
    # Remove the waters and write a pdb file
    fixer.removeHeterogens(False)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(file_pathway,
                            '%s_fixed_ph%s_apo_nowater.pdb' % (pdbid, ph)), 'w'), keepIds=keepNumbers)


#############
# FILENAMES #
#############

exception_filename = os.path.join(output_path, 'exceptions.out') # to store exceptions
name = pdbfilename


#################################
#          MAIN                 #
#################################

# Create output directory.
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Create temporary directory.
tmp_path = tempfile.mkdtemp()
print("Working in temporary directory: %s" % tmp_path)

# Open file to write all exceptions that occur during execution.
exception_outfile = open(exception_filename, 'a')

simulation=None

# Create directory to store files in.
workdir = os.path.join(tmp_path, name)
if not os.path.exists(workdir):
    os.makedirs(workdir)
    print("Creating path %s" % workdir)

# retrieve PDB from code
if verbose: print("Downloading PDB file...")
download_pdb(pdbfilename, output_path)

# Write PDB file for fixed.
if verbose: print("Writing ProteinPrep output...")
if fix:
    if use_schrodinger:
        pdb_fix_schrodinger(pdbfilename, output_path, pH)
    else:
        pdb_fix_pdbfixer(pdbfilename, output_path, pH)

pdb = app.PDBFile(os.path.join(output_path, '%s-fixed.pdb' % pdbfilename))

# Solvate
if verbose: print("Loading forcefield...")
forcefield = app.ForceField('gaff.xml', ff_name+'.xml', water_name+'.xml')
forcefield.registerTemplateGenerator(gaffTemplateGenerator)

if verbose: print("Creating Modeller object...")
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller = discard_organic(modeller)
num_chains = modeller.topology.getNumChains()

if solvate:
    # Solvate and create system.
    if verbose: print("Solvating with %s..." % water_name)
    modeller.addSolvent(forcefield, padding=padding, model=water_name, ionicStrength=ionicStrength)

    # Separate waters into a separate 'W' chain with renumbered residues.
    # This is kind of a hack, as waters are numbered 1-9999 and then we repeat
    print("Renumbering waters...")
    for chain in modeller.topology.chains():
        if chain.id == str(num_chains + 1): chain.id = 'W'
    nwaters_and_ions = 0
    for residue in modeller.topology.residues():
        if residue.chain.id == 'W':
            residue.id = (nwaters_and_ions % 9999) + 1
            nwaters_and_ions += 1
    print("System contains %d waters and ions." % nwaters_and_ions)

# Create OpenMM system.
if verbose: print("Creating OpenMM system...")
system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
if verbose: print("Adding barostat...")
system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

# Create simulation.
if verbose: print("Creating simulation...")
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Write modeller positions.
if verbose: print("Writing modeller output...")
filename = os.path.join(workdir, 'modeller.pdb')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

# Minimize energy.
if verbose: print("Minimizing energy...")
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN before minimization.")
if verbose: print("Initial potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))
simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN after minimization.")
if verbose: print("Final potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))

# Write initial positions.
filename = os.path.join(workdir, 'minimized.pdb')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

# Assign temperature
simulation.context.setVelocitiesToTemperature(temperature)

# Take a few steps to relax structure.
if verbose: print("Taking a few steps...")
simulation.step(nsteps)

# Write initial positions.
if verbose: print("Writing positions...")
filename = os.path.join(workdir, 'system.pdb')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)


# Serialize to XML files.
if verbose: print("Serializing to XML...")
system_filename = os.path.join(workdir, 'system.xml')
integrator_filename = os.path.join(workdir, 'integrator.xml')
write_file(system_filename, openmm.XmlSerializer.serialize(system))
write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
simulation.context.setVelocitiesToTemperature(temperature)
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
state_filename = os.path.join(workdir, 'state.xml')
serialized = openmm.XmlSerializer.serialize(state)
write_file(state_filename, serialized)

# If everything worked, add this RUN.
run_name = 'RUN%s' % run
run_dir = os.path.join(output_path, run_name)
shutil.move(workdir, run_dir)

# Clean up.
del simulation.context
del simulation
del system
del positions

exception_outfile.close()
