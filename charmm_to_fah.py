"""
Set up AurKA spin probe simulations from charmm-gui output files (Solvator on the website)
@author John D. Chodera, Steven K. Albanese
@date 20 Mar 2017
"""

#
# IMPORTS
#
from __future__ import division, print_function

import os, os.path
import copy
import numpy
import shutil
import tempfile
import numpy as np
from math import pi

import pdbfixer
from simtk import openmm, unit
from simtk.openmm import app
import mdtraj as md
import sys
import argparse

# ParmEd Imports
#from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet
from simtk.openmm.app import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet

from parmed import unit as u ### ???

##########
# Parser #
##########
parser = argparse.ArgumentParser(description="Script to setup FAH projects from openMM ")
parser.add_argument('--input', required=True, dest='pdb_file',
                    help='THe pdb file that you would like to run through this program')
parser.add_argument('--output', required=True, dest='output',
                    help='the name of the output file')
parser.add_argument('--run', dest='run_number', action='store', required=False, default=0,
                    help='number of run being set up. Default is 0')
parser.add_argument('--id', dest='content', action='store', required=False, default='empty',
                    help='text to be written out in the file in run')
args = parser.parse_args()

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)


def fix_charmm_impropers(system):
    print('Repairing CHARMM impropers...')
    for (force_index, force) in enumerate(system.getForces()):
        if (force.__class__.__name__ in ['CustomTorsionForce']):
            print(force_index, force.__class__.__name__)
            energy_function  = 'k*dtheta_torus^2;'
            energy_function += 'dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi);'
            energy_function += 'dtheta = theta - theta0;'
            energy_function += 'pi = %f;' % pi
            force.setEnergyFunction(energy_function)

if __name__ == '__main__':

    print("Input PDB structure: %s" % args.pdb_file)
    runs = args.run_number

    output_path = args.output

    # Source PDB
    pdbfilename = 'step2_solvator.pdb'

    print("Source PDB filename: %s" % pdbfilename)
    print("Output directory for mutations: %s" % output_path)

    #
    # PARAMETERS: these should be included in the file dump from Charmm gui
    #

    universal_parameter_files = [
        'toppar/top_all36_prot.rtf',
        'toppar/par_all36_prot.prm',
        'toppar/toppar_all36_prot_retinol.str',
        'toppar/toppar_all36_na_rna_modified.str',
        'toppar/toppar_all36_prot_fluoro_alkanes.str',
        'toppar/toppar_all36_prot_na_combined.str',
        'toppar/toppar_water_ions.str',
        'toppar/toppar_all36_na_nad_ppi.str',
        'toppar/toppar_dum_noble_gases.str',
        'toppar/toppar_all36_label_spin.str',
        'toppar/par_all36_na.prm',
        'toppar/top_all36_na.rtf'
    ]

    charmm_parameter_files = [
        'step2_solvator.str',
        'step2.2_ions.prm',
        'step2.1_waterbox.prm',
        'step1_pdbreader.str',
        'step1_labelrot.str',
        'step1_pdbreader.str',
    ]

    #psf_file = 'step2_solvator.xplor.psf'
    psf_file = 'step2_solvator.psf'

    #
    # Minimization and low-temperature NVT relaxation (50 ps)
    #

    nonbonded_cutoff = 9.0 * unit.angstroms
    nonbonded_method = app.PME
    max_minimization_iterations = 10000
    temperature = 50.0 * unit.kelvin
    pressure = 1.0 * unit.atmospheres
    collision_rate = 90.0 / unit.picoseconds
    barostat_frequency = 50
    timestep = 1.0 * unit.femtoseconds
    nsteps = 50000
    ionicStrength = 20 * unit.millimolar

    # Verbosity level
    verbose = True

    def write_file(filename, contents):
        with open(filename, 'w') as outfile:
            outfile.write(contents)

    exception_filename = os.path.join(output_path, 'exceptions.out') # to store exceptions
    run_index_filename = os.path.join(output_path, 'run-index.txt') # to store index of which mutants are which

    # Create temporary directory.
    tmp_path = tempfile.mkdtemp()
    print("Working in temporary directory: %s" % tmp_path)
    tmp_path = output_path

    # Open file to write all exceptions that occur during execution.
    exception_outfile = open(exception_filename, 'a')
    run_index_outfile = open(run_index_filename, 'a')

    name = str(runs)
    simulation = None

    # Load the CHARMM files
    print('Loading CHARMM files...')
    param_files = universal_parameter_files + charmm_parameter_files
    params = CharmmParameterSet(*param_files)
    psf = CharmmPsfFile(psf_file)
    pdb = app.PDBFile(pdbfilename)
    system_coords = CharmmCrdFile('step2_solvator.crd')

    coords = system_coords.positions
    min_crds = [coords[0][0], coords[0][1], coords[0][2]]
    max_crds = [coords[0][0], coords[0][1], coords[0][2]]

    for coord in coords:
        min_crds[0] = min(min_crds[0], coord[0])
        min_crds[1] = min(min_crds[1], coord[1])
        min_crds[2] = min(min_crds[2], coord[2])
        max_crds[0] = max(max_crds[0], coord[0])
        max_crds[1] = max(max_crds[1], coord[1])
        max_crds[2] = max(max_crds[2], coord[2])

    psf.setBox(max_crds[0]-min_crds[0] + (3 * unit.angstrom),
               max_crds[1]-min_crds[1] + (3 * unit.angstrom),
               max_crds[2]-min_crds[2] + (3 * unit.angstrom),
               )

    # Create PDBFixer, retrieving PDB template
    print("creating Modeller...")
    modeller = app.Modeller(psf.topology, system_coords.positions)

    # Create directory to store files in.
    workdir = os.path.join(tmp_path, 'RUN'+name)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
        print("Creating path %s" % workdir)

    # Write PDB file for solute only.
    if verbose: print("Writing initial output...")
    pdb_filename = os.path.join(workdir, 'initial.pdb')
    outfile = open(pdb_filename, 'w')
    app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    outfile.close()

    # Create OpenMM system.
    if verbose: print("Creating OpenMM system for initial minimization")
    system = psf.createSystem(params, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds, flexibleConstraints=False, verbose=True, removeCMMotion=False)

    # Fix CHARMM impropers
    fix_charmm_impropers(system)

    # Set properties for default platform
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator)
    platform = context.getPlatform()
    del context, integrator
    print('Using platform %s' % platform.getName())
    try:
        platform.setPropertyDefaultValue('Precision', 'mixed') # use double precision
        print('Using mixed precision')
    except:
        pass

    # Create simulation for NVT minimization and equilibration.
    if verbose: print("Creating NVT system for equillibration")
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    simulation = app.Simulation(modeller.topology, system, integrator)
    try:
        simulation.context.setPositions(modeller.positions)
    except Exception as e:
        print(len(modeller.positions))
        print(simulation.context.getSystem().getNumParticles())
        raise(e)

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

    # Write modeller positions.
    if verbose: print("Writing modeller output...")
    filename = os.path.join(workdir, 'modeller.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    print(abs(positions / unit.nanometers).max())
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

    # Write initial positions.
    filename = os.path.join(workdir, 'minimized.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

    # Serialize System
    system_filename = os.path.join(workdir, 'system.xml')
    write_file(system_filename, openmm.XmlSerializer.serialize(system))

    # Take steps at NVT to equilibrated
    if verbose: print("Performing NVT equilibration")
    simulation.step(nsteps)
    
    # Write initial positions.
    if verbose: print("Writing positions...")
    filename = os.path.join(workdir, 'equilibrated-NVT.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

    # Update positions
    positions = simulation.context.getState(getPositions=True).getPositions()

    del (integrator)
    del (simulation.context)
    del (simulation)

    #
    # NPT box size equilibration (1 ns)
    #

    # Create NPT system and equilibrate
    if verbose: print("Creating NPT equillibration system now...")
    temperature = 300.0 * unit.kelvin
    pressure = 1.0 * unit.atmospheres
    collision_rate = 90.0 / unit.picoseconds
    barostat_frequency = 50
    timestep = 2 * unit.femtoseconds
    nsteps = 500000  # number of steps to take for testing

    # Change parameters in the integrator
    if verbose: print("Changing to NPT integrator ")
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency)
    barostat_index = system.addForce(barostat)

    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.step(nsteps)

    if verbose: print("Writing positions...")
    filename = os.path.join(workdir, 'equilibrated-NPT.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

    # Retrieve the periodic box vectors and positions
    v1, v2, v3 = simulation.context.getState().getPeriodicBoxVectors()
    system.setDefaultPeriodicBoxVectors(v1, v2, v3)
    positions = simulation.context.getState(getPositions=True).getPositions()

    del (integrator)
    del (simulation.context)
    del (simulation)

    #
    # NPT production test and FAH packaging (1 ns)
    #

    # Create production system and test
    if verbose: print("Creating production system now...")
    temperature = 300.0 * unit.kelvin
    collision_rate = 5.0 / unit.picoseconds
    timestep = 2 * unit.femtoseconds
    nsteps = 500000  # number of steps to take for testing

    # Change parameters in the integrator
    if verbose: print("Changing to production integrator ")
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(positions)
    if verbose: print("Taking a few steps to test the structure at production conditions")
    simulation.step(nsteps)

    # Serialize to XML files.
    if verbose: print("This system passed the test! Serializing to XML...")
    system_filename = os.path.join(workdir, 'system.xml')
    integrator_filename = os.path.join(workdir, 'integrator.xml')
    write_file(system_filename, openmm.XmlSerializer.serialize(system))
    write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    state_filename = os.path.join(workdir, 'state.xml')
    serialized = openmm.XmlSerializer.serialize(state)
    write_file(state_filename, serialized)

    # Write txt file
    text_filename = os.path.join(workdir, 'run-info.txt')
    write_file(text_filename, args.content)

    # If everything worked, add this RUN.
    run_name = 'RUN%s' % runs
    run_dir = os.path.join(output_path, run_name)
    shutil.move(workdir, run_dir)
    run_index_outfile.write('%s %s\n' % (run_name, name))
    run_index_outfile.flush()


    # Clean up.
    del simulation.context
    del simulation
    del system
    del positions


    exception_outfile.close()
    run_index_outfile.close()
