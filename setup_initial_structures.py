"""
Set up T4 lysozyme L99A apo simulations for Folding@home

Contributors:
* Julie M. Behr <julie.behr@choderalab.org>
* John D. Chodera <john.chodera@choderalab.org>
"""

#===============================================================================
# IMPORTS
#===============================================================================

import os, os.path
import copy
import numpy
import shutil
import tempfile
import numpy as np

import pdbfixer
from simtk import openmm, unit
from simtk.openmm import app
import mdtraj as md


#===============================================================================
# TEMPLATES
#===============================================================================

# All templates must have identical SEQRES entries!
template_pdbids = [
    '4W51', # L99A mutant without ligands present
    '4W53', # L99A mutant with a "closed" ligand present: toluene
    '4W56', # L99A mutant with an "intermediate" ligand present: sec-butylbenzene
    '4W59'  # L99A mutant with an "open" ligand present: n-hexylbenzene
    ]

# Path to put all output in
output_basepath = 'output'

#
# PARAMETERS
#

chain_id_to_mutate = 'A' # chain to mutate
pH = 7.4 # pH to model
keep_crystallographic_water = False # keep crystallographic waters?

# Single point mutants
point_mutants = ['Q185C', 'Q185L','Q185M','Q185N']

# Forcefield
ff_name = 'amber99sbildn'
water_name = 'tip3p'

solvate = True # if True, will add water molecules using simtk.openm.app.modeller
padding = 11.0 * unit.angstroms
nonbonded_cutoff = 9.0 * unit.angstroms
nonbonded_method = app.PME
max_minimization_iterations = 5000
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
collision_rate = 5.0 / unit.picoseconds
barostat_frequency = 50
timestep = 2.0 * unit.femtoseconds
nsteps = 5000 # number of steps to take for testing
ionicStrength = 300 * unit.millimolar

# Verbosity level
verbose = True

#===============================================================================
# DATA
#===============================================================================

three_letter_code = {
    'A' : 'ALA',
    'C' : 'CYS',
    'D' : 'ASP',
    'E' : 'GLU',
    'F' : 'PHE',
    'G' : 'GLY',
    'H' : 'HIS',
    'I' : 'ILE',
    'K' : 'LYS',
    'L' : 'LEU',
    'M' : 'MET',
    'N' : 'ASN',
    'P' : 'PRO',
    'Q' : 'GLN',
    'R' : 'ARG',
    'S' : 'SER',
    'T' : 'THR',
    'V' : 'VAL',
    'W' : 'TRP',
    'Y' : 'TYR'
}

one_letter_code = dict()
for one_letter in three_letter_code.keys():
    three_letter = three_letter_code[one_letter]
    one_letter_code[three_letter] = one_letter

def decompose_mutation(mutation):
    import re
    match = re.match('(\D)(\d+)(\D)', mutation)
    original_residue_name = three_letter_code[match.group(1)]
    residue_index = int(match.group(2))
    mutated_residue_name = three_letter_code[match.group(3)]
    return (original_residue_name, residue_index, mutated_residue_name)

def generate_pdbfixer_mutation_code(original_residue_name, residue_index, mutated_residue_name):
    return '%s-%d-%s' % (original_residue_name, residue_index, mutated_residue_name)

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

#
# Read reference PDB file to create a list of possible alterations.
#

nwaters_and_ions = None
for (run, pdbid) in enumerate(template_pdbids):
    print('Processing %s (RUN%d)...' % (pdbid, run))

    # Create output directory.
    output_path = os.path.join(output_basepath, 'RUN%d' % run)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    workdir = output_path

    # Note template
    write_file(os.path.join(workdir, 'template.txt'), pdbid + '\n')

    # Create PDBFixer for this pdbid
    print "Creating pdbfixer..."
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    #fixer.topology.createStandardBonds()
    print "    findMissingResidues..."
    fixer.missingResidues = {}
    print "    findNonstandardResidues..."
    fixer.findNonstandardResidues()
    print "    replaceNonstandardResidues..."
    fixer.replaceNonstandardResidues()
    print "    findMissingAtoms..."
    fixer.findMissingAtoms()
    print "    addMissingAtoms..."
    fixer.addMissingAtoms()
    print "    Removing all heterogens, including water..."
    fixer.removeHeterogens(False)
    print "    addingmissinghydrogens..."
    fixer.addMissingHydrogens(pH)

    # Create forcefield
    if verbose: print "Loading forcefield..."
    forcefield = app.ForceField(ff_name+'.xml',water_name+'.xml')

    if verbose: print "Creating Modeller object..."
    modeller = app.Modeller(fixer.topology, fixer.positions)

    # Solvate and create system.
    if verbose: print "Solvating with %s..." % water_name
    if nwaters_and_ions is None:
        # For first template, add solvent using padding argument
        modeller.addSolvent(forcefield, padding=padding, model=water_name, ionicStrength=ionicStrength)

        # Count waters.
        print "Renumbering and counting waters..."
        for chain in modeller.topology.chains(): water_id = chain.id
        nwaters_and_ions = 0
        for residue in modeller.topology.residues():
            if residue.chain.id == water_id:
                residue.id = (nwaters_and_ions % 9999) + 1
                nwaters_and_ions += 1
        print "System contains %d waters and ions." % nwaters_and_ions
    else:
        # For later templates, add solvent using specified number
        modeller.addSolvent(forcefield, numAdded=nwaters_and_ions, model=water_name, ionicStrength=ionicStrength)

        # Count waters.
        print "Renumbering and counting waters..."
        for chain in modeller.topology.chains(): water_id = chain.id
        nwaters_and_ions = 0
        for residue in modeller.topology.residues():
            if residue.chain.id == water_id:
                residue.id = (nwaters_and_ions % 9999) + 1
                nwaters_and_ions += 1
        print "System contains %d waters and ions." % nwaters_and_ions


    # Convert positions to numpy format and remove offset
    if verbose: print "Subtracting offset..."
    modeller.positions = unit.Quantity(np.array(modeller.positions / unit.nanometer), unit.nanometer)
    offset = modeller.positions[0,:]
    if verbose: print "Shifting model by %s" % str(offset)
    modeller.positions -= offset

    # Write PDB file for solute only.
    if verbose: print "Writing pdbfixer output..."
    pdb_filename = os.path.join(workdir, 'pdbfixer.pdb')
    outfile = open(pdb_filename, 'w')
    app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile, keepIds=True)
    outfile.close()

    # Create OpenMM system.
    if verbose: print "Creating OpenMM system..."
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
    if verbose: print "Adding barostat..."
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

    # Create simulation.
    if verbose: print "Creating solvated simulation..."
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    #platform = openmm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Write modeller positions.
    if verbose: print "Writing modeller output..."
    filename = os.path.join(workdir, 'modeller_solvent.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

    # Minimize energy.
    if verbose: print "Minimizing energy..."
    potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
        raise Exception("Potential energy is NaN before minimization.")
    if verbose: print "Initial solvated potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)
    simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
    potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
        raise Exception("Potential energy is NaN after minimization.")
    if verbose: print "Final solvated potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)

    # Write initial positions.
    filename = os.path.join(workdir, 'minimized.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

    # Assign temperature
    simulation.context.setVelocitiesToTemperature(temperature)

    # Take a few steps to relax structure.
    if verbose: print "Taking a few steps..."
    simulation.step(nsteps)

    # Write initial positions.
    if verbose: print "Writing positions..."
    filename = os.path.join(workdir, 'system.pdb')
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

    # Serialize to XML files.
    if verbose: print "Serializing to XML..."
    system_filename = os.path.join(workdir, 'system.xml')
    integrator_filename = os.path.join(workdir, 'integrator.xml')
    write_file(system_filename, openmm.XmlSerializer.serialize(system))
    write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
    simulation.context.setVelocitiesToTemperature(temperature)
    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    state_filename = os.path.join(workdir, 'state.xml')
    serialized = openmm.XmlSerializer.serialize(state)
    write_file(state_filename, serialized)

    # Clean up.
    del simulation.context
    del simulation
    del system
    del positions
