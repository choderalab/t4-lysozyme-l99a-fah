"""
Test that simulations can be resumed from serialized XML files.

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
import glob
import numpy as np

from simtk import openmm, unit
from simtk.openmm import app

#===============================================================================
# TEMPLATES
#===============================================================================

# Path to put all output in
output_basepath = 'output'

#
# PARAMETERS
#

nsteps = 500 # number of steps to simulate for testing

# Verbosity level
verbose = True

#===============================================================================
# SUBROUTINES
#===============================================================================

def read_file(filename):
    infile = open(filename, 'r')
    contents = infile.read()
    infile.close()
    return contents

#===============================================================================
# MAIN
#===============================================================================

#
# Attempt to resume from each RUN
#

runs = glob.glob(os.path.join(output_basepath, 'RUN*'))

for run in runs:
    workdir = run
    print('Testing %s: %s' % (run, workdir))

    # Deserialize from XML files.
    system_filename = os.path.join(workdir, 'system.xml')
    integrator_filename = os.path.join(workdir, 'integrator.xml')
    state_filename = os.path.join(workdir, 'state.xml')

    print('  Deserializing System...')
    system = openmm.XmlSerializer.deserialize(read_file(system_filename))
    print('  Deserializing Integrator...')
    integrator = openmm.XmlSerializer.deserialize(read_file(integrator_filename))
    print('  Deserializing State...')
    state = openmm.XmlSerializer.deserialize(read_file(state_filename))

    # Simulate
    context = openmm.Context(system, integrator)
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    context.setVelocities(state.getVelocities())
    potential = context.getState(getEnergy=True).getPotentialEnergy()
    print('  Initial potential energy is %.3f kcal/mol' % (potential / unit.kilocalories_per_mole))
    print('  Simulating for %d steps' % nsteps)
    integrator.step(nsteps)
    potential = context.getState(getEnergy=True).getPotentialEnergy()
    print('  Final potential energy is %.3f kcal/mol' % (potential / unit.kilocalories_per_mole))

    # Clean up.
    del context
