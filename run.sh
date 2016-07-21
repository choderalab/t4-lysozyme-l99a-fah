#!/bin/bash

echo "Setting up initial structures..."
conda create -n model openmm>=7.0.1 mdtraj pdbfixer
source activate model
python setup_initial_structures.py
deactivate

echo "Packaging for FAH..."
conda create -n package openmm==6.3
surce activate package
python package_for_fah.py
deactivate
