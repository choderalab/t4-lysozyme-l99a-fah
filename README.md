# T4 lysozyme L99A apo simulations for Folding@home

## Contributors
* Mehtap Isik `<mehtap.isik@choderalab.org>`
* John D. Chodera `<john.chodera@choderalab.org>`
* Samuel Gill `<sgill2@uci.edu>`
* David L. Mobley `<dmobley@uci.edu>`

## Description

This project sets up initial structures of T4 lysozyme L99A apo protein for simulation on Folding@home to construct Markov state models.

Initial structures:
* `4W51`                                    (L99A mutant without ligands present)
* `4W53`     Ligand Chain name: `MBN`       (L99A mutant with a "closed" ligand present: toluene)
* `4W56`     Ligand Chain name: `3GY`       (L99A mutant with an "intermediate" ligand present: sec-butylbenzene)
* `4W59`     Ligand Chain name: `3GZ`       (L99A mutant with an "open" ligand present: n-hexylbenzene)

There are a handful of missing heavy atoms that slightly vary between structures; otherwise the set of pdb structures is well defined and vary by the bound ligand and F helix conformation. The SEQRES blocks are the same. Also, residues 165-172 don't have atom data, but those residues correspond to a HIS tag so they shouldn't be included.

Forcefield: AMBER99SB-ildn

## Usage

### Setting up models

Modeling of templates must use openmm >=7.0.1
```bash
conda install --yes openmm>=7.0.1 mdtraj pdbfixer
python setup_initial_structures.py
```
Test that this works with OpenMM 6.2 (which Folding@home core 21 v0.0.17 is built on)
```bash
conda install --yes openmm==6.2
python test_resume.py
```
Don't forget to go back to a recent version of `openmm` afterwards!

You may want to use `conda env` to create conda environments for OpenMM 6.2 and 7.0.1 to simplify things!
http://conda.pydata.org/docs/using/envs.html
