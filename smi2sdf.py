"""
Conversion of SMILES 2D representation to 3D Coordinate generation

This program uses RDKits EmbedMultipleConfs from AllChem to generate 3D
structures of molecules based on an input SMILES string.

The code is based on GB-GA (https://github.com/cstein/GB-GA) which is
a fork of Jan Jensen's GB-GA (https://github.com/jensengroup/GB-GA).
"""
import errno
import multiprocessing as mp
import os
import random
import string
import subprocess
import sys
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem


def safe_create_dir(path: str):
    """ Creates a directory safely by not raising an error if it already exists

        :param path: the pathname to create
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def read_smi_file(filename: str, i_from: int, i_to: int):
    """ Reads a file with SMILES between two indices

        :param filename: the filename to read SMILES from
        :param i_from: where to search from in the file
        :param i_to: where to search to (but not included)
    """
    mol_list = []
    with open(filename, 'r') as smiles_file:
        for i, line in enumerate(smiles_file):
            if i_from <= i < i_to:
                tokens = line.split()
                smiles = tokens[0]
                mol_list.append(Chem.MolFromSmiles(smiles))

    return mol_list


def shell(cmd: str):
    """ Executes the command given on input through the current shell

        :param cmd: the command to execute
    """
    try:
        p = subprocess.run(cmd, capture_output=True, shell=True)
    except AttributeError:
        # p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # output, err = p.communicate()
        pass
    else:
        if p.returncode > 0:
            print("shell:", p)
            raise ValueError("Error with Docking.")


def get_structure(mol, num_conformations):
    """ Converts an RDKit molecule (2D representation) to a 3D representation

    :param Chem.Mol mol: the RDKit molecule
    :param int num_conformations:
    :return: an RDKit molecule with 3D structure information
    """
    try:
        s_mol = Chem.MolToSmiles(mol)
    except ValueError:
        print("get_structure: could not convert molecule to SMILES")
        return None

    try:
        mol = Chem.AddHs(mol)
    except ValueError as e:
        print("get_structure: could not kekulize the molecule '{}'".format(s_mol))
        return None

    new_mol = Chem.Mol(mol)

    try:
        if num_conformations > 0:
            AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            conformer_energies = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000, nonBondedThresh=100.0)
            energies = [e[1] for e in conformer_energies]
            min_energy_index = energies.index(min(energies))
            new_mol.AddConformer(mol.GetConformer(min_energy_index))
        else:
            AllChem.EmbedMolecule(new_mol)
            AllChem.MMFFOptimizeMolecule(new_mol)
    except ValueError:
        print("get_structure: '{}' could not converted to 3D".format(s_mol))
        new_mol = None
    finally:
        return new_mol


def choices(sin, nin=6):
    result = []
    try:
        result = random.choices(sin, k=nin)
    except AttributeError:
        for i in range(nin):
            result.append( random.choice(sin) )
    finally:
        return result


def molecules_to_structure(population, num_conformations, num_cpus):
    """ Converts RDKit molecules to structures

    """

    with mp.Pool(num_cpus) as pool:
        args = [(p, num_conformations) for p in population]
        generated_molecules = pool.starmap(get_structure, args)
   
        molecules = [mol for mol in generated_molecules if mol is not None]
        names = [''.join(choices(string.ascii_uppercase + string.digits, 6)) for pop in molecules]
        updated_population = [p for (p, m) in zip(population, generated_molecules) if m is not None]

        return molecules, names, updated_population


def molecule_to_sdf(mol: Chem.Mol, output_filename: str, name: Optional[str] = None):
    """ Saves an RDKit molecule to an SDF file

        Optionally writes an internal name to the file as well.

        :param mol: the RDKit molecule to save to file
        :param output_filename: the filename to save to
        :param name: optional internal name to use

    """
    if name is not None:
        mol.SetProp("_Name", name)
    Chem.SDWriter("{}".format(output_filename)).write(mol)


if __name__ == '__main__':
    filename = sys.argv[1]
    index = int(sys.argv[2])
    width = int(sys.argv[3])
    i_from = (index-1) * width
    i_to = i_from + width
    num_cpus = 1
    population = read_smi_file(filename, i_from, i_to)
    pop_names = ["{0:05d}".format(i) for i in range(i_from+1, i_to+1)]
    basename, _ = os.path.splitext(filename)
    wrk_dir = basename
    safe_create_dir(wrk_dir)

    # change to work directory
    os.chdir(wrk_dir)
    molecules, _, _ = molecules_to_structure(population, 5, num_cpus)
    filenames = ["{0:s}.sdf".format(s) for s in pop_names]
    for molecule, filename, mol_name in zip(molecules, filenames, pop_names):
        molecule_to_sdf(molecule, filename, name=mol_name)

    # go back from work directory
    os.chdir("..")
