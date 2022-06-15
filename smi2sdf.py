"""
Conversion of SMILES 2D representation to 3D Coordinate generation

This program uses RDKits EmbedMultipleConfs from AllChem to generate 3D
structures of molecules based on an input SMILES string.

The code is based on GB-GA (https://github.com/cstein/GB-GA) which is
a fork of Jan Jensen's GB-GA (https://github.com/jensengroup/GB-GA).
"""
import argparse
import errno
import multiprocessing as mp
import os
import random
import string
import subprocess
import sys
from typing import Optional, List

from rdkit import Chem
from rdkit.Chem import AllChem

has_protonator = False
try:
    from protonator import protonator
except ImportError:
    print("disabling protonation states for smi2sdf")
    protonator = None
else:
    has_protonator = True


def safe_create_dir(path: str):
    """ Creates a directory safely by not raising an error if it already exists

        :param path: the pathname to create
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def read_smi_file(filename: str, i_from: int, i_to: int) -> List[Chem.Mol]:
    """ Reads a file with SMILES between two indices

        :param filename: the filename to read SMILES from
        :param i_from: where to search from in the file
        :param i_to: where to search to (but not included)
    """
    mol_list: List[Chem.Mol] = []
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
        pass
    else:
        if p.returncode > 0:
            print("shell:", p)
            raise ValueError("Error with structure generation.")


def get_structure(mol: Chem.Mol, num_conformations: int, index: int) -> Optional[Chem.Mol]:
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

    if has_protonator:
        mol = protonator(mol)

    try:
        mol = Chem.AddHs(mol)
    except ValueError as e:
        print("get_structure: could not kekulize the molecule '{}'".format(s_mol))
        return None

    new_mol = Chem.Mol(mol)

    conformer_energies = []
    try:
        AllChem.EmbedMultipleConfs(mol, numConfs=num_conformations, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        conformer_energies = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000, nonBondedThresh=100.0)
    except ValueError:
        print("get_structure: '{}' could not converted to 3D".format(s_mol))
        return None
    # print(conformer_energies)
    # print("index on input:", index)
    if index == 0:
        try:
            i = conformer_energies.index(min(conformer_energies))
        except ValueError:
            print("get_structure: '{}' no energies for conformers".format(s_mol))
            return None
    elif index > 0:
        i = index - 1
    else:
        raise ValueError("index cannot be less that zero.")
    new_mol.AddConformer(mol.GetConformer(i))
    return new_mol


def molecules_to_structure(population: List[Chem.Mol],
                           num_conformations: int,
                           index: int,
                           num_cpus: int):
    """ Converts RDKit molecules to structures

        :param population: the entire set of molecules to obtain 3D structures for
        :param num_conformations: the number of conformations to generate for each molecule
        :param index: the index to return for each structure. 0 is minimum. > 0 is an index
        :param num_cpus: allocate a number of CPUs to parallelize the task
    """

    with mp.Pool(num_cpus) as pool:
        args = [(p, num_conformations, index) for p in population]
        generated_molecules = pool.starmap(get_structure, args)
   
        # molecules = [mol for mol in generated_molecules if mol is not None]
        names = [''.join(random.choices(string.ascii_uppercase + string.digits, k=6)) for pop in generated_molecules]
        updated_population = [p for (p, m) in zip(population, generated_molecules) if m is not None]

        return generated_molecules, names, updated_population


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
    ap = argparse.ArgumentParser()
    ap.add_argument("filename", type=str, metavar="file")
    ap.add_argument("-i", "--index", dest="index", type=int, metavar="number", default=1,
                    help="index to start from. Default is %(default)s.")
    ap.add_argument("-w", "--width", dest="width", type=int, metavar="number", default=1,
                    help="number of entries to parse. Default is %(default)s.")
    ap.add_argument("--cpus", dest="num_cpus", type=int, default=1, metavar="number",
                    help="number of CPUs to use in parallel. Default is %(default)s.")
    ap.add_argument("--confs", dest="num_confs", type=int, default=5, metavar="number",
                    help="number of conformations to generate. Sorted by energy. Default is %(default)s.")
    ap.add_argument("-c", dest="index_conformer", type=int, default=0, metavar="number",
                    help="index of generated conformer to write to file. A value of 0 means minimum energy. Default is %(default)s.")
    if has_protonator:
        ap.add_argument("-p", dest="protonator", default=False, action="store_true",
                        help="add this flag to give the most probably protonation states for each molecule")
    args = ap.parse_args()
    filename = args.filename
    index = args.index
    width = args.width
    num_cpus = args.num_cpus
    num_confs = args.num_confs
    idx_conformer = args.index_conformer
    print(args)
    i_from = index-1
    i_to = i_from + width
    initial_population = read_smi_file(filename, i_from, i_to)
    pop_names = ["{0:06d}".format(i) for i in range(i_from+1, i_to+1)]
    basename, _ = os.path.splitext(filename)
    wrk_dir = basename
    safe_create_dir(wrk_dir)

    # change to work directory
    os.chdir(wrk_dir)
    molecules, _, _ = molecules_to_structure(initial_population, num_confs, idx_conformer, num_cpus)
    filenames = ["{0:s}.sdf".format(s) for s in pop_names]
    for molecule, filename, mol_name in zip(molecules, filenames, pop_names):
        if molecule is not None:
            molecule_to_sdf(molecule, filename, name=mol_name)

    # go back from work directory
    os.chdir("..")

