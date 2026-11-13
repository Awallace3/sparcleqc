import re
from typing import Tuple

import numpy as np
import pandas as pd
import parmed as pmd
from sparcleqc.combine_data import mol2_atom_metadata, mol2_atom_records


def id_to_coords(pdb_file1: str, atom_id1: str) -> Tuple[str, str, str]:
    """
    extracts the x,y,z coordinates of a specified atom in a specified pdb

    Parameters
    ----------
    pdb_file: str
        path to pdb
    
    Returns
    -------
    x_coord, y_coord, z_coord: Tuple[str, str, str]
        coordinates of specified atom
    """

    with open(pdb_file1, 'r') as file:
        lines = file.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                if l[6:11].strip() == atom_id1:
                    x_coord = l[29:38].strip()
                    y_coord = l[38:46].strip()
                    z_coord = l[46:54].strip()
                    return x_coord, y_coord, z_coord

def coords_to_new_id(x:str, y:str, z:str, pdb_file2:str) -> str:
    """
    returns the atom number that corresponds to the specified x,y,z
    coordinates in the specified pdb_file

    Parameters
    ----------
    pdb_file: str
        path to uncapped complex pdb

    Returns
    -------
    atom_id: str
        atom id that corresponds to the specified coordinates in the specified pdb file
    """
    with open(pdb_file2, 'r') as file:
        lines = file.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
               x_coord = l[29:38].strip()
               y_coord = l[38:46].strip()
               z_coord = l[46:54].strip()
               if x_coord == x and y_coord == y and z_coord == z:
                   return l[6:11].strip()

def convert_seed(seed: str, seed_file:str, new_pdb:str ='cx_autocap_fixed.pdb') -> str:
    """
    finds the atom the corresponds to atom number = seed in seed_file
    pdb and then maps this to the newly created complex pdb

    Parameters
    ----------
    seed: str
        atom number in seed_file
    seed_file: str
        file to original pdb that contains seed
    new_pdb: str
        path to new pdb

    Returns
    -------
    cx_id: str
        atom in new_pdb that corresponds to the atom seed in seed_file
    coords: list
        coordinates of seed atom
    """
    x1,y1,z1 = id_to_coords(seed_file, seed)
    cx_id = coords_to_new_id(x1, y1, z1, new_pdb)
    coords = [float(x1),float(y1),float(z1)]
    return cx_id, coords

def check_mol2_charges(mol2_file: str) -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the
    mol2 and ensures that it sums to an integer

    Parameters
    ----------
    mol2_file: str
        path to mol2 file

    Returns
    -------
    error_message: Tuple[int, str]
        tuple where the first entry is an integer. 
        0 = at least one residue is not integer charge. 
        1 = all residues are integer charge.
        The second entry is a string containing information about which 
        residue is not integer charge and the corresponding fractional charge.
    """
    #tolerance is how far away a residue can be from an integer without throwing an error
    #this was determined because amber lists charge to 3 decimal places, 
    #so anything different to the fourth decimal place or further is insignificant
    tolerance = 0.0001
    # Residue-charge validation should see the same logical atom records as the dataframe builder.
    charges = {}
    for fields in mol2_atom_records(mol2_file):
        _, residue_number, residue_name, charge = mol2_atom_metadata(fields)
        residue_id = f"{residue_name}_{residue_number}"
        charges[residue_id] = charges.get(residue_id, 0.0) + charge
    if not charges:
        return (0, 'No atom records found in mol2. Check Amber topology generation and restart')
    # Include the last residue (and single-residue systems), not only boundaries.
    for residue_id, charge in charges.items():
        if not np.isfinite(charge) or np.abs(np.round(charge, 0) - charge) > tolerance:
            return (0, f'Residue {residue_id} charge: {charge}. Fix pdb/mol2 such that this residue has integer charge and restart')
    return (1, 'passed')

def check_df_charges() -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the
    dataframe and ensures that it sums to an integer

    Parameters
    ----------
    None

    Returns
    -------
    return_message: Tuple[int, str]
        tuple where the first entry is an integer. 
        0 = at least one residue is not integer charge. 
        1 = all residues are integer charge.
        The second entry is a string containing information about which 
        residue is not integer charge and the corresponding fractional charge.
    """
    
    tolerance = 0.0001
    df = pd.read_csv('dataframe.csv', index_col='PDB_ID')
    # Blank PDB labels and inserted water extra points round-trip through CSV
    # as NaN. Normalize them before string operations, without writing "nan".
    df['AT_LABEL'] = df['AT_LABEL'].astype('string').fillna('')
    if df['q'].isna().any() or not np.isfinite(df['q']).all():
        return (0, 'Missing or non-finite partial charges in dataframe.csv. Check PDB/topology atom mapping and restart')

    for residue, atoms in df.groupby('PDB_RES', sort=False):
        charge = atoms['q'].sum()
        if np.abs(np.round(charge) - charge) > tolerance:
            return (0, f'Residue {residue} charge: {charge}. Fix pdb/mol2 such that this residue has integer charge and restart')
        if ('WAT' in residue or 'HOH' in residue) and abs(charge) > tolerance:
            return (0, f'Water charge ({residue}) does not equal zero. Please choose a new water model or add custom charges in the input. See dataframe.csv for current charges.')
        int_charge = int(np.round(charge))
        if int_charge == 0:
            continue
        suffix = f"{abs(int_charge)}{'+' if int_charge > 0 else '-'}"
        labelled = [
            (idx, re.search(r'\d+[+-]$', label))
            for idx, label in atoms['AT_LABEL'].items()
        ]
        labelled = [(idx, match) for idx, match in labelled if match is not None]
        if len(labelled) == 1 and labelled[0][1].group() == suffix:
            continue
        # Reconcile conflicting/multiple labels once per residue. The old
        # df_charge variable leaked across residues and left later ones unlabelled.
        if labelled:
            print(
                f'PDB formal charge labels for {residue} do not match the sum '
                f'of partial charges ({suffix}). Replacing them with one residue charge label.'
            )
        for idx, match in labelled:
            label = df.at[idx, 'AT_LABEL']
            df.at[idx, 'AT_LABEL'] = label[:match.start()]
        cb_indices = atoms.index[atoms['PDB_AT'] == 'CB']
        target = cb_indices[0] if len(cb_indices) else atoms.index[0]
        df.at[target, 'AT_LABEL'] += suffix

    df.to_csv('dataframe.csv')
    return (1, 'worked')

def closest_contact(pdb1, pdb2):
    """
    Calculate the distance of the closest contact between two PDBs, or one PDB and one coordinate.
    
    Args:
    pdb1 (np.ndarray): A PDB representing the first set of XYZ coordinates.
    pdb2 (np.ndarray): A PDB representing the second set of XYZ coordinates or list of coords of one atom (seed).
    
    Returns:
    float: The minimum distance between any two points from pdb1 and pdb2.
    """
    pmdpdb1 = pmd.load_file(pdb1)
    coords1 = pmdpdb1.coordinates
    if 'pdb' in pdb2:
        pmdpdb2 = pmd.load_file(pdb2)
        coords2 = pmdpdb2.coordinates
    else:
        coords2 = pdb2

    # Ensure the input arrays are NumPy arrays
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # Compute the distance matrix
    distances = np.linalg.norm(coords1[:, np.newaxis] - coords2, axis=2)

    # Find the minimum distance
    min_distance = np.min(distances)

    return min_distance
