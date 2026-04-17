import pandas as pd
from typing import Dict


def _token_is_float(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


# Amber can wrap long mol2 atom records across two physical lines.
# Parse the ATOM block into logical records before downstream code indexes fields.
def mol2_atom_records(MOL2_PATH: str) -> list[list[str]]:
    records = []
    in_atom_block = False
    current_record = None
    with open(MOL2_PATH, 'r', encoding='latin1', errors='ignore') as molfile:
        for raw_line in molfile:
            if raw_line.startswith('@<TRIPOS>ATOM'):
                in_atom_block = True
                current_record = None
                continue
            if raw_line.startswith('@<TRIPOS>') and in_atom_block:
                break
            if not in_atom_block:
                continue
            fields = raw_line.split()
            if not fields:
                continue
            if len(fields) >= 8 and fields[0].isdigit() and _token_is_float(fields[2]) and _token_is_float(fields[3]) and _token_is_float(fields[4]):
                if current_record is not None:
                    records.append(current_record)
                current_record = fields
            elif current_record is not None:
                current_record.extend(fields)
        if current_record is not None:
            records.append(current_record)
    return records


# Wrapped mol2 records can carry status tokens after the charge field.
# Recover the charge from the trailing numeric token instead of fixed negative indexing.
def mol2_atom_charge(fields: list[str]) -> float:
    for token in reversed(fields[5:]):
        if _token_is_float(token):
            return float(token)
    raise ValueError(f"Could not locate atom charge in mol2 record: {' '.join(fields)}")


# Amber occasionally omits the atom-type token on capped atoms and can also spill status bits across lines.
# Infer the semantic mol2 fields from the logical record instead of relying on one fixed column layout.
def mol2_atom_metadata(fields: list[str]) -> tuple[str, str, str, float]:
    tail = fields[5:]
    atom_type = ''
    if tail and not tail[0].isdigit() and not _token_is_float(tail[0]):
        atom_type = tail[0]
        tail = tail[1:]
    subst_id = None
    subst_name = None
    for index, token in enumerate(tail):
        if token.isdigit() or (token.startswith('-') and token[1:].isdigit()):
            subst_id = token
            if index + 1 < len(tail):
                subst_name = tail[index + 1]
            break
    if subst_id is None or subst_name is None:
        raise ValueError(f"Could not locate residue information in mol2 record: {' '.join(fields)}")
    charge = None
    for token in reversed(tail):
        if _token_is_float(token):
            charge = float(token)
            break
    if charge is None:
        raise ValueError(f"Could not locate atom charge in mol2 record: {' '.join(fields)}")
    return atom_type, subst_id, subst_name, charge

def prot_pdb_to_df(PDB_PATH: str, d:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the pdb
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    PDB_PATH: str
        path to protein pdb to populate the dataframe
    d: Dict
        Dictionary that has one key, PDB_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the PDB

    """
    with open(PDB_PATH, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                d['PDB_ID'].append(l[6:11].strip())
        df = pd.DataFrame(d)
        df = df.set_index('PDB_ID')
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                df.loc[l[6:11].strip(), 'PDB_AT'] = l[11:16].strip()

                x_coord = l[29:38].strip()
                y_coord = l[38:46].strip()
                z_coord = l[46:54].strip()

                df.loc[l[6:11].strip(), 'PDB_RES'] = l[16:20].strip()+'_'+l[22:26].strip()
                df.loc[l[6:11].strip(), 'X'] = float(x_coord)
                df.loc[l[6:11].strip(), 'Y'] = float(y_coord)
                df.loc[l[6:11].strip(), 'Z'] = float(z_coord)
                df.loc[l[6:11].strip(), 'AT_LABEL'] = l[66:87].strip() 
    return df

def mol2_to_df(MOL2_PATH:str, m:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the mol2
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    MOL2_PATH: str
        path to protein mol2 to populate the dataframe
    d: Dict
        Dictionary that has one key, MOL2_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the MOL2

    """
    records = mol2_atom_records(MOL2_PATH)
    for fields in records:
        m['MOL2_ID'].append(fields[0])
    df = pd.DataFrame(m)
    df = df.set_index('MOL2_ID')
    for fields in records:
        atom_type, _, residue_name, charge = mol2_atom_metadata(fields)
        df.loc[fields[0], 'MOL2_AT'] = atom_type
        df.loc[fields[0], 'MOL2_RES'] = residue_name
        x_mol = float("{:.3f}".format(float(fields[2])))
        y_mol = float("{:.3f}".format(float(fields[3])))
        z_mol = float("{:.3f}".format(float(fields[4])))
        df.loc[fields[0], 'X'] = x_mol
        df.loc[fields[0], 'Y'] = y_mol
        df.loc[fields[0], 'Z'] = z_mol
        df.loc[fields[0], 'q'] = charge
    return df        

def combine_prot_dfs(pdb_df:pd.DataFrame, mol2_df:pd.DataFrame) -> pd.DataFrame:
    """ 
    combines information from a pdb dataframe with a
    mol2 dataframe based on matching the coordinates

    Parameters
    ----------
    pdb_df: pd.DataFrame
        dataframe containing information from the protein pdb
    mol2_df: pd.DataFrame
        dataframe containing information from the protein mol2

    Returns
    -------
    df: pd.DataFrame
        dataframe containing information combined from the pdb_df and the mol2_df with one entry per atom

    """
    for idx in pdb_df.index:
        mol2_idx_x = mol2_df.loc[mol2_df['X'] == pdb_df.loc[idx, 'X']].index.tolist()
        mol2_idx_y = mol2_df.loc[mol2_df['Y'] == pdb_df.loc[idx, 'Y']].index.tolist()
        mol2_idx_z = mol2_df.loc[mol2_df['Z'] == pdb_df.loc[idx, 'Z']].index.tolist()
        mol2_idx_xyz = mol2_idx_x + mol2_idx_y + mol2_idx_z
        for x in mol2_idx_xyz:
            if mol2_idx_xyz.count(x) == 3:
                mol2_idx = x
                pdb_df.loc[idx,'MOL2_ID'] = mol2_idx
                pdb_df.loc[idx,'MOL2_AT'] = mol2_df.loc[mol2_idx,'MOL2_AT']
                pdb_df.loc[idx,'MOL2_RES'] = mol2_df.loc[mol2_idx,'MOL2_RES']
                pdb_df.loc[idx,'q'] = mol2_df.loc[mol2_idx,'q']
        # insert fourth point of water, if present in mol2
        if 'HOH' in pdb_df.loc[idx, 'PDB_RES'] and pdb_df.loc[idx, 'PDB_AT'] == 'O':
            EPW_idx = str(int(mol2_idx)+3)
            if EPW_idx in mol2_df.index and 'EP' in mol2_df.loc[EPW_idx, 'MOL2_AT']:
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_AT'] = mol2_df.loc[EPW_idx, 'MOL2_AT']
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_RES'] = mol2_df.loc[EPW_idx, 'MOL2_RES']
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_ID'] = EPW_idx
                pdb_df.loc[str(float(EPW_idx)+.5), 'q'] = mol2_df.loc[EPW_idx, 'q']
                pdb_df.loc[str(float(EPW_idx)+.5), 'X'] = mol2_df.loc[EPW_idx, 'X']
                pdb_df.loc[str(float(EPW_idx)+.5), 'Y'] = mol2_df.loc[EPW_idx, 'Y']
                pdb_df.loc[str(float(EPW_idx)+.5), 'Z'] = mol2_df.loc[EPW_idx, 'Z']
                pdb_df.loc[str(float(EPW_idx)+.5), 'PDB_RES'] = pdb_df.loc[idx, 'PDB_RES']
    return pdb_df

def cx_pdf_to_df(CX_PDB_PATH:str, d:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the pdb
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    CX_PDB_PATH: str
        path to complex pdb to populate the dataframe
    d: Dict
        Dictionary that has one key, CX_PDB_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the PDB

    """
    with open(CX_PDB_PATH, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                d['CX_PDB_ID'].append(l[6:11].strip())
        df = pd.DataFrame(d)
        df = df.set_index('CX_PDB_ID')
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                x_coord = l[29:38].strip()
                y_coord = l[38:46].strip()
                z_coord = l[46:54].strip()
                df.loc[l[6:11].strip(), 'X'] = float(x_coord)
                df.loc[l[6:11].strip(), 'Y'] = float(y_coord)
                df.loc[l[6:11].strip(), 'Z'] = float(z_coord)
    return df

def combine_all_dfs(pdb_df:pd.DataFrame, cx_pdb_df:pd.DataFrame) -> pd.DataFrame:
    """
    combines information from a dataframe from the protein pdb with a
    dataframe from the cx pdb based on matching the coordinates

    Parameters
    ----------
    pdb_df: pd.DataFrame
        dataframe containing information from the protein pdb
    cx_pdb_df: pd.DataFrame
        dataframe containing information from the cx pdb

    Returns
    -------
    df: pd.DataFrame
        dataframe containing information combined from the pdb_df and the cx_df with one entry per atom

    """
    for idx in pdb_df.index:
        cx_idx_x = cx_pdb_df.loc[cx_pdb_df['X'] == pdb_df.loc[idx, 'X']].index.tolist()
        cx_idx_y = cx_pdb_df.loc[cx_pdb_df['Y'] == pdb_df.loc[idx, 'Y']].index.tolist()
        cx_idx_z = cx_pdb_df.loc[cx_pdb_df['Z'] == pdb_df.loc[idx, 'Z']].index.tolist()
        cx_idx_xyz = cx_idx_x + cx_idx_y + cx_idx_z
        for x in cx_idx_xyz:
            if cx_idx_xyz.count(x) == 3:
                cx_idx = x
                pdb_df.loc[idx,'CX_PDB_ID'] = cx_idx
        if pdb_df.loc[idx, 'MOL2_AT'] == 'EP':
            pdb_df.loc[idx,'CX_PDB_ID'] = idx
    return pdb_df

def change_water_charges(df: pd.DataFrame, o: str, h: str, ep:str = None) -> pd.DataFrame:
    """
    Changes charges of the water atoms (and possibly EP) depending on the user specified inputs

    Parameters
    ----------
    df: pd.DataFrame
        dataframe containing information from the mol2s and pdbs
    o: str
        desired oxygen charge
    h: str
        desired hydrogen charge
    ep: str
        desired extra point charge

    Returns
    -------
    df: pd.DataFrame
        original dataframe with water charges updated accordingly 

    """
    for idx in df.index:
        resi = df.loc[idx, 'PDB_RES']
        if 'HOH' in resi or 'WAT' in resi:
            atom = df.loc[idx, 'MOL2_AT']
            if 'O' in atom:
                df.loc[idx, 'q'] = float(o)
            elif 'H' in atom:
                df.loc[idx, 'q'] = float(h)
            elif 'EP' in atom:
                if ep != None:
                    df.loc[idx, 'q'] = float(ep)
    return df

def create_csv(o_charge:str = None, h_charge:str = None, ep_charge:str = None) -> None:
    """
    Creates a dataframe from the provided protein pdb, complex pdb,
    and mol2 information and then updates it with the desired water
    atom charges

    Parameters
    ----------
    o: str
        desired oxygen charge
    h: str
        desired hydrogen charge
    ep: str
        desired extra point charge

    Returns
    -------

    """
    PDB_PATH = 'prot_autocap_fixed.pdb'
    MOL2_PATH = 'prot_autocap_fixed.mol2'
    CX_PDB_PATH = 'cx_autocap_fixed.pdb'
    p = {'PDB_ID':[]}
    m = {'MOL2_ID':[]}
    c = {'CX_PDB_ID':[]}
    pdb_info = prot_pdb_to_df(PDB_PATH, p)
    mol2_info = mol2_to_df(MOL2_PATH, m)
    cx_pdb_info = cx_pdf_to_df(CX_PDB_PATH, c)
    combined = combine_prot_dfs(pdb_info, mol2_info)
    combined2 = combine_all_dfs(combined, cx_pdb_info)
    if o_charge!= None and h_charge != None and ep_charge !=None:
        final = change_water_charges(combined2, o_charge, h_charge, ep_charge)
    elif o_charge != None and h_charge != None:
        final = change_water_charges(combined2, o_charge , h_charge)
    else:
        final = combined2
    final.to_csv('dataframe.csv')
