from pymol.cgo import *
from pymol import cmd, editor

# These residue names are treated as solvent so blank chain IDs can be
# normalized without mixing waters into inferred protein chains.
WATER_RESNAMES = {'WAT', 'HOH', 'TIP', 'TIP3', 'TIP3P', 'OPC', 'SPC', 'SOL'}

# Some PL-REX structures carry residue name SEM for atoms that are exactly a
# standard serine sidechain. Normalize those local copies before tleap so Amber
# does not create an untyped placeholder residue that later breaks capping.
AMBER_PREP_RESIDUE_ALIASES = {
    ('SEM', frozenset({'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'C', 'O'})): 'SER',
}


def _is_atom_record(line: str) -> bool:
    return line[0:6].strip() in {'ATOM', 'HETATM'}


def _residue_block_key(line: str):
    return (line[16:20].strip(), line[21], line[22:26], line[26])


def _parse_resseq(line: str):
    try:
        return int(line[22:26].strip())
    except ValueError:
        return None


def normalize_residue_names_for_amber_prep(pdb_file: str) -> bool:
    """
    Rename known residue aliases on the local copied PDB before Amber prep.

    This is only for cases where the atom-name signature already matches a
    standard Amber residue and the original label is what prevents typing.
    """
    with open(pdb_file, 'r') as handle:
        lines = handle.readlines()

    residues = {}
    for idx, line in enumerate(lines):
        if not _is_atom_record(line):
            continue
        key = _residue_block_key(line)
        block = residues.setdefault(key, {'line_indices': [], 'atom_names': set()})
        block['line_indices'].append(idx)
        block['atom_names'].add(line[11:16].strip())

    changed = False
    for key, block in residues.items():
        resname = key[0]
        alias = AMBER_PREP_RESIDUE_ALIASES.get((resname, frozenset(block['atom_names'])))
        if alias is None:
            continue
        for line_index in block['line_indices']:
            line = lines[line_index]
            lines[line_index] = line[:17] + f'{alias:>3}' + line[20:]
        changed = True

    if changed:
        with open(pdb_file, 'w') as handle:
            handle.writelines(lines)
    return changed


def _starts_new_chain(previous_block: dict, current_block: dict) -> bool:
    if 'OXT' in previous_block['atom_names']:
        return True
    previous_resseq = previous_block['resseq']
    current_resseq = current_block['resseq']
    if previous_resseq is not None and current_resseq is not None and current_resseq < previous_resseq:
        return True
    return False


# Some input PDBs mark separate protein chains only through terminal atom
# patterns such as OXT on the previous residue and H2/H3 on the next N-terminus.
# This repair fills in missing chain IDs before Amber prep so tleap does not try
# to reconnect distinct chains and then mis-handle legitimate terminal atoms.
def normalize_chain_ids_for_amber_prep(pdb_file: str) -> bool:
    """
    Infer missing chain IDs before Amber preparation.

    This is aimed at structures where terminal atoms such as OXT/H2/H3
    indicate separate protein chains, but the input PDB leaves the chain
    column blank. In that case downstream prep can incorrectly reconnect
    chains and mis-handle legitimate termini.
    """
    with open(pdb_file, 'r') as handle:
        lines = handle.readlines()

    blocks = []
    current_block = None
    for idx, line in enumerate(lines):
        if not _is_atom_record(line):
            continue
        block_key = _residue_block_key(line)
        if current_block is None or current_block['key'] != block_key:
            current_block = {
                'key': block_key,
                'line_indices': [],
                'resname': block_key[0],
                'chain': block_key[1].strip(),
                'resseq': _parse_resseq(line),
                'atom_names': set(),
                'record_types': set(),
            }
            blocks.append(current_block)
        current_block['line_indices'].append(idx)
        current_block['atom_names'].add(line[11:16].strip())
        current_block['record_types'].add(line[0:6].strip())

    if not blocks:
        return False

    existing_chains = {block['chain'] for block in blocks if block['chain']}
    available_protein_chains = [
        chain_id for chain_id in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
        if chain_id not in existing_chains and chain_id != 'W'
    ]
    if not available_protein_chains:
        return False

    water_chain = 'W' if 'W' not in existing_chains else None
    current_chain = None
    previous_protein_block = None
    changed = False

    for block in blocks:
        is_water = block['resname'] in WATER_RESNAMES
        is_protein = 'ATOM' in block['record_types'] and not is_water

        if is_water and not block['chain'] and water_chain is not None:
            assigned_chain = water_chain
        elif is_protein and not block['chain']:
            if previous_protein_block is None or _starts_new_chain(previous_protein_block, block):
                if not available_protein_chains:
                    return changed
                current_chain = available_protein_chains.pop(0)
            assigned_chain = current_chain
        else:
            assigned_chain = None

        if assigned_chain is not None:
            for line_index in block['line_indices']:
                line = lines[line_index]
                lines[line_index] = line[:21] + assigned_chain + line[22:]
            block['chain'] = assigned_chain
            changed = True

        if is_protein:
            previous_protein_block = block
            if block['chain']:
                current_chain = block['chain']

    if changed:
        with open(pdb_file, 'w') as handle:
            handle.writelines(lines)
    return changed


# Wrapper-generated structures often place the user-supplied ligand in residue
# name LIG and chain L while cofactors remain as other non-polymer residues.
# Prefer that specific residue when it exists so bound cofactors are not
# accidentally exported as the main ligand.
def get_ligand_selection() -> str:
    preferred_selections = [
        'resn LIG and chain L',
        'resn LIG',
        'chain L and not polymer and not metals and not solvent and not resn nme and not resn ace',
    ]
    fallback_selection = 'all and not polymer and not metals and not solvent and not resn nme and not resn ace'
    for selection in preferred_selections:
        if cmd.count_atoms(selection) > 0:
            return selection
    return fallback_selection


def reorder_atoms_amber(pdb_file: str) -> None:
    """
    When given a pdb, creates a new copy {pdb_file}_fixed.pdb that
    has the protein residues followed by waters and then the ligand.
    Corrects for any mistakes in atom or residue numbering that may
    have been caused by manipulation of the system in pymol.
    Ensures that the ligand atoms are labeled as HETATM.

    Parameters
    ----------
    pdb_file: str
        path to pdb

    Returns
    -------
    None
    """

    with open('ligand.pdb') as lig:
        lig_lines = lig.readlines()
    for line in lig_lines:
        if line[0:6].strip() =='ATOM' or line[0:6].strip() =='HETATM':
            lig_name = line[16:20].strip()
            break
    
    out = open(f'{pdb_file[:-4]}_fixed.pdb', 'w') 
    with open(pdb_file) as w:
        lines = w.readlines()
    resnum =0
    atomnum = 0
    ligand_lines = []
    HOH_lines = []
    oldres = ''
    for line in lines:
        if 'WAT' not in line and 'HOH' not in line and 'TIP' not in line and len(line)>70 and line[16:20].strip() !=lig_name and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}           {line[66:len(line)].strip():<3}\n')
        elif 'WAT' in line or 'HOH' in line or 'TIP' in line and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            HOH_lines.append(line)
        elif lig_name in line and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            ligand_lines.append(line)
        else:
            pass
    
    for line in HOH_lines:
        if len(line)>70:
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:len(line)].strip():>12}\n')
    for line in ligand_lines:
        if len(line)>70 and line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM':
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            out.write(f'HETATM{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:len(line)].strip():>12}\n')
    if 'cx' in pdb_file:
        out.write('CONECT\n')
    out.write('END')
    out.close()
def autocap(pdb_file: str) -> None:
    """ 
    When given an amber complex pdb, caps the ends of each chain
    with a capping group (NME on C-terminal and ACE on N-terminal)
    and then renames the default pymol atom names to amber atom names
    for the new capping residues. Saves the entire new, capped complex
    as cx_autocap.pdb, just the capped protein as prot_autocap.pdb,
    and the ligand as ligand.pdb

    Parameters
    ----------
    pdb_file: str
        path to uncapped complex pdb

    Returns
    -------
    None
    """
    cmd.reinitialize()
    
    # Load PDB
    cmd.load(pdb_file,"pdb")
    cmd.show("sticks", "all")
    cmd.label("all", "name")
    
    cmd.select("sidechains", "sidechain")
    chains = []
    chains = cmd.get_chains("sidechains")
    for chain in chains:
        #clean termini
        cmd.select('capped_n', 'not name H1 and element H and bound_to (first name N and chain %s)' % chain)
        cmd.remove('capped_n')
        cmd.select("capped_c", "last name OXT and chain %s" % chain)
        cmd.remove("capped_c")
        cmd.alter('first name H1 and chain %s' % chain, "name='H'")
        #cmd.set('retain_order', 0)
        #cap termini
        editor.attach_amino_acid("last name C and chain %s" % chain, 'nme' )
        editor.attach_amino_acid("first name N and chain %s" % chain, 'ace')
    
    
    
    cmd.alter('resname NME and name 1HH3', "name='H1'")
    cmd.alter('resname NME and name 2HH3', "name='H2'")
    cmd.alter('resname NME and name 3HH3', "name='H3'")
    cmd.alter('resname NME and name CH3', "name='C'")
    
    
    cmd.alter('resname ACE and name 1HH3', "name='H1'")
    cmd.alter('resname ACE and name 2HH3', "name='H2'")
    cmd.alter('resname ACE and name 3HH3', "name='H3'")
    
    cmd.save(f"cx_autocap.pdb", "pdb")
    
    cmd.select("ligand", get_ligand_selection())
    cmd.save("ligand.pdb", "ligand")
    cmd.remove("ligand")
    
    cmd.save(f"prot_autocap.pdb", "pdb")

def skip_autocap(pdb_file: str) -> None:
    """
    When given an amber complex pdb that already has capping residues at
    each termini, this function is called instead of autocap(pdb_file)
    Saves just the capped protein as prot_autocap.pdb and the ligand
    as ligand.pdb

    Parameters
    ----------
    pdb_file: str
        path to capped complex pdb

    Returns
    -------
    None
    """
    cmd.reinitialize()
    # Load PDB
    
    cmd.load(pdb_file,"pdb")
    cmd.show("sticks", "all")
    cmd.label("all", "name")
    
    cmd.select("ligand", get_ligand_selection())
    cmd.save("ligand.pdb", "ligand")
    cmd.remove("ligand")
    
    cmd.save(f"prot_autocap.pdb", "pdb")

def write_cpptraj(pdb_file: str)-> None:
    """
    Writes a cpptraj input file for the provided (un-capped) pdb

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout uncapped.pdb nosugar\n')

def write_cpptraj_skip_autocap(pdb_file: str) -> None:
    """
    Writes a cpptraj input file for the provided (capped) pdb

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout cx_autocap.pdb nosugar\n')

def write_tleap(forcefield: str, water_model: str, other_ffs: list) -> None:
    """
    Writes a tleap input file that loads the given forcefield and water model

    Parameters
    ----------
    forcefield: str
        name of amber forcefield (ex. ff19SB)
    water_model: str
        name of the desired water model (ex. OPC)
    other_ffs: list
        list of names of non-water, non-protein model (ex. ['leaprc.DNA.OL15', 'leaprc.gaff2'])

    Returns
    -------
    None
    """
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.protein.{forcefield}\n')
        f.write(f'source leaprc.water.{water_model}\n')
        for ff in other_ffs:
            f.write(f'source {ff}\n')
        f.write(f'mol = loadPdb "prot_autocap_fixed.pdb"\n')
        f.write('check mol\n')
        f.write(f'savemol2 mol prot_autocap_fixed.mol2 1\n')
        f.write('quit')
