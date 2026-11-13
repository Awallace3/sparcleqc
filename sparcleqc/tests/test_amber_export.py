"""Regression coverage for topology-backed MOL2 export and charge labels."""

import json
import shutil
import subprocess

import numpy as np
import pandas as pd
import parmed as pmd
import pytest

from sparcleqc import amber_prep, combine_data, complex_tools, create_est_inp


@pytest.mark.validation
def test_tleap_writes_topology_instead_of_mol2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    amber_prep.write_tleap("ff19SB", "opc", [])
    script = (tmp_path / "tleap.in").read_text()
    assert "saveAmberParm mol prot_autocap_fixed.prmtop prot_autocap_fixed.inpcrd" in script
    assert "savemol2" not in script.lower()


@pytest.mark.amber
@pytest.mark.skipif(shutil.which("tleap") is None, reason="AmberTools is required")
@pytest.mark.parametrize("water_name", ["HOH", "WAT"])
def test_parmed_export_preserves_large_ids_charges_and_extra_points(tmp_path, monkeypatch, water_name):
    """Cross both the atom-999 and residue-99 boundaries with real OPC topology."""
    monkeypatch.chdir(tmp_path)
    pdb = tmp_path / "prot_autocap_fixed.pdb"
    pdb.write_text(
        "".join(
            f"HETATM{i:5d}  O   {water_name} A{i:4d}    "
            f"{i * 3.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
            "  1.00  0.00           O  \n"
            for i in range(1, 252)
        )
        + "END\n"
    )
    amber_prep.write_tleap("ff19SB", "opc", [])
    result = subprocess.run(["tleap", "-f", "tleap.in"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    amber_prep.amber_topology_to_mol2()

    topology = pmd.load_file("prot_autocap_fixed.prmtop", xyz="prot_autocap_fixed.inpcrd")
    records = combine_data.mol2_atom_records("prot_autocap_fixed.mol2")
    assert len(records) == len(topology.atoms) == 1004
    assert [int(record[0]) for record in records] == list(range(1, 1005))
    metadata = [combine_data.mol2_atom_metadata(record) for record in records]
    assert [int(item[1]) for item in metadata] == [
        atom.residue.idx + 1 for atom in topology.atoms
    ]
    assert [item[0] for item in metadata] == [atom.type for atom in topology.atoms]
    assert [item[3] for item in metadata] == pytest.approx(
        [atom.charge for atom in topology.atoms], abs=5.1e-7
    )
    np.testing.assert_allclose(
        np.array([record[2:5] for record in records], dtype=float),
        topology.coordinates,
        atol=5.1e-5,
        rtol=0,
    )
    exported = pmd.load_file("prot_autocap_fixed.mol2", structure=True)
    assert [(bond.atom1.idx, bond.atom2.idx) for bond in exported.bonds] == [
        (bond.atom1.idx, bond.atom2.idx) for bond in topology.bonds
    ]
    assert sum(item[0] == "EP" for item in metadata) == 251
    assert complex_tools.check_mol2_charges("prot_autocap_fixed.mol2")[0] == 1

    pdb_df = combine_data.prot_pdb_to_df(str(pdb), {"PDB_ID": []})
    mol2_df = combine_data.mol2_to_df("prot_autocap_fixed.mol2", {"MOL2_ID": []})
    merged = combine_data.combine_prot_dfs(pdb_df, mol2_df)
    assert len(merged) == 502  # original oxygen atoms plus their OPC extra points
    assert merged["q"].notna().all()
    assert (merged["MOL2_AT"] == "EP").sum() == 251


@pytest.mark.validation
def test_parmed_export_requires_valid_topology(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(RuntimeError, match="(?i)topology.*prepar"):
        amber_prep.amber_topology_to_mol2()


@pytest.mark.validation
def test_parmed_export_rejects_invalid_topology_without_fallback(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "prot_autocap_fixed.prmtop").write_text("invalid topology\n")
    (tmp_path / "prot_autocap_fixed.inpcrd").write_text("invalid coordinates\n")
    with pytest.raises(RuntimeError, match="(?i)topology.*prepar"):
        amber_prep.amber_topology_to_mol2()
    assert not (tmp_path / "prot_autocap_fixed.mol2").exists()


@pytest.mark.validation
def test_mol2_dataframe_rejects_duplicate_ids(tmp_path):
    path = tmp_path / "duplicate.mol2"
    path.write_text(
        "@<TRIPOS>MOLECULE\nbad\n2 0 1 0 0\nSMALL\nUSER_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        "100 O 0.0 0.0 0.0 O 1 ASP -0.5\n"
        "100 O 1.0 0.0 0.0 O 1 ASP -0.5\n"
        "@<TRIPOS>BOND\n"
    )
    with pytest.raises(ValueError, match="Duplicate MOL2 atom ID"):
        combine_data.mol2_to_df(str(path), {"MOL2_ID": []})


@pytest.mark.validation
def test_mol2_without_charge_column_is_neutral(tmp_path):
    path = tmp_path / "neutral.mol2"
    path.write_text(
        "@<TRIPOS>MOLECULE\nneutral\n1 0 1 0 0\nSMALL\nNO_CHARGES\n"
        "@<TRIPOS>ATOM\n1 C 0.0 0.0 0.0 CT 1 LIG\n@<TRIPOS>BOND\n"
    )
    df = combine_data.mol2_to_df(str(path), {"MOL2_ID": []})
    assert df.at["1", "q"] == 0.0
    assert complex_tools.check_mol2_charges(str(path))[0] == 1


def _write_charge_dataframe(tmp_path, labels):
    df = pd.DataFrame(
        {
            "PDB_ID": range(1, 6),
            "PDB_AT": ["N", "CB", "N", "CB", "O"],
            "PDB_RES": ["ASP_1", "ASP_1", "LYS_2", "LYS_2", "HOH_3"],
            "X": range(5),
            "Y": [0.0] * 5,
            "Z": [0.0] * 5,
            "AT_LABEL": labels,
            "MOL2_ID": range(1, 6),
            "MOL2_AT": ["N", "CT", "N", "CT", "OW"],
            "MOL2_RES": ["ASP", "ASP", "LYS", "LYS", "WAT"],
            "q": [-0.5, -0.5, 0.5, 0.5, 0.0],
            "CX_PDB_ID": range(1, 6),
        }
    ).set_index("PDB_ID")
    df.to_csv(tmp_path / "dataframe.csv")


@pytest.mark.validation
@pytest.mark.parametrize(
    "labels",
    [
        [None, None, None, None, None],
        ["N1+", "C", "N", "C", "O"],  # conflicting label in first residue
        ["N1-", "C", "N", "C", ""],  # second residue must get its own label
    ],
)
def test_charge_labels_handle_missing_values_and_each_residue(tmp_path, monkeypatch, labels):
    monkeypatch.chdir(tmp_path)
    _write_charge_dataframe(tmp_path, labels)
    assert complex_tools.check_df_charges() == (1, "worked")
    path = tmp_path / "dataframe.csv"
    df = pd.read_csv(path, index_col="PDB_ID", keep_default_na=False)
    for residue, suffix in [("ASP_1", "1-"), ("LYS_2", "1+")]:
        residue_labels = df.loc[df["PDB_RES"] == residue, "AT_LABEL"]
        assert sum(label.endswith(suffix) for label in residue_labels) == 1
        assert not any("nan" in label for label in residue_labels)
    assert not any(label.endswith("1+") for label in df.loc[df["PDB_RES"] == "ASP_1", "AT_LABEL"])
    first_pass = path.read_text()
    assert complex_tools.check_df_charges() == (1, "worked")
    assert path.read_text() == first_pass


@pytest.mark.validation
def test_charge_label_without_beta_carbon(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pd.DataFrame(
        {"PDB_ID": [1], "PDB_AT": ["ZN"], "PDB_RES": ["ZN_1"], "AT_LABEL": [None], "q": [2.0]}
    ).set_index("PDB_ID").to_csv("dataframe.csv")
    assert complex_tools.check_df_charges() == (1, "worked")
    df = pd.read_csv("dataframe.csv", index_col="PDB_ID", keep_default_na=False)
    assert df.at[1, "AT_LABEL"] == "2+"


@pytest.mark.validation
def test_qm_charge_count_handles_blank_labels_after_csv_reload(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pd.DataFrame(
        {"CX_PDB_ID": [1, 2, 3], "AT_LABEL": [None, "C1-", None]}
    ).to_csv("dataframe.csv", index=False)
    (tmp_path / "with_HL.dat").write_text(json.dumps({"QM": [1, 2, 3]}))
    # Isolate charge accounting from the later geometry/boundary construction.
    for name in ("atoms_to_pdb_lines", "write_capped_qm_pdb", "pdb_to_xyz", "ligand_pdb_lines", "Z1"):
        monkeypatch.setattr(create_est_inp, name, lambda *args: [])
    _, charge, _, _ = create_est_inp.make_monomers("Z1")
    assert charge == "-1"


@pytest.mark.validation
@pytest.mark.parametrize(
    "atom_id,charge,message",
    [(1, float("nan"), "Missing"), (1, -0.4, "integer charge"), (5, 1.0, "Water charge")],
)
def test_charge_label_fix_does_not_hide_invalid_charges(
    tmp_path, monkeypatch, atom_id, charge, message
):
    monkeypatch.chdir(tmp_path)
    _write_charge_dataframe(tmp_path, [None] * 5)
    df = pd.read_csv("dataframe.csv", index_col="PDB_ID")
    df.at[atom_id, "q"] = charge
    df.to_csv("dataframe.csv")
    status, diagnostic = complex_tools.check_df_charges()
    assert status == 0
    assert message in diagnostic


@pytest.mark.validation
def test_charge_check_rejects_fractional_last_mol2_residue(tmp_path):
    path = tmp_path / "fractional.mol2"
    path.write_text(
        "@<TRIPOS>MOLECULE\nbad\n2 0 2 0 0\nSMALL\nUSER_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        "1 O 0.0 0.0 0.0 O 1 ASP -1.0\n"
        "2 N 1.0 0.0 0.0 N 2 LYS 0.5\n"
        "@<TRIPOS>BOND\n"
    )
    status, message = complex_tools.check_mol2_charges(str(path))
    assert status == 0
    assert "LYS_2" in message
