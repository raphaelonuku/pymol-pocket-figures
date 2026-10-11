"""Integration check using a synthetic structure; run with PyMOL's Python."""

from pathlib import Path
import sys

import pymol

pymol.finish_launching(["pymol", "-cq"])
from pymol import cmd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pocketfigure import pocket_figure  # noqa: E402


def atom_line(record, number, name, resn, chain, resi, x, y, z, element):
    return (
        f"{record:<6}{number:5d} {name:<4s} {resn:>3s} {chain}{resi:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}"
    )


atoms = [
    ("ATOM", 1, "N", "ALA", "A", 1, -2.2, 0.0, 0.0, "N"),
    ("ATOM", 2, "CA", "ALA", "A", 1, -1.0, 0.0, 0.0, "C"),
    ("ATOM", 3, "C", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
    ("ATOM", 4, "O", "ALA", "A", 1, 0.7, 0.9, 0.0, "O"),
    ("ATOM", 5, "CB", "ALA", "A", 1, -1.0, 1.5, 0.0, "C"),
    ("ATOM", 6, "N", "SER", "A", 2, 1.4, -0.2, 0.0, "N"),
    ("ATOM", 7, "CA", "SER", "A", 2, 2.8, 0.0, 0.0, "C"),
    ("ATOM", 8, "C", "SER", "A", 2, 3.8, 0.0, 0.0, "C"),
    ("ATOM", 9, "O", "SER", "A", 2, 4.5, 0.9, 0.0, "O"),
    ("ATOM", 10, "CB", "SER", "A", 2, 2.8, 1.5, 0.0, "C"),
    ("ATOM", 11, "OG", "SER", "A", 2, 2.8, 2.7, 0.0, "O"),
    ("HETATM", 12, "C1", "LIG", "X", 1, 0.0, 0.0, 3.0, "C"),
    ("HETATM", 13, "C2", "LIG", "X", 1, 1.4, 0.0, 3.0, "C"),
    ("HETATM", 14, "O1", "LIG", "X", 1, 2.2, 1.0, 3.0, "O"),
]
pdb = "\n".join(atom_line(*atom) for atom in atoms) + "\nEND\n"
cmd.read_pdbstr(pdb, "example")
before = [(atom.index, tuple(atom.coord)) for atom in cmd.get_model("example").atom]

prefix = pocket_figure(
    "example",
    "resn LIG and chain X",
    residue_selection="chain A and resi 1+2",
    ligand_label="Test ligand",
)
fragments = sorted(name for name in cmd.get_names("objects") if name.startswith(prefix + "_res"))
assert len(fragments) == 2, fragments
assert all(
    len({(atom.chain, atom.resi, atom.resn) for atom in cmd.get_model(name).atom}) == 1
    for name in fragments
)

labels = []
for selection in fragments + [prefix + "_ligand"]:
    cmd.iterate(selection, "labels.append(label)", space={"labels": labels})
assert sorted(label for label in labels if label) == ["ALA1", "SER2", "Test ligand"]

after = [(atom.index, tuple(atom.coord)) for atom in cmd.get_model("example").atom]
assert before == after, "The source object's atom coordinates changed"

auto_prefix = pocket_figure(
    "example", "resn LIG and chain X", residue_radius=4.5,
    polar_contacts="yes", isolate="no",
)
auto_fragments = [name for name in cmd.get_names("objects") if name.startswith(auto_prefix + "_res")]
assert len(auto_fragments) == 2, auto_fragments
assert auto_prefix + "_polar" in cmd.get_names("objects")
assert before == [(atom.index, tuple(atom.coord)) for atom in cmd.get_model("example").atom]
print("PASS: explicit and automatic pockets, individualized residues, labels, polar contacts, unchanged source")
