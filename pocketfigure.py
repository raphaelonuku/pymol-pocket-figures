"""Make labeled, publication-ready PyMOL views of a protein–ligand pocket.

Load this file in PyMOL with ``run pocketfigure.py``. The module registers the
``pocket_figure`` and ``pocket_save`` commands. It does not load or alter a
structure file; the user supplies a molecular object already open in PyMOL.
"""

from __future__ import annotations

from collections import defaultdict
from math import dist
from pathlib import Path
import re

from pymol import CmdException, cmd


_PREFIX = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_ELEMENT_COLORS = {"N": "blue", "O": "red", "S": "yellow", "P": "orange"}


def _positive_number(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise CmdException(f"{name} must be a positive number") from error
    if not 0 < result < float("inf"):
        raise CmdException(f"{name} must be a positive finite number")
    return result


def _yes(value: object) -> bool:
    text = str(value).strip().lower()
    if text in {"1", "yes", "true", "on"}:
        return True
    if text in {"0", "no", "false", "off"}:
        return False
    raise CmdException("Expected yes/no, true/false, or 1/0")


def _prefix(requested: str) -> str:
    names = set(cmd.get_names("all"))
    if requested:
        if not _PREFIX.fullmatch(requested):
            raise CmdException("prefix must start with a letter and contain only letters, digits or underscores")
        if any(name == requested or name.startswith(requested + "_") for name in names):
            raise CmdException(f"prefix '{requested}' already exists; choose another prefix")
        return requested
    for number in range(1, 10000):
        candidate = f"pf{number:03d}"
        if not any(name == candidate or name.startswith(candidate + "_") for name in names):
            return candidate
    raise CmdException("No unused figure prefix is available")


def _residue_label(atom: object, duplicated: bool) -> str:
    base = f"{atom.resn}{atom.resi}"
    chain = atom.chain.strip()
    return f"{chain}:{base}" if duplicated and chain else base


def _style_atoms(selection: str, carbon_color: str) -> None:
    cmd.color(carbon_color, f"({selection}) and elem C")
    for element, color in _ELEMENT_COLORS.items():
        cmd.color(color, f"({selection}) and elem {element}")


def _label_atom(selection: str, text: str) -> None:
    atoms = [atom for atom in cmd.get_model(selection).atom if atom.symbol.upper() != "H"]
    if not atoms:
        raise CmdException(f"No heavy atom is available for label '{text}'")
    anchor = next((atom for atom in atoms if atom.name == "CA" and not atom.alt), None)
    if anchor is None:
        anchor = next((atom for atom in atoms if atom.name == "CA"), atoms[0])
    cmd.label(f"({selection}) and index {anchor.index}", repr(text))


def pocket_figure(
    object_name: str,
    ligand: str,
    residue_selection: str = "auto",
    residue_radius: object = 4.5,
    cartoon_radius: object = 8.0,
    ligand_label: str = "",
    polar_contacts: object = "no",
    contact_cutoff: object = 3.6,
    isolate: object = "yes",
    prefix: str = "",
) -> str:
    """Style a loaded complex and label individual nearby protein residues.

    ``ligand`` and ``residue_selection`` are PyMOL selections within
    ``object_name``. With ``residue_selection=auto``, protein residues within
    ``residue_radius`` Å of the ligand are displayed. The function creates one
    visual object per residue so that inter-residue peptide bonds are not drawn
    in the stick layer; it never edits the input object's coordinates or bonds.
    ``polar_contacts=yes`` draws PyMOL mode-2 geometric polar contacts.
    """

    if object_name not in cmd.get_names("objects"):
        raise CmdException(f"Molecular object '{object_name}' is not loaded")
    if not str(ligand).strip():
        raise CmdException("Provide a ligand selection, for example 'resn LIG and chain A'")
    residue_radius = _positive_number(residue_radius, "residue_radius")
    cartoon_radius = _positive_number(cartoon_radius, "cartoon_radius")
    contact_cutoff = _positive_number(contact_cutoff, "contact_cutoff")
    contacts_enabled = _yes(polar_contacts)
    isolate_enabled = _yes(isolate)
    tag = _prefix(str(prefix).strip())
    ligand_sel = f"{tag}_ligand"
    pocket_sel = f"{tag}_pocket"
    cartoon_sel = f"{tag}_cartoon"

    cmd.select(ligand_sel, f"({object_name}) and ({ligand})")
    if cmd.count_atoms(ligand_sel) == 0:
        cmd.delete(ligand_sel)
        raise CmdException("The ligand selection contains no atoms in the specified object")

    if str(residue_selection).strip().lower() == "auto":
        pocket_expr = (
            f"byres ((({object_name}) and polymer.protein) "
            f"within {residue_radius:g} of ({ligand_sel}))"
        )
    else:
        pocket_expr = (
            f"byres ((({object_name}) and polymer.protein) "
            f"and ({residue_selection}))"
        )
    cmd.select(pocket_sel, pocket_expr)
    atoms = cmd.get_model(pocket_sel).atom
    if not atoms:
        cmd.delete(ligand_sel)
        cmd.delete(pocket_sel)
        raise CmdException("No protein residues were selected; increase the radius or supply residue_selection")

    groups: dict[tuple[str, str, str, str], list[object]] = defaultdict(list)
    for atom in atoms:
        groups[(atom.segi, atom.chain, atom.resi, atom.resn)].append(atom)
    if len(groups) > 30:
        cmd.delete(ligand_sel)
        cmd.delete(pocket_sel)
        raise CmdException("More than 30 residues selected; use a smaller radius or an explicit residue_selection")

    cmd.select(
        cartoon_sel,
        f"byres ((({object_name}) and polymer.protein) within {cartoon_radius:g} of ({ligand_sel}))",
    )
    if isolate_enabled:
        cmd.disable("all")
        cmd.enable(object_name)
    cmd.hide("everything", object_name)
    cmd.set_color(f"{tag}_protein_green", [0.42, 0.77, 0.28])
    cmd.set_color(f"{tag}_ligand_cream", [0.992, 0.894, 0.812])
    cmd.set_color(f"{tag}_contact_amber", [0.80, 0.45, 0.02])

    cmd.bg_color("white")
    for setting, value in (
        ("opaque_background", 0),
        ("ray_opaque_background", 0),
        ("orthoscopic", 1),
        ("ray_trace_mode", 1),
        ("ray_shadows", 0),
        ("ambient", 0.5),
        ("antialias", 3),
        ("sphere_quality", 4),
        ("stick_quality", 24),
        ("float_labels", 1),
        ("label_size", 16),
        ("label_shadow_mode", 0),
        ("label_color", "black"),
        ("label_outline_color", "white"),
        ("label_bg_color", "white"),
        ("label_bg_transparency", 1.0),
    ):
        cmd.set(setting, value)

    cmd.show("cartoon", cartoon_sel)
    cmd.color(f"{tag}_protein_green", cartoon_sel)
    cmd.set("cartoon_transparency", 0.15, cartoon_sel)
    cmd.hide("cartoon", f"byres (({cartoon_sel}) within 3.5 of ({ligand_sel}))")

    duplicates: dict[tuple[str, str], int] = defaultdict(int)
    for _, chain, resi, resn in groups:
        duplicates[(resn, resi)] += 1
    fragments = []
    for number, (key, residue_atoms) in enumerate(sorted(groups.items()), start=1):
        fragment = f"{tag}_res{number:02d}"
        indices = "+".join(str(atom.index) for atom in residue_atoms)
        cmd.create(fragment, f"({pocket_sel}) and index {indices}")
        fragments.append(fragment)
        cmd.hide("everything", fragment)
        cmd.show("sticks", fragment)
        cmd.set("stick_radius", 0.17, fragment)
        _style_atoms(fragment, f"{tag}_protein_green")
        label = _residue_label(residue_atoms[0], duplicates[(key[3], key[2])] > 1)
        _label_atom(fragment, label)
        cmd.set("label_position", [0.8, 0.6, 8.0], fragment)

    cmd.show("sticks", ligand_sel)
    cmd.show("spheres", ligand_sel)
    cmd.set("stick_radius", 0.17, ligand_sel)
    cmd.set("sphere_scale", 0.18, ligand_sel)
    cmd.set("valence", 1, ligand_sel)
    _style_atoms(ligand_sel, f"{tag}_ligand_cream")
    cmd.hide("everything", f"({ligand_sel}) and elem H")

    ligand_atoms = [atom for atom in cmd.get_model(ligand_sel).atom if atom.symbol.upper() != "H"]
    centroid = tuple(sum(atom.coord[axis] for atom in ligand_atoms) / len(ligand_atoms) for axis in range(3))
    anchor_candidates = [atom for atom in ligand_atoms if atom.symbol.upper() == "C"] or ligand_atoms
    anchor = min(anchor_candidates, key=lambda atom: dist(atom.coord, centroid))
    ligand_text = str(ligand_label).strip() or ligand_atoms[0].resn or "Ligand"
    cmd.label(f"({ligand_sel}) and index {anchor.index}", repr(ligand_text))
    cmd.set("label_position", [0.8, 0.6, 8.0], ligand_sel)

    if contacts_enabled:
        contact_name = f"{tag}_polar"
        cmd.distance(
            contact_name,
            ligand_sel,
            " or ".join(f"({fragment})" for fragment in fragments),
            contact_cutoff,
            2,
        )
        cmd.hide("labels", contact_name)
        cmd.set("dash_color", f"{tag}_contact_amber", contact_name)
        cmd.set("dash_width", 3.0, contact_name)

    focus = f"({ligand_sel}) or ({pocket_sel})"
    cmd.orient(focus)
    cmd.zoom(focus, buffer=2.5)
    print(f"Created {tag}: {len(fragments)} labeled residues; rotate the view, then use pocket_save")
    return tag


def pocket_save(filename: str, width: object = 3600, height: object = 2700, dpi: object = 600) -> None:
    """Ray-render the current camera view directly to a transparent PNG."""

    path = Path(str(filename)).expanduser()
    if path.suffix.lower() != ".png":
        raise CmdException("pocket_save requires a .png filename")
    width = int(_positive_number(width, "width"))
    height = int(_positive_number(height, "height"))
    dpi = _positive_number(dpi, "dpi")
    if width < 100 or height < 100:
        raise CmdException("width and height must each be at least 100 pixels")
    path.parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(path), width, height, dpi, 1)
    print(f"Saved {path} ({width} × {height} pixels, {dpi:g} dpi)")


cmd.extend("pocket_figure", pocket_figure)
cmd.extend("pocket_save", pocket_save)
