# PyMOL Pocket Figures

A small, original PyMOL extension for making labeled protein–ligand pocket figures from a **structure you load yourself**. It separates displayed residues into individual stick objects, colors a ligand in a thick ball-and-stick style, optionally draws geometric polar contacts, and saves the camera view directly to a high-resolution transparent PNG.

The code contains no bundled structures, ligand coordinates, project results or fixed residue numbers. It was tested with PyMOL 3.1.1; other PyMOL versions may differ in available rendering settings.

## Quick start

1. Load your prepared complex in PyMOL. Check that ligand bond orders are correct; a PDB alone may not preserve aromaticity or double-bond order, so use a chemically typed ligand file or verify the loaded connectivity where necessary.
2. In the PyMOL command line, run:

   ```text
   run /path/to/pocketfigure.py
   pocket_figure complex, resn LIG and chain A, ligand_label=Compound-1, polar_contacts=yes
   ```

   Replace `complex` with your PyMOL object name and `resn LIG and chain A` with a selection that identifies **one ligand** in that object. Nearby protein residues are selected automatically within 4.5 Å. The command prints its generated prefix, such as `pf001`.

3. Rotate and zoom the scene to your preferred angle. Save the current view in one command:

   ```text
   pocket_save pocket_figure.png, 3600, 2700, 600
   ```

   `pocket_save` performs the ray render itself. There is no separate `ray` step. The default PNG has a transparent background; the PyMOL viewport can still appear white.

For selected residues rather than an automatic distance shell:

```text
pocket_figure complex, resn LIG and chain A, residue_selection=chain A and resi 42+57+103, ligand_label=Compound-1
```

For a second complex in the same session, run `pocket_figure` again with its object and ligand selection. Each run receives a fresh prefix (`pf001`, `pf002`, ...). By default, other loaded objects are disabled so that the chosen complex is isolated in the figure; set `isolate=no` to preserve their visibility.

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `object_name` | required | Name of a loaded PyMOL molecular object. |
| `ligand` | required | PyMOL selection for the ligand inside that object. |
| `residue_selection` | `auto` | PyMOL selection for protein residues; `auto` uses the ligand-distance shell. |
| `residue_radius` | `4.5` Å | Distance used to select residues in automatic mode. |
| `cartoon_radius` | `8.0` Å | Distance used for the surrounding protein cartoon. |
| `ligand_label` | ligand residue name | Text placed near a central ligand atom. |
| `polar_contacts` | `no` | `yes` asks PyMOL to draw mode-2 polar-contact dashes. |
| `contact_cutoff` | `3.6` Å | Maximum distance for the optional polar-contact search. |
| `isolate` | `yes` | Disable other loaded objects in the view. |
| `prefix` | automatic | Prefix for generated selections and residue objects. |

`pocket_save filename.png, width, height, dpi` defaults to 3600 × 2700 pixels at 600 dpi. The command exports the **current** camera view.

## What the figure means

One temporary PyMOL object is created for each displayed amino-acid residue. This suppresses peptide links **between displayed residue sticks** so individual side chains are easier to read. The original molecular object, coordinates and bonds are left intact. Residue labels use three-letter residue names plus numbers; a chain identifier is added automatically if two displayed chains have the same residue name and number. Labels are attached to the residue Cα atom where present.

Amber dashes from `polar_contacts=yes` are PyMOL's geometric polar-contact suggestions. They are **not covalent bonds** and should be checked against atom types, hydrogen orientation and distances before a manuscript assigns specific hydrogen bonds. The tool does not calculate binding affinity or validate a docking pose.

For accurate chemical depictions, inspect the structure before exporting: correct ligand identity and bond order, biological assembly, alternate conformers, atom naming, relevant cofactors and waters, and residue numbering. The script styles the supplied model; it does not repair it.

## Testing

Run the integration smoke test with a Python interpreter from a working PyMOL installation:

```text
/path/to/pymol/bin/python tests/smoke_test.py
```

The test creates a small synthetic protein–ligand structure in memory, checks residue separation and labels, and verifies that the input object's coordinates remain unchanged. No project structure is needed.

## License

MIT; see [LICENSE](LICENSE). This implementation was written specifically as a reusable PyMOL extension and does not include code from third-party rendering gists.
