#!/usr/bin/env python3
"""
constrained_smina_docking.py
============================

End-to-end constrained SMINA docking pipeline, consolidating what used to be
four separate scripts:

    constrained_docking_prep_thorough.py    -> stage "prep"
    run_smina_docking.sh                    -> stage "dock"
    generate_conformers.py                  -> stage "solution"
    conformer_information_addition_to_smina.py + create_pymol_session.py
                                            -> stage "enrich" (+ visualisation)

Pipeline
--------
  prep      Find the maximum common substructure (MCS) between each query
            compound and the reference core (center.pdb), then embed conformers
            with the MCS atoms pinned to the reference coordinates.  Non-chair
            six-membered saturated rings are rejected via Cremer-Pople analysis.

  dock      Run SMINA on every constrained conformer.  Default mode is
            --minimize (local relaxation), which is the correct mode for
            pre-positioned constrained poses.  Both the affinity and the
            pose-drift RMSD reported by SMINA are captured; the drift tells you
            whether the constraint actually survived minimisation.

  solution  Generate an unconstrained (solution-phase) conformer ensemble,
            align each member to the best docked pose, and compute MMFF
            energies, Boltzmann populations and RMSD-to-pose.

  enrich    Merge the solution ensemble summary into the docking summary and
            build PyMOL sessions.

Protonation
-----------
RDKit's Chem.AddHs() fills *neutral* valences; it performs no protonation-state
assignment.  MMFF94 is an all-atom force field, so the solution-phase energies
(and hence the Boltzmann populations) are only meaningful for the species that
actually exists in solution.  That is where the protomer earns its keep.

It is NOT a meaningful lever on the docking score.  Measured on identical
coordinates, the neutral and +2 forms of the demo ligand score -5.685 and
-5.723 kcal/mol -- a 0.04 kcal/mol difference.  Vina-family scoring has no
electrostatic term, so charge state barely registers; the H-bond term is
non-directional and keys off heavy-atom donor/acceptor typing.  For the same
reason, adding hydrogens to the *receptor* changes the score by exactly zero
(verified: 2365 added hydrogens, identical affinity to five decimals), so
receptor protonation is not required and is not performed.

The pipeline assigns a dominant protomer with Open Babel (``obabel -p <pH>``)
and uses that single species consistently across every stage -- not because
docking demands it, but so the solution ensemble and the docked pose describe
the same molecule.  Open Babel's -p is a coarse SMARTS transform table rather
than a real pKa predictor, so --protomer-col lets you pin a protomer per
compound in the input CSV when you disagree with it.

Environment
-----------
Designed to run inside a conda environment that provides both SMINA and RDKit,
e.g. ``conda activate smina``.  Kept Python 3.8 compatible on purpose (no PEP
604 unions, no builtin-generic annotations, no argparse.BooleanOptionalAction).

Example
-------
    python constrained_smina_docking.py \\
        -i target_compounds.csv \\
        -r receptor.pdb \\
        -C center.pdb \\
        -o results
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

STAGES = ("prep", "dock", "solution", "enrich")

# Accepted spellings for the two required input columns, matched case-insensitively.
ID_COL_CANDIDATES = ("Compound ID", "Compound #", "compound_id", "ID", "Name", "Compound")
SMILES_COL_CANDIDATES = ("SMILES", "Canonical SMILES", "canonical_smiles", "smi")

# Columns appended to the docking summary by the enrich stage, in order.
CONFORMER_COLUMNS = ["Min E", "Min E RMSD", "Nearest E", "Nearest RMSD"]

GAS_CONSTANT = 0.001987  # kcal/(mol K)

# PyMOL object colours, cycled across compounds. All are real PyMOL colour names.
PYMOL_COLORS = (
    "cyan", "magenta", "orange", "green", "purple", "pink", "lime", "marine",
    "salmon", "lightblue", "wheat", "palegreen", "lightpink", "deepteal",
    "yelloworange", "slate", "olive", "raspberry",
)

MCS_ATOM_COMPARE = {
    "elements": rdFMCS.AtomCompare.CompareElements,
    "isotopes": rdFMCS.AtomCompare.CompareIsotopes,
    "any": rdFMCS.AtomCompare.CompareAny,
}
MCS_BOND_COMPARE = {
    "any": rdFMCS.BondCompare.CompareAny,
    "order": rdFMCS.BondCompare.CompareOrder,
    "orderexact": rdFMCS.BondCompare.CompareOrderExact,
}


# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------

class Log(object):
    """Minimal leveled logger; avoids a logging.config dependency."""

    def __init__(self, verbose=False, quiet=False):
        self.verbose = verbose
        self.quiet = quiet
        self.warnings = []

    def rule(self, title):
        if self.quiet:
            return
        print("=" * 78)
        print(title)
        print("=" * 78)

    def info(self, msg=""):
        if not self.quiet:
            print(msg)

    def debug(self, msg):
        if self.verbose and not self.quiet:
            print("    " + msg)

    def warn(self, msg):
        self.warnings.append(msg)
        print("  WARNING: " + msg, file=sys.stderr)

    def error(self, msg):
        print("  ERROR: " + msg, file=sys.stderr)


LOG = Log()


# --------------------------------------------------------------------------
# External tool discovery
# --------------------------------------------------------------------------

def _conda_env_candidates(tool):
    """Plausible locations for `tool` inside sibling conda environments."""
    roots = []
    for var in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX"):
        value = os.environ.get(var)
        if value:
            roots.append(value)
            roots.append(os.path.dirname(os.path.dirname(value)))  # .../envs/<x> -> root
    roots.append(sys.prefix)
    roots.extend([
        os.path.expanduser("~/miniconda3"), os.path.expanduser("~/anaconda3"),
        os.path.expanduser("~/miniforge3"), "/opt/miniconda3", "/opt/anaconda3",
        "/opt/homebrew", "/usr/local",
    ])

    candidates = []
    seen = set()
    for root in roots:
        if not root or root in seen:
            continue
        seen.add(root)
        candidates.append(os.path.join(root, "bin", tool))
        candidates.extend(sorted(glob.glob(os.path.join(root, "envs", "*", "bin", tool))))
    return candidates


def find_tool(tool, explicit=None, required=True):
    """Resolve an external executable, searching PATH then nearby conda envs."""
    if explicit:
        resolved = shutil.which(explicit) or (explicit if os.path.isfile(explicit) else None)
        if resolved and os.access(resolved, os.X_OK):
            return resolved
        if required:
            raise SystemExit("Error: --{0}-path '{1}' is not an executable file".format(tool, explicit))
        return None

    found = shutil.which(tool)
    if found:
        return found

    for candidate in _conda_env_candidates(tool):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    if required:
        raise SystemExit(
            "Error: could not find '{0}' on PATH or in a nearby conda environment.\n"
            "       Activate the environment that provides it, or pass --{0}-path.".format(tool)
        )
    return None


# --------------------------------------------------------------------------
# Input CSV
# --------------------------------------------------------------------------

def resolve_column(df, explicit, candidates, what):
    """Pick a column by explicit name, else by case-insensitive candidate match."""
    if explicit:
        if explicit not in df.columns:
            raise SystemExit(
                "Error: column '{0}' not found. Available: {1}".format(explicit, ", ".join(df.columns))
            )
        return explicit

    lowered = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        hit = lowered.get(candidate.lower())
        if hit is not None:
            return hit

    raise SystemExit(
        "Error: could not auto-detect the {0} column.\n"
        "       Looked for: {1}\n"
        "       Found: {2}\n"
        "       Pass it explicitly.".format(what, ", ".join(candidates), ", ".join(df.columns))
    )


def load_compounds(args):
    """Read the input CSV into a list of {id, smiles, protomer_override} records."""
    if not os.path.exists(args.input):
        raise SystemExit("Error: input CSV not found: {0}".format(args.input))

    df = pd.read_csv(args.input)
    id_col = resolve_column(df, args.id_col, ID_COL_CANDIDATES, "compound ID")
    smiles_col = resolve_column(df, args.smiles_col, SMILES_COL_CANDIDATES, "SMILES")

    protomer_col = None
    if args.protomer_col:
        protomer_col = resolve_column(df, args.protomer_col, (), "protomer")

    LOG.info("Input CSV: {0}".format(args.input))
    LOG.info("  ID column     : {0!r}".format(id_col))
    LOG.info("  SMILES column : {0!r}".format(smiles_col))
    if protomer_col:
        LOG.info("  Protomer column: {0!r}".format(protomer_col))

    compounds = []
    for _, row in df.iterrows():
        compound_id = str(row[id_col]).strip()
        smiles = str(row[smiles_col]).strip()
        if not compound_id or compound_id.lower() == "nan" or not smiles or smiles.lower() == "nan":
            LOG.warn("skipping row with empty ID or SMILES")
            continue
        override = None
        if protomer_col:
            raw = row[protomer_col]
            if isinstance(raw, str) and raw.strip() and raw.strip().lower() != "nan":
                override = raw.strip()
        compounds.append({"id": compound_id, "smiles": smiles, "protomer_override": override})

    if not compounds:
        raise SystemExit("Error: no usable rows in {0}".format(args.input))
    return compounds


def safe_name(text):
    """Sanitise an identifier for use as a filename and a PyMOL object name."""
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


# --------------------------------------------------------------------------
# Protonation
# --------------------------------------------------------------------------

def obabel_protomer(smiles, ph, obabel_path):
    """
    Return the dominant protomer SMILES at `ph` according to Open Babel.

    Open Babel's -p applies a SMARTS transform table (phmodel.txt), not a
    trained pKa model.  It reliably handles carboxylic acids and aliphatic
    amines; it is less trustworthy for heteroaryl nitrogens and for basic
    centres whose pKa is perturbed by nearby electron-withdrawing groups.
    Returns None on any failure so the caller can fall back to the input.
    """
    cmd = [obabel_path, "-:{0}".format(smiles), "-osmi", "-p", str(ph)]
    try:
        proc = subprocess.run(cmd, capture_output=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        LOG.warn("obabel protonation failed ({0})".format(exc))
        return None

    text = proc.stdout.decode("utf-8", "replace").strip()
    if not text:
        return None
    candidate = text.splitlines()[0].split()[0].strip()
    return candidate or None


def prepare_species(compound, args, obabel_path):
    """
    Decide the chemical species to carry through every stage.

    Precedence: explicit per-compound override > obabel at --ph > input SMILES.
    Returns (mol, species_smiles, formal_charge, source).  `mol` has no explicit
    hydrogens; each stage adds them itself.
    """
    original = compound["smiles"]

    if compound["protomer_override"]:
        mol = Chem.MolFromSmiles(compound["protomer_override"])
        if mol is not None:
            return mol, Chem.MolToSmiles(mol), Chem.GetFormalCharge(mol), "override"
        LOG.warn("{0}: protomer override did not parse; falling back".format(compound["id"]))

    if args.protonate == "obabel" and obabel_path:
        candidate = obabel_protomer(original, args.ph, obabel_path)
        if candidate:
            mol = Chem.MolFromSmiles(candidate)
            if mol is not None:
                return mol, Chem.MolToSmiles(mol), Chem.GetFormalCharge(mol), "obabel"
            LOG.warn("{0}: obabel protomer did not parse in RDKit; using input SMILES".format(compound["id"]))

    mol = Chem.MolFromSmiles(original)
    if mol is None:
        return None, original, 0, "unparsed"
    return mol, Chem.MolToSmiles(mol), Chem.GetFormalCharge(mol), "input"


# --------------------------------------------------------------------------
# Ring-pucker (Cremer-Pople) chair filter
# --------------------------------------------------------------------------

def get_saturated_6_rings(mol):
    """
    Six-membered rings in which every ring bond is a single bond.

    Benzo-fused rings are deliberately excluded: two of their ring bonds are
    aromatic, so they are half-chairs rather than chairs and the chair/boat
    dichotomy does not apply.
    """
    rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        all_single = True
        for i in range(6):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % 6])
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                all_single = False
                break
        if all_single:
            rings.append(list(ring))
    return rings


def cremer_pople_theta(conf, ring_atoms):
    """
    Cremer-Pople polar puckering angle theta, in degrees, for a 6-ring.

    Reference: Cremer & Pople, JACS 1975, 97, 1354.
        theta ~   0 deg -> 4C1 chair
        theta ~ 180 deg -> 1C4 chair
        theta ~  90 deg -> boat / twist-boat / half-chair / envelope

    Returns 90.0 (i.e. "not a chair") on any numerical degeneracy, so that
    failures are conservative.
    """
    n = 6
    positions = np.array([
        [conf.GetAtomPosition(idx).x, conf.GetAtomPosition(idx).y, conf.GetAtomPosition(idx).z]
        for idx in ring_atoms
    ])
    centered = positions - positions.mean(axis=0)

    # Mean-plane basis vectors R' and R'' (Eqs. 3-4).
    r_prime = np.zeros(3)
    r_dprime = np.zeros(3)
    for j in range(n):
        angle = 2.0 * np.pi * j / n
        r_prime += centered[j] * np.sin(angle)
        r_dprime += centered[j] * np.cos(angle)

    normal = np.cross(r_prime, r_dprime)
    norm = np.linalg.norm(normal)
    if norm < 1e-10:
        return 90.0
    normal /= norm

    # Signed out-of-plane displacements (Eq. 6).
    z = centered.dot(normal)

    # Puckering coordinates for n = 6 (Eqs. 9-11).
    q2_cos = 0.0
    q2_sin = 0.0
    q3 = 0.0
    for j in range(n):
        angle = 4.0 * np.pi * j / n
        q2_cos += z[j] * np.cos(angle)
        q2_sin -= z[j] * np.sin(angle)
        q3 += z[j] * ((-1) ** j)
    q2_cos *= np.sqrt(2.0 / n)
    q2_sin *= np.sqrt(2.0 / n)
    q3 *= np.sqrt(1.0 / n)

    q2 = math.sqrt(q2_cos ** 2 + q2_sin ** 2)
    amplitude = math.sqrt(q2 ** 2 + q3 ** 2)
    if amplitude < 1e-10:
        return 90.0  # planar ring

    return math.degrees(math.acos(max(-1.0, min(1.0, q3 / amplitude))))


def all_rings_are_chair(mol, conf_id, theta_max, sat_rings=None):
    """True if every saturated 6-ring in this conformer is a chair."""
    if sat_rings is None:
        sat_rings = get_saturated_6_rings(mol)
    if not sat_rings:
        return True
    conf = mol.GetConformer(conf_id)
    for ring_atoms in sat_rings:
        theta = cremer_pople_theta(conf, ring_atoms)
        if not (theta <= theta_max or theta >= 180.0 - theta_max):
            return False
    return True


def apply_chair_filter(mol, conf_ids, theta_max, label):
    """Keep only chair conformers; fall back to the unfiltered set if none survive."""
    sat_rings = get_saturated_6_rings(mol)
    if not sat_rings:
        LOG.debug("{0}: no saturated 6-rings, chair filter is a no-op".format(label))
        return conf_ids, 0

    kept = [cid for cid in conf_ids if all_rings_are_chair(mol, cid, theta_max, sat_rings)]
    rejected = len(conf_ids) - len(kept)
    if rejected:
        LOG.info("  Chair filter ({0} saturated 6-ring(s)): rejected {1}, kept {2}".format(
            len(sat_rings), rejected, len(kept)))
    if not kept and conf_ids:
        LOG.warn("{0}: every conformer had a non-chair ring; keeping unfiltered set".format(label))
        return conf_ids, rejected
    return kept, rejected


# --------------------------------------------------------------------------
# Molecule I/O helpers
# --------------------------------------------------------------------------

def write_conformer(mol, conf_id, path, name, fmt):
    """Write a single conformer to SDF (preferred) or PDB."""
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", name)
    if fmt == "sdf":
        writer = Chem.SDWriter(str(path))
        writer.write(mol, confId=conf_id)
        writer.close()
    else:
        writer = Chem.PDBWriter(str(path))
        writer.write(mol, confId=conf_id)
        writer.close()


def read_molecule(path, remove_hs=False):
    """Read a single molecule from SDF or PDB, returning None on failure."""
    path = str(path)
    if not os.path.exists(path):
        return None
    if path.lower().endswith((".sdf", ".mol")):
        supplier = Chem.SDMolSupplier(path, removeHs=remove_hs, sanitize=True)
        for mol in supplier:
            if mol is not None:
                return mol
        # Retry without sanitisation for poses SMINA wrote with odd valences.
        supplier = Chem.SDMolSupplier(path, removeHs=remove_hs, sanitize=False)
        for mol in supplier:
            if mol is not None:
                return mol
        return None
    return Chem.MolFromPDBFile(path, removeHs=remove_hs)


def heavy_atom_indices(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]


def direct_rmsd(mol, cid_a, cid_b, indices):
    """
    In-place RMSD between two conformers over `indices`, with no superposition.

    Constrained conformers all live in the receptor frame, so superposing them
    before comparison would hide exactly the differences we are pruning on.
    """
    conf_a = mol.GetConformer(cid_a)
    conf_b = mol.GetConformer(cid_b)
    total = 0.0
    for idx in indices:
        pa = conf_a.GetAtomPosition(idx)
        pb = conf_b.GetAtomPosition(idx)
        total += (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
    return math.sqrt(total / max(1, len(indices)))


def prune_by_rmsd(mol, conf_ids, threshold, indices):
    """Greedily drop conformers within `threshold` of an already-kept one."""
    if not threshold or threshold <= 0:
        return conf_ids
    kept = []
    for cid in conf_ids:
        if all(direct_rmsd(mol, cid, other, indices) >= threshold for other in kept):
            kept.append(cid)
    if len(kept) < len(conf_ids):
        LOG.info("  RMS pruning at {0} A: {1} -> {2} conformers".format(
            threshold, len(conf_ids), len(kept)))
    return kept


def build_forcefield(mol, conf_id, which):
    """MMFF94 force field for a conformer, falling back to UFF."""
    if which == "mmff":
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
            if ff is not None:
                return ff
    return AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)


# --------------------------------------------------------------------------
# Reference core and MCS
# --------------------------------------------------------------------------

def load_reference(path, reference_smiles=None):
    """
    Load the reference core from PDB.

    Two deliberate departures from the original prep script:

    1. No Chem.AddHs().  A PDB without CONECT records is read with every bond
       single, so aromatic carbons look sp3 and AddHs inflates the molecule with
       chemically wrong hydrogens (26 -> 67 atoms for center.pdb).  The
       reference is only ever used for heavy-atom MCS matching and for
       coordinates, so the hydrogens were pure downside.

    2. Optional bond-order repair from a template SMILES, which restores
       aromaticity and makes strict bond comparison modes usable.
    """
    LOG.info("Reference core: {0}".format(path))
    ref = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
    if ref is None:
        ref = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
    if ref is None:
        raise SystemExit("Error: could not parse reference PDB: {0}".format(path))

    if reference_smiles:
        template = Chem.MolFromSmiles(reference_smiles)
        if template is None:
            LOG.warn("--reference-smiles did not parse; keeping PDB-perceived bond orders")
        else:
            try:
                ref = AllChem.AssignBondOrdersFromTemplate(template, ref)
                LOG.info("  Bond orders assigned from --reference-smiles")
            except Exception as exc:
                LOG.warn("could not apply --reference-smiles template ({0}); "
                         "keeping PDB-perceived bond orders".format(exc))

    aromatic = sum(1 for b in ref.GetBonds() if b.GetIsAromatic())
    LOG.info("  {0} atoms, {1} heavy, {2} aromatic bonds".format(
        ref.GetNumAtoms(), ref.GetNumHeavyAtoms(), aromatic))
    if aromatic == 0 and not reference_smiles:
        LOG.info("  (PDB has no bond orders; pass --reference-smiles to restore aromaticity)")
    return ref


def find_core(ref_mol, query_mol, args):
    """
    Locate the constrained core.

    Returns (core_mol, description) where core_mol is a query-able pattern, or
    (None, reason) if no usable core exists.
    """
    if args.mcs_smarts:
        core = Chem.MolFromSmarts(args.mcs_smarts)
        if core is None:
            return None, "--mcs-smarts did not parse"
        return core, "user SMARTS ({0} atoms)".format(core.GetNumAtoms())

    result = rdFMCS.FindMCS(
        [ref_mol, query_mol],
        timeout=args.mcs_timeout,
        atomCompare=MCS_ATOM_COMPARE[args.mcs_atom_compare],
        bondCompare=MCS_BOND_COMPARE[args.mcs_bond_compare],
        ringMatchesRingOnly=args.mcs_ring_matches_ring_only,
        completeRingsOnly=args.mcs_complete_rings_only,
        matchValences=False,
    )
    if result.numAtoms == 0:
        return None, "no common substructure found"
    if result.canceled:
        LOG.warn("MCS search hit the {0}s timeout; result may be suboptimal".format(args.mcs_timeout))

    core = Chem.MolFromSmarts(result.smartsString)
    if core is None:
        return None, "MCS SMARTS did not parse"
    return core, "MCS ({0} atoms)".format(result.numAtoms)


def choose_mapping(ref_mol, query_mol, core, args):
    """
    Pick the reference/query atom correspondence for the core.

    A symmetric core matches its own coordinates in several orientations, and
    GetSubstructMatch() returns an arbitrary one -- which can pin the query onto
    a flipped copy of the reference.  With --all-symmetry-matches we trial-embed
    one conformer per candidate mapping and keep whichever reproduces the
    reference coordinates most closely.
    """
    ref_matches = ref_mol.GetSubstructMatches(core, uniquify=False, maxMatches=64)
    query_matches = query_mol.GetSubstructMatches(core, uniquify=False, maxMatches=64)
    if not ref_matches or not query_matches:
        return None, None, "core did not match both molecules"

    if not args.all_symmetry_matches or len(ref_matches) == 1:
        return ref_matches[0], query_matches[0], "first match"

    ref_conf = ref_mol.GetConformer()
    query_h = Chem.AddHs(query_mol)
    best = None

    for ref_match in ref_matches[: args.max_symmetry_trials]:
        coord_map = {}
        for core_idx, query_idx in enumerate(query_matches[0]):
            coord_map[query_idx] = ref_conf.GetAtomPosition(ref_match[core_idx])
        trial = Chem.Mol(query_h)
        trial_map = [(query_idx, ref_match[core_idx])
                     for core_idx, query_idx in enumerate(query_matches[0])]
        cids = embed_with_core(trial, coord_map, 1, args.constrained_seed, args,
                               ref_mol, trial_map)
        if not cids:
            continue

        conf = trial.GetConformer(cids[0])
        total = 0.0
        for core_idx, query_idx in enumerate(query_matches[0]):
            pos = conf.GetAtomPosition(query_idx)
            target = ref_conf.GetAtomPosition(ref_match[core_idx])
            total += (pos.x - target.x) ** 2 + (pos.y - target.y) ** 2 + (pos.z - target.z) ** 2
        score = math.sqrt(total / len(query_matches[0]))
        if best is None or score < best[0]:
            best = (score, ref_match)

    if best is None:
        return ref_matches[0], query_matches[0], "first match (all trials failed)"
    return best[1], query_matches[0], "best of {0} symmetry mappings (core RMSD {1:.3f} A)".format(
        min(len(ref_matches), args.max_symmetry_trials), best[0])


# --------------------------------------------------------------------------
# Stage 1 -- constrained conformer preparation
# --------------------------------------------------------------------------

def embed_with_core(mol, coord_map, num_confs, seed, args, ref_mol, atom_map):
    """
    Embed conformers with the core pinned, then place them in the reference frame.

    Two RDKit behaviours make this less obvious than it looks, and they mask
    each other:

    1. `coordMap` constrains the *internal* geometry only. The embedded
       conformer comes out in an arbitrary frame -- in practice centred on the
       origin, ~230 A away from a receptor-frame reference. RDKit's own
       AllChem.ConstrainedEmbed has exactly this problem and solves it exactly
       this way, by aligning onto the core afterwards.

    2. `useRandomCoords=True` together with a coordMap returns zero conformers
       on RDKit >= 2024. On 2023.09 it returned conformers that were already in
       the reference frame, which hid problem 1 completely.

    So: embed without random coordinates, which works on every version tested,
    then align explicitly. Measured pin deviation after alignment is 0.09 A on
    RDKit 2025.09 and 0.16 A on 2023.09; without the alignment step it is 231 A
    on both, and nothing raises to tell you.

    `atom_map` is a list of (query_atom_idx, reference_atom_idx) pairs.
    Returns a list of conformer IDs, empty if every attempt failed.
    """
    conf_ids = []
    for random_coords in (False, True):
        try:
            conf_ids = list(AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                coordMap=coord_map,
                randomSeed=seed,
                useRandomCoords=random_coords,
                numThreads=args.cpu,
                enforceChirality=True,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
                ETversion=2,
            ))
        except Exception as exc:
            LOG.debug("embedding raised (useRandomCoords={0}): {1}".format(random_coords, exc))
            continue
        if conf_ids:
            break

    if not conf_ids:
        return []

    # Move every conformer onto the reference core. Without this the ligand is
    # embedded correctly but positioned nowhere near the binding site.
    for conf_id in conf_ids:
        try:
            rdMolAlign.AlignMol(mol, ref_mol, prbCid=conf_id, refCid=0, atomMap=atom_map)
        except Exception as exc:
            LOG.debug("core alignment failed for conformer {0}: {1}".format(conf_id, exc))
    return conf_ids


def generate_constrained_conformers(query_mol, ref_mol, ref_match, query_match, args, label):
    """
    Embed conformers with the core atoms pinned to the reference coordinates.

    Oversamples, minimises with the core held fixed, then applies the chair
    filter and optional RMS pruning.  Returns (mol_with_hs, conformer_ids).
    """
    mol = Chem.AddHs(query_mol)

    ref_conf = ref_mol.GetConformer()
    coord_map = {}
    for core_idx, query_idx in enumerate(query_match):
        coord_map[query_idx] = ref_conf.GetAtomPosition(ref_match[core_idx])
    LOG.info("  Constraining {0} atoms to the reference".format(len(coord_map)))

    # EmbedMultipleConfs clears existing conformers by default, so each attempt
    # replaces the previous set rather than extending it.
    atom_map = [(query_idx, ref_match[core_idx])
                for core_idx, query_idx in enumerate(query_match)]

    conf_ids = []
    for attempt in range(args.max_embed_attempts):
        seed = args.constrained_seed + attempt * 1000
        conf_ids = embed_with_core(
            mol, coord_map, args.n_constrained * args.oversample, seed, args,
            ref_mol, atom_map)
        if len(conf_ids) >= args.n_constrained:
            break

    if not conf_ids:
        return mol, []
    LOG.info("  Embedded {0} raw conformers".format(len(conf_ids)))

    # Relax everything except the pinned core. Two passes: the force field is
    # rebuilt between them so ring geometry converges from a clean state.
    for conf_id in conf_ids:
        for _ in range(2):
            ff = build_forcefield(mol, conf_id, args.ff)
            if ff is None:
                LOG.warn("{0}: no force field available for conformer {1}".format(label, conf_id))
                break
            for query_idx in query_match:
                ff.AddFixedPoint(query_idx)
            ff.Minimize(maxIts=args.minimize_iters)

    if args.chair_filter:
        conf_ids, _ = apply_chair_filter(mol, conf_ids, args.chair_theta_max, label)

    # Prune over the *unconstrained* heavy atoms only. The pinned core is
    # identical across conformers by construction, so including it just dilutes
    # the metric -- with 26 of 28 atoms pinned, a full rotation of the free
    # ethyl registers as barely 0.5 A when averaged over everything.
    constrained = set(query_match)
    free_heavy = [i for i in heavy_atom_indices(mol) if i not in constrained]
    if not free_heavy:
        free_heavy = heavy_atom_indices(mol)
    conf_ids = prune_by_rmsd(mol, conf_ids, args.prune_rms, free_heavy)
    return mol, conf_ids[: args.n_constrained]


def stage_prep(compounds, args, paths, obabel_path):
    """Build constrained conformers for every compound."""
    LOG.rule("STAGE 1/4  PREP -- constrained conformer generation")

    ref_mol = load_reference(args.center, args.reference_smiles)
    protonate = args.protonate == "obabel" and "prep" in args.protonate_stages
    LOG.info("Protonation: {0}".format(
        "obabel at pH {0}".format(args.ph) if protonate else "none (input SMILES as given)"))
    LOG.info("Chair filter: {0} (theta_max {1} deg)".format(
        "on" if args.chair_filter else "off", args.chair_theta_max))
    LOG.info("")

    records = []
    for compound in compounds:
        compound_id = compound["id"]
        LOG.info("{0}".format(compound_id))

        record = {
            "compound_id": compound_id,
            "input_smiles": compound["smiles"],
            "protomer_smiles": "",
            "protomer_source": "",
            "formal_charge": 0,
            "core_atoms": 0,
            "n_conformers": 0,
            "success": False,
            "error": "",
        }

        species_args = args if protonate else argparse.Namespace(**vars(args))
        if not protonate:
            species_args.protonate = "none"
        mol, species_smiles, charge, source = prepare_species(compound, species_args, obabel_path)

        record["protomer_smiles"] = species_smiles
        record["protomer_source"] = source
        record["formal_charge"] = charge
        compound["species_smiles"] = species_smiles
        compound["formal_charge"] = charge

        if mol is None:
            record["error"] = "SMILES did not parse"
            LOG.error("{0}: SMILES did not parse".format(compound_id))
            records.append(record)
            continue

        LOG.info("  Species: {0}  (charge {1:+d}, from {2})".format(species_smiles, charge, source))

        core, description = find_core(ref_mol, mol, args)
        if core is None:
            record["error"] = description
            LOG.error("{0}: {1}".format(compound_id, description))
            records.append(record)
            continue

        ref_match, query_match, how = choose_mapping(ref_mol, mol, core, args)
        if ref_match is None:
            record["error"] = how
            LOG.error("{0}: {1}".format(compound_id, how))
            records.append(record)
            continue

        LOG.info("  Core: {0} via {1}".format(description, how))
        record["core_atoms"] = len(query_match)

        if args.min_mcs_atoms and len(query_match) < args.min_mcs_atoms:
            record["error"] = "core has {0} atoms, below --min-mcs-atoms {1}".format(
                len(query_match), args.min_mcs_atoms)
            LOG.error("{0}: {1}".format(compound_id, record["error"]))
            records.append(record)
            continue

        mol_3d, conf_ids = generate_constrained_conformers(
            mol, ref_mol, ref_match, query_match, args, compound_id)
        if not conf_ids:
            record["error"] = "no conformers could be embedded"
            LOG.error("{0}: no conformers could be embedded".format(compound_id))
            records.append(record)
            continue

        out_dir = paths["conformers"] / safe_name(compound_id)
        if out_dir.exists():
            for stale in list(out_dir.glob("*.sdf")) + list(out_dir.glob("*.pdb")):
                stale.unlink()
        out_dir.mkdir(parents=True, exist_ok=True)

        for index, conf_id in enumerate(conf_ids, start=1):
            name = "{0}_conf_{1:02d}".format(safe_name(compound_id), index)
            write_conformer(mol_3d, conf_id, out_dir / "{0}.{1}".format(name, args.ligand_format),
                            name, args.ligand_format)

        record["n_conformers"] = len(conf_ids)
        record["success"] = True
        LOG.info("  Wrote {0} conformers to {1}/".format(len(conf_ids), out_dir))
        records.append(record)

    summary = pd.DataFrame(records)
    summary.to_csv(paths["prep_summary"], index=False)

    successful = int(summary["success"].sum())
    LOG.info("")
    LOG.info("Prep complete: {0}/{1} compounds, summary -> {2}".format(
        successful, len(records), paths["prep_summary"]))
    return summary


# --------------------------------------------------------------------------
# Stage 2 -- SMINA docking
# --------------------------------------------------------------------------

AFFINITY_RE = re.compile(r"^Affinity:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
RMSD_RE = re.compile(r"^RMSD:\s*(-?\d+(?:\.\d+)?)")
INTRA_RE = re.compile(r"^Intramolecular energy:\s*(-?\d+(?:\.\d+)?)")
MODE_ROW_RE = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")


def parse_smina_log(path):
    """
    Extract affinity / pose RMSD / intramolecular energy from a SMINA log.

    Output shape differs by mode, so all four are handled:
        --score_only   'Affinity: X (kcal/mol)' + 'Intramolecular energy: Y', no RMSD
        --minimize     'Affinity: X Y (kcal/mol)' + 'RMSD: Z'
        --local_only   as --minimize
        (docking)      a mode table; the first data row is the best mode
    """
    result = {"affinity": None, "pose_rmsd": None, "intramolecular": None}
    if not os.path.exists(path):
        return result

    with open(path, "r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            match = AFFINITY_RE.match(line)
            if match and result["affinity"] is None:
                result["affinity"] = float(match.group(1))
                continue
            match = RMSD_RE.match(line)
            if match and result["pose_rmsd"] is None:
                result["pose_rmsd"] = float(match.group(1))
                continue
            match = INTRA_RE.match(line)
            if match and result["intramolecular"] is None:
                result["intramolecular"] = float(match.group(1))
                continue
            if result["affinity"] is None:
                match = MODE_ROW_RE.match(line)
                if match:
                    result["affinity"] = float(match.group(1))
    return result


def build_smina_command(smina_path, receptor, ligand, out_path, log_path, args):
    """Assemble the SMINA argument vector for the configured mode."""
    cmd = [smina_path, "-r", str(receptor), "-l", str(ligand),
           "-o", str(out_path), "--log", str(log_path)]

    if args.box_center and args.box_size:
        cx, cy, cz = args.box_center
        sx, sy, sz = args.box_size
        cmd += ["--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz)]
    else:
        cmd += ["--autobox_ligand", str(args.center), "--autobox_add", str(args.autobox_add)]

    if args.mode == "minimize":
        cmd += ["--minimize", "--minimize_iters", str(args.smina_minimize_iters)]
    elif args.mode == "local_only":
        cmd += ["--local_only"]
    elif args.mode == "score_only":
        cmd += ["--score_only"]
    else:  # full docking search
        cmd += ["--exhaustiveness", str(args.exhaustiveness),
                "--num_modes", str(args.num_modes),
                "--energy_range", str(args.energy_range),
                "--min_rmsd_filter", str(args.min_rmsd_filter)]

    if args.scoring:
        cmd += ["--scoring", args.scoring]
    if args.smina_seed is not None:
        cmd += ["--seed", str(args.smina_seed)]
    if args.cpu:
        cmd += ["--cpu", str(args.cpu)]
    if not args.addh:
        cmd += ["--addH", "0"]
    if args.smina_extra:
        cmd += shlex.split(args.smina_extra)
    return cmd


def stage_dock(compounds, args, paths, smina_path):
    """Dock every constrained conformer and select the best pose per compound."""
    LOG.rule("STAGE 2/4  DOCK -- SMINA ({0})".format(args.mode))
    LOG.info("SMINA    : {0}".format(smina_path))
    LOG.info("Receptor : {0}".format(args.receptor))
    LOG.info("Box      : {0}".format(
        "explicit centre/size" if (args.box_center and args.box_size)
        else "autobox on {0} + {1} A".format(args.center, args.autobox_add)))
    if args.mode == "minimize":
        LOG.info("Iters    : --minimize_iters {0}".format(args.smina_minimize_iters))
    if args.max_pose_rmsd:
        LOG.info("Drift    : poses above {0} A from their constrained input are rejected".format(
            args.max_pose_rmsd))
    LOG.info("")

    conformer_files = []
    for compound in compounds:
        compound_dir = paths["conformers"] / safe_name(compound["id"])
        if not compound_dir.is_dir():
            continue
        for path in sorted(compound_dir.glob("*.{0}".format(args.ligand_format))):
            conformer_files.append((compound["id"], path))

    if not conformer_files:
        LOG.warn("no conformers found in {0}; did the prep stage run?".format(paths["conformers"]))
        return pd.DataFrame()

    rows = []
    total = len(conformer_files)
    for index, (compound_id, ligand_path) in enumerate(conformer_files, start=1):
        stem = ligand_path.stem
        out_path = paths["dock_all"] / "{0}_min.{1}".format(stem, args.ligand_format)
        log_path = paths["dock_all"] / "{0}.log".format(stem)
        cmd = build_smina_command(smina_path, args.receptor, ligand_path, out_path, log_path, args)

        if args.dry_run:
            LOG.info("[{0}/{1}] {2}".format(index, total, " ".join(shlex.quote(c) for c in cmd)))
            continue

        try:
            subprocess.run(cmd, capture_output=True, timeout=args.smina_timeout)
        except subprocess.TimeoutExpired:
            LOG.warn("{0}: SMINA timed out after {1}s".format(stem, args.smina_timeout))
            continue
        except OSError as exc:
            raise SystemExit("Error: could not run SMINA ({0})".format(exc))

        parsed = parse_smina_log(log_path)
        if parsed["affinity"] is None:
            LOG.warn("{0}: no affinity in the SMINA log".format(stem))
            continue

        drift = parsed["pose_rmsd"]
        rejected = bool(args.max_pose_rmsd and drift is not None and drift > args.max_pose_rmsd)
        rows.append({
            "compound_id": compound_id,
            "conformer_id": stem,
            "binding_affinity": parsed["affinity"],
            "pose_rmsd": drift,
            "intramolecular": parsed["intramolecular"],
            "pose_file": str(out_path),
            "log_file": str(log_path),
            "drift_rejected": rejected,
        })

        LOG.info("[{0}/{1}] {2}: {3:.4f} kcal/mol{4}{5}".format(
            index, total, stem, parsed["affinity"],
            "  drift {0:.3f} A".format(drift) if drift is not None else "",
            "  [REJECTED: drift]" if rejected else ""))

    if args.dry_run:
        LOG.info("\nDry run: {0} SMINA invocations planned.".format(total))
        return pd.DataFrame()

    results = pd.DataFrame(rows)
    if results.empty:
        LOG.warn("SMINA produced no usable results")
        return results
    results.to_csv(paths["dock_all_csv"], index=False)

    # Pick the best pose per compound, preferring poses that kept the constraint.
    LOG.info("")
    LOG.info("Selecting best pose per compound")
    best_rows = []
    for compound_id, group in results.groupby("compound_id", sort=False):
        eligible = group[~group["drift_rejected"]]
        if eligible.empty:
            LOG.warn("{0}: every pose exceeded --max-pose-rmsd; selecting from all poses".format(
                compound_id))
            eligible = group

        ranked = eligible.sort_values("binding_affinity", ascending=True)
        keep = ranked.head(args.top_n) if args.top_n and args.top_n > 0 else ranked.head(1)

        for rank, (_, row) in enumerate(keep.iterrows(), start=1):
            suffix = "best" if rank == 1 else "rank{0}".format(rank)
            pose_dst = paths["poses"] / "{0}_{1}.{2}".format(
                safe_name(compound_id), suffix, args.ligand_format)
            log_dst = paths["logs"] / "{0}_{1}.log".format(safe_name(compound_id), suffix)
            if os.path.exists(row["pose_file"]):
                shutil.copy2(row["pose_file"], pose_dst)
            if os.path.exists(row["log_file"]):
                shutil.copy2(row["log_file"], log_dst)

            entry = {
                "compound_id": compound_id,
                "conformer_id": row["conformer_id"],
                "binding_affinity": row["binding_affinity"],
                "pose_rmsd": row["pose_rmsd"],
                "rank": rank,
                "pose_file": str(pose_dst),
            }
            best_rows.append(entry)
            if rank == 1:
                LOG.info("  {0}: {1} at {2:.4f} kcal/mol{3}".format(
                    compound_id, row["conformer_id"], row["binding_affinity"],
                    "  (drift {0:.3f} A)".format(row["pose_rmsd"])
                    if row["pose_rmsd"] is not None and not pd.isna(row["pose_rmsd"]) else ""))

    summary = pd.DataFrame(best_rows)
    summary.to_csv(paths["docking_summary"], index=False)
    LOG.info("")
    LOG.info("Docking summary -> {0}".format(paths["docking_summary"]))

    if not args.keep_all_poses:
        for path in paths["dock_all"].glob("*.{0}".format(args.ligand_format)):
            path.unlink()

    return summary


# --------------------------------------------------------------------------
# Stage 3 -- solution-phase conformer ensemble
# --------------------------------------------------------------------------

def boltzmann_populations(energies, temperature):
    """
    Boltzmann populations (as percentages) from relative MMFF energies.

        p_i = exp(-dE_i / RT) / sum_j exp(-dE_j / RT),   dE_i = E_i - E_min
    """
    energies = np.asarray(energies, dtype=float)
    relative = energies - energies.min()
    factors = np.exp(-relative / (GAS_CONSTANT * temperature))
    return (factors / factors.sum()) * 100.0


def map_to_pose(mol, pose_mol):
    """
    Atom correspondence between the generated molecule and the docked pose.

    Both come from the same protomer, but SMINA round-trips the ligand through
    Open Babel, which may reorder atoms -- so match on the graph rather than
    trusting index order.  Returns a list of (mol_idx, pose_idx) heavy-atom
    pairs, or None if no consistent mapping exists.
    """
    mol_heavy = Chem.RemoveHs(Chem.Mol(mol))
    pose_heavy = Chem.RemoveHs(Chem.Mol(pose_mol))
    if mol_heavy.GetNumAtoms() != pose_heavy.GetNumAtoms():
        return None

    match = mol_heavy.GetSubstructMatch(pose_heavy)
    if not match:
        match = pose_heavy.GetSubstructMatch(mol_heavy)
        if not match:
            return None
        # Invert: the match maps pose atoms onto mol atoms.
        pairs = [(match[i], i) for i in range(len(match))]
    else:
        # match[i] is the mol atom corresponding to pose atom i.
        pairs = [(match[i], i) for i in range(len(match))]

    # Translate heavy-atom indices back to full-molecule indices.
    mol_heavy_ids = heavy_atom_indices(mol)
    pose_heavy_ids = heavy_atom_indices(pose_mol)
    if len(mol_heavy_ids) != len(pose_heavy_ids):
        return None
    try:
        return [(mol_heavy_ids[m], pose_heavy_ids[p]) for m, p in pairs]
    except IndexError:
        return None


def generate_solution_conformers(compound, args, paths, pose_path):
    """
    Build the unconstrained ensemble for one compound and score it.

    Explicit hydrogens are mandatory here: MMFF94 is an all-atom force field, so
    both the energies and the Boltzmann populations derived from them are
    meaningless without them.
    """
    compound_id = compound["id"]
    smiles = compound.get("species_smiles") or compound["smiles"]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        LOG.error("{0}: solution SMILES did not parse".format(compound_id))
        return None

    pose_mol = read_molecule(pose_path, remove_hs=False) if pose_path else None
    if pose_path and pose_mol is None:
        LOG.warn("{0}: could not read docked pose {1}".format(compound_id, pose_path))

    mol = Chem.AddHs(mol)

    params = getattr(AllChem, "ETKDGv{0}".format(args.etkdg_version), AllChem.ETKDGv3)()
    params.randomSeed = args.solution_seed
    params.useRandomCoords = True
    params.numThreads = args.cpu
    try:
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=args.n_solution, params=params))
    except Exception as exc:
        LOG.warn("{0}: solution embedding failed ({1})".format(compound_id, exc))
        return None
    if not conf_ids:
        LOG.warn("{0}: no solution conformers generated".format(compound_id))
        return None

    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=args.cpu, maxIters=args.minimize_iters)

    if args.solution_chair_filter:
        conf_ids, _ = apply_chair_filter(mol, conf_ids, args.chair_theta_max, compound_id)

    atom_map = map_to_pose(mol, pose_mol) if pose_mol is not None else None
    if pose_mol is not None and atom_map is None:
        LOG.warn("{0}: could not map conformers onto the docked pose; RMSD unavailable".format(
            compound_id))

    out_dir = paths["solution"] / safe_name(compound_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, conf_id in enumerate(conf_ids, start=1):
        ff = build_forcefield(mol, conf_id, "mmff")
        energy = ff.CalcEnergy() if ff is not None else float("nan")

        rmsd = None
        if atom_map:
            try:
                # AlignMol moves the probe conformer onto the pose, so this must
                # happen before the conformer is written out.
                rmsd = rdMolAlign.AlignMol(mol, pose_mol, prbCid=conf_id, refCid=0, atomMap=atom_map)
            except Exception as exc:
                LOG.debug("{0} conf {1}: alignment failed ({2})".format(compound_id, index, exc))

        name = "{0}_conf_{1}".format(safe_name(compound_id), index)
        write_conformer(mol, conf_id, out_dir / "{0}.{1}".format(name, args.ligand_format),
                        name, args.ligand_format)

        row = {"Conformer": index, "MMFF_Energy_kcal_mol": energy}
        if rmsd is not None:
            row["RMSD_to_Docked_Pose"] = rmsd
        rows.append(row)

    populations = boltzmann_populations([r["MMFF_Energy_kcal_mol"] for r in rows], args.temperature)
    for row, population in zip(rows, populations):
        row["Boltzmann_Population_%"] = population

    energies_csv = out_dir / "{0}_energies.csv".format(safe_name(compound_id))
    pd.DataFrame(rows).to_csv(energies_csv, index=False)

    LOG.info("  {0}: {1} conformers -> {2}".format(compound_id, len(rows), energies_csv.name))
    return energies_csv


def stage_solution(compounds, args, paths):
    """Generate and score the solution ensemble for every compound."""
    LOG.rule("STAGE 3/4  SOLUTION -- unconstrained conformer ensemble")
    LOG.info("Conformers  : {0} per compound (ETKDGv{1}, seed {2})".format(
        args.n_solution, args.etkdg_version, args.solution_seed))
    LOG.info("Temperature : {0} K".format(args.temperature))
    LOG.info("Chair filter: {0}".format("on" if args.solution_chair_filter else "off"))
    LOG.info("Hydrogens   : explicit (required by MMFF94)")
    LOG.info("")

    for compound in compounds:
        pose_path = paths["poses"] / "{0}_best.{1}".format(
            safe_name(compound["id"]), args.ligand_format)
        generate_solution_conformers(
            compound, args, paths, pose_path if pose_path.exists() else None)


# --------------------------------------------------------------------------
# Stage 4 -- enrichment
# --------------------------------------------------------------------------

def summarize_ensemble(energies_csv):
    """
    Reduce one compound's ensemble to the four summary values.

    Returns None when the file lacks the RMSD column, which happens whenever the
    solution stage ran without a docked pose to align against.
    """
    df = pd.read_csv(energies_csv)
    if "MMFF_Energy_kcal_mol" not in df.columns or "RMSD_to_Docked_Pose" not in df.columns:
        return None
    df = df.dropna(subset=["MMFF_Energy_kcal_mol", "RMSD_to_Docked_Pose"])
    if df.empty:
        return None

    lowest_energy = df.loc[df["MMFF_Energy_kcal_mol"].idxmin()]
    closest_pose = df.loc[df["RMSD_to_Docked_Pose"].idxmin()]
    return {
        "Min E": float(lowest_energy["MMFF_Energy_kcal_mol"]),
        "Min E RMSD": float(lowest_energy["RMSD_to_Docked_Pose"]),
        "Nearest E": float(closest_pose["MMFF_Energy_kcal_mol"]),
        "Nearest RMSD": float(closest_pose["RMSD_to_Docked_Pose"]),
    }


def stage_enrich(args, paths):
    """Append the solution-ensemble columns to the docking summary."""
    LOG.rule("STAGE 4/4  ENRICH -- merge solution ensemble into the docking summary")

    if not paths["docking_summary"].exists():
        LOG.warn("no docking summary at {0}; nothing to enrich".format(paths["docking_summary"]))
        return None

    summary = pd.read_csv(paths["docking_summary"])
    if "compound_id" not in summary.columns:
        LOG.error("docking summary has no 'compound_id' column")
        return None

    values = {}
    missing = []
    for compound_id in summary["compound_id"].astype(str):
        energies_csv = paths["solution"] / safe_name(compound_id) / "{0}_energies.csv".format(
            safe_name(compound_id))
        stats = summarize_ensemble(energies_csv) if energies_csv.exists() else None
        if stats is None:
            missing.append(compound_id)
            stats = dict((column, float("nan")) for column in CONFORMER_COLUMNS)
        values[compound_id] = stats

    for column in CONFORMER_COLUMNS:
        summary[column] = [
            round(values[str(cid)][column], args.precision) for cid in summary["compound_id"]
        ]

    summary.to_csv(paths["docking_summary"], index=False)
    LOG.info("Enriched {0}/{1} compounds with: {2}".format(
        len(summary) - len(missing), len(summary), ", ".join(CONFORMER_COLUMNS)))
    if missing:
        LOG.warn("no usable solution ensemble for: {0}".format(", ".join(missing)))
    LOG.info("Summary -> {0}".format(paths["docking_summary"]))
    return summary


# --------------------------------------------------------------------------
# Visualisation
# --------------------------------------------------------------------------

def run_pymol(pml_path, pymol_path):
    """Execute a .pml headlessly. Never imports pymol -- the binary may live in another env."""
    if not pymol_path:
        LOG.info("  PyMOL not found; run manually: pymol -cq {0}".format(pml_path))
        return False
    try:
        proc = subprocess.run([pymol_path, "-cq", str(pml_path)], capture_output=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        LOG.warn("PyMOL failed on {0} ({1})".format(pml_path, exc))
        return False
    if proc.returncode != 0:
        LOG.warn("PyMOL exited {0} on {1}".format(proc.returncode, pml_path))
        return False
    return True


def build_overview_session(args, paths, pymol_path):
    """Receptor + reference core + every best pose, in one session."""
    poses = sorted(paths["poses"].glob("*_best.{0}".format(args.ligand_format)))
    if not poses:
        LOG.info("  No best poses found; skipping the overview session.")
        return

    affinities = {}
    if paths["docking_summary"].exists():
        summary = pd.read_csv(paths["docking_summary"])
        if {"compound_id", "binding_affinity"}.issubset(summary.columns):
            for _, row in summary.iterrows():
                affinities[safe_name(row["compound_id"])] = row["binding_affinity"]

    lines = [
        "# Docking overview session (generated by constrained_smina_docking.py)",
        "delete all",
        "bg_color white",
        "load {0}, receptor".format(os.path.abspath(args.receptor)),
        "hide everything, receptor",
        "show cartoon, receptor",
        "color gray70, receptor",
        "set cartoon_transparency, 0.3, receptor",
        "load {0}, reference".format(os.path.abspath(args.center)),
        "hide everything, reference",
        "show sticks, reference",
        "color yellow, reference",
        "",
    ]

    names = []
    for index, pose in enumerate(poses):
        name = pose.stem.replace("_best", "")
        names.append(name)
        affinity = affinities.get(name)
        label = "{0:.4f} kcal/mol".format(affinity) if affinity is not None else "n/a"
        lines += [
            "# {0} ({1})".format(name, label),
            "load {0}, {1}".format(os.path.abspath(str(pose)), name),
            "hide everything, {0}".format(name),
            "show sticks, {0}".format(name),
            "color {0}, {1}".format(PYMOL_COLORS[index % len(PYMOL_COLORS)], name),
            "set stick_radius, 0.12, {0}".format(name),
            "",
        ]

    lines += [
        "group docked_ligands, {0}".format(" ".join(names)),
        "zoom reference, 10",
        "set antialias, 2",
        "set ray_shadows, 0",
        "disable reference",
        "save {0}".format(os.path.abspath(str(paths["overview_pse"]))),
    ]

    pml_path = paths["pymol"] / "overview.pml"
    pml_path.write_text("\n".join(lines) + "\n")
    if run_pymol(pml_path, pymol_path):
        LOG.info("  Overview session -> {0}".format(paths["overview_pse"]))
    else:
        LOG.info("  Overview script -> {0}".format(pml_path))


def build_compound_sessions(compounds, args, paths, pymol_path):
    """Per-compound session: docked pose plus its aligned solution ensemble."""
    for compound in compounds:
        name = safe_name(compound["id"])
        solution_dir = paths["solution"] / name
        pose = paths["poses"] / "{0}_best.{1}".format(name, args.ligand_format)
        if not solution_dir.is_dir() or not pose.exists():
            continue

        conformers = sorted(solution_dir.glob("{0}_conf_*.{1}".format(name, args.ligand_format)))
        if not conformers:
            continue

        lines = [
            "# {0}: docked pose (green) + solution ensemble (cyan)".format(name),
            "delete all",
            "bg_color white",
            "load {0}, {1}_docked".format(os.path.abspath(str(pose)), name),
            "color green, {0}_docked".format(name),
            "show sticks, {0}_docked".format(name),
            "",
        ]
        for conformer in conformers:
            lines.append("load {0}, {1}".format(os.path.abspath(str(conformer)), conformer.stem))
        lines += [
            "",
            "color cyan, {0}_conf_*".format(name),
            "show sticks, {0}_conf_*".format(name),
            "set stick_radius, 0.15",
            "set stick_transparency, 0.3, {0}_conf_*".format(name),
            "disable {0}_conf_*".format(name),
            "zoom {0}_docked".format(name),
            "save {0}".format(os.path.abspath(str(solution_dir / "{0}_session.pse".format(name)))),
        ]

        pml_path = solution_dir / "{0}_session.pml".format(name)
        pml_path.write_text("\n".join(lines) + "\n")
        run_pymol(pml_path, pymol_path)


def stage_visualise(compounds, args, paths, pymol_path):
    LOG.rule("VISUALISATION -- PyMOL sessions")
    build_overview_session(args, paths, pymol_path)
    if args.pymol_per_compound:
        build_compound_sessions(compounds, args, paths, pymol_path)
        LOG.info("  Per-compound sessions -> {0}/<compound>/".format(paths["solution"]))


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def add_bool(group, name, dest, default, help_text):
    """argparse.BooleanOptionalAction is 3.9+, so pair the flags by hand."""
    exclusive = group.add_mutually_exclusive_group()
    exclusive.add_argument("--" + name, dest=dest, action="store_true", default=default,
                           help=help_text + (" (default)" if default else ""))
    exclusive.add_argument("--no-" + name, dest=dest, action="store_false",
                           help=argparse.SUPPRESS if default else "disable --" + name)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="constrained_smina_docking.py",
        description="Constrained SMINA docking: prep -> dock -> solution -> enrich.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # full pipeline
  python constrained_smina_docking.py -i target_compounds.csv -r receptor.pdb -C center.pdb -o results

  # prep only, more conformers, tighter RMS pruning
  python constrained_smina_docking.py -i in.csv -r receptor.pdb -C center.pdb \\
      --stages prep --n-constrained 30 --prune-rms 0.5

  # re-dock existing conformers with a full search instead of minimisation
  python constrained_smina_docking.py -i in.csv -r receptor.pdb -C center.pdb \\
      --stages dock,enrich --mode dock --exhaustiveness 16
""")

    core = parser.add_argument_group("core input/output")
    core.add_argument("-i", "--input", required=True, help="input CSV of compound IDs and SMILES")
    core.add_argument("-r", "--receptor", required=True, help="receptor PDB (rigid)")
    core.add_argument("-C", "--center", required=True,
                      help="reference core PDB; defines both the MCS target and the autobox")
    core.add_argument("-o", "--output-dir", default="results", help="output root (default: results)")
    core.add_argument("--id-col", default=None,
                      help="compound ID column; auto-detected from {0}".format(
                          ", ".join(repr(c) for c in ID_COL_CANDIDATES)))
    core.add_argument("--smiles-col", default=None, help="SMILES column (auto-detected)")
    core.add_argument("--ligand-format", choices=("sdf", "pdb"), default="sdf",
                      help="ligand file format; SDF preserves bond orders (default: sdf)")

    flow = parser.add_argument_group("pipeline control")
    flow.add_argument("--stages", default=",".join(STAGES),
                      help="comma-separated subset of {0} (default: all)".format(",".join(STAGES)))
    flow.add_argument("--resume", action="store_true",
                      help="skip any stage whose primary output already exists")
    flow.add_argument("--dry-run", action="store_true",
                      help="print the SMINA commands that would run, without running them")
    flow.add_argument("--seed", type=int, default=42,
                      help="master random seed; per-stage seeds derive from it (default: 42)")
    flow.add_argument("-j", "--cpu", type=int, default=0,
                      help="worker threads; 0 lets RDKit/SMINA decide (default: 0)")
    flow.add_argument("-v", "--verbose", action="store_true")
    flow.add_argument("--quiet", action="store_true")

    tools = parser.add_argument_group("external tools (auto-detected on PATH and in conda envs)")
    tools.add_argument("--smina-path", default=None)
    tools.add_argument("--obabel-path", default=None)
    tools.add_argument("--pymol-path", default=None)

    prot = parser.add_argument_group("protonation (applied consistently across stages)")
    prot.add_argument("--protonate", choices=("obabel", "none"), default="obabel",
                      help="assign a dominant protomer with obabel (default: obabel)")
    prot.add_argument("--ph", type=float, default=7.4, help="pH for protonation (default: 7.4)")
    prot.add_argument("--protomer-col", default=None,
                      help="CSV column holding a per-compound protomer SMILES that overrides obabel")
    prot.add_argument("--protonate-stages", default="prep,solution",
                      help="stages the protomer applies to (default: prep,solution)")

    prep = parser.add_argument_group("stage: prep")
    prep.add_argument("--n-constrained", type=int, default=10,
                      help="constrained conformers kept per compound (default: 10)")
    prep.add_argument("--constrained-seed", type=int, default=None,
                      help="override the derived prep seed")
    prep.add_argument("--mcs-timeout", type=int, default=5, help="MCS search timeout, s (default: 5)")
    prep.add_argument("--mcs-atom-compare", choices=tuple(MCS_ATOM_COMPARE), default="elements")
    prep.add_argument("--mcs-bond-compare", choices=tuple(MCS_BOND_COMPARE), default="order")
    prep.add_argument("--mcs-smarts", default=None,
                      help="use this SMARTS as the constrained core instead of running MCS")
    prep.add_argument("--min-mcs-atoms", type=int, default=0,
                      help="fail a compound whose core is smaller than this (default: 0, off)")
    prep.add_argument("--reference-smiles", default=None,
                      help="template SMILES restoring bond orders on the reference PDB")
    prep.add_argument("--max-symmetry-trials", type=int, default=8,
                      help="candidate mappings tried by --all-symmetry-matches (default: 8)")
    prep.add_argument("--oversample", type=int, default=3,
                      help="embed this multiple of --n-constrained before filtering (default: 3)")
    prep.add_argument("--max-embed-attempts", type=int, default=3)
    prep.add_argument("--ff", choices=("mmff", "uff"), default="mmff",
                      help="force field for constrained minimisation (default: mmff)")
    prep.add_argument("--minimize-iters", type=int, default=1000,
                      help="RDKit force-field iterations (default: 1000)")
    prep.add_argument("--prune-rms", type=float, default=0.0,
                      help="drop conformers within this RMSD of a kept one, in A. Measured over "
                           "the unconstrained heavy atoms only, since the pinned core is identical "
                           "by construction (default: 0, off)")
    add_bool(prep, "chair-filter", "chair_filter", True,
             "reject non-chair saturated 6-rings via Cremer-Pople analysis")
    add_bool(prep, "all-symmetry-matches", "all_symmetry_matches", False,
             "trial every symmetry-equivalent core mapping and keep the best")
    prep.add_argument("--chair-theta-max", type=float, default=30.0,
                      help="max Cremer-Pople theta deviation from an ideal chair (default: 30)")
    add_bool(prep, "mcs-ring-matches-ring-only", "mcs_ring_matches_ring_only", True,
             "require ring atoms to match ring atoms")
    add_bool(prep, "mcs-complete-rings-only", "mcs_complete_rings_only", True,
             "only allow whole rings in the MCS")

    dock = parser.add_argument_group("stage: dock")
    dock.add_argument("--mode", choices=("minimize", "local_only", "score_only", "dock"),
                      default="minimize",
                      help="SMINA mode; minimize is correct for pre-positioned poses (default)")
    dock.add_argument("--smina-minimize-iters", type=int, default=10,
                      help="SMINA --minimize_iters. SMINA's own default of 0 drifts furthest from "
                           "the constrained input; 10 trades ~1 kcal/mol for ~1 A less drift "
                           "(default: 10)")
    dock.add_argument("--max-pose-rmsd", type=float, default=0.0,
                      help="reject poses that drifted more than this from their constrained input, "
                           "in A (default: 0, off)")
    dock.add_argument("--autobox-add", type=float, default=4.0,
                      help="autobox padding in A (default: 4)")
    dock.add_argument("--box-center", type=float, nargs=3, metavar=("X", "Y", "Z"),
                      help="explicit box centre, replacing the autobox")
    dock.add_argument("--box-size", type=float, nargs=3, metavar=("X", "Y", "Z"),
                      help="explicit box size, replacing the autobox")
    dock.add_argument("--scoring", default=None, help="alternative SMINA scoring function")
    dock.add_argument("--exhaustiveness", type=int, default=8, help="dock mode only (default: 8)")
    dock.add_argument("--num-modes", type=int, default=9, help="dock mode only (default: 9)")
    dock.add_argument("--energy-range", type=float, default=3.0, help="dock mode only (default: 3)")
    dock.add_argument("--min-rmsd-filter", type=float, default=1.0, help="dock mode only (default: 1)")
    dock.add_argument("--smina-seed", type=int, default=None)
    dock.add_argument("--smina-extra", default=None, help="raw extra arguments passed to SMINA")
    dock.add_argument("--smina-timeout", type=int, default=3600,
                      help="per-conformer SMINA timeout, s (default: 3600)")
    dock.add_argument("--top-n", type=int, default=1,
                      help="poses retained per compound (default: 1)")
    add_bool(dock, "keep-all-poses", "keep_all_poses", False,
             "keep every minimised pose, not just the retained ones")
    add_bool(dock, "addh", "addh", True, "let SMINA add hydrogens to the ligand")

    sol = parser.add_argument_group("stage: solution")
    sol.add_argument("--n-solution", type=int, default=10,
                     help="solution conformers per compound (default: 10)")
    sol.add_argument("--solution-seed", type=int, default=None,
                     help="override the derived solution seed")
    sol.add_argument("--temperature", type=float, default=298.15,
                     help="temperature for Boltzmann populations, K (default: 298.15)")
    sol.add_argument("--etkdg-version", type=int, choices=(1, 2, 3), default=3,
                     help="ETKDG version for solution embedding (default: 3)")
    add_bool(sol, "solution-chair-filter", "solution_chair_filter", False,
             "apply the chair filter to the solution ensemble too. Off by default: it truncates "
             "the ensemble the Boltzmann populations normalise over, and inflates Nearest RMSD "
             "when the docked pose sits on a non-chair ring")

    post = parser.add_argument_group("stage: enrich and visualisation")
    post.add_argument("--precision", type=int, default=5,
                      help="decimal places for appended columns (default: 5)")
    add_bool(post, "pymol", "pymol", True, "build PyMOL sessions")
    add_bool(post, "pymol-per-compound", "pymol_per_compound", True,
             "also build one session per compound")

    return parser


def resolve_args(parser):
    args = parser.parse_args()

    requested = [s.strip() for s in args.stages.split(",") if s.strip()]
    unknown = [s for s in requested if s not in STAGES]
    if unknown:
        parser.error("unknown stage(s): {0}. Choose from {1}".format(
            ", ".join(unknown), ", ".join(STAGES)))
    args.stage_set = requested

    protonate_stages = [s.strip() for s in args.protonate_stages.split(",") if s.strip()]
    args.protonate_stages = protonate_stages if protonate_stages != ["all"] else list(STAGES)

    # Derive per-stage seeds from the master seed unless explicitly overridden.
    if args.constrained_seed is None:
        args.constrained_seed = args.seed
    if args.solution_seed is None:
        args.solution_seed = args.seed + 513

    if bool(args.box_center) != bool(args.box_size):
        parser.error("--box-center and --box-size must be given together")

    for label, path in (("--receptor", args.receptor), ("--center", args.center)):
        if not os.path.exists(path):
            parser.error("{0} not found: {1}".format(label, path))

    return args


def make_paths(args):
    root = Path(args.output_dir)
    paths = {
        "root": root,
        "conformers": root / "conformers",
        "docking": root / "docking",
        "dock_all": root / "docking" / "all",
        "poses": root / "docking" / "poses",
        "logs": root / "docking" / "logs",
        "solution": root / "solution",
        "pymol": root / "pymol",
        "prep_summary": root / "preparation_summary.csv",
        "docking_summary": root / "docking_summary.csv",
        "dock_all_csv": root / "docking" / "all_results.csv",
        "overview_pse": root / "pymol" / "docking_overview.pse",
    }
    for key in ("conformers", "dock_all", "poses", "logs", "solution", "pymol"):
        paths[key].mkdir(parents=True, exist_ok=True)
    return paths


def main():
    global LOG
    parser = build_parser()
    args = resolve_args(parser)
    LOG = Log(verbose=args.verbose, quiet=args.quiet)

    LOG.rule("CONSTRAINED SMINA DOCKING")
    LOG.info("Stages     : {0}".format(", ".join(args.stage_set)))
    LOG.info("Output root: {0}".format(args.output_dir))
    LOG.info("")

    compounds = load_compounds(args)
    LOG.info("Loaded {0} compound(s)".format(len(compounds)))
    LOG.info("")

    paths = make_paths(args)

    obabel_path = None
    if args.protonate == "obabel":
        obabel_path = find_tool("obabel", args.obabel_path, required=False)
        if obabel_path is None:
            LOG.warn("obabel not found; falling back to --protonate none")
            args.protonate = "none"

    smina_path = None
    if "dock" in args.stage_set:
        smina_path = find_tool("smina", args.smina_path, required=True)

    pymol_path = None
    if args.pymol:
        pymol_path = find_tool("pymol", args.pymol_path, required=False)

    if "prep" in args.stage_set:
        if args.resume and paths["prep_summary"].exists():
            LOG.info("Resume: prep output exists, skipping.\n")
        else:
            stage_prep(compounds, args, paths, obabel_path)
            LOG.info("")

    # Later stages need the protomer even when prep was skipped.
    for compound in compounds:
        if "species_smiles" not in compound:
            protonate = args.protonate == "obabel" and "solution" in args.protonate_stages
            species_args = argparse.Namespace(**vars(args))
            if not protonate:
                species_args.protonate = "none"
            _, species_smiles, charge, _ = prepare_species(compound, species_args, obabel_path)
            compound["species_smiles"] = species_smiles
            compound["formal_charge"] = charge

    if "dock" in args.stage_set:
        if args.resume and paths["docking_summary"].exists():
            LOG.info("Resume: docking summary exists, skipping.\n")
        else:
            stage_dock(compounds, args, paths, smina_path)
            LOG.info("")

    if "solution" in args.stage_set and not args.dry_run:
        stage_solution(compounds, args, paths)
        LOG.info("")

    if "enrich" in args.stage_set and not args.dry_run:
        stage_enrich(args, paths)
        LOG.info("")

    if args.pymol and not args.dry_run:
        stage_visualise(compounds, args, paths, pymol_path)
        LOG.info("")

    LOG.rule("DONE")
    LOG.info("Results: {0}".format(paths["root"]))
    if paths["docking_summary"].exists():
        LOG.info("Summary: {0}".format(paths["docking_summary"]))
    if LOG.warnings:
        LOG.info("")
        LOG.info("{0} warning(s) were emitted; see stderr above.".format(len(LOG.warnings)))


if __name__ == "__main__":
    main()
