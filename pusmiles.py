"""
PUSMILES: Polyurethane SMILES for Materials Informatics
=======================================================

A notation for encoding molecular structures with weight percent composition.

Syntax:
    SMILES@wt%                     Pure compound (100%)
    SMILES@wt%{key=val}            With metadata
    SMILES1@wt1|SMILES2@wt2|...    Mixture/formulation

Examples:
    CCO@100                        Pure ethanol
    CCO@50|O@50                    50/50 ethanol/water
    *CC*@100{Mn=50000}             Polyethylene with MW
    *OCC(C)*@55{role=soft}|...     PU formulation (soft/hard segments)

Key design principles:
    - Minimal output: only include metadata that was explicitly provided
    - Polymer-agnostic: works for any polymer, not just polyurethanes
    - Composition-focused: direct wt% for ML models
    - RDKit-compatible: uses * for attachment points (dummy atoms)

Isocyanate hard-segment SMILES convention:
    * represents the polyol/chain-extender oxygen attachment point.
    Pattern: *OC(=O)N-[isocyanate core]-NC(=O)O*
    N remains bonded to the isocyanate backbone; O faces the polyol/chain-extender.

Chain-extender segment SMILES convention:
    * represents the isocyanate nitrogen attachment point.
    Pattern: *NC(=O)O-[diol core]-OC(=O)N*
    These two conventions are complementary: the O on the isocyanate segment
    and the N on the chain-extender segment meet at the urethane bond.
"""

import re
import json
import warnings
import urllib.request
import urllib.error
import urllib.parse
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# =============================================================================
# EXCEPTIONS
# =============================================================================

class PUSMILESValidationError(ValueError):
    """Raised when a PUSMILES object has one or more validation errors."""
    def __init__(self, message: str, result: "ValidationResult"):
        super().__init__(message)
        self.result = result


# =============================================================================
# VALIDATION TYPES
# =============================================================================

@dataclass
class ValidationIssue:
    """A single validation issue found during PUSMILES validation."""
    severity: str                           # "error" or "warning"
    message: str
    component_index: Optional[int] = None  # None = formulation-level issue
    component_smiles: Optional[str] = None

    def __str__(self):
        if self.component_index is not None:
            loc = f"component[{self.component_index}] ({self.component_smiles})"
        else:
            loc = "formulation"
        return f"[{self.severity.upper()}] {loc}: {self.message}"


@dataclass
class ValidationResult:
    """
    Result of calling PUSMILES.validate().

    Attributes:
        is_valid:  True if there are no errors. Warnings do not affect validity.
        issues:    All issues found (errors and warnings combined).
    """
    is_valid: bool
    issues: List[ValidationIssue]

    def __bool__(self) -> bool:
        return self.is_valid

    @property
    def errors(self) -> List[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> List[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    def print_report(self) -> None:
        """Print a formatted validation report to stdout."""
        if self.is_valid and not self.warnings:
            print("PUSMILES validation passed: no issues found.")
            return
        status = "PASSED" if self.is_valid else "FAILED"
        print(f"PUSMILES validation {status}  "
              f"({len(self.errors)} error(s), {len(self.warnings)} warning(s))")
        print("-" * 60)
        for issue in self.issues:
            print(f"  {issue}")
        if not self.is_valid:
            print("-" * 60)
            print("Fix all errors before using this PUSMILES in a pipeline.")


# =============================================================================
# POLYMER LIBRARY
# =============================================================================
#
# Each entry: name -> SMILES string
#
# Isocyanate hard-segment pattern:  *OC(=O)N-[core]-NC(=O)O*
# Chain-extender segment pattern:   *NC(=O)O-[diol]-OC(=O)N*
# Polymer repeat units:             *-[repeat unit]-*
# Small molecules (CE/crosslinker): full SMILES, no *
#
# NOTE on graft/polymer polyols:
#   Graft polyols (polymer polyols) are colloidal dispersions of SAN
#   (styrene-acrylonitrile) particles grafted into a polyether matrix.
#   They cannot be represented by a single repeat-unit SMILES. The entries
#   below encode the polyether backbone only. Use name= metadata to record
#   the commercial grade (e.g. name=Arcol_HS-100) and fn= for functionality.

POLYMER_REPEAT_UNITS: Dict[str, str] = {

    # -------------------------------------------------------------------------
    # COMMODITY POLYMERS
    # -------------------------------------------------------------------------
    "PE":    "*CC*",
    "PP":    "*CC(C)*",
    "PS":    "*CC(c1ccccc1)*",
    "PVC":   "*CC(Cl)*",
    "PMMA":  "*CC(C)(C(=O)OC)*",
    "PAN":   "*CC(C#N)*",
    "PVA":   "*CC(O)*",
    "PVDF":  "*CC(F)(F)*",
    "PVAC":  "*CC(OC(C)=O)*",

    # -------------------------------------------------------------------------
    # POLYETHER POLYOLS  (soft segments)
    # -------------------------------------------------------------------------
    "PEG":   "*OCC*",           # Polyethylene glycol (EO-based)
    "PPG":   "*OCC(C)*",        # Polypropylene glycol (PO-based); dominant flexible-foam polyol
    "PTMEG": "*OCCCC*",         # Polytetramethylene ether glycol (THF-based); premier TPU soft segment
    "PBO":   "*OCC(CC)*",       # Poly(1,2-butylene oxide); hydrophobic alternative to PPG

    # Graft/polymer polyol backbones (SAN phase not representable in SMILES)
    "graft_polyol_PPG": "*OCC(C)*",   # SAN-grafted PPG backbone (high-load flexible foam)
    "graft_polyol_PEG": "*OCC*",      # EO-capped graft polyol backbone

    # -------------------------------------------------------------------------
    # POLYESTER POLYOLS  (soft segments)
    # -------------------------------------------------------------------------
    # Adipoyl core = C(=O)CCCCC(=O): 2 carbonyls flanking 4 methylenes (adipic acid)
    # Succinyl core = C(=O)CCC(=O):  2 carbonyls flanking 2 methylenes (succinic acid)

    "PBA":  "*OCCCCOC(=O)CCCCC(=O)*",      # Poly(butylene adipate)       BDO + adipic acid
    "PEA":  "*OCCOC(=O)CCCCC(=O)*",        # Poly(ethylene adipate)       EG + adipic acid
    "PHA":  "*OCCCCCCOC(=O)CCCCC(=O)*",    # Poly(hexamethylene adipate)  HDO + adipic acid
    "PCL":  "*OCCCCCC(=O)*",               # Polycaprolactone             ring-opening of epsilon-caprolactone
    "PES":  "*OCCOC(=O)CCC(=O)*",          # Poly(ethylene succinate)     EG + succinic acid
    "PET":  "*OC(=O)c1ccc(C(=O)OCCO)cc1*", # Polyethylene terephthalate
    "PLA":  "*OC(C)C(=O)*",                # Polylactic acid

    # -------------------------------------------------------------------------
    # POLYCARBONATE DIOLS  (soft segments - best hydrolysis/oxidation resistance)
    # -------------------------------------------------------------------------
    # Carbonate linkage = -O-C(=O)-O-  (more hydrolysis-resistant than ester)
    "PCDL": "*OCCCCCCOC(=O)*",    # Poly(1,6-hexanediol carbonate) diol; preferred for medical/auto TPU

    # -------------------------------------------------------------------------
    # POLYAMIDES
    # -------------------------------------------------------------------------
    "nylon6":  "*NCCCCC(=O)*",
    "nylon66": "*NCCCCCCNC(=O)CCCCC(=O)*",

    # -------------------------------------------------------------------------
    # CELLULOSICS
    # -------------------------------------------------------------------------
    "cellulose": "*OC[C@H]1O[C@@H](O*)[C@H](O)[C@@H](O)[C@@H]1O",

    # =========================================================================
    # POLYURETHANE HARD SEGMENTS  (reacted isocyanate repeat units)
    # =========================================================================
    #
    # Convention: * = attachment to polyol/chain-extender oxygen
    # Pattern:    *OC(=O)N-[isocyanate core]-NC(=O)O*
    #
    # Urethane formation: R-NCO + HO-R' -> R-NH-C(=O)-O-R'
    # N stays bonded to the isocyanate core; O connects to the polyol/CE.
    # Therefore * (representing the polyol/CE side) is on oxygen.

    # Aromatic diisocyanates
    "MDI_urethane":
        "*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1",
        # 4,4'-MDI; methylene-bridged diphenyl; CAS 101-68-8

    "TDI_urethane":
        "*OC(=O)Nc1ccc(C)c(NC(=O)O*)c1",
        # 2,4-TDI (dominant commercial isomer); CAS 584-84-9

    "TDI_24_urethane":
        "*OC(=O)Nc1ccc(C)c(NC(=O)O*)c1",
        # 2,4-TDI explicit alias

    "TDI_26_urethane":
        "*OC(=O)Nc1c(C)cccc1NC(=O)O*",
        # 2,6-TDI; CH3 flanked by both urethane nitrogens; CAS 91-08-7

    "NDI_urethane":
        "*OC(=O)Nc1cccc2c(NC(=O)O*)cccc12",
        # 1,5-Naphthalene diisocyanate; rigid fused-ring; CAS 3173-72-6

    "PPDI_urethane":
        "*OC(=O)Nc1ccc(NC(=O)O*)cc1",
        # p-Phenylene diisocyanate; smallest symmetric aromatic; CAS 104-49-4

    "TODI_urethane":
        "*OC(=O)Nc1ccc(-c2ccc(NC(=O)O*)c(C)c2)cc1C",
        # 3,3'-Dimethylbiphenyl-4,4'-diisocyanate; biphenyl core; CAS 91-97-4

    # Araliphatic diisocyanates (N on benzylic CH2, not ring -> non-yellowing)
    "mXDI_urethane":
        "*OC(=O)NCc1cccc(CNC(=O)O*)c1",
        # m-Xylylene diisocyanate; CAS 3634-83-1

    "pXDI_urethane":
        "*OC(=O)NCc1ccc(CNC(=O)O*)cc1",
        # p-Xylylene diisocyanate; CAS 4538-44-5

    # Cycloaliphatic diisocyanates (non-yellowing)
    "HDI_urethane":
        "*OC(=O)NCCCCCCNC(=O)O*",
        # Hexamethylene diisocyanate; 6-carbon aliphatic chain; CAS 822-06-0

    "IPDI_urethane":
        "CC1(C)CC(CC(C)(CNC(=O)O*)C1)NC(=O)O*",
        # Isophorone diisocyanate; 3,5,5-trimethylcyclohexane; asymmetric reactivity; CAS 4098-71-9

    "H12MDI_urethane":
        "*OC(=O)NC1CCC(CC2CCC(NC(=O)O*)CC2)CC1",
        # 4,4'-Methylenebis(cyclohexyl isocyanate); hydrogenated MDI; CAS 5124-30-1

    # --- Polymeric and modified MDI variants ---
    #
    # pMDI (polymeric MDI / crude MDI):
    #   Commercial pMDI is ~50% 4,4'-MDI monomer + ~35% 3-ring oligomer + ~15%
    #   higher oligomers; average NCO functionality 2.6-2.8. The entry below is
    #   the 4,4'-MDI repeat unit (the dominant structural motif in all fractions).
    #   For accurate formulation encoding, use the BLEND_LIBRARY entry "pMDI"
    #   which mixes the 2-ring and 3-ring fractions. Set fn=2.7 in metadata.
    "pMDI_urethane":
        "*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1",
        # Polymeric MDI - same aromatic urethane SMILES as MDI (IR-identical);
        # encode branching via fn=2.7 metadata; CAS 9016-87-9

    # pMDI 3-ring oligomer hard segment (trifunctional):
    #   Central ring is 1,4-bridged by two CH2 groups with an ortho urethane arm.
    #   fn=3; represents the dominant higher-molecular-weight fraction of pMDI.
    "pMDI_3ring_urethane":
        "*OC(=O)Nc1ccc(Cc2c(NC(=O)O*)ccc(Cc3ccc(NC(=O)O*)cc3)c2)cc1",
        # 3-ring polymeric MDI hard segment; fn=3

    # Modified MDI (carbodiimide-modified / liquid MDI):
    #   Produced by partial carbodiimidization of 4,4'-MDI. Contains ~80% MDI
    #   monomer + ~20% carbodiimide/uretonimine species. NCO content 28-31%
    #   vs 33.6% for pure MDI; liquid at room temperature. In the reacted polymer
    #   the uretonimine ring opens to give urethane/allophanate linkages; the IR
    #   signature remains MDI-aromatic-urethane. Set fn=2.1 and NCO_pct in metadata.
    #   For blend-level encoding use BLEND_LIBRARY entry "MDI_modified".
    "MDI_modified_urethane":
        "*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1",
        # Carbodiimide-modified MDI; same core SMILES as MDI_urethane

    # =========================================================================
    # CHAIN-EXTENDER SEGMENTS  (reacted diol between two isocyanate groups)
    # =========================================================================
    #
    # Convention: * = attachment to isocyanate core nitrogen
    # Pattern:    *NC(=O)O-[diol core]-OC(=O)N*
    #
    # Complements isocyanate repeat units: O on the isocyanate segment bonds
    # to N on the chain-extender segment at the urethane linkage.

    "EG_urethane":   "*NC(=O)OCCOC(=O)N*",         # Ethylene glycol CE segment
    "BDO_urethane":  "*NC(=O)OCCCCOC(=O)N*",        # 1,4-Butanediol CE segment (dominant TPU CE with MDI)
    "HDO_urethane":  "*NC(=O)OCCCCCCOC(=O)N*",      # 1,6-Hexanediol CE segment
    "HQEE_urethane": "*NC(=O)OCCOc1ccc(OCCOC(=O)N*)cc1",  # HQEE CE segment; aromatic diol

    # =========================================================================
    # CHAIN EXTENDERS  (small molecules - no * attachment points)
    # =========================================================================

    # Diol chain extenders
    "EG":   "OCCO",                   # Ethylene glycol;               CAS 107-21-1
    "PDO":  "OCCCO",                  # 1,3-Propanediol;               CAS 504-63-2
    "BDO":  "OCCCCO",                 # 1,4-Butanediol;                CAS 110-63-4
    "HDO":  "OCCCCCCO",               # 1,6-Hexanediol;                CAS 629-11-8
    "DEG":  "OCCOCCO",                # Diethylene glycol;             CAS 111-46-6
    "NPG":  "OCC(C)(C)CO",            # Neopentyl glycol;              CAS 126-30-7
    "HQEE": "OCCOc1ccc(OCCO)cc1",     # Hydroquinone bis(2-hydroxyethyl) ether; CAS 104-38-1
    "DMPA": "CC(CO)(CO)C(=O)O",       # Dimethylolpropionic acid; waterborne PU; CAS 4767-03-7

    # Amine chain extenders (form urea linkages, not urethane)
    "MOCA":  "Nc1ccc(Cc2ccc(N)c(Cl)c2)cc1Cl",  # 4,4'-Methylenebis(2-chloroaniline); CAS 101-14-4
    "DETDA": "CCc1cc(C)c(N)c(CC)c1N",            # 3,5-Diethyl-2,4-toluenediamine; CAS 68479-98-1
    "IPDA":  "CC1(C)CC(N)CC(C)(CN)C1",           # Isophoronediamine; CAS 2855-13-2
    "EDA":   "NCCN",                              # Ethylenediamine; CAS 107-15-3

    # =========================================================================
    # CROSSLINKERS
    # =========================================================================

    "glycerol": "OCC(O)CO",          # Glycerol;          f=3; CAS 56-81-5
    "TMP":      "CCC(CO)(CO)CO",     # Trimethylolpropane; f=3; CAS 77-99-6
    "DEA":      "OCCNCCO",           # Diethanolamine;    f=3 (1 NH + 2 OH); CAS 111-42-2
    "TEA":      "OCCN(CCO)CCO",      # Triethanolamine;   f=3 + autocatalytic tertiary N; CAS 102-71-6
    "sorbitol": "OCC(O)C(O)C(O)C(O)CO",  # Sorbitol;     f=6; CAS 50-70-4
}


def get_repeat_unit(name: str) -> Optional[str]:
    """Return the SMILES for a named entry in the polymer library, or None."""
    return POLYMER_REPEAT_UNITS.get(name)


def list_polymers() -> List[str]:
    """Return all names currently in the polymer library."""
    return list(POLYMER_REPEAT_UNITS.keys())


def resolve_smiles(name_or_smiles: str) -> str:
    """
    Resolve a library name or raw SMILES to a SMILES string.
    Checks the polymer library first; returns the input unchanged if not found.
    Does NOT expand blend names - use the builder for that.
    """
    resolved = POLYMER_REPEAT_UNITS.get(name_or_smiles)
    return resolved if resolved is not None else name_or_smiles


# =============================================================================
# BLEND LIBRARY
# =============================================================================
#
# Named multi-component blends. Each entry maps a blend name to a list of
# (component_name, fraction) tuples where fraction is 0-1 and fractions sum to 1.
#
# When a blend name is passed to PUSMILESBuilder.add() or add_weight(), it is
# automatically expanded into its constituent components, each weighted
# proportionally. Metadata (role, Mn, etc.) is applied to all sub-components.
#
# TDI grades:
#   TDI-80 (CAS 26471-62-5): 80 wt% 2,4-TDI + 20 wt% 2,6-TDI
#     Most common grade globally; used in flexible slabstock and moulded foam.
#   TDI-65 (CAS 26471-62-5): 65 wt% 2,4-TDI + 35 wt% 2,6-TDI
#     Higher 2,6-content; sometimes used in flexible foam for altered reactivity.
#
# pMDI blend:
#   Approximates commercial polymeric MDI: ~50% 4,4'-MDI monomer (fn=2) +
#   ~50% higher oligomers represented by the 3-ring species (fn=3), giving a
#   blend-average fn ~2.5. Adjust fractions to match a specific grade
#   (e.g. Lupranate M20, Mondur MR, Rubinate M).
#
# MDI_modified blend:
#   Approximates carbodiimide-modified liquid MDI: ~80% MDI monomer +
#   ~20% branched fraction (represented by pMDI_3ring_urethane). Add
#   NCO_pct=29 metadata manually to reflect lower isocyanate content.

BLEND_LIBRARY: Dict[str, List[Tuple[str, float]]] = {
    "TDI_80":       [("TDI_24_urethane",     0.80), ("TDI_26_urethane",     0.20)],
    "TDI_65":       [("TDI_24_urethane",     0.65), ("TDI_26_urethane",     0.35)],
    "pMDI":         [("pMDI_urethane",       0.50), ("pMDI_3ring_urethane", 0.50)],
    "MDI_modified": [("MDI_modified_urethane", 0.80), ("pMDI_3ring_urethane", 0.20)],
}


def get_blend(name: str) -> Optional[List[Tuple[str, float]]]:
    """Return the component list for a named blend, or None if not found."""
    return BLEND_LIBRARY.get(name)


def list_blends() -> List[str]:
    """Return all names currently in the blend library."""
    return list(BLEND_LIBRARY.keys())


def is_blend_name(name: str) -> bool:
    """Return True if the name refers to a blend in BLEND_LIBRARY."""
    return name in BLEND_LIBRARY


# =============================================================================
# ADDITIVE ROLE VOCABULARY
# =============================================================================
#
# Controlled set of valid role strings for additives. Using anything outside
# this set won't break anything, but consistent role labels are essential for
# ML models to align features across datasets.
#
# Polymer/formulation roles (already in use):
#   "soft"          - soft segment (polyol)
#   "hard"          - hard segment (reacted isocyanate)
#   "chain_ext"     - chain extender
#   "ionic_diol"    - ionic diol (e.g. DMPA in waterborne PU)
#
# Additive roles:
#   "filler"          - inorganic particle fillers (CaCO3, talc, BaSO4, silica)
#   "flame_retardant" - FR additives (ATH, MEL, phosphate esters, halogenated)
#   "plasticizer"     - plasticizers (phthalates, adipates, phosphates)
#   "pigment"         - colorants and opacifiers (TiO2, carbon black)
#   "blowing_agent"   - physical or chemical blowing agents (water, pentane, CO2)
#   "catalyst"        - urethane/urea reaction catalysts (DABCO, DBTDL)
#   "surfactant"      - foam stabilizers and emulsifiers (silicone surfactants)
#   "antioxidant"     - thermal/oxidative stabilizers (hindered phenols, phosphites)
#   "UV_stabilizer"   - light stabilizers (HALS, UV absorbers)
#   "adhesion_promoter" - adhesion promoters (silanes, titanates)
#   "additive"        - catch-all for unclassified additives

ADDITIVE_ROLES: frozenset = frozenset({
    # Formulation roles
    "soft", "hard", "chain_ext", "ionic_diol",
    # Additive roles
    "filler", "flame_retardant", "plasticizer", "pigment",
    "blowing_agent", "catalyst", "surfactant", "antioxidant",
    "UV_stabilizer", "adhesion_promoter", "additive",
})

# =============================================================================
# ADDITIVE LIBRARY
# =============================================================================
#
# Common PU additives. Each entry maps a name to a dict with keys:
#   smiles      - SMILES string (or "[*]" for non-encodable inorganics)
#   role        - controlled role string from ADDITIVE_ROLES
#   name        - human-readable name
#   cas         - CAS registry number (str)
#   encodable   - False if the material has no discrete molecular SMILES
#   formula     - molecular formula (for non-encodable entries)
#
# Non-encodable materials (minerals, amorphous carbons, polymeric surfactants)
# use "[*]" as a placeholder SMILES. The validator skips attachment-point checks
# when encodable=False. Store composition in the formula= metadata field.
# Use add_by_name() or add_by_cas() for anything not in this library.

ADDITIVE_LIBRARY: Dict[str, Dict[str, Any]] = {

    # -------------------------------------------------------------------------
    # INORGANIC FILLERS
    # -------------------------------------------------------------------------
    "CaCO3": {
        "smiles": "[Ca+2].[O-]C([O-])=O",
        "role": "filler", "name": "calcium_carbonate",
        "cas": "471-34-1",
    },
    "BaSO4": {
        "smiles": "[Ba+2].[O-]S([O-])(=O)=O",
        "role": "filler", "name": "barium_sulfate",
        "cas": "7727-43-7",
    },
    "silica": {
        "smiles": "O=[Si]=O",
        "role": "filler", "name": "silicon_dioxide",
        "cas": "7631-86-9",
    },
    "talc": {
        "smiles": "[*]", "encodable": False,
        "formula": "Mg3Si4O10(OH)2",
        "role": "filler", "name": "talc",
        "cas": "14807-96-6",
    },
    "kaolin": {
        "smiles": "[*]", "encodable": False,
        "formula": "Al2Si2O5(OH)4",
        "role": "filler", "name": "kaolin",
        "cas": "1318-74-7",
    },
    "glass_fiber": {
        "smiles": "[*]", "encodable": False,
        "formula": "SiO2_Al2O3_mixed",
        "role": "filler", "name": "glass_fiber",
    },
    "carbon_black": {
        "smiles": "[*]", "encodable": False,
        "formula": "C",
        "role": "pigment", "name": "carbon_black",
        "cas": "1333-86-4",
    },

    # -------------------------------------------------------------------------
    # FLAME RETARDANTS
    # -------------------------------------------------------------------------
    # Inorganic
    "ATH": {
        "smiles": "[Al+3].[OH-].[OH-].[OH-]",
        "role": "flame_retardant", "name": "aluminum_trihydroxide",
        "cas": "21645-51-2",
    },
    "MEL": {
        "smiles": "Nc1nc(N)nc(N)n1",
        "role": "flame_retardant", "name": "melamine",
        "cas": "108-78-1",
    },
    # Phosphate ester FRs (halogenated)
    "TCPP": {
        "smiles": "CC(CCl)OP(=O)(OC(C)CCl)OC(C)CCl",
        "role": "flame_retardant", "name": "tris_1_chloro_2_propyl_phosphate",
        "cas": "13674-84-5",
    },
    "TCEP": {
        "smiles": "ClCCOP(=O)(OCCCl)OCCCl",
        "role": "flame_retardant", "name": "tris_2_chloroethyl_phosphate",
        "cas": "115-96-8",
    },
    # Phosphate ester FRs (non-halogenated)
    "TEP": {
        "smiles": "CCOP(=O)(OCC)OCC",
        "role": "flame_retardant", "name": "triethyl_phosphate",
        "cas": "78-40-0",
    },
    "DMMP": {
        "smiles": "COP(=O)(OC)C",
        "role": "flame_retardant", "name": "dimethyl_methylphosphonate",
        "cas": "756-79-6",
    },
    "DOPO": {
        "smiles": "O=P1Oc2ccccc2-c2ccccc21",
        "role": "flame_retardant", "name": "9_10_dihydro_9_oxa_10_phosphaphenanthrene_10_oxide",
        "cas": "35948-25-5",
    },

    # -------------------------------------------------------------------------
    # PIGMENTS / OPACIFIERS
    # -------------------------------------------------------------------------
    "TiO2": {
        "smiles": "[*]", "encodable": False,
        "formula": "TiO2",
        "role": "pigment", "name": "titanium_dioxide",
        "cas": "13463-67-7",
    },

    # -------------------------------------------------------------------------
    # PLASTICIZERS
    # -------------------------------------------------------------------------
    "DOA": {
        "smiles": "CCCCC(CC)COC(=O)CCCCC(=O)OCC(CC)CCCC",
        "role": "plasticizer", "name": "dioctyl_adipate",
        "cas": "103-23-1",
    },
    "DEP": {
        "smiles": "CCOC(=O)c1ccccc1C(=O)OCC",
        "role": "plasticizer", "name": "diethyl_phthalate",
        "cas": "84-66-2",
    },
    "DBP": {
        "smiles": "CCCCOC(=O)c1ccccc1C(=O)OCCCC",
        "role": "plasticizer", "name": "dibutyl_phthalate",
        "cas": "84-74-2",
    },

    # -------------------------------------------------------------------------
    # BLOWING AGENTS
    # -------------------------------------------------------------------------
    "water": {
        "smiles": "O",
        "role": "blowing_agent", "name": "water",
        "cas": "7732-18-5",
    },
    "n_pentane": {
        "smiles": "CCCCC",
        "role": "blowing_agent", "name": "n_pentane",
        "cas": "109-66-0",
    },
    "cyclopentane": {
        "smiles": "C1CCCC1",
        "role": "blowing_agent", "name": "cyclopentane",
        "cas": "287-92-3",
    },
    "isopentane": {
        "smiles": "CCC(C)C",
        "role": "blowing_agent", "name": "isopentane",
        "cas": "78-78-4",
    },
    "HCBA": {
        "smiles": "FC(F)(Cl)C(F)(F)Cl",
        "role": "blowing_agent", "name": "HCFC_141b",
        "cas": "1717-00-6",
    },

    # -------------------------------------------------------------------------
    # CATALYSTS
    # -------------------------------------------------------------------------
    "DABCO": {
        "smiles": "C1CN2CCN1CC2",
        "role": "catalyst", "name": "triethylenediamine",
        "cas": "280-57-9",
    },
    "DBTDL": {
        "smiles": "[*]", "encodable": False,
        "formula": "C32H64O4Sn",
        "role": "catalyst", "name": "dibutyltin_dilaurate",
        "cas": "77-58-7",
    },
    "TEA_cat": {
        "smiles": "OCCN(CCO)CCO",
        "role": "catalyst", "name": "triethanolamine_catalyst",
        "cas": "102-71-6",
    },

    # -------------------------------------------------------------------------
    # SURFACTANTS / FOAM STABILIZERS
    # -------------------------------------------------------------------------
    "silicone_surfactant": {
        "smiles": "[*]", "encodable": False,
        "formula": "polysiloxane_polyether_copolymer",
        "role": "surfactant", "name": "silicone_surfactant",
    },

    # -------------------------------------------------------------------------
    # ANTIOXIDANTS
    # -------------------------------------------------------------------------
    "BHT": {
        "smiles": "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1",
        "role": "antioxidant", "name": "butylated_hydroxytoluene",
        "cas": "128-37-0",
    },

    # -------------------------------------------------------------------------
    # UV STABILIZERS
    # -------------------------------------------------------------------------
    "benzophenone": {
        "smiles": "O=C(c1ccccc1)c1ccccc1",
        "role": "UV_stabilizer", "name": "benzophenone",
        "cas": "119-61-9",
    },
}


def get_additive(name: str) -> Optional[Dict[str, Any]]:
    """Return the additive entry dict for a named additive, or None."""
    return ADDITIVE_LIBRARY.get(name)


def list_additives(role: Optional[str] = None) -> List[str]:
    """
    Return additive library names, optionally filtered by role.

        list_additives()                    # all additives
        list_additives("flame_retardant")   # only FRs
    """
    if role is None:
        return list(ADDITIVE_LIBRARY.keys())
    return [k for k, v in ADDITIVE_LIBRARY.items() if v.get("role") == role]


def is_additive_name(name: str) -> bool:
    """Return True if the name is in ADDITIVE_LIBRARY."""
    return name in ADDITIVE_LIBRARY


# =============================================================================
# PUBCHEM LOOKUP
# =============================================================================

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
_PUBCHEM_PROPS = "IsomericSMILES,IUPACName,MolecularFormula,MolecularWeight"
_PUBCHEM_TIMEOUT = 8  # seconds


def _pubchem_fetch(identifier: str, id_type: str = "name") -> Dict[str, Any]:
    """
    Internal: fetch compound properties from PubChem by name or CAS number.

    Args:
        identifier: compound name or CAS number string.
        id_type:    "name" (works for both common names and CAS numbers).

    Returns:
        Dict with keys: smiles, iupac_name, formula, mw (float), cid (int).

    Raises:
        ValueError:  compound not found, ambiguous, or request failed.
        RuntimeError: network error.
    """
    encoded = urllib.parse.quote(identifier.strip())
    url = f"{_PUBCHEM_BASE}/{id_type}/{encoded}/property/{_PUBCHEM_PROPS}/JSON"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PUSMILES/1.0"})
        with urllib.request.urlopen(req, timeout=_PUBCHEM_TIMEOUT) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(
                f"PubChem: '{identifier}' not found. "
                "Try a different name, the IUPAC name, or a CAS number."
            ) from e
        raise RuntimeError(f"PubChem HTTP error {e.code} for '{identifier}'.") from e
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"PubChem network error for '{identifier}': {e.reason}. "
            "Check your internet connection."
        ) from e

    props = data.get("PropertyTable", {}).get("Properties", [])
    if not props:
        raise ValueError(f"PubChem returned no properties for '{identifier}'.")

    p = props[0]
    return {
        "smiles":     p.get("IsomericSMILES", ""),
        "iupac_name": p.get("IUPACName", ""),
        "formula":    p.get("MolecularFormula", ""),
        "mw":         float(p.get("MolecularWeight", 0)),
        "cid":        int(p.get("CID", 0)),
    }


def pubchem_smiles_from_name(name: str) -> str:
    """
    Look up an isomeric SMILES from PubChem by compound name.
    Accepts common names, synonyms, and IUPAC names.

    Requires an internet connection. Returns the SMILES string.

        pubchem_smiles_from_name("calcium carbonate")  -> "[Ca+2].[O-]C([O-])=O"
        pubchem_smiles_from_name("melamine")           -> "Nc1nc(N)nc(N)n1"

    Raises ValueError if the name is not found.
    """
    return _pubchem_fetch(name)["smiles"]


def pubchem_smiles_from_cas(cas: str) -> str:
    """
    Look up an isomeric SMILES from PubChem by CAS registry number.

    Requires an internet connection. Returns the SMILES string.

        pubchem_smiles_from_cas("471-34-1")  -> "[Ca+2].[O-]C([O-])=O"

    Raises ValueError if the CAS number is not found.
    """
    return _pubchem_fetch(cas)["smiles"]


def pubchem_info(identifier: str) -> Dict[str, Any]:
    """
    Return full compound info from PubChem (smiles, iupac_name, formula, mw, cid).
    Accepts names, synonyms, or CAS numbers.
    """
    return _pubchem_fetch(identifier)
# =============================================================================

def _validate_smiles_string(smiles: str) -> Tuple[bool, str]:
    """
    Check whether a SMILES string is chemically valid.
    Returns (is_valid, reason) where reason is empty on success.
    """
    if not smiles or not smiles.strip():
        return False, "SMILES string is empty."
    if not HAS_RDKIT:
        return True, ""   # Cannot validate without RDKit; pass through.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, (
            f"RDKit could not parse '{smiles}'. "
            "Check for typos, invalid valence, or mismatched ring closures."
        )
    return True, ""


def _run_validation(components: List["QComponent"]) -> "ValidationResult":
    """Run all checks and return a ValidationResult."""
    issues: List[ValidationIssue] = []

    if not components:
        issues.append(ValidationIssue(
            severity="error",
            message="Formulation has no components.",
        ))
        return ValidationResult(is_valid=False, issues=issues)

    seen_smiles: Dict[str, int] = {}

    for i, c in enumerate(components):

        # 1. SMILES chemical validity (error)
        # Skip for non-encodable placeholder entries (minerals, amorphous materials)
        if c.metadata.get("encodable") is not False:
            ok, reason = _validate_smiles_string(c.smiles)
            if not ok:
                issues.append(ValidationIssue(
                    severity="error",
                    message=reason,
                    component_index=i,
                    component_smiles=c.smiles,
                ))

        # 2. Zero weight percent (warning)
        if c.weight_pct == 0:
            issues.append(ValidationIssue(
                severity="warning",
                message="weight_pct is 0. This component contributes nothing to the formulation.",
                component_index=i,
                component_smiles=c.smiles,
            ))

        # 3. Attachment point count for polymer SMILES
        # Skip this check for non-encodable placeholder [*] (minerals, amorphous materials)
        is_placeholder = c.smiles == "[*]" or c.metadata.get("encodable") is False
        if "*" in c.smiles and not is_placeholder:
            n_stars = c.smiles.count("*")
            if n_stars == 1:
                issues.append(ValidationIssue(
                    severity="warning",
                    message=(
                        "Polymer SMILES has only 1 attachment point (*). "
                        "Linear repeat units require exactly 2. "
                        "If this is a deliberate chain-end or dendrimer arm, add fn=1 to metadata."
                    ),
                    component_index=i,
                    component_smiles=c.smiles,
                ))
            elif n_stars > 2:
                issues.append(ValidationIssue(
                    severity="warning",
                    message=(
                        f"Polymer SMILES has {n_stars} attachment points (*), "
                        "indicating a branch or crosslink point. "
                        f"Add fn={n_stars} to metadata if intentional."
                    ),
                    component_index=i,
                    component_smiles=c.smiles,
                ))

        # 4. Duplicate SMILES (warning)
        if c.smiles in seen_smiles:
            issues.append(ValidationIssue(
                severity="warning",
                message=(
                    f"Duplicate SMILES: same structure as component[{seen_smiles[c.smiles]}]. "
                    "If intentional (e.g. same backbone, different Mn), "
                    "add a 'name' or 'grade' metadata field to distinguish them."
                ),
                component_index=i,
                component_smiles=c.smiles,
            ))
        else:
            seen_smiles[c.smiles] = i

        # 5. Metadata type checks for numeric fields
        for key in ("MW", "Mn", "Mw", "PDI", "fn"):
            val = c.metadata.get(key)
            if val is not None and not isinstance(val, (int, float)):
                issues.append(ValidationIssue(
                    severity="warning",
                    message=(
                        f"Metadata '{key}' should be numeric, "
                        f"got {type(val).__name__} ({val!r}). "
                        "This may cause errors during feature extraction."
                    ),
                    component_index=i,
                    component_smiles=c.smiles,
                ))

        # 6. Non-positive molecular weight
        for key in ("MW", "Mn", "Mw"):
            val = c.metadata.get(key)
            if isinstance(val, (int, float)) and val <= 0:
                issues.append(ValidationIssue(
                    severity="warning",
                    message=f"Metadata '{key}' = {val}. Molecular weights must be > 0.",
                    component_index=i,
                    component_smiles=c.smiles,
                ))

    # 7. Weight sum (formulation-level warning)
    total = sum(c.weight_pct for c in components)
    if not (99.0 <= total <= 101.0):
        issues.append(ValidationIssue(
            severity="warning",
            message=(
                f"Component weights sum to {total:.4f} wt%, expected 100 +/- 1. "
                "Call .normalize() to rescale, or check for missing components."
            ),
        ))

    is_valid = all(i.severity != "error" for i in issues)
    return ValidationResult(is_valid=is_valid, issues=issues)


# =============================================================================
# CORE CLASSES
# =============================================================================

@dataclass
class QComponent:
    """A single component in a PUSMILES string."""
    smiles: str
    weight_pct: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.weight_pct < 0 or self.weight_pct > 100:
            raise ValueError(f"weight_pct must be 0-100, got {self.weight_pct}.")

    @property
    def mw(self) -> Optional[float]:
        """Molecular weight from metadata (prefers MW, falls back to Mn)."""
        return self.metadata.get("MW") or self.metadata.get("Mn")

    @property
    def name(self) -> Optional[str]:
        return self.metadata.get("name")

    @property
    def role(self) -> Optional[str]:
        return self.metadata.get("role")

    @property
    def is_polymer(self) -> bool:
        return "*" in self.smiles

    @property
    def attachment_count(self) -> int:
        return self.smiles.count("*")

    def to_string(self) -> str:
        wt_str = (str(int(self.weight_pct)) if self.weight_pct == int(self.weight_pct)
                  else str(self.weight_pct))
        s = f"{self.smiles}@{wt_str}"
        if self.metadata:
            pairs = []
            for k, v in self.metadata.items():
                if isinstance(v, str) and any(ch in v for ch in (" ", ",", "|")):
                    pairs.append(f'{k}="{v}"')
                elif isinstance(v, float) and v == int(v):
                    pairs.append(f"{k}={int(v)}")
                else:
                    pairs.append(f"{k}={v}")
            s += "{" + ",".join(pairs) + "}"
        return s

    def __repr__(self):
        return self.to_string()

    def canonicalize(self) -> "QComponent":
        """Return a copy with a canonicalized SMILES (requires RDKit)."""
        if not HAS_RDKIT:
            return self
        try:
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                return self
            return QComponent(Chem.MolToSmiles(mol, canonical=True),
                              self.weight_pct, self.metadata.copy())
        except Exception:
            return self


@dataclass
class PUSMILES:
    """A complete PUSMILES formulation (one or more components)."""
    components: List[QComponent]

    # =========================================================================
    # Constructors
    # =========================================================================

    @classmethod
    def parse(cls, pusmiles_string: str) -> "PUSMILES":
        """
        Parse a PUSMILES string.

            PUSMILES.parse("CCO@100")
            PUSMILES.parse("PPG@55{role=soft,Mn=2000}|MDI_urethane@45{role=hard}")
        """
        parts: List[str] = []
        current = ""
        brace_depth = 0

        for ch in pusmiles_string:
            if ch == "{":
                brace_depth += 1
                current += ch
            elif ch == "}":
                brace_depth -= 1
                current += ch
            elif ch == "|" and brace_depth == 0:
                parts.append(current.strip())
                current = ""
            else:
                current += ch
        if current.strip():
            parts.append(current.strip())

        components = []
        for part in parts:
            if not part:
                continue
            metadata: Dict[str, Any] = {}
            brace_match = re.search(r"\{([^}]+)\}$", part)
            if brace_match:
                meta_str = brace_match.group(1)
                part = part[: brace_match.start()].strip()
                for m in re.finditer(r'(\w+)=("([^"]*)"|([^,}]+))', meta_str):
                    key = m.group(1)
                    value: Any = m.group(3) if m.group(3) is not None else m.group(4)
                    value = value.strip()
                    try:
                        value = float(value) if "." in str(value) else int(value)
                    except (ValueError, TypeError):
                        pass
                    metadata[key] = value

            if "@" in part:
                smiles, wt_str = part.rsplit("@", 1)
                weight_pct = float(wt_str.strip())
            else:
                smiles, weight_pct = part, 100.0

            components.append(QComponent(smiles=smiles.strip(),
                                         weight_pct=weight_pct,
                                         metadata=metadata))
        return cls(components=components)

    @classmethod
    def from_smiles(cls, smiles: str, weight_pct: float = 100.0, **metadata) -> "PUSMILES":
        """
        Single-component PUSMILES. Accepts a library name or raw SMILES.

            PUSMILES.from_smiles("CCO")
            PUSMILES.from_smiles("PPG", weight_pct=100, Mn=2000)
        """
        smiles = resolve_smiles(smiles)
        return cls([QComponent(smiles, weight_pct, metadata or {})])

    @classmethod
    def from_mixture(cls, *components) -> "PUSMILES":
        """
        Create PUSMILES from explicit weight-percent components.

        Each element:  (name_or_smiles, weight_pct)
                   or  (name_or_smiles, weight_pct, metadata_dict)

            PUSMILES.from_mixture(
                ("PPG",          55, {"role": "soft", "Mn": 2000}),
                ("MDI_urethane", 45, {"role": "hard"}),
            )
        """
        parsed = []
        for c in components:
            smiles = resolve_smiles(c[0])
            weight_pct = c[1]
            meta = c[2] if len(c) >= 3 else {}
            parsed.append(QComponent(smiles, weight_pct, meta or {}))
        return cls(components=parsed)

    @classmethod
    def from_weights(cls, *components) -> "PUSMILES":
        """
        Create PUSMILES from absolute weights; auto-normalized to 100 wt%.

        Each element:  (name_or_smiles, weight)
                   or  (name_or_smiles, weight, metadata_dict)

            PUSMILES.from_weights(
                ("PPG",          110, {"role": "soft", "Mn": 2000}),
                ("MDI_urethane",  90, {"role": "hard"}),
            )
        """
        raw = []
        for c in components:
            smiles = resolve_smiles(c[0])
            weight = float(c[1])
            meta = c[2] if len(c) >= 3 else {}
            raw.append((smiles, weight, meta or {}))
        total = sum(w for _, w, _ in raw)
        if total == 0:
            raise ValueError("Total weight is zero.")
        return cls([
            QComponent(s, round(w / total * 100, 4), m)
            for s, w, m in raw
        ])

    @classmethod
    def from_name(cls, name: str, weight_pct: float = 100.0, **metadata) -> "PUSMILES":
        """
        Create a single-component PUSMILES by looking up a compound name on PubChem.

        Checks the additive library first; falls back to a live PubChem query.
        Requires an internet connection for compounds not in the additive library.

            PUSMILES.from_name("calcium carbonate", weight_pct=15, role="filler")
            PUSMILES.from_name("melamine", weight_pct=10, role="flame_retardant")
            PUSMILES.from_name("triethyl phosphate", weight_pct=8, role="flame_retardant")

        Raises:
            ValueError: if PubChem cannot find the compound.
            RuntimeError: on network error.
        """
        # Try additive library first (no network needed)
        entry = ADDITIVE_LIBRARY.get(name)
        if entry:
            meta = {k: v for k, v in entry.items()
                    if k not in ("smiles", "encodable")}
            meta.update(metadata)
            if entry.get("encodable") is False:
                meta["encodable"] = False
            return cls([QComponent(entry["smiles"], weight_pct, meta)])

        # Fall back to PubChem
        info = _pubchem_fetch(name)
        meta = {
            "name": name,
            "iupac": info["iupac_name"],
            "formula": info["formula"],
            "MW": info["mw"],
            "cid": info["cid"],
        }
        meta.update(metadata)
        return cls([QComponent(info["smiles"], weight_pct, meta)])

    @classmethod
    def from_cas(cls, cas: str, weight_pct: float = 100.0, **metadata) -> "PUSMILES":
        """
        Create a single-component PUSMILES by CAS registry number via PubChem.

        Checks the additive library first; falls back to a live PubChem query.

            PUSMILES.from_cas("471-34-1", weight_pct=15, role="filler")

        Raises:
            ValueError: if the CAS number is not found.
            RuntimeError: on network error.
        """
        # Try additive library by CAS match first
        for lib_name, entry in ADDITIVE_LIBRARY.items():
            if entry.get("cas") == cas.strip():
                meta = {k: v for k, v in entry.items()
                        if k not in ("smiles", "encodable")}
                meta.update(metadata)
                if entry.get("encodable") is False:
                    meta["encodable"] = False
                return cls([QComponent(entry["smiles"], weight_pct, meta)])

        # Fall back to PubChem
        info = _pubchem_fetch(cas)
        meta = {
            "cas": cas,
            "iupac": info["iupac_name"],
            "formula": info["formula"],
            "MW": info["mw"],
            "cid": info["cid"],
        }
        meta.update(metadata)
        return cls([QComponent(info["smiles"], weight_pct, meta)])

    # =========================================================================
    # String output
    # =========================================================================

    def to_string(self) -> str:
        return "|".join(c.to_string() for c in self.components)

    def __repr__(self):
        return self.to_string()

    def __str__(self):
        return self.to_string()

    # =========================================================================
    # Validation
    # =========================================================================

    def validate(self) -> "ValidationResult":
        """
        Run full validation and return a ValidationResult.

        Checks performed:
            ERROR:   invalid SMILES (RDKit parse failure, empty string)
            ERROR:   empty formulation
            WARNING: weight_pct = 0
            WARNING: polymer SMILES with != 2 attachment points
            WARNING: duplicate SMILES within the formulation
            WARNING: non-numeric MW/Mn/Mw/PDI/fn metadata
            WARNING: non-positive molecular weight metadata
            WARNING: weights sum outside 99-101 wt%

        Example:
            result = q.validate()
            result.print_report()
            if not result:
                raise RuntimeError("Invalid PUSMILES")
        """
        return _run_validation(self.components)

    # =========================================================================
    # Analysis
    # =========================================================================

    def total_weight(self) -> float:
        return sum(c.weight_pct for c in self.components)

    def is_valid(self) -> bool:
        """True if no validation errors and weights sum to 100 +/- 1 wt%."""
        return _run_validation(self.components).is_valid and 99.0 <= self.total_weight() <= 101.0

    def normalize(self) -> "PUSMILES":
        """Return a copy with all weights rescaled to sum to exactly 100 wt%."""
        total = self.total_weight()
        if total == 0:
            return self
        factor = 100.0 / total
        return PUSMILES([
            QComponent(c.smiles, round(c.weight_pct * factor, 4), c.metadata.copy())
            for c in self.components
        ])

    def by_role(self, role: str) -> List[QComponent]:
        """Return all components whose 'role' metadata matches the given string."""
        return [c for c in self.components if c.metadata.get("role") == role]

    def sum_by_role(self, role: str) -> float:
        return sum(c.weight_pct for c in self.by_role(role))

    def get_smiles_list(self) -> List[str]:
        return [c.smiles for c in self.components]

    def get_polymer_components(self) -> List[QComponent]:
        """Return components that are polymer repeat units (contain *)."""
        return [c for c in self.components if c.is_polymer]

    def canonicalize(self) -> "PUSMILES":
        """Return a copy with all SMILES canonicalized (requires RDKit)."""
        return PUSMILES([c.canonicalize() for c in self.components])

    # =========================================================================
    # ML feature extraction
    # =========================================================================

    def to_features(self) -> Dict[str, Any]:
        """
        Flat dict of ML-ready features, suitable for a pandas DataFrame row.

        Always present:   n_components, total_weight, n_polymers, n_additives
        Per-role:         {role}_pct for every role found in the formulation
                          Covers both polymer roles (soft, hard, chain_ext) and
                          additive roles (filler, flame_retardant, plasticizer, etc.)
        Per-component:    c{i}_wt, c{i}_mw (if set), c{i}_fn (if set)
        """
        additive_roles = ADDITIVE_ROLES - {"soft", "hard", "chain_ext", "ionic_diol"}
        n_additives = sum(
            1 for c in self.components
            if c.role in additive_roles or c.metadata.get("encodable") is False
        )

        features: Dict[str, Any] = {
            "n_components": len(self.components),
            "total_weight": self.total_weight(),
            "n_polymers": len(self.get_polymer_components()),
            "n_additives": n_additives,
        }
        for role in {c.role for c in self.components if c.role}:
            features[f"{role}_pct"] = self.sum_by_role(role)
        for i, c in enumerate(self.components):
            features[f"c{i}_wt"] = c.weight_pct
            if c.mw:
                features[f"c{i}_mw"] = c.mw
            if c.metadata.get("fn"):
                features[f"c{i}_fn"] = c.metadata["fn"]
        return features

    # =========================================================================
    # Serialization
    # =========================================================================

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pusmiles": self.to_string(),
            "components": [
                {"smiles": c.smiles, "weight_pct": c.weight_pct, "metadata": c.metadata}
                for c in self.components
            ],
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


# =============================================================================
# FLUENT BUILDER
# =============================================================================

class PUSMILESBuilder:
    """
    Fluent builder for PUSMILES formulations.

    Resolves polymer library names, blend names, and additive library names
    automatically. Validates on build(): errors raise PUSMILESValidationError,
    warnings are issued via Python's warnings module.

    Usage (weight percent):
        q = (PUSMILESBuilder()
             .add("PPG",          55, role="soft", Mn=2000)
             .add("MDI_urethane", 45, role="hard")
             .build())

    With additives:
        q = (PUSMILESBuilder()
             .add("PPG",          80, role="soft", Mn=3000)
             .add("MDI_urethane", 15, role="hard")
             .add_additive("CaCO3", 4)
             .add_additive("BHT",   1)
             .build())

    Usage (absolute weights, auto-normalized):
        q = (PUSMILESBuilder()
             .add_weight("PPG",          110, role="soft", Mn=2000)
             .add_weight("MDI_urethane",  90, role="hard")
             .build())

    Mixing add() and add_weight() raises ValueError immediately.
    """

    def __init__(self):
        self._components: List[Dict[str, Any]] = []
        self._use_weights: bool = False

    def add(self, name_or_smiles: str, weight_pct: float, **metadata) -> "PUSMILESBuilder":
        """
        Add a component by weight percent (0-100).
        All components should sum to 100.

        Args:
            name_or_smiles: Library name (e.g. "PPG") or raw SMILES.
            weight_pct:     Weight percent.
            **metadata:     role, Mn, Mw, MW, PDI, fn, name, cas, etc.
        """
        if self._use_weights:
            raise ValueError(
                "Cannot mix add() and add_weight(). "
                "Use add_weight() for all components or add() for all."
            )
        if is_blend_name(name_or_smiles):
            for comp_name, fraction in BLEND_LIBRARY[name_or_smiles]:
                self._components.append({
                    "smiles": resolve_smiles(comp_name),
                    "value": round(float(weight_pct) * fraction, 6),
                    "metadata": dict(metadata),
                })
        else:
            self._components.append({
                "smiles": resolve_smiles(name_or_smiles),
                "value": float(weight_pct),
                "metadata": metadata,
            })
        return self

    def add_weight(self, name_or_smiles: str, weight: float, **metadata) -> "PUSMILESBuilder":
        """
        Add a component by absolute weight (grams, parts, or any consistent unit).
        Weights are normalized to wt% on build().

        Args:
            name_or_smiles: Library name or raw SMILES.
            weight:         Absolute weight in any consistent unit.
            **metadata:     role, Mn, Mw, MW, PDI, fn, name, cas, etc.
        """
        if self._components and not self._use_weights:
            raise ValueError(
                "Cannot mix add() and add_weight(). "
                "Use add_weight() for all components or add() for all."
            )
        self._use_weights = True
        if is_blend_name(name_or_smiles):
            for comp_name, fraction in BLEND_LIBRARY[name_or_smiles]:
                self._components.append({
                    "smiles": resolve_smiles(comp_name),
                    "value": round(float(weight) * fraction, 6),
                    "metadata": dict(metadata),
                })
        else:
            self._components.append({
                "smiles": resolve_smiles(name_or_smiles),
                "value": float(weight),
                "metadata": metadata,
            })
        return self

    def build(self) -> PUSMILES:
        """
        Build, validate, and return the PUSMILES object.

        Errors (invalid SMILES, empty formulation) raise PUSMILESValidationError.
        Warnings (weight-sum drift, duplicates, metadata type issues) are emitted
        via warnings.warn() but do not block construction.

        Raises:
            ValueError: if no components were added or weights sum to zero.
            PUSMILESValidationError: if any validation errors are present.
        """
        if not self._components:
            raise ValueError("No components added. Call .add() or .add_weight() first.")

        if self._use_weights:
            total = sum(c["value"] for c in self._components)
            if total == 0:
                raise ValueError("Total weight is zero.")
            qcomponents = [
                QComponent(c["smiles"], round(c["value"] / total * 100, 4), c["metadata"])
                for c in self._components
            ]
        else:
            qcomponents = [
                QComponent(c["smiles"], c["value"], c["metadata"])
                for c in self._components
            ]

        result = _run_validation(qcomponents)

        if result.errors:
            raise PUSMILESValidationError(
                f"PUSMILES construction failed with {len(result.errors)} error(s). "
                "Call e.result.print_report() for details.",
                result,
            )

        for w in result.warnings:
            warnings.warn(str(w), stacklevel=2)

        return PUSMILES(components=qcomponents)

    def add_additive(self, name: str, weight_pct: float, **metadata) -> "PUSMILESBuilder":
        """
        Add an additive from the ADDITIVE_LIBRARY by name, with weight percent.

        The role, name, cas, and encodable flag are pulled from the library
        automatically and can be overridden by keyword arguments.

        For non-encodable materials (talc, glass fiber, TiO2, etc.) the
        placeholder SMILES "[*]" is stored and the validator skips SMILES
        checks for that component.

            builder.add_additive("CaCO3",   15, role="filler")
            builder.add_additive("TCPP",    8,  role="flame_retardant")
            builder.add_additive("BHT",     0.5)

        For additives not in the library, use add() with a raw SMILES or use
        add_additive_by_name() for a live PubChem lookup.

        Raises:
            KeyError: if the name is not in ADDITIVE_LIBRARY.
        """
        if name not in ADDITIVE_LIBRARY:
            raise KeyError(
                f"'{name}' is not in ADDITIVE_LIBRARY. "
                f"Use add() with a SMILES string, or add_additive_by_name() "
                f"for a PubChem lookup. Available: {list_additives()}"
            )
        if self._use_weights:
            raise ValueError(
                "Cannot mix add_additive() and add_weight(). "
                "Use add_additive_weight() instead."
            )
        entry = ADDITIVE_LIBRARY[name]
        meta: Dict[str, Any] = {k: v for k, v in entry.items()
                                 if k not in ("smiles",)}
        meta.update(metadata)
        self._components.append({
            "smiles": entry["smiles"],
            "value": float(weight_pct),
            "metadata": meta,
        })
        return self

    def add_additive_weight(self, name: str, weight: float, **metadata) -> "PUSMILESBuilder":
        """
        Add an additive from ADDITIVE_LIBRARY by absolute weight (auto-normalized).

        Mirrors add_additive() but for use with add_weight() workflows.

            builder.add_additive_weight("CaCO3",  30)   # 30 g
            builder.add_additive_weight("TCPP",   15)   # 15 g

        Raises:
            KeyError: if the name is not in ADDITIVE_LIBRARY.
        """
        if name not in ADDITIVE_LIBRARY:
            raise KeyError(
                f"'{name}' is not in ADDITIVE_LIBRARY. "
                f"Use add_weight() with a SMILES string, or look up via pubchem_info()."
            )
        if self._components and not self._use_weights:
            raise ValueError(
                "Cannot mix add_additive_weight() and add(). "
                "Use add_additive() instead."
            )
        self._use_weights = True
        entry = ADDITIVE_LIBRARY[name]
        meta: Dict[str, Any] = {k: v for k, v in entry.items()
                                 if k not in ("smiles",)}
        meta.update(metadata)
        self._components.append({
            "smiles": entry["smiles"],
            "value": float(weight),
            "metadata": meta,
        })
        return self

    def add_additive_by_name(self, compound_name: str, weight_pct: float,
                             **metadata) -> "PUSMILESBuilder":
        """
        Look up a compound on PubChem by name and add it as an additive.

        Checks ADDITIVE_LIBRARY first (no network needed); falls back to a live
        PubChem query. Requires an internet connection for unlisted compounds.

            builder.add_additive_by_name("triethyl phosphate", 8, role="flame_retardant")
            builder.add_additive_by_name("calcium carbonate",  15, role="filler")

        Raises:
            ValueError: if PubChem cannot find the compound.
            RuntimeError: on network error.
        """
        if self._use_weights:
            raise ValueError(
                "Cannot mix add_additive_by_name() and add_weight(). "
                "Use add_additive_weight() or a dedicated weight-based workflow."
            )
        # Try additive library first
        entry = ADDITIVE_LIBRARY.get(compound_name)
        if entry:
            meta: Dict[str, Any] = {k: v for k, v in entry.items()
                                     if k not in ("smiles",)}
            meta.update(metadata)
            self._components.append({
                "smiles": entry["smiles"],
                "value": float(weight_pct),
                "metadata": meta,
            })
            return self

        # Live PubChem lookup
        info = _pubchem_fetch(compound_name)
        meta = {
            "name": compound_name,
            "iupac": info["iupac_name"],
            "formula": info["formula"],
            "MW": info["mw"],
            "cid": info["cid"],
            "role": metadata.pop("role", "additive"),
        }
        meta.update(metadata)
        self._components.append({
            "smiles": info["smiles"],
            "value": float(weight_pct),
            "metadata": meta,
        })
        return self

    def reset(self) -> "PUSMILESBuilder":
        """Clear all components and return self for reuse."""
        self._components.clear()
        self._use_weights = False
        return self

    def __repr__(self):
        mode = "weights" if self._use_weights else "wt%"
        return f"PUSMILESBuilder({len(self._components)} components, mode={mode})"


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize a SMILES string using RDKit. Returns input unchanged if RDKit unavailable."""
    if not HAS_RDKIT:
        return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return smiles


def validate_smiles(smiles: str) -> bool:
    """Return True if the SMILES string is chemically valid (requires RDKit)."""
    ok, _ = _validate_smiles_string(smiles)
    return ok


def is_polymer_smiles(smiles: str) -> bool:
    """Return True if the SMILES contains * attachment points."""
    return "*" in smiles


# =============================================================================
# MAIN (examples and validation demo)
# =============================================================================

if __name__ == "__main__":
    print("PUSMILES Module - Examples")
    print("=" * 60)

    # 1. Pure compound
    print("\n1. Pure compound:")
    q1 = PUSMILES.from_smiles("CCO")
    print(f"   {q1}")

    # 2. Flexible PU foam (builder, wt%)
    print("\n2. Flexible PU foam (builder, wt%):")
    q2 = (PUSMILESBuilder()
          .add("PPG",          45, role="soft",      Mn=3000)
          .add("PEG",          10, role="soft",      Mn=1000)
          .add("MDI_urethane", 40, role="hard")
          .add("EG_urethane",   5, role="chain_ext")
          .build())
    print(f"   {q2}")
    print(f"   Soft: {q2.sum_by_role('soft'):.1f}%  Hard: {q2.sum_by_role('hard'):.1f}%")

    # 3. TPU from bench weights (PTMEG/MDI/BDO system)
    print("\n3. TPU from bench weights (add_weight):")
    q3 = (PUSMILESBuilder()
          .add_weight("PTMEG",         100, role="soft",      Mn=2000)
          .add_weight("MDI_urethane",   72, role="hard")
          .add_weight("BDO_urethane",   18, role="chain_ext")
          .build())
    print(f"   {q3}")

    # 4. Filled flexible foam with additives (add_additive)
    print("\n4. Filled flexible foam with additives:")
    import warnings as _w
    with _w.catch_warnings(record=True):
        _w.simplefilter("always")
        q4 = (PUSMILESBuilder()
              .add("PPG",          72, role="soft",  Mn=3000)
              .add("MDI_urethane", 18, role="hard")
              .add_additive("CaCO3",  5)           # auto-gets role=filler from library
              .add_additive("DABCO",  2)           # auto-gets role=catalyst
              .add_additive("water",  1)           # auto-gets role=blowing_agent
              .add_additive("BHT",    1)           # auto-gets role=antioxidant
              .add_additive("silicone_surfactant", 1)  # non-encodable, placeholder [*]
              .build())
    print(f"   {q4}")
    print(f"   Soft: {q4.sum_by_role('soft'):.1f}%  "
          f"Filler: {q4.sum_by_role('filler'):.1f}%  "
          f"Blowing: {q4.sum_by_role('blowing_agent'):.1f}%")

    # 5. Flame-retarded rigid foam
    print("\n5. Flame-retarded rigid foam (non-halogenated FRs):")
    with _w.catch_warnings(record=True):
        _w.simplefilter("always")
        q5 = (PUSMILESBuilder()
              .add("PPG",          35, role="soft",  Mn=400, fn=3)
              .add("pMDI",         45, role="hard",  fn=2.7)
              .add_additive("TEP",  10, role="flame_retardant")
              .add_additive("MEL",   5, role="flame_retardant")
              .add_additive("water", 3, role="blowing_agent")
              .add_additive("DABCO", 2)
              .build())
    print(f"   {q5}")
    print(f"   FR total: {q5.sum_by_role('flame_retardant'):.1f}%")

    # 6. Absolute-weight workflow with additives
    print("\n6. Bench-weight foam with CaCO3 filler:")
    with _w.catch_warnings(record=True):
        _w.simplefilter("always")
        q6 = (PUSMILESBuilder()
              .add_weight("PPG",          100, role="soft", Mn=3000)
              .add_weight("MDI_urethane",  45, role="hard")
              .add_additive_weight("CaCO3", 20)
              .add_additive_weight("DABCO",  1)
              .add_additive_weight("water",  2)
              .build())
    print(f"   {q6}")
    print(f"   Filler: {q6.sum_by_role('filler'):.1f}%  "
          f"Blowing: {q6.sum_by_role('blowing_agent'):.1f}%")

    # 7. TDI-80 blend with halogenated FR
    print("\n7. Flexible foam, TDI-80, halogenated FR:")
    with _w.catch_warnings(record=True):
        _w.simplefilter("always")
        q7 = (PUSMILESBuilder()
              .add("PPG",    60, role="soft", Mn=3000)
              .add("TDI_80", 30, role="hard")
              .add_additive("TCPP",  8, role="flame_retardant")
              .add_additive("water", 2, role="blowing_agent")
              .build())
    print(f"   {q7}")

    # 8. PubChem lookup demo (requires internet)
    print("\n8. PubChem lookup (from_name):")
    try:
        q8 = PUSMILES.from_name("calcium carbonate", weight_pct=100, role="filler")
        print(f"   {q8}")
    except (ValueError, RuntimeError) as e:
        print(f"   (offline) Would produce: [Ca+2].[O-]C([O-])=O@100{{role=filler}}")
        print(f"   Error: {e}")

    # 9. PubChem via builder (add_additive_by_name)
    print("\n9. PubChem via builder (add_additive_by_name):")
    try:
        q9 = (PUSMILESBuilder()
              .add("PPG",          85, role="soft", Mn=3000)
              .add("MDI_urethane", 10, role="hard")
              .add_additive_by_name("triethyl phosphate", 5, role="flame_retardant")
              .build())
        print(f"   {q9}")
    except (ValueError, RuntimeError) as e:
        print(f"   (offline) Would resolve triethyl phosphate SMILES from PubChem.")
        print(f"   Error: {e}")

    # 10. ML features showing additive roles
    print("\n10. ML features (q4 - filled foam):")
    for k, v in q4.to_features().items():
        print(f"    {k}: {v}")

    # 11. Validation - valid formulation
    print("\n11. Validation (valid formulation):")
    result = q2.validate()
    result.print_report()

    # 12. Validation - invalid SMILES
    print("\n12. Validation (invalid SMILES - raises PUSMILESValidationError):")
    try:
        q_bad = (PUSMILESBuilder()
                 .add("TOTALLY_BOGUS###", 50)
                 .add("PPG", 50)
                 .build())
    except PUSMILESValidationError as e:
        print("   Caught PUSMILESValidationError as expected:")
        e.result.print_report()

    # 13. Library overview
    print(f"\n13. Polymer library: {len(list_polymers())} entries")
    categories = {
        "Commodity":         ["PE", "PP", "PS", "PVC", "PMMA"],
        "Polyether polyols": ["PEG", "PPG", "PTMEG", "PBO",
                              "graft_polyol_PPG", "graft_polyol_PEG"],
        "Polyester polyols": ["PBA", "PEA", "PHA", "PCL", "PES", "PCDL"],
        "Hard segments":     ["MDI_urethane", "pMDI_urethane", "pMDI_3ring_urethane",
                              "MDI_modified_urethane",
                              "TDI_urethane", "TDI_24_urethane", "TDI_26_urethane",
                              "HDI_urethane", "IPDI_urethane", "H12MDI_urethane",
                              "NDI_urethane"],
        "CE segments":       ["EG_urethane", "BDO_urethane", "HDO_urethane"],
        "Diol CEs":          ["EG", "BDO", "PDO", "HDO", "NPG", "DMPA"],
        "Amine CEs":         ["MOCA", "DETDA", "IPDA"],
        "Crosslinkers":      ["glycerol", "TMP", "DEA", "TEA", "sorbitol"],
    }
    for cat, names in categories.items():
        present = [n for n in names if n in POLYMER_REPEAT_UNITS]
        print(f"   {cat}: {', '.join(present)}")

    print(f"\n   Blend library: {', '.join(list_blends())}")

    print(f"\n   Additive library: {len(list_additives())} entries")
    for role in sorted({v['role'] for v in ADDITIVE_LIBRARY.values()}):
        names = list_additives(role)
        print(f"   {role}: {', '.join(names)}")

    print("\n" + "=" * 60)
    print(f"RDKit available: {HAS_RDKIT}")
    if HAS_RDKIT:
        test = canonicalize_smiles("C(CO)O")
        print(f"Canonicalization test: C(CO)O -> {test}")

