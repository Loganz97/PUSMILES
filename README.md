# PUSMILES

**Polyurethane SMILES**: a quantitative, machine-readable notation for encoding polyurethane formulations as strings suitable for machine learning pipelines.

```
*OCC(C)*@55{role=soft,Mn=3000}|*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1@40{role=hard}|OCCCCO@5{role=chain_ext}
```

Each component is its **SMILES repeat unit** at its **weight percent**, with optional metadata. The full formulation is a single pipe-delimited string, analogous to SMILES for individual molecules, but for multi-component polymer systems.

---

## The problem PUSMILES solves

Standard approaches to representing polymer formulations for ML (one-hot component name vectors, raw wt% arrays) discard molecular structure. A model that sees `"PPG": 55` learns nothing about polypropylene glycol's backbone, chain topology, or how it differs from PTMEG. PUSMILES encodes structure and composition together:

- **Structure**: every component is a validated SMILES string, not just a name
- **Composition**: weight percents are embedded directly in the string
- **Metadata**: molecular weight, functionality, and role travel with the component

---

## Repository contents

```
pusmiles/
├── pusmiles.py                   Core PUSMILES library
├── pusmiles_encoder.ipynb       Colab notebook: CSV → PUSMILES batch encoder
├── examples/
│   └── pusmiles_template.csv    Example formulations CSV (three PU systems)
├── README.md
└── LICENSE                      Apache 2.0
```

---

## Quick start

### Google Colab batch encoder

Open the notebook in Colab. Run all cells top to bottom. The notebook downloads `pusmiles.py` automatically; the only file you supply is your formulation CSV.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Loganz97/pusmiles/blob/main/pusmiles_encoder.ipynb)

See [examples/pusmiles_template.csv](examples/pusmiles_template.csv) for the expected CSV format.

### Python

```python
from pusmiles import PUSMILESBuilder

q = (PUSMILESBuilder()
     .add("PPG",          55, role="soft", Mn=3000, fn=3)
     .add("MDI_urethane", 40, role="hard")
     .add("BDO",           5, role="chain_ext")
     .build())

print(q)
# *OCC(C)*@55{role=soft,Mn=3000,fn=3}|*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1@40{role=hard}|OCCCCO@5{role=chain_ext}
```

### Dependencies

RDKit is optional but strongly recommended for SMILES validation on `build()`. Everything else is standard library.

```bash
pip install rdkit
```

Place `pusmiles.py` on your Python path or working directory. No package installation required.

---

## Notation

### Syntax

```
SMILES@wt%{key=val,key=val}
```

A formulation is multiple components joined by `|`:

```
SMILES1@wt1{...}|SMILES2@wt2{...}|...
```

### Metadata keys

| Key | Type | Description |
|---|---|---|
| `role` | string | Component role (see vocabulary below) |
| `Mn` | number | Number-average molecular weight (g/mol) |
| `Mw` | number | Weight-average molecular weight (g/mol) |
| `fn` | number | Chain-end functionality |
| `PDI` | number | Polydispersity index |
| `cas` | string | CAS registry number |
| `name` | string | Common or trade name |
| `encodable` | bool | `false` for non-molecular materials (talc, glass fiber, etc.) |

### Role vocabulary

| Role | Description |
|---|---|
| `soft` | Soft segment (polyol) |
| `hard` | Hard segment (reacted isocyanate) |
| `chain_ext` | Chain extender |
| `ionic_diol` | Ionic diol (e.g. DMPA in waterborne PU dispersions) |
| `filler` | Inorganic particle fillers |
| `flame_retardant` | Flame retardant additives |
| `plasticizer` | Plasticizers |
| `pigment` | Colorants and opacifiers |
| `blowing_agent` | Physical or chemical blowing agents |
| `catalyst` | Urethane/urea reaction catalysts |
| `surfactant` | Foam stabilizers and surfactants |
| `antioxidant` | Thermal/oxidative stabilizers |
| `UV_stabilizer` | Light stabilizers (HALS, UV absorbers) |
| `adhesion_promoter` | Adhesion promoters |
| `additive` | Catch-all for unclassified additives |

### SMILES conventions for polymer repeat units

Linear repeat units carry exactly two `*` attachment points:

```
*OCC(C)*              PPG repeat unit
*OC(=O)CCCCC(=O)*    PBA repeat unit
```

Isocyanate hard segments follow `*OC(=O)N-[core]-NC(=O)O*`. The `*` sits on oxygen because the polyol hydroxyl oxygen attacks the isocyanate carbonyl during urethane formation. Nitrogen stays bonded to the aromatic or aliphatic core:

```
*OC(=O)Nc1ccc(Cc2ccc(NC(=O)O*)cc2)cc1    MDI hard segment
*OC(=O)Nc1ccc(C)c(NC(=O)O*)c1            2,4-TDI hard segment
```

Chain-extender segments use the complementary pattern `*NC(=O)O-[diol]-OC(=O)N*`:

```
*NC(=O)OCCCCOC(=O)N*    BDO chain-extender segment
*NC(=O)OCCOC(=O)N*      EG chain-extender segment
```

Non-molecular materials with no discrete SMILES (talc, glass fiber, silicone surfactants, carbon black) use `[*]` as a placeholder with `encodable=false`:

```
[*]@2{role=surfactant,name=silicone_surfactant,encodable=false}
```

---

## Component library

`pusmiles.py` ships with three built-in libraries. Use these names as column headers in the batch encoder CSV or as arguments to `PUSMILESBuilder.add()`.

### Polymer repeat-unit library (62 entries)

#### Polyether polyols
| Name | Description |
|---|---|
| `PPG` | Polypropylene glycol |
| `PEG` | Polyethylene glycol |
| `PTMEG` | Poly(tetramethylene ether) glycol |
| `PBO` | Poly(butylene oxide) |
| `graft_polyol_PPG` | SAN-grafted PPG polymer polyol |
| `graft_polyol_PEG` | SAN-grafted PEG polymer polyol |

#### Polyester polyols
| Name | Description |
|---|---|
| `PBA` | Poly(butylene adipate) |
| `PEA` | Poly(ethylene adipate) |
| `PHA` | Poly(hexamethylene adipate) |
| `PCL` | Polycaprolactone diol |
| `PES` | Poly(ethylene succinate) |
| `PET` | Poly(ethylene terephthalate) |
| `PLA` | Polylactic acid |
| `PCDL` | Polycarbonate diol |

#### Isocyanate hard segments
| Name | Description | CAS |
|---|---|---|
| `MDI_urethane` | 4,4'-Diphenylmethane diisocyanate | 101-68-8 |
| `pMDI_urethane` | Polymeric MDI monomer fraction (fn=2) | 9016-87-9 |
| `pMDI_3ring_urethane` | Polymeric MDI 3-ring oligomer (fn=3) | N/A |
| `MDI_modified_urethane` | Carbodiimide-modified liquid MDI | N/A |
| `TDI_urethane` | 2,4-Toluene diisocyanate (alias) | 584-84-9 |
| `TDI_24_urethane` | 2,4-Toluene diisocyanate | 584-84-9 |
| `TDI_26_urethane` | 2,6-Toluene diisocyanate | 91-08-7 |
| `HDI_urethane` | Hexamethylene diisocyanate | 822-06-0 |
| `IPDI_urethane` | Isophorone diisocyanate | 4098-71-9 |
| `H12MDI_urethane` | 4,4'-Methylenebis(cyclohexyl isocyanate) | 5124-30-1 |
| `NDI_urethane` | 1,5-Naphthalene diisocyanate | 3173-72-6 |
| `PPDI_urethane` | p-Phenylene diisocyanate | 104-49-4 |
| `TODI_urethane` | 3,3'-Dimethyl-4,4'-biphenylene diisocyanate | 91-97-4 |
| `mXDI_urethane` | m-Xylylene diisocyanate | 3634-83-1 |
| `pXDI_urethane` | p-Xylylene diisocyanate | 4538-44-5 |

#### Chain-extender segments (reacted form, between isocyanate groups)
| Name | Description |
|---|---|
| `EG_urethane` | Ethylene glycol segment |
| `BDO_urethane` | 1,4-Butanediol segment |
| `HDO_urethane` | 1,6-Hexanediol segment |
| `HQEE_urethane` | Hydroquinone bis(2-hydroxyethyl ether) segment |

#### Diol chain extenders (small molecule, unreacted form)
`EG`, `BDO`, `PDO`, `HDO`, `DEG`, `NPG`, `HQEE`, `DMPA`

#### Amine chain extenders
`MOCA`, `DETDA`, `IPDA`, `EDA`

#### Crosslinkers
`glycerol`, `TMP`, `DEA`, `TEA`, `sorbitol`

#### Commodity polymers
`PE`, `PP`, `PS`, `PVC`, `PMMA`, `PAN`, `PVA`, `PVDF`, `PVAC`

### Blend library (4 entries)

Commercial multi-isomer grades that automatically expand to constituent fractions when used in the builder:

| Name | Composition | Notes |
|---|---|---|
| `TDI_80` | 80% 2,4-TDI + 20% 2,6-TDI | CAS 26471-62-5, standard flexible foam grade |
| `TDI_65` | 65% 2,4-TDI + 35% 2,6-TDI | Higher 2,6-content grade |
| `pMDI` | 50% MDI monomer (fn=2) + 50% 3-ring oligomer (fn=3) | Avg fn ≈ 2.5, approximates commercial polymeric MDI |
| `MDI_modified` | 80% MDI monomer + 20% branched fraction | Approximates carbodiimide-modified liquid MDI |

### Additive library (29 entries)

| Category | Names |
|---|---|
| Fillers | `CaCO3`, `BaSO4`, `silica`, `talc`, `kaolin`, `glass_fiber` |
| Flame retardants | `ATH`, `MEL`, `TCPP`, `TCEP`, `TEP`, `DMMP`, `DOPO` |
| Plasticizers | `DOA`, `DEP`, `DBP` |
| Blowing agents | `water`, `n_pentane`, `cyclopentane`, `isopentane`, `HCBA` |
| Catalysts | `DABCO`, `DBTDL`, `TEA_cat` |
| Surfactants | `silicone_surfactant` |
| Antioxidants | `BHT` |
| UV stabilizers | `benzophenone` |
| Pigments | `carbon_black`, `TiO2` |

---

## Python API reference

### `PUSMILESBuilder`

Fluent builder: the primary interface for constructing PUSMILES strings.

```python
from pusmiles import PUSMILESBuilder

# By weight percent
q = (PUSMILESBuilder()
     .add("PPG",          55, role="soft", Mn=3000, fn=3)
     .add("MDI_urethane", 40, role="hard")
     .add("BDO",           5, role="chain_ext")
     .build())

# By absolute weight, auto-normalized to 100 wt%
q = (PUSMILESBuilder()
     .add_weight("PTMEG",         100, role="soft", Mn=2000)
     .add_weight("MDI_urethane",   72, role="hard")
     .add_weight("BDO_urethane",   18, role="chain_ext")
     .build())

# With additives: role and CAS pulled from library automatically
q = (PUSMILESBuilder()
     .add("PPG",          70, role="soft", Mn=3000)
     .add("TDI_80",       25, role="hard")   # expands to 2,4 and 2,6 fractions
     .add_additive("CaCO3",   3)             # role=filler set automatically
     .add_additive("DABCO",   1)
     .add_additive("water",   1, role="blowing_agent")
     .build())

# PubChem lookup for names not in the library
q = (PUSMILESBuilder()
     .add("PPG",          85, role="soft", Mn=3000)
     .add_additive_by_name("triethyl phosphate", 15, role="flame_retardant")
     .build())
```

### `PUSMILES` object

```python
from pusmiles import PUSMILES

# Parse an existing PUSMILES string
q = PUSMILES.parse("*OCC(C)*@55{role=soft,Mn=3000}|OCCCCO@5{role=chain_ext}")

# Inspect components
for c in q.components:
    print(c.smiles, c.weight_pct, c.role, c.metadata)

# Composition queries
soft_pct = q.sum_by_role("soft")
total    = q.total_weight()

# Validation
result = q.validate()
result.print_report()

# Normalize to exactly 100 wt%
q_norm = q.normalize()

# ML feature dict: n_components, soft_pct, hard_pct, c0_wt, c0_mw, ...
features = q.to_features()

# Serialization
print(str(q))        # PUSMILES string
print(q.to_json())   # JSON with component breakdown
```

### PubChem utilities

```python
from pusmiles import pubchem_smiles_from_name, pubchem_smiles_from_cas, pubchem_info

smiles = pubchem_smiles_from_name("calcium carbonate")
smiles = pubchem_smiles_from_cas("471-34-1")
info   = pubchem_info("melamine")
# {"smiles": "...", "iupac_name": "...", "formula": "C3H6N6", "mw": 126.12, "cid": 7955}
```

### Library inspection

```python
from pusmiles import list_polymers, list_additives, list_blends

list_polymers()                       # all 62 polymer repeat-unit names
list_additives()                      # all 29 additive names
list_additives("flame_retardant")     # filtered by role
list_blends()                         # ['TDI_80', 'TDI_65', 'pMDI', 'MDI_modified']
```

---

## Batch encoder CSV format

| Column | Example header | Value |
|---|---|---|
| Sample ID | `Sample_ID` | Any text label; must be the first column |
| Component | `PPG`, `MDI_urethane`, `CaCO3` | Weight percent (0–100) |
| Molecular weight | `PPG_Mn` | Mn for the component named before `_Mn` |
| Functionality | `PPG_fn` | Chain-end functionality for the named component |

Zero or blank cells are skipped per row. Rows auto-normalize if total is outside 99–101%.
See [examples/pusmiles_template.csv](examples/pusmiles_template.csv) for a ready-to-use example.

---

## Design notes

**Non-encodable materials.** Talc, glass fiber, carbon black, and silicone surfactants have no discrete molecular SMILES. These are stored as `[*]` with `encodable=false`. ML models should handle this case explicitly, typically with a zero molecular feature vector, using `role` and `wt%` as the signal.

**Blend expansion.** Commercial isocyanate grades are mixtures of isomers. When a blend name is used, it expands to constituent components each carrying its proportional wt%. A TDI-80 formulation has two separate hard-segment entries in the output, which is structurally accurate and allows models to learn grade-dependent behavior.

**Urethane bond orientation.** The `*` placement in hard-segment SMILES follows reaction chemistry. In `*OC(=O)N-[core]-NC(=O)O*`, `*` is on oxygen because the polyol hydroxyl oxygen attacks the isocyanate carbon. Hard-segment and chain-extender SMILES are complementary: the O on the hard segment and the N on the chain-extender segment meet at the urethane bond.

---

## Citation

If you use PUSMILES in published research, please cite:

```
Hessefort, L. (2025). PUSMILES: Polyurethane SMILES for Materials Informatics.
https://github.com/Loganz97/pusmiles
```

---

## License

Apache License 2.0. See [LICENSE](LICENSE).
