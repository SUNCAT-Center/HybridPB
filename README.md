# HybridPB: Hybrid Pourbaix Diagram Generation Tool

![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)

HybridPB is a Python tool for generating Pourbaix diagrams (potential–pH diagrams) for electrochemical systems. It combines DFT surface energies with experimental thermodynamic data, supporting surface-only, hybrid surface–bulk, and Grand Canonical DFT (GC-DFT) calculations.

The main entry point is `HybridPourbaix.py`.

## Key Features

- **Hybrid calculations**: Combine DFT surface slabs with bulk/solution species (ions, solids, gases, liquids)
- **Grand Canonical DFT**: GC-DFT corrections with potential-dependent energy terms (`A·U² + B·U + C`)
- **Thermodynamic integration**: Element-wise thermodynamic database in JSONC format
- **Selectable references**: Per-functional reference energies (`--functional`) and per-element reference models such as metal vs. oxide (`--ref-model`), regenerated from a structure database
- **Flexible activity control**: Per-species, per-element, or global ion/gas activity via `conditions.jsonc`
- **Dual visualization**: 2D stability map and 1D energy profile at a fixed pH
- **Two renderings**: Every 2D and bulk diagram is drawn twice — once from a (pH, U) grid scan, once by pymatgen's `PourbaixDiagram` with analytic phase boundaries (`--no-pymatgen` to skip)
- **Customizable plots**: Separate colormaps for bulk, 2D, and 1D plots; explicit color lists; legend placement
- **Electrochemical references**: HER/OER lines, custom potential lines, and shaded regions
- **Debug tools**: Thermodynamic data inspection, element counts, minimum-coordination diagnostics

## Requirements

### System Requirements
- Python 3.7 or higher
- Linux, macOS, or Windows

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | >=1.19.0 | Grid calculations and numerical arrays |
| `pandas` | >=1.2.0 | `label.csv` processing |
| `matplotlib` | >=3.3.0 | Pourbaix diagram plotting |
| `ase` | >=3.21.0 | Read/write JSON structure files |
| `mendeleev` | >=0.9.0 | Element sorting by atomic number |
| `pymatgen` | >=2022.0.0 | Ion formula parsing |

### Installation
```bash
git clone https://github.com/SUNCAT-Center/HybridPB.git
cd HybridPB
pip install -r requirements.txt
```

### Verification
```bash
python HybridPourbaix.py --help
```

## Quick Start

### Basic Usage
```bash
# Surface-only Pourbaix diagram
python HybridPourbaix.py --json-dir ./structures --HER --OER --show-fig

# Custom pH and potential ranges
python HybridPourbaix.py --json-dir ./structures --pHmin 0 --pHmax 14 --Umin -2 --Umax 2
```

### Hybrid Mode
```bash
# Hybrid surface + bulk with default thermodynamic data
python HybridPourbaix.py --hybrid --json-dir ./structures --HER --OER

# Hybrid without per-element bulk diagrams
python HybridPourbaix.py --hybrid --no-bulk --legend-in --pH 7
```

### Advanced Usage
```bash
# GC-DFT with explicit Gibbs corrections from label.csv
python HybridPourbaix.py --hybrid --gc --gibbs --csv-dir ./labels

# Custom thermodynamic data and reference energies
python HybridPourbaix.py --thermo-data ./custom_thermo.jsonc --ref-energies ./custom_ref.jsonc

# Use the reference energies of another functional
python HybridPourbaix.py --hybrid --functional BEEF

# Reference the metals to their oxides instead of the elemental metals
python HybridPourbaix.py --hybrid --ref-model oxide

# Per-species activity overrides
python HybridPourbaix.py --hybrid --conditions ./conditions.jsonc --concentration 1e-6
```

## Input Files

### 1. Structure Files (Required)
- **Format**: ASE-readable JSON files with atomic structures and DFT energies
- **Location**: `--json-dir` (default: current directory)

### 2. Label File (Optional)
- **Filename**: `label.csv`
- **Location**: `--csv-dir` or `--label-csv`
- **Columns**: `json_name, label, #OH, G_corr, A, B, C`

| Column | Description |
|--------|-------------|
| `json_name` | JSON filename |
| `label` | Display name in the legend |
| `#OH` | Number of OH groups (used for default Gibbs correction) |
| `G_corr` | Explicit Gibbs correction (eV); applied only with `--gibbs` |
| `A`, `B`, `C` | GC-DFT parameters; used with `--gc` (`ΔGh = C`, energy = `A·U² + B·U + C`) |

**Example `label.csv`:**
```csv
structure1.json,Fe-OH,1,0.05,0.0,0.0,-123.45
structure2.json,Fe-O,0,0.03,0.1,-0.02,-124.67
structure3.json,Fe,0,0.0,0.0,0.0,-125.89
```

If `label.csv` is absent, chemical formulas from the JSON files are used as labels.

### 3. Thermodynamic Data (Optional)
- **Filename**: `thermodynamic_data.jsonc` (default in package root)
- **Format**: Element-keyed JSONC with phase categories; energies in **kcal/mol**
- **Phases**: `ions`, `solids`, `gases`, `liquids`
- **Comments**: Single-line `//` comments are supported (JSONC)

**Example:**
```jsonc
{
  "N": {
    "gases": {
      "N2": 0,
      "NH3": -3.976,
      "NO2": 12.390
    }
  },
  "Fe": {
    "solids": {
      "Fe": 0,
      "FeO": -58.880,
      "Fe2O3": -177.100,
      "Fe(OH)3": -166.000
    },
    "ions": {
      "Fe++": -20.300,
      "Fe+++": -2.530,
      "FeO4--": -111.685
    }
  }
}
```

Ion formulas are parsed with pymatgen (`Fe++`, `MnO4-`, etc.). Condensed phases (solids, liquids) always have activity 1. Element blocks are ordered by atomic number.

A species may be given the value `null` to reserve a slot whose formation energy is not known yet. Such entries are skipped with a warning instead of breaking the run, and `build_reference_energies.py` reports exactly which placeholder is still empty.

### 4. Reference Energies (Optional)
- **Filename**: `reference_energies.jsonc` (default in package root, or `--ref-energies`)
- **Format**: Functional name → `gases` (H2/H2O DFT total energies) + `elements` (reference energy per element atom, keyed by the formula of the reference compound), in eV
- **Selection**: `--functional NAME` picks the functional block (default: `PBE`, matched case-insensitively); `--ref-model` picks the reference within it (default: `metal`)

**Example:**
```jsonc
{
  "PBE": {
    "gases": {
      "H2": -6.82099190,
      "H2O": -14.24491949
    },
    "elements": {
      "N": -8.56951971,
      "P": { "P": -5.40833708, "P2O5": -5.1334630293750045 },
      "Fe": { "Fe": -8.2408819, "Fe2O3": -5.51563177 },  // Fe2O3: Ueff = 4.3 eV
      "Ni": { "Ni": -5.47060251 }
    }
  },
  "BEEF": {
    "gases": {
      "H2": -7.17103764,
      "H2O": -12.82990184
    },
    "elements": {
      "P": { "P2O5": -4.1600280112500005 }
    }
  }
}
```

H2 and H2O are required — they set the H and O chemical potentials. Elements are ordered by atomic number. When a reference was calculated with DFT+U (`ldau: true` and a non-zero `U − J` for that element), the effective U is noted in a trailing comment.

An element maps either to a plain number (one model-independent reference) or to a `{formula: energy}` dict. A `--ref-model` token is matched against the keys directly, then as a category resolved from each key's composition — `metal` (element only), `oxide`, `hydride`, `hydroxide`:

```bash
--ref-model metal            # every element referenced to its elemental phase
--ref-model oxide            # every element referenced to its oxide
--ref-model metal Mn=oxide   # metal by default, oxide for Mn
--ref-model metal Fe=Fe3O4   # an explicit formula when several oxides exist
```

If the requested model is missing but the element has exactly one entry, that one is used and a note is printed; if a category matches several keys, the run aborts and lists them. Missing H2/H2O energies or an element absent from the block also abort with an explicit message. Zero-point, heat-capacity, and entropy corrections stay in `HybridPourbaix.py`; only DFT total energies and derived references live in this file.

Example directories ship their own pinned `reference_energies.jsonc` and pass it via `--ref-energies`, so they
reproduce even if the package-level file changes.

### 5. Reference Structure Database (Optional)
- **Location**: `reference_energies/<functional>/<formula>.json`
- **Purpose**: The DFT structures the numbers above are derived from, so `reference_energies.jsonc` can be regenerated instead of hand-edited

Each file is an ASE-readable structure with an attached calculator (e.g. a VASP `final_with_calculator.json`), named after its reduced formula. The element it references is the single constituent that is not O or H, so no per-element folders are needed:

```
reference_energies/
  PBE/
    Fe.json        # cell Fe2   -> reference for Fe
    Fe2O3.json     # cell Fe4O6 -> reference for Fe
    Ru.json
    RuO2.json
```

A file whose formula contains two non-O/H elements is skipped (the referenced element would be ambiguous), and a filename that disagrees with the structure's composition is reported with the expected name.

Rebuild the JSONC from it with:

```bash
python build_reference_energies.py                 # rebuild every functional found
python build_reference_energies.py --dry-run       # print without writing
python build_reference_energies.py --show-source   # report the structure behind each value
python build_reference_energies.py --functional PBE
```

Each structure is converted to a reference energy per element atom:

```
mu_M = E_total/n_M − ΔGf/n_M − (n_O/n_M)·g_O − (n_H/n_M)·g_H
g_O  = g_H2O − g_H2 + ΔGf(H2O)      # O referenced to ½O2
g_H  = g_H2 / 2
```

`ΔGf` is the experimental formation energy of the reduced compound, looked up in `thermodynamic_data.jsonc` (Gibbs corrections of the compound are assumed zero). A pure element has `ΔGf = 0` and no O/H, so it reduces to `E_total/n_M`. Structures whose reduced formula has no entry — or only a `null` placeholder — in `thermodynamic_data.jsonc` are reported as `SKIPPED` and left out. Elements absent from the database keep their current values.

### 6. Conditions File (Optional)
- **Filename**: `conditions.jsonc` (default in package root, or `--conditions`)
- **Purpose**: Override ion/gas activity on a per-species or per-element basis

**Example:**
```jsonc
{
  "defaults": {
    "ions": 1e-6,
    "gases": 1e-6
  },
  "elements": {
    "Fe": { "ions": 1e-3 }
  },
  "species": {
    "MnO4-": 1e-4
  }
}
```

**Priority**: `species` > `elements.{el}.{phase}` > `defaults.{phase}` > CLI (`--concentration` / `--pressure`)

## Command Line Options

Options are grouped in `HybridPourbaix.py --help`. Summary below.

### Input Paths
```bash
--json-dir PATH         # Directory with JSON structure files (default: .)
--csv-dir PATH          # Directory with label.csv (default: .)
--label-csv PATH        # Explicit path to label.csv
--thermo-data PATH      # Path to thermodynamic_data.jsonc
--ref-energies PATH     # Path to reference_energies.jsonc
--conditions PATH       # Path to conditions.jsonc
```

### Calculation Modes
```bash
--hybrid                # Enable hybrid surface–bulk mode
--no-bulk               # Skip per-element bulk Pourbaix diagrams in hybrid mode
--gc                    # Apply Grand Canonical DFT (A, B, C from label.csv)
--gibbs                 # Use G_corr from label.csv instead of #OH-based correction
--ref-json FILE         # Reference surface JSON (default: auto-detect pure-metal slab)
--suffix STRING         # Output filename suffix
```

### Thermodynamic Conditions
```bash
--concentration FLOAT   # Default ion activity in M (default: 1e-6)
--pressure FLOAT        # Default gas activity in atm (default: 1e-6)
--functional NAME       # Functional block in reference_energies.jsonc (default: PBE)
--ref-model MODEL ...   # Reference model + per-element exceptions, e.g. metal Mn=oxide
```

### Axis Range
```bash
--pHmin, --pHmax FLOAT  # pH range (default: 0–14)
--Umin, --Umax FLOAT    # Potential vs. SHE in V (default: -1–3)
--tick FLOAT            # Grid resolution (default: 0.01)
--pH INT                # Fixed pH for 1D energy plot (default: 0)
--Gmin, --Gmax FLOAT    # Y-axis limits for 1D plot
```

### Figure
```bash
--figx, --figy FLOAT    # Figure size in inches (default: 4 × 4)
--HER, --OER            # Draw HER/OER reference lines
--line FLOAT            # Custom reference line (V vs. SHE at pH 0)
--fill LOW HIGH         # Shaded region between two reference lines
--legend-in             # Legend inside plot
--legend-out            # Legend outside plot (right)
--legend-up             # Legend above plot
--label-fontsize FLOAT  # Size of the pymatgen domain labels (default: 10)
--label-color COLOR     # Color of those labels; name or hex (default: black)
```

The pymatgen figure names each domain in place rather than building a legend
box, so any of the three `--legend-*` flags turns those labels on and they all
produce the same picture.

### Colormaps
Separate settings for bulk/combination, 2D original surfaces, and 1D plots:

```bash
# Bulk / hybrid combination phases
--cmap, --cmin, --cmax, --cgap STRING/FLOAT   # default: Greys, 0.1, 0.7, 0.0
--colors-bulk COLOR [COLOR ...]                 # explicit colors (overrides --cmap)

# Original DFT surfaces (2D plot)
--cmap-2d, --cmin-2d, --cmax-2d, --cgap-2d    # default: RdBu, 0.0, 1.0, 0.2
--colors-2d COLOR [COLOR ...]

# 1D energy plot
--cmap-1d, --cmin-1d, --cmax-1d, --cgap-1d    # default: Spectral, 0.0, 1.0, 0.0
--colors-1d COLOR [COLOR ...]
```

`--cgap` skips the center of diverging colormaps (useful for neutral reference states).

### Display / Debug
```bash
--show-fig              # Display matplotlib window
--show-thermo           # Print parsed thermodynamic species
--show-element          # Print element lists
--show-count            # Print minimum atom counts per element
--show-label            # Print structure labels
--show-min-coord        # Print lowest (pH, U) coordinate per stable phase
```

### Output
```bash
--png                   # Export structure PNGs from JSON files
--png-rotation STRING   # ASE view rotation (default: '-90x, -90y, 0z')
--no-pymatgen           # Skip the pymatgen rendering, drawn by default
```

## Examples

Eight examples are provided under `examples/`, grouped by system:

| Example | System | Focus |
|---------|--------|-------|
| [1_MNC/1_CoNC](examples/1_MNC/1_CoNC/) | Co–N₄–C SAC | Surface adsorbates, concentration effects, 1D pH scans |
| [1_MNC/2_TiNC](examples/1_MNC/2_TiNC/) | Ti–N₄–C catalyst | NO₃RR pathways, 43 surface species |
| [1_MNC/3_FeNC](examples/1_MNC/3_FeNC/) | Fe–N₄–C + GC-DFT | Spin states, GC-DFT vs. standard DFT |
| [2_MnO2/1_MnO2_100](examples/2_MnO2/1_MnO2_100/) | MnO₂ (100) surface | Oxide surface, K co-adsorption, custom thermo data |
| [2_MnO2/2_MnO2_110](examples/2_MnO2/2_MnO2_110/) | MnO₂ (110) surface | Facet comparison, OH coverage effects |
| [3_RuO2/1_RuO2](examples/3_RuO2/1_RuO2/) | RuO₂ (110) surface | Bridge/cus site occupancy, bridge-row vacancies |
| [3_RuO2/2_ReRuO2](examples/3_RuO2/2_ReRuO2/) | Re-doped RuO₂ (110) | Ru leaving the bridge row |
| [3_RuO2/3_ReRuO2](examples/3_RuO2/3_ReRuO2/) | Re-doped RuO₂ (110) | Re leaving the bridge row |

Each example directory contains its own `reference_energies.jsonc`, and every command in its `command*.sh` passes it with `--ref-energies`. The examples therefore stay reproducible no matter how the package-level reference energies change.

### Example Commands

```bash
# CoNC: hybrid without bulk phases
cd examples/1_MNC/1_CoNC
python ../../../HybridPourbaix.py --hybrid --no-bulk --legend-in \
  --figx 4 --figy 4 --cmap-2d Purples --pH 7

# TiNC: nitrate reduction with custom concentration
cd examples/1_MNC/2_TiNC
python ../../../HybridPourbaix.py --suffix NO3RR_mono \
  --figx 6 --figy 4 --cmap-2d Reds --concentration 1e-3

# FeNC: GC-DFT analysis
cd examples/1_MNC/3_FeNC
python ../../../HybridPourbaix.py --gc --legend-in \
  --figx 6 --figy 4 --cmap-2d RdBu --suffix gc_analysis

# MnO2 (100): hybrid with custom bulk colors
cd examples/2_MnO2/1_MnO2_100
python ../../../HybridPourbaix.py --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
  --colors-2d dodgerblue lightskyblue --figx 6 --figy 6

# RuO2 (110): bridge/cus site occupancy
cd examples/3_RuO2/1_RuO2
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk \
  --Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --legend-out
```

Each example directory includes `command.sh` and `command-simple.sh` scripts with full reproduction workflows.

## Methodology

### Surface Energy Correction
DFT total energies are converted to relative formation energies using:
1. A **reference surface** — the highest-energy pure-metal slab (H/O only), or the structure specified by `--ref-json`
2. **Element and H2/H2O reference energies** from the `--functional` / `--ref-model` selection in `reference_energies.jsonc`
3. **Gibbs corrections** — either explicit `G_corr` (`--gibbs`) or estimated from `#OH` counts

In hybrid mode, a new reference is selected from surface+bulk combinations (preferring neutral phases with basic solid references).

### Hybrid Combinations
For each surface slab, missing elements (relative to the maximum count across all slabs) are compensated by adding thermodynamic species (ions, solids, gases, liquids). Combined states appear as `clean+Fe(s)` in the legend.

### Relative Gibbs Energy
At each (pH, U) grid point, the most stable phase is determined by:

```
ΔG = (A·U² + B·U + ΔGh)_surf − (A·U² + B·U + ΔGh)_ref
     + (H − 2O − e)·U + kT·ln(10)·(H − 2O)·pH
```

where `A = B = 0` unless `--gc` is set.

### Pourbaix Diagram Output
- **2D map**: Stable phase regions over pH and potential
- **1D profile**: Relative energies vs. potential at `--pH`, with dashed lines for second-lowest phases
- **Bulk diagrams** (`--hybrid` only, unless `--no-bulk`): Per-element bulk Pourbaix from thermodynamic data alone

## Output Files

| File | Description |
|------|-------------|
| `pourbaix_surface.pdf` | 2D diagram (surface-only mode) |
| `pourbaix_hybrid.pdf` | 2D diagram (hybrid mode) |
| `pourbaix_hybrid_pH{N}.pdf` | 1D energy profile at fixed pH |
| `pourbaix_bulk_{El}.pdf` | Bulk diagram per element (hybrid, no `--no-bulk`) |
| `pourbaix_surface_pymatgen.pdf` | The 2D diagram again, from pymatgen |
| `pourbaix_hybrid_pymatgen.pdf` | ditto, hybrid mode |
| `pourbaix_bulk_{El}_pymatgen.pdf` | ditto, bulk |
| `{structure}.png` | Structure images (`--png`) |

Suffixes `_gc`, `_legend_in`, `_legend_out`, `_legend_up`, and `--suffix` are appended automatically.
The pymatgen names carry the same suffixes except that the three legend flags collapse to a plain
`_legend`, since they all give the same figure. There is no pymatgen counterpart to the 1D profile,
and none is written under `--gc`, whose `A·U² + B·U` term `PourbaixEntry` cannot express — bulk
species never carry that term, so their diagram still appears.

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Citation

Jung, H.; Carlson, E. Z.; Hossain, M. D.; Bajdich, M. Bridging Bulk and Surface Thermodynamics: A Hybrid Pourbaix Framework for Electrocatalyst Stability. *ChemRxiv* **2026**. DOI: 10.26434/chemrxiv.15005305/v1
