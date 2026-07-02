# HybridPB: Hybrid Pourbaix Diagram Generation Tool

![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)

HybridPB is a Python tool for generating Pourbaix diagrams (potential–pH diagrams) for electrochemical systems. It combines DFT surface energies with experimental thermodynamic data, supporting surface-only, hybrid surface–bulk, and Grand Canonical DFT (GC-DFT) calculations.

The main entry point is `HybridPourbaix.py`.

## Key Features

- **Hybrid calculations**: Combine DFT surface slabs with bulk/solution species (ions, solids, gases, liquids)
- **Grand Canonical DFT**: GC-DFT corrections with potential-dependent energy terms (`A·U² + B·U + C`)
- **Thermodynamic integration**: Element-wise thermodynamic database in JSONC format
- **Flexible activity control**: Per-species, per-element, or global ion/gas activity via `conditions.jsonc`
- **Dual visualization**: 2D stability map and 1D energy profile at a fixed pH
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
git clone https://github.com/hailey-suncat/HybridPB.git
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
  "Fe": {
    "solids": {
      "Fe": 0,
      "Fe2O3": -197.3,
      "Fe3O4": -310.1
    },
    "ions": {
      "Fe++": -20.3,
      "Fe+++": -2.5
    },
    "gases": {
      "Fe(CO)5": -151.0
    }
  }
}
```

Ion formulas are parsed with pymatgen (`Fe++`, `MnO4-`, etc.). Condensed phases (solids, liquids) always have activity 1.

### 4. Reference Energies (Optional)
- **Filename**: `reference_energies.jsonc` (default in package root)
- **Format**: Element symbol → DFT reference energy (eV)

**Example:**
```json
{
  "Fe": -5.041720865,
  "Co": -4.55701137,
  "Ni": -1.6850032
}
```

### 5. Conditions File (Optional)
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
```

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
```

## Examples

Five examples are provided under `examples/`:

| Example | System | Focus |
|---------|--------|-------|
| [1_CoNC](examples/1_CoNC/) | Co–N₄–C SAC | Surface adsorbates, concentration effects, 1D pH scans |
| [2_TiNC](examples/2_TiNC/) | Ti–N₄–C catalyst | NO₃RR pathways, 43 surface species |
| [3_FeNC](examples/3_FeNC/) | Fe–N₄–C + GC-DFT | Spin states, GC-DFT vs. standard DFT |
| [4_MnO2_100](examples/4_MnO2_100/) | MnO₂ (100) surface | Oxide surface, K co-adsorption, custom thermo data |
| [5_MnO2_110](examples/5_MnO2_110/) | MnO₂ (110) surface | Facet comparison, OH coverage effects |

### Example Commands

```bash
# CoNC: hybrid without bulk phases
cd examples/1_CoNC
python ../../HybridPourbaix.py --hybrid --no-bulk --legend-in \
  --figx 4 --figy 4 --cmap-2d Purples --pH 7

# TiNC: nitrate reduction with custom concentration
cd examples/2_TiNC
python ../../HybridPourbaix.py --suffix NO3RR_mono \
  --figx 6 --figy 4 --cmap-2d Reds --concentration 1e-3

# FeNC: GC-DFT analysis
cd examples/3_FeNC
python ../../HybridPourbaix.py --gc --legend-in \
  --figx 6 --figy 4 --cmap-2d RdBu --suffix gc_analysis

# MnO2 (100): hybrid with custom bulk colors
cd examples/4_MnO2_100
python ../../HybridPourbaix.py --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
  --colors-2d dodgerblue lightskyblue --figx 6 --figy 6
```

Each example directory includes `command.sh` and `command-simple.sh` scripts with full reproduction workflows.

## Methodology

### Surface Energy Correction
DFT total energies are converted to relative formation energies using:
1. A **reference surface** — the highest-energy pure-metal slab (H/O only), or the structure specified by `--ref-json`
2. **Element reference energies** from `reference_energies.jsonc`
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
| `pourbaix_surface.png` | 2D diagram (surface-only mode) |
| `pourbaix_hybrid.png` | 2D diagram (hybrid mode) |
| `pourbaix_hybrid_pH{N}.png` | 1D energy profile at fixed pH |
| `pourbaix_bulk_{El}.png` | Bulk diagram per element (hybrid, no `--no-bulk`) |
| `{structure}.png` | Structure images (`--png`) |

Suffixes `_gc`, `_legend_in`, `_legend_out`, `_legend_up`, and `--suffix` are appended automatically.

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Citation

```bibtex
@article{jung2025hybridpb,
  title={Bridging Bulk and Surface Thermodynamics: A Hybrid Pourbaix Framework for Electrocatalyst Stability},
  author={Hyeonjung Jung and Evan Z. Carlson and Md Delowar Hossain and Michal Bajdich},
  year={2025},
  doi={10.26434/chemrxiv.15005305/v1}
}
```
