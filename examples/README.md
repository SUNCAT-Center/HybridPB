# HybridPB Examples

This directory contains eight worked examples, grouped into three systems, covering different electrochemical
systems and calculation modes. Each example includes DFT structure files (`*.json`), a `label.csv`, and shell scripts that reproduce the published figures.

Run all commands from the example directory using `HybridPourbaix.py` in the package root:

```bash
cd examples/1_MNC/1_CoNC
python ../../../HybridPourbaix.py --help
```

## Examples Overview

| Example | System | Key Features | Applications |
|---------|--------|--------------|--------------|
| [1_MNC/1_CoNC](1_MNC/1_CoNC/) | Co–N₄–C SAC | Surface adsorbates, concentration sweeps, 1D pH profiles | ORR, HER |
| [1_MNC/2_TiNC](1_MNC/2_TiNC/) | Ti–N₄–C catalyst | 43 surface species, NOxHy vs. OxHy subsets | NO₃RR, selectivity |
| [1_MNC/3_FeNC](1_MNC/3_FeNC/) | Fe–N₄–C + GC-DFT | High/intermediate/low spin states, `--gc`, custom thermo data | Advanced ORR |
| [2_MnO2/1_MnO2_100](2_MnO2/1_MnO2_100/) | MnO₂ (100) surface | Oxide slab, K co-adsorption, custom `thermodynamic_data.jsonc` | OER, oxide electrochemistry |
| [2_MnO2/2_MnO2_110](2_MnO2/2_MnO2_110/) | MnO₂ (110) surface | Facet comparison with (100), OH coverage variants | Surface facet effects |
| [3_RuO2/1_RuO2](3_RuO2/1_RuO2/) | RuO₂ (110) surface | Bridge/cus site occupancy, bridge-row vacancies | OER, oxide dissolution |
| [3_RuO2/2_ReRuO2](3_RuO2/2_ReRuO2/) | Re-doped RuO₂ (110) | Same site model, Ru leaving the bridge row | Dopant effect on dissolution |
| [3_RuO2/3_ReRuO2](3_RuO2/3_ReRuO2/) | Re-doped RuO₂ (110) | Re rather than Ru leaving the bridge row | Selective dopant dissolution |

### Files in Each Example

| File | Description |
|------|-------------|
| `*.json` | ASE-readable DFT structures with total energies |
| `label.csv` | Structure labels, `#OH` counts, optional `G_corr` and GC-DFT `A,B,C` |
| `reference_energies.jsonc` | Gas and element reference energies pinned per example, passed via `--ref-energies` |
| `thermodynamic_data.jsonc` | Element-specific bulk/solution data, where the package defaults need overriding |
| `command.sh` | Full reproduction workflow with legend and color settings |
| `command-simple.sh` | Simplified workflow without legend placement flags |

Some scripts temporarily move JSON files into subdirectories (`vac/`, `OH/`, `no3rr/`) to select structure subsets. Read the script before running, or run individual `python` lines by hand.

## Learning Path

### Beginner — Example 1 (CoNC)

Learn surface-only and hybrid modes, ion concentration effects, and 1D energy profiles at fixed pH.

```bash
cd 1_MNC/1_CoNC
python ../../../HybridPourbaix.py --HER --OER --png --figx 4 --figy 4
python ../../../HybridPourbaix.py --hybrid --no-bulk --legend-in --pH 7 --Gmin -4 --Gmax 12
```

Concepts: reference surface selection (`--ref-json`), `--no-bulk`, `--concentration`, `--colors-2d`

### Intermediate — Example 2 (TiNC)

Multi-species nitrate reduction with competing HER. Compare OxHy-only vs. full NOxHy adsorbate sets.

```bash
cd 1_MNC/2_TiNC
python ../../../HybridPourbaix.py --line -0.7 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --suffix OxHy
python ../../../HybridPourbaix.py --hybrid --no-bulk --line -0.7 \
  --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --concentration 1e-3
```

Concepts: custom potential lines (`--line`), wide pH windows, subset filtering via file moves

### Advanced — Example 3 (FeNC)

Grand Canonical DFT with multiple spin states and bulk Fe oxide phases.

```bash
cd 1_MNC/3_FeNC
python ../../../HybridPourbaix.py --gc --fill 0.72 0.92 --Umin -1.5 --Umax 2.0 --legend-in
python ../../../HybridPourbaix.py --hybrid --no-bulk --gc --thermo-data ./thermodynamic_data.jsonc \
  --fill 0.72 0.92 --Gmin -6 --Gmax 8
```

Concepts: `--gc`, `--fill`, custom thermodynamic data, 1D y-axis limits (`--Gmin`/`--Gmax`)

Additional script: `command-toc.sh` (table-of-contents figures)

### Oxide Surfaces — Examples 4 & 5 (MnO₂)

MnO₂ (100) and (110) facet calculations with K co-adsorption and OH coverage variants.

```bash
cd 2_MnO2/1_MnO2_100
python ../../../HybridPourbaix.py --hybrid --HER --OER --Umin -0.5 --Umax 2.0 --suffix K
python ../../../HybridPourbaix.py --hybrid --no-bulk --thermo-data ./thermodynamic_data.jsonc \
  --Umin -0.5 --Umax 2.0 --suffix MnO4

cd ../../2_MnO2/2_MnO2_110
python ../../../HybridPourbaix.py --hybrid --no-bulk --Umin -0.5 --Umax 2.0
```

Concepts: per-example `thermodynamic_data.jsonc`, explicit bulk colors (`--colors-bulk`), facet-to-facet comparison

### Oxide Dissolution — Examples 6–8 (RuO₂)

RuO₂ (110) bridge/cus site occupancy, and how Re doping changes which cation leaves the bridge row.
`2_ReRuO2` lets Ru go, `3_ReRuO2` lets Re go.

```bash
cd 3_RuO2/1_RuO2
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk \
  --Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --legend-out \
  --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0
```

Concepts: two-site occupancy labels, colormap windowing (`--cmin-2d`/`--cmax-2d`/`--cgap-2d`), wide 1D energy ranges

## Getting Started

### Prerequisites

```bash
cd HybridPB
pip install -r requirements.txt
python HybridPourbaix.py --help
```

### Run a Full Example Script

```bash
cd examples/1_MNC/1_CoNC
bash command-simple.sh
```

### Quick Smoke Test (All Examples)

Runs a minimal surface calculation in each example directory:

```bash
cd HybridPB
for example in 1_MNC/1_CoNC 1_MNC/2_TiNC 1_MNC/3_FeNC \
               2_MnO2/1_MnO2_100 2_MnO2/2_MnO2_110 \
               3_RuO2/1_RuO2 3_RuO2/2_ReRuO2 3_RuO2/3_ReRuO2; do
    cd examples/$example
    python ../../../HybridPourbaix.py --figx 4 --figy 4 --suffix test
    cd ../../..
done
```

## Contributing Examples

We welcome new examples. Each submission should include:

- **Scientific relevance** — a real electrochemical system or methodology
- **Complete inputs** — all `*.json`, `label.csv`, and thermo data as needed
- **Reproducibility** — a `command.sh` or `command-simple.sh` that regenerates the figures
- **Thermodynamic consistency** — verified species energies and reference states

To submit: create a numbered directory under a system group (e.g. `5_MySystem/1_MyCase/`), add inputs and
scripts, then open a pull request.

## Further Reading

- Main documentation: [../README.md](../README.md)
- Repository: [https://github.com/SUNCAT-Center/HybridPB](https://github.com/SUNCAT-Center/HybridPB)

## Troubleshooting

| Problem | Suggestion |
|---------|------------|
| Missing thermo species | Check warnings; add entries to `thermodynamic_data.jsonc` or use `--show-thermo` |
| Wrong reference surface | Pass `--ref-json clean.json` (or your clean-slab file) |
| Too many legend entries | Use `--no-bulk` or filter structures by moving JSON files into subdirectories |
| GC-DFT not applied | Ensure `--gc` is set and `A,B,C` columns are filled in `label.csv` |

- **Issues**: [GitHub Issues](https://github.com/SUNCAT-Center/HybridPB/issues)
- **Questions**: [GitHub Discussions](https://github.com/SUNCAT-Center/HybridPB/discussions)
