"""Generate Pourbaix diagrams from DFT surface energies and thermodynamic data.

Same as HybridPourbaix.py, except that every run also draws the phase diagram
through pymatgen's PourbaixDiagram, which solves the boundaries analytically
rather than scanning a (pH, U) grid. Pass --no-pymatgen to skip it.
"""

import argparse
import glob
import json
import os
import re
from collections import Counter
from math import log10

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from ase.io import read, write
from mendeleev import element
from pymatgen.analysis.pourbaix_diagram import (
    MU_H2O,
    IonEntry,
    PourbaixDiagram,
    PourbaixEntry,
    PourbaixPlotter,
)
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SUBSCRIPT_NUMS = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}
SUPERSCRIPT_NUMS = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}
DIVERGING_COLORMAPS = ['RdBu', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

# Thermodynamic constants (set by init_thermo_constants)
const = water_formation_energy = gh = go = goh = dgh = dgo = dgoh = None

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_jsonc(file_path):
    """Load JSON file with single-line comments (JSONC format)."""
    with open(file_path, 'r') as f:
        content = re.sub(r'//.*', '', f.read())
    return json.loads(content)


def format_name(formula):
    """Convert ion formula to display name with sub/superscripts."""
    plus_count = formula.count('+')
    minus_count = formula.count('-')
    base_formula = formula.replace('+', '').replace('-', '')

    formatted = ''.join(
        SUBSCRIPT_NUMS[char] if char.isdigit() else char
        for char in base_formula
    )

    if plus_count > 0:
        charge = '' if plus_count == 1 else SUPERSCRIPT_NUMS.get(str(plus_count), str(plus_count))
        return formatted + charge + '⁺'
    if minus_count > 0:
        charge = '' if minus_count == 1 else SUPERSCRIPT_NUMS.get(str(minus_count), str(minus_count))
        return formatted + charge + '⁻'
    return formatted


def format_df_for_display(df):
    """Format numeric columns for console display."""
    df_copy = df.copy()
    for col in df_copy.select_dtypes(include=[np.number]).columns:
        if col == 'conc':
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.0e}")
        else:
            df_copy[col] = df_copy[col].apply(lambda x: f"{x:.2f}")
    return df_copy


# ---------------------------------------------------------------------------
# Thermodynamic calculations
# ---------------------------------------------------------------------------

def dg(surf, pH, U, ref_surf):
    """Gibbs free energy relative to reference surface."""
    surface_term = (
        (surf['A'] * U ** 2 + surf['B'] * U + surf['ΔGh'])
        - (ref_surf['A'] * U ** 2 + ref_surf['B'] * U + ref_surf['ΔGh'])
    )
    U_coeff = surf['H'] - 2 * surf['O'] - surf['e']
    pH_coeff = surf['H'] - 2 * surf['O']
    return surface_term + U_coeff * U + const * pH_coeff * pH


def init_thermo_constants(gas_energies):
    """Initialize module-level thermodynamic constants.

    gas_energies holds the H2/H2O DFT total energies (eV) of the selected
    functional, as read from reference_energies.jsonc.
    """
    global const, water_formation_energy, gh, go, goh, dgh, dgo, dgoh

    calmol = 23.061
    kb = 8.617e-5  # eV/K
    T = 298.15  # K
    const = kb * T * np.log(10)  # 0.0592 eV
    water_formation_energy = 56.690 / calmol # 2.46 eV at 300 K

    if T != 298.15:
        print(f"Warning: Temperature is not 300 K. The water formation energy is not valid.")

    h2 = gas_energies['H2']
    h2o = gas_energies['H2O']
    gh2o = h2o + 0.558 - 0.675 + 0.103
    gh2 = h2 + 0.268 - 0.408 + 0.0905

    gh = gh2 / 2
    go = gh2o - gh2
    goh = gh2o - gh2 / 2

    dgo = 0.064 + 0.034 - 0.060
    dgoh = 0.376 + 0.042 - 0.066
    dgh = dgoh - dgo


def parse_thermo_entry(formula, energy, el, remaining_elements,
                       conc, phase_suffix):
    """Parse one thermodynamic entry into a row dict, or None on failure."""
    calmol = 23.061
    try:
        # reduced_dict = Ion.from_formula(formula).to_reduced_dict
        full_dict = Ion.from_formula(formula)
    except Exception:
        return None, None
    
    o_count = full_dict['O']
    row = {
        # 'dG': energy / calmol + o_count * 4.92/2 + const * log10(conc), # water formation energy at 300 K
        'ΔGf': energy / calmol,
        'ΔGh': energy / calmol + o_count * water_formation_energy + const * log10(conc),
        'e': int(full_dict.charge),
        'conc': conc,
        'name': format_name(formula) + phase_suffix,
        'A': 0.0,
        'B': 0.0,
    }

    for elem in remaining_elements:
        row[elem] = int(full_dict.get(elem, 0))

    el_count = full_dict.get(el, 1)
    if el_count > 1:
        row['ΔGf'] /= el_count
        row['ΔGh'] /= el_count
        row['e'] /= el_count
        for elem in remaining_elements:
            row[elem] /= el_count
        subscript = SUBSCRIPT_NUMS.get(str(int(el_count)), str(int(el_count)))
        # row['name'] = f'¹⁄{subscript}' + row['name']

    return row, full_dict


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def make_color_values(n_colors, cmap_name, cmin, cmax, cgap):
    """Generate colormap sample values, optionally skipping the center gap."""
    if n_colors == 0:
        return np.array([])

    if cgap > 0 and cmap_name in DIVERGING_COLORMAPS:
        left_end = 0.5 - cgap / 2
        right_start = 0.5 + cgap / 2
        n_left = int(np.ceil(n_colors / 2))
        n_right = n_colors - n_left
        parts = []
        if n_left > 0:
            parts.append(np.linspace(cmin, left_end, n_left, endpoint=True))
        if n_right > 0:
            parts.append(np.linspace(right_start, cmax, n_right, endpoint=True))
        return np.concatenate(parts) if parts else np.linspace(cmin, cmax, n_colors)

    return np.linspace(cmin, cmax, n_colors)


def parse_colors(color_names):
    """Convert color name strings to matplotlib RGBA tuples."""
    if color_names is None:
        return None

    colors = []
    for color_name in color_names:
        try:
            colors.append(mcolors.to_rgba(color_name))
        except ValueError:
            try:
                hex_name = color_name if color_name.startswith('#') else f'#{color_name}'
                colors.append(mcolors.to_rgba(hex_name))
            except ValueError:
                print(f"Warning: Invalid color name '{color_name}', using black instead.")
                colors.append(mcolors.to_rgba('black'))
    return colors


def resolve_plot_colors(n_colors, explicit_color_names, cmap_name, cmin, cmax, cgap, default_cmap):
    """Return plot colors from explicit names or a colormap."""
    if explicit_color_names:
        parsed = parse_colors(explicit_color_names)
        return [parsed[i % len(parsed)] for i in range(n_colors)]

    colormap = getattr(plt.cm, cmap_name, default_cmap)
    color_values = make_color_values(n_colors, cmap_name, cmin, cmax, cgap)
    return colormap(color_values)


def show_used_colors(title, colors, labels):
    """Print hex codes of colors actually used in a plot."""
    print(f"\n{title} colors:")
    for color, label in zip(colors, labels):
        print(f"  {mcolors.to_hex(color)}  {label}")


def apply_legend(args, fontsize=12):
    """Apply legend placement based on CLI flags."""
    legend_kwargs = dict(fontsize=fontsize, handlelength=3, edgecolor='black')
    if args.legend_in:
        plt.legend(ncol=1, loc='upper right', **legend_kwargs)
    elif args.legend_out:
        plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., ncol=1, **legend_kwargs)
    elif args.legend_up:
        plt.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', borderaxespad=0., ncol=3, **legend_kwargs)


def plot_water_stability(pHrange, args):
    """Draw HER/OER stability lines."""
    if args.OER:
        plt.plot(pHrange, 1.23 - pHrange * const, '--', lw=1, color='mediumblue')
    if args.HER:
        plt.plot(pHrange, 0 - pHrange * const, '--', lw=1, color='mediumblue')
    if args.line:
        plt.plot(pHrange, args.line - pHrange * const, '--', lw=1, color='mediumblue')
    if args.fill:
        plt.fill_between(pHrange, args.fill[0] - pHrange * const, args.fill[1] - pHrange * const, alpha=0.2, color='red')

    # plt.plot(pHrange, -0.768925499892501+pHrange*0, '--', lw=1, color='white') #oh(hs)
    # plt.plot(pHrange, -1.10346332406801+pHrange*0, '--', lw=1, color='white') #clean(hs)
    # plt.plot(pHrange, -1.01280687446284+pHrange*0, '--', lw=1, color='white') #clean(is)
    # plt.plot(pHrange, -1.01280687446284+pHrange*0, '--', lw=1, color='green') #clean(is)
    # plt.plot(pHrange, -0.092131432+pHrange*0, '--', lw=1, color='green') #vac(*)
    # plt.plot(pHrange, -0.934679247+pHrange*0, '--', lw=1, color='green') #vac(*H2)


def find_lowest_indices(surfs, pHrange, Urange, ref_surf):
    """Return 2D array of lowest-energy surface indices over (U, pH) grid."""
    lowest = np.full((len(Urange), len(pHrange)), np.nan)
    for j, pH in enumerate(pHrange):
        for i, U in enumerate(Urange):
            values = [dg(s, pH, U, ref_surf) for s in surfs]
            lowest[i, j] = min(range(len(values)), key=lambda k: values[k])
    return lowest


def build_output_names(args):
    """Build output PNG basename and suffix from CLI args."""
    png_name = 'pourbaix_hybrid' if args.hybrid else 'pourbaix_surface'
    suffix = ''
    if args.gc:
        suffix += '_gc'
    if args.legend_in:
        suffix += '_legend_in'
    elif args.legend_out:
        suffix += '_legend_out'
    elif args.legend_up:
        suffix += '_legend_up'
    if args.suffix:
        suffix += '_' + args.suffix
    return png_name, suffix


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate Pourbaix diagram')

    paths = parser.add_argument_group('input paths')
    paths.add_argument('--json-dir', type=str, default='.', help='Folder with JSON structure files')
    paths.add_argument('--csv-dir', type=str, default='.', help='Folder with label.csv')
    paths.add_argument('--label-csv', type=str, help='Path to label.csv (overrides --csv-dir)')
    paths.add_argument('--thermo-data', type=str, help='Path to thermodynamic_data.jsonc')
    paths.add_argument('--ref-energies', type=str, help='Path to reference_energies.jsonc')

    mode = parser.add_argument_group('mode')
    mode.add_argument('--hybrid', action='store_true', help='Hybrid bulk + surface mode')
    mode.add_argument('--gc', action='store_true', help='Grand Canonical DFT (A, B, C columns)')
    mode.add_argument('--no-bulk', action='store_true', help='Skip bulk phases in hybrid mode')
    mode.add_argument('--suffix', type=str, default='', help='Output filename suffix')

    thermo = parser.add_argument_group('thermodynamics')
    thermo.add_argument('--concentration', type=float, default=1e-6,
                        help='Default ion activity (M); used when conditions.jsonc has no override')
    thermo.add_argument('--pressure', type=float, default=1e-6,
                        help='Default gas activity (atm); used when conditions.jsonc has no override')
    thermo.add_argument('--conditions', type=str,
                        help='Path to conditions.jsonc (per-species activity/pressure overrides)')
    thermo.add_argument('--functional', type=str, default='PBE',
                        help='Functional block to read from reference_energies.jsonc (default: PBE)')
    thermo.add_argument('--ref-model', type=str, nargs='+', default=['metal'], metavar='MODEL',
                        help='Reference model for element energies: a global default plus optional '
                             'per-element exceptions, e.g. --ref-model metal Mn=oxide. Accepts a '
                             'category (metal, oxide, hydride, hydroxide) or an explicit formula '
                             'such as Fe=Fe3O4 (default: metal)')
    thermo.add_argument('--gibbs', action='store_true', help='Apply G_corr from label.csv')
    thermo.add_argument('--ref-json', type=str, help='JSON filename for reference surface (default: auto-detect)')

    axis = parser.add_argument_group('axis range')
    axis.add_argument('--tick', type=float, default=0.01, help='Grid tick size')
    axis.add_argument('--pHmin', type=float, default=0, help='Minimum pH')
    axis.add_argument('--pHmax', type=float, default=14, help='Maximum pH')
    axis.add_argument('--Umin', type=float, default=-1, help='Minimum electrode potential (V)')
    axis.add_argument('--Umax', type=float, default=3, help='Maximum electrode potential (V)')
    axis.add_argument('--pH', type=int, default=0, help='pH for 1D plot')
    axis.add_argument('--Gmin', type=float, help='Y-axis min for 1D plot')
    axis.add_argument('--Gmax', type=float, help='Y-axis max for 1D plot')

    figure = parser.add_argument_group('figure')
    figure.add_argument('--figx', type=float, default=4, help='Figure width (inches)')
    figure.add_argument('--figy', type=float, default=4, help='Figure height (inches)')
    figure.add_argument('--HER', action='store_true', help='Draw hydrogen evolution reaction line')
    figure.add_argument('--OER', action='store_true', help='Draw oxygen evolution reaction line')
    figure.add_argument('--line', type=float, help='Draw custom line')
    figure.add_argument('--fill', type=float, nargs=2, help='Draw fill between two lines')
    figure.add_argument('--legend-in', action='store_true', help='Place legend inside plot')
    figure.add_argument('--legend-out', action='store_true', help='Place legend outside plot')
    figure.add_argument('--legend-up', action='store_true', help='Place legend above plot')
    figure.add_argument('--label-fontsize', type=float, default=10,
                        help='Font size of the domain labels pymatgen writes inside each region (default: 10)')
    figure.add_argument('--label-color', type=str, default='black',
                        help='Color of those domain labels; any matplotlib color name '
                             'or hex code (default: black)')

    cmap_bulk = parser.add_argument_group('colormap (bulk / combination)')
    cmap_bulk.add_argument('--cmap', type=str, default='Greys', help='Colormap name')
    cmap_bulk.add_argument('--cmin', type=float, default=0.1, help='Colormap value min')
    cmap_bulk.add_argument('--cmax', type=float, default=0.7, help='Colormap value max')
    cmap_bulk.add_argument('--cgap', type=float, default=0.0, help='Gap at colormap center')
    cmap_bulk.add_argument('--colors-bulk', type=str, nargs='+',
                           help='Explicit color names (overrides --cmap)')

    cmap_2d = parser.add_argument_group('colormap (2D original surfaces)')
    cmap_2d.add_argument('--cmap-2d', type=str, default='RdBu', help='Colormap name')
    cmap_2d.add_argument('--cmin-2d', type=float, default=0.0, help='Colormap value min')
    cmap_2d.add_argument('--cmax-2d', type=float, default=1.0, help='Colormap value max')
    cmap_2d.add_argument('--cgap-2d', type=float, default=0.2, help='Gap at colormap center')
    cmap_2d.add_argument('--colors-2d', type=str, nargs='+',
                         help='Explicit color names (overrides --cmap-2d)')

    cmap_1d = parser.add_argument_group('colormap (1D plot)')
    cmap_1d.add_argument('--cmap-1d', type=str, default='Spectral', help='Colormap name')
    cmap_1d.add_argument('--cmin-1d', type=float, default=0.0, help='Colormap value min')
    cmap_1d.add_argument('--cmax-1d', type=float, default=1.0, help='Colormap value max')
    cmap_1d.add_argument('--cgap-1d', type=float, default=0.0, help='Gap at colormap center')
    cmap_1d.add_argument('--colors-1d', type=str, nargs='+',
                         help='Explicit color names (overrides --cmap-1d)')

    display = parser.add_argument_group('display / debug')
    display.add_argument('--show-fig', action='store_true', help='Show matplotlib figure interactively')
    display.add_argument('--show-thermo', action='store_true', help='Print thermodynamic corrections')
    display.add_argument('--show-element', action='store_true', help='Print element list')
    display.add_argument('--show-count', action='store_true', help='Print atom counts per structure')
    display.add_argument('--show-label', action='store_true', help='Print structure labels')
    display.add_argument('--show-min-coord', action='store_true', help='Print minimum coordination info')
    display.add_argument('--show-colors', action='store_true', help='Print hex codes of colors used in each plot')

    output = parser.add_argument_group('output')
    output.add_argument('--png', action='store_true', help='Export structure PNGs from JSON files')
    output.add_argument('--png-rotation', type=str, default='-90x, -90y, 0z', help='ASE view rotation for structure PNGs')
    output.add_argument('--no-pymatgen', action='store_true',
                        help='Skip the pymatgen PourbaixDiagram rendering, which is '
                             'drawn alongside the grid one by default (it is skipped '
                             'under --gc either way)')

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_labels(args, json_files):
    """Load labels and optional corrections from label.csv."""
    label_csv_path = args.label_csv or os.path.join(args.csv_dir, 'label.csv')
    file_labels = {}
    file_gibbs_corrections = {}
    file_gc_params = {}
    file_oh_counts = {}

    if os.path.exists(label_csv_path):
        label_df = pd.read_csv(
            label_csv_path, header=None,
            names=['json_name', 'label', '#OH', 'G_corr', 'A', 'B', 'C'],
        )
        for _, row in label_df.iterrows():
            json_name = row['json_name']
            file_labels[json_name] = row['label']
            if pd.notna(row.get('#OH')):
                file_oh_counts[json_name] = int(row['#OH'])
            if args.gibbs and pd.notna(row.get('G_corr')):
                file_gibbs_corrections[json_name] = float(row['G_corr'])
            if args.gc and all(pd.notna(row.get(col)) for col in ('A', 'B', 'C')):
                file_gc_params[json_name] = {
                    'A': float(row['A']), 'B': float(row['B']), 'C': float(row['C']),
                }
    else:
        for json_file in json_files:
            file_labels[os.path.basename(json_file)] = read(json_file).get_chemical_formula()

    return file_labels, file_gibbs_corrections, file_gc_params, file_oh_counts


def build_surfs(json_files, file_labels, file_gibbs_corrections, file_gc_params, args):
    """Build surface rows from JSON structure files."""
    elements = set()
    for json_file in json_files:
        elements.update(read(json_file).get_chemical_symbols())
    sorted_elements = sorted(elements, key=lambda x: element(x).atomic_number)
    if 'H' not in sorted_elements:
        sorted_elements.append('H')
    if 'O' not in sorted_elements:
        sorted_elements.append('O')
    sorted_elements = sorted(sorted_elements, key=lambda x: element(x).atomic_number)
    min_counts = {el: None for el in sorted_elements}

    surfs = []
    for json_file in json_files:
        atoms = read(json_file)
        symbol_count = Counter(atoms.get_chemical_symbols())
        row = {'ΔGh': atoms.get_potential_energy(), 'e': 0, 'conc': 1}
        for el in sorted_elements:
            count = symbol_count[el]
            row[el] = float(count)
            min_counts[el] = count if min_counts[el] is None else min(min_counts[el], count)

        json_basename = os.path.basename(json_file)
        row['name'] = file_labels.get(json_basename, json_file)
        row['_json'] = json_basename
        row['gibbs_corr'] = file_gibbs_corrections.get(json_basename, 0.0) if args.gibbs else 0.0

        if args.gc and json_basename in file_gc_params:
            gc = file_gc_params[json_basename]
            row['A'], row['B'] = gc['A'], gc['B']
            row['ΔGh'] = gc['C']
        else:
            row['A'], row['B'] = 0.0, 0.0
        surfs.append(row)

    for row in surfs:
        for el in sorted_elements:
            row[el] -= min_counts[el]

    zero_cols = [el for el in sorted_elements if all(row[el] == 0 for row in surfs) and el not in ('H', 'O')]
    remaining_elements = [el for el in sorted_elements if el not in zero_cols]
    for row in surfs:
        for col in zero_cols:
            del row[col]

    unique_elements = [el for el in remaining_elements if el not in ('O', 'H')]
    return surfs, sorted_elements, remaining_elements, unique_elements, min_counts


def find_ref_surface(surfs, unique_elements, args):
    """Find reference surface (highest energy among pure-metal surfaces)."""
    if args.ref_json:
        ref_name = os.path.basename(args.ref_json)
        if not ref_name.endswith('.json'):
            ref_name += '.json'
        for i, row in enumerate(surfs):
            if row.get('_json') == ref_name:
                return i, surfs[i]
        raise SystemExit(f"Reference surface not found: {ref_name}")

    ref_candidates = [
        (i, float(row['ΔGh']))
        for i, row in enumerate(surfs)
        if all(row.get(elem, 0.0) == 0.0 for elem in unique_elements)
    ]
    if not ref_candidates:
        return None, None

    ref_surf_idx, _ = max(ref_candidates, key=lambda x: x[1])
    return ref_surf_idx, surfs[ref_surf_idx]


def _is_ref_surf_combination(surf, ref_surf):
    """Return True if surf is a compound built from ref_surf."""
    return surf['name'].startswith(ref_surf['name'] + '+')


def _has_basic_solid_elements(surf, unique_elements):
    """Return True if every unique element is present as a pure (s) phase (e.g. clean+P(s)+N(s))."""
    if not all(surf.get(el, 0.0) > 0.0 for el in unique_elements):
        return False
    solid_parts = [part for part in surf['name'].split('+')[1:] if part.endswith('(s)')]
    for el in unique_elements:
        el_solid = format_name(el) + '(s)'
        if not any(part == el_solid for part in solid_parts):
            return False
    return True


def find_new_ref_surface(surfs, ref_surf, unique_elements):
    """Pick a new reference surface from ref_surf-based combinations.

    Priority:
    1. All unique_elements present as basic solids (e.g. clean+P(s)+N(s))
    2. charge=0 with lowest dG among ref_surf combinations
    3. Lowest dG among all surfaces
    """
    ref_combos = [
        (i, surf) for i, surf in enumerate(surfs)
        if _is_ref_surf_combination(surf, ref_surf)
    ]

    def pick_lowest(candidates):
        idx, surf = min(candidates, key=lambda item: float(item[1]['ΔGh']))
        return surf, idx

    basic_solid = [
        item for item in ref_combos
        if _has_basic_solid_elements(item[1], unique_elements)
    ]
    if basic_solid:
        return pick_lowest(basic_solid)

    neutral = [
        item for item in ref_combos
        if int(item[1].get('e', 0)) == 0
    ]
    if neutral:
        return pick_lowest(neutral)

    return pick_lowest(list(enumerate(surfs)))


def surface_formation_corrections(surfs, ref_surf, file_oh_counts, args, unique_elements, ref_energies):
    """Apply formation energy and optional Gibbs corrections."""
    for row in surfs:
        if args.gibbs:
            gibbs_correction = row['gibbs_corr']
        else:
            oh_count = file_oh_counts.get(row.get('_json', ''), 0)
            gibbs_correction = (row['H'] - oh_count) * dgh + (row['O'] - oh_count) * dgo + oh_count * dgoh

        formation_correction = -(row['H'] * gh + row['O'] * go)
        for el in unique_elements:
            if el not in ref_energies:
                raise SystemExit(
                    f"Element '{el}' has no reference energy for functional "
                    f"'{args.functional}'. Add it to reference_energies.jsonc."
                )
            formation_correction -= row[el] * ref_energies[el]

        row['ΔGh'] = row['ΔGh'] + gibbs_correction + formation_correction

    ref_surf_energy = ref_surf['ΔGh']
    for row in surfs:
        row['ΔGh'] = row['ΔGh'] - ref_surf_energy
    
    for row in surfs:
        row['ΔGf'] = row['ΔGh'] - row['O'] * water_formation_energy

def load_conditions(args):
    """Load activity/pressure overrides from conditions.jsonc."""
    script_dir = os.path.dirname(__file__)
    cond_path = args.conditions or os.path.join(script_dir, 'conditions.jsonc')
    return load_jsonc(cond_path) if os.path.exists(cond_path) else {}


def _cli_phase_default(phase_key, args):
    """CLI fallback activity for a thermodynamic phase."""
    if phase_key == 'ions':
        return args.concentration
    if phase_key == 'gases':
        return args.pressure
    return 1.0


def resolve_activity(formula, el, phase_key, conditions, args):
    """Resolve species activity: species > element+phase > defaults > CLI.

    Solids and liquids are always activity 1 (pure condensed phases).
    """
    if phase_key in ('solids', 'liquids'):
        return 1.0

    species_overrides = conditions.get('species', {})
    if formula in species_overrides:
        return float(species_overrides[formula])

    element_overrides = conditions.get('elements', {}).get(el, {})
    if phase_key in element_overrides:
        return float(element_overrides[phase_key])

    defaults = conditions.get('defaults', {})
    if phase_key in defaults:
        return float(defaults[phase_key])

    return float(_cli_phase_default(phase_key, args))


def process_thermo_phase(phase_data, phase_suffix, phase_key, el, remaining_elements,
                         conditions, args):
    """Process one phase (ions/solids/gases/liquids) from thermodynamic data."""
    rows = []
    for formula, energy in phase_data.items():
        if energy is None:
            print(f"WARNING: '{formula}' has no formation energy yet in thermodynamic data; skipped.")
            continue
        conc = resolve_activity(formula, el, phase_key, conditions, args)
        row, full_dict = parse_thermo_entry(
            formula, energy, el, remaining_elements, conc, phase_suffix,
        )
        if args.show_thermo:
            if full_dict is not None:
                print(f"    {formula}: {full_dict}, energy: {energy} kcal/mol, "
                      f"activity: {conc:.0e}")
            else:
                print(f"    {formula}: parsing failed, energy: {energy} kcal/mol")
        if row is not None:
            rows.append(row)
    return rows


def load_ref_energies(args):
    """Load gas and element reference energies for the selected functional.

    Returns (gas_energies, element_energies) from the --functional block of
    reference_energies.jsonc.
    """
    script_dir = os.path.dirname(__file__)
    ref_path = args.ref_energies or os.path.join(script_dir, 'reference_energies.jsonc')
    if not os.path.exists(ref_path):
        raise SystemExit(f"Reference energy file not found: {ref_path}")

    data = load_jsonc(ref_path)

    functional = args.functional
    if functional not in data:
        matches = [key for key in data if key.lower() == functional.lower()]
        if not matches:
            available = ', '.join(data) or '(none)'
            raise SystemExit(
                f"Functional '{functional}' not found in {ref_path}. Available: {available}"
            )
        functional = matches[0]

    block = data[functional] or {}
    gas_energies = block.get('gases') or {}
    element_energies = block.get('elements') or {}

    missing_gases = [gas for gas in ('H2', 'H2O') if gas not in gas_energies]
    if missing_gases:
        raise SystemExit(
            f"Functional '{functional}' in {ref_path} is missing gas energies: "
            f"{', '.join(missing_gases)}"
        )
    if not element_energies:
        print(f"WARNING: Functional '{functional}' has no element reference energies in "
              f"{ref_path}! Energy corrections may be inaccurate.")

    return gas_energies, element_energies


REF_MODEL_ALIASES = {
    'metal': frozenset(),
    'element': frozenset(),
    'oxide': frozenset({'O'}),
    'hydride': frozenset({'H'}),
    'hydroxide': frozenset({'O', 'H'}),
}


def classify_ref_key(key, el):
    """Category ('metal', 'oxide', ...) of a reference key from its formula."""
    try:
        composition = Ion.from_formula(key)
        others = frozenset(str(symbol) for symbol in composition.elements) - {el}
    except Exception:
        return None

    for alias, elements in REF_MODEL_ALIASES.items():
        if others == elements:
            return alias
    return None


def parse_ref_model(tokens):
    """Split --ref-model tokens into (global model, {element: model})."""
    global_model = 'metal'
    per_element = {}
    for token in tokens:
        if '=' in token:
            el, _, model = token.partition('=')
            if not el or not model:
                raise SystemExit(f"Invalid --ref-model entry '{token}'. Use the form El=model.")
            per_element[el] = model
        else:
            global_model = token
    return global_model, per_element


def resolve_element_energies(element_energies, unique_elements, args):
    """Pick one reference energy per element according to --ref-model.

    An element maps either to a plain number (single, model-independent
    reference) or to a {formula: energy} dict. A --ref-model token is matched
    against the keys directly, then as a category alias ('metal', 'oxide',
    'hydride', 'hydroxide') resolved from each key's composition. When neither
    matches but the element has exactly one entry, that one is used with a
    notice.
    """
    global_model, per_element = parse_ref_model(args.ref_model)
    resolved = {}
    used_models = {}

    for el in unique_elements:
        entry = element_energies.get(el)
        if entry is None:
            continue
        if not isinstance(entry, dict):
            resolved[el] = float(entry)
            used_models[el] = 'single'
            continue

        model = per_element.get(el, global_model)
        if model in entry:
            resolved[el] = float(entry[model])
            used_models[el] = model
            continue

        matches = [key for key in entry if classify_ref_key(key, el) == model]
        if len(matches) == 1:
            resolved[el] = float(entry[matches[0]])
            used_models[el] = matches[0]
            continue
        if len(matches) > 1:
            raise SystemExit(
                f"Element '{el}' has several '{model}' references for functional "
                f"'{args.functional}': {', '.join(matches)}. Pick one with "
                f"--ref-model {el}={matches[0]}"
            )

        if len(entry) == 1:
            only_model, only_energy = next(iter(entry.items()))
            print(f"NOTE: Element '{el}' has no '{model}' reference energy for functional "
                  f"'{args.functional}'; using the only available one ('{only_model}').")
            resolved[el] = float(only_energy)
            used_models[el] = only_model
            continue

        available = ', '.join(entry) or '(none)'
        raise SystemExit(
            f"Element '{el}' has no '{model}' reference energy for functional "
            f"'{args.functional}'. Available references: {available}"
        )

    if args.show_thermo and resolved:
        print(f"\nElement reference energies (functional: {args.functional}, "
              f"default model: {global_model}):")
        for el, energy in resolved.items():
            print(f"    {el}: {energy} eV ({used_models[el]})")

    return resolved


def load_thermo_data(args):
    """Load thermodynamic species data from JSONC file."""
    script_dir = os.path.dirname(__file__)
    thermo_path = args.thermo_data or os.path.join(script_dir, 'thermodynamic_data.jsonc')
    return load_jsonc(thermo_path) if os.path.exists(thermo_path) else {}


def load_thermo_species(unique_elements, remaining_elements, thermo_data, ref_energies, args):
    """Load ions, solids, gases, liquids from thermodynamic data."""
    conditions = load_conditions(args)
    ions, solids, gases, liquids = [], [], [], []
    phase_config = [
        ('ions', '(aq)', ions),
        ('solids', '(s)', solids),
        ('gases', '(g)', gases),
        ('liquids', '(l)', liquids),
    ]

    for el in unique_elements:
        if el not in ref_energies:
            print(f"WARNING: Element '{el}' not found in the '{args.functional}' block of "
                  "reference_energies.jsonc! Energy corrections may be inaccurate.")

        if el not in thermo_data:
            print(f"WARNING: Element '{el}' not found in thermodynamic data! "
                  "This element will be ignored in Pourbaix diagram calculations.")
            if args.show_thermo:
                print(f"\n{el}: No thermodynamic data found")
            continue

        if args.show_thermo:
            print(f"\n{el} thermodynamic data:")

        el_ions, el_solids, el_gases, el_liquids = [], [], [], []
        el_lists = {'ions': el_ions, 'solids': el_solids, 'gases': el_gases, 'liquids': el_liquids}

        for phase_key, suffix, target_list in phase_config:
            phase_data = thermo_data[el].get(phase_key, {})
            if not phase_data:
                continue
            if args.show_thermo:
                print(f"  {phase_key.capitalize()} reduced dict:")
            rows = process_thermo_phase(
                phase_data, suffix, phase_key, el, remaining_elements, conditions, args,
            )
            target_list.extend(rows)
            el_lists[phase_key].extend(rows)

        yield el, el_ions, el_solids, el_gases, el_liquids


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_bulk_diagram(el, bulks, pHrange, Urange, pHmin, pHmax, Umin, Umax, args, suffix):
    """Plot bulk Pourbaix diagram for one element."""
    nbulks = len(bulks)
    if nbulks == 0:
        return

    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    ax.axis([pHmin, pHmax, Umin, Umax])
    ax.set_xlabel('pH', fontsize=12)
    ax.set_ylabel('Potential (V vs. SHE)', fontsize=12)
    plt.xticks(np.arange(pHmin, pHmax + 1, 2))

    lowest_bulks = find_lowest_indices(bulks, pHrange, Urange, bulks[0])
    unique_ids = []
    for i in range(len(Urange)):
        for j in range(len(pHrange)):
            bulk_id = int(lowest_bulks[i, j])
            if bulk_id not in unique_ids:
                unique_ids.append(bulk_id)

    colors = resolve_plot_colors(
        len(unique_ids), args.colors_bulk, args.cmap, args.cmin, args.cmax, args.cgap, plt.cm.Greys,
    )
    if args.show_colors:
        show_used_colors(f'Bulk ({el})', colors, [bulks[int(b)]['name'] for b in unique_ids])
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(np.arange(len(unique_ids) + 1) - 0.5, cmap.N)
    id_map = {val: idx for idx, val in enumerate(unique_ids)}
    mapped_bulks = np.vectorize(id_map.get)(lowest_bulks)

    for idx, bulk_id in enumerate(reversed(unique_ids)):
        plt.plot([], [], color=colors[len(unique_ids) - 1 - idx], linewidth=5,
                 label=bulks[int(bulk_id)]['name'])

    pH_grid, U = np.meshgrid(pHrange, Urange)
    plt.pcolormesh(pH_grid, U, mapped_bulks, cmap=cmap, norm=norm, rasterized=True)
    plot_water_stability(pHrange, args)
    apply_legend(args)

    out_name = f'pourbaix_bulk_{el}{suffix}.pdf'
    plt.savefig(out_name, dpi=300, bbox_inches='tight')
    print(f"Bulk Pourbaix diagram saved as {out_name}")
    if args.show_fig:
        plt.tight_layout()
        plt.show()
    plt.close(fig)

    return {bulks[bid]['name']: colors[i] for i, bid in enumerate(unique_ids)}


def build_hybrid_combinations(surfs, nsurfs, unique_elements, ions, solids, gases, liquids):
    """Build surface+bulk hybrid combinations."""
    new_surfs = []
    max_counts = {unique_elem: max(surfs[k][unique_elem] for k in range(nsurfs)) for unique_elem in unique_elements}
    for k in range(nsurfs):
        for unique_elem in unique_elements:
            if surfs[k][unique_elem] == max_counts[unique_elem]:
                continue
            else:
                diff = max_counts[unique_elem] - surfs[k][unique_elem]
            for compound in ions + solids + gases + liquids:
                if compound[unique_elem] == 1:
                    new_surf = {}
                    for key in surfs[k]:
                        if key == 'name':
                            new_surf[key] = surfs[k][key] + '+' + compound[key]
                        elif key == 'conc':
                            new_surf[key] = surfs[k][key] * compound[key] ** diff
                        elif key in ('gibbs_corr', 'A', 'B', '_json'):
                            new_surf[key] = surfs[k][key]
                        else:
                            new_surf[key] = surfs[k][key] + compound[key] * diff
                    new_surfs.append(new_surf)
    return max_counts, new_surfs


def surf_df_columns(remaining_elements, is_gc):
    """Return DataFrame column list for surface rows."""
    cols = ['ΔGf', 'ΔGh', 'e'] + remaining_elements + ['conc', 'name']
    if is_gc:
        cols += ['A', 'B']
    return cols


def plot_2d_diagram(surfs, save_surfs, lowest_surfaces, pHrange, Urange,
                    pHmin, pHmax, Umin, Umax, ref_surf_idx, args, png_name, suffix):
    """Plot 2D Pourbaix diagram."""
    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    ax.axis([pHmin, pHmax, Umin, Umax])
    ax.set_xlabel('pH', fontsize=12)
    ax.set_ylabel('Potential (V vs. SHE)', fontsize=12)
    plt.xticks(np.arange(pHmin, pHmax + 1, 2))

    unique_ids = []
    for i in range(len(Urange)):
        for j in range(len(pHrange)):
            surf_id = int(lowest_surfaces[i, j])
            if surf_id not in unique_ids:
                unique_ids.append(surf_id)

    save_names = {s['name'] for s in save_surfs}
    save_ids = [sid for sid in unique_ids if surfs[sid]['name'] in save_names]
    new_ids = [sid for sid in unique_ids if surfs[sid]['name'] not in save_names]

    colors_save = resolve_plot_colors(
        len(save_ids), args.colors_2d, args.cmap_2d, args.cmin_2d, args.cmax_2d, args.cgap_2d,
        plt.cm.RdBu,
    ) if save_ids else []
    colors_new = resolve_plot_colors(
        len(new_ids), args.colors_bulk, args.cmap, args.cmin, args.cmax, args.cgap, plt.cm.Greys,
    ) if new_ids else []

    if args.show_colors:
        if len(save_ids) > 0:
            show_used_colors('2D original surfaces', colors_save, [surfs[int(i)]['name'] for i in save_ids])
        if len(new_ids) > 0:
            show_used_colors('2D combinations', colors_new, [surfs[int(i)]['name'] for i in new_ids])

    all_colors, id_map = [], {}
    for i, surf_id in enumerate(save_ids):
        all_colors.append(colors_save[i])
        id_map[surf_id] = len(all_colors) - 1
    for i, surf_id in enumerate(new_ids):
        all_colors.append(colors_new[i])
        id_map[surf_id] = len(all_colors) - 1

    cmap = mcolors.ListedColormap(all_colors)
    norm = mcolors.BoundaryNorm(np.arange(len(all_colors) + 1) - 0.5, cmap.N)
    mapped = np.vectorize(id_map.get)(lowest_surfaces)

    for i, surf_id in enumerate(reversed(save_ids)):
        plt.plot([], [], color=colors_save[len(save_ids) - 1 - i], linewidth=5,
                 label=surfs[int(surf_id)]['name'])
    for i, surf_id in enumerate(reversed(new_ids)):
        plt.plot([], [], color=colors_new[len(new_ids) - 1 - i], linewidth=5,
                 label=surfs[int(surf_id)]['name'])

    pH_grid, U = np.meshgrid(pHrange, Urange)
    plt.pcolormesh(pH_grid, U, mapped, cmap=cmap, norm=norm, rasterized=True)
    plot_water_stability(pHrange, args)
    apply_legend(args)

    out_name = f'{png_name}{suffix}.pdf'
    plt.savefig(out_name, dpi=300, bbox_inches='tight')
    print(f"Pourbaix diagram saved as {out_name}")
    if args.show_fig:
        plt.tight_layout()
        plt.show()

    # Hand the exact pairing back so the pymatgen figure paints each surface the
    # same color rather than re-deriving it and risking a drift.
    color_map = {surfs[sid]['name']: colors_save[i] for i, sid in enumerate(save_ids)}
    color_map.update({surfs[sid]['name']: colors_new[i] for i, sid in enumerate(new_ids)})
    return color_map


def plot_1d_diagram(surfs, nsurfs, Urange, Umin, Umax, target_pH, ref_surf, unique_ids,
                    args, png_name, suffix, tick):
    """Plot 1D energy diagram at fixed pH."""
    lowest_pH = np.zeros(len(Urange))
    second_lowest_pH = np.zeros(len(Urange))
    for Uindex, U in enumerate(Urange):
        values = [dg(s, target_pH, U, ref_surf) for s in surfs]
        sorted_idx = sorted(range(len(values)), key=lambda k: values[k])
        lowest_pH[Uindex] = sorted_idx[0]
        if len(sorted_idx) > 1:
            second_lowest_pH[Uindex] = sorted_idx[1]

    unique_ids_set = set(unique_ids)
    unique_second_set = set(np.unique(second_lowest_pH.astype(int)))
    all_relevant = unique_ids_set | unique_second_set

    energies_at_U0 = [
        (k, surfs[k]['e'] - surfs[k]['H'] + 2 * surfs[k]['O'],
         dg(surfs[k], pH=target_pH, U=0, ref_surf=ref_surf))
        for k in all_relevant
    ]
    sorted_ids = [k for k, _, _ in sorted(energies_at_U0, key=lambda x: (x[1], x[2]))]

    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    ax.set_xlabel('Potential (V vs. SHE)', fontsize=12)
    ax.set_ylabel('Relative Energy (ΔG, eV)', fontsize=12)
    ax.tick_params()

    colors = resolve_plot_colors(
        len(sorted_ids), args.colors_1d, args.cmap_1d, args.cmin_1d, args.cmax_1d, args.cgap_1d,
        plt.cm.RdBu,
    )
    if args.show_colors:
        show_used_colors(f'1D (pH {target_pH})', colors, [surfs[k]['name'] for k in sorted_ids])

    Urange_1d = np.arange(Umin, Umax, tick)
    all_energies = []

    for k in range(nsurfs):
        energies = np.array([dg(surfs[k], target_pH, U, ref_surf) for U in Urange_1d])
        all_energies.extend(energies)
        if k not in sorted_ids:
            ax.plot(Urange_1d, energies, color='lightgray', alpha=0.3, lw=0.5)

    for j, k in enumerate(sorted_ids):
        energies = np.array([dg(surfs[k], target_pH, U, ref_surf) for U in Urange_1d])
        if k in unique_ids_set:
            ax.plot(Urange_1d, energies, label=surfs[k]['name'], lw=1, color=colors[j])
        elif k in unique_second_set:
            ax.plot(Urange_1d, energies, '--', label=surfs[k]['name'], lw=1, color=colors[j])

    Gmin = args.Gmin if args.Gmin is not None else min(all_energies)
    Gmax = args.Gmax if args.Gmax is not None else max(all_energies)
    ax.set_ylim(Gmin, Gmax)
    ax.set_xlim(Umin, Umax)
    apply_legend(args)

    out_name = f'{png_name}_pH{target_pH}{suffix}.pdf'
    plt.savefig(out_name, bbox_inches='tight')
    print(f"Pourbaix diagram saved as {out_name}")
    if args.show_fig:
        plt.tight_layout()
        plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# pymatgen normalizes every Pourbaix energy per non-H/O atom, while dg() does
# not normalize at all and never reads the metal count: its slopes are
# H - 2O for pH and H - 2O - charge for potential. Composing each entry as one
# dummy atom plus its H and O therefore reproduces dg() exactly, keeps the
# normalization factor at 1 for every entry, and works for adsorbate-only runs
# that carry no metal at all.
PMG_DUMMY_ELEMENT = 'X'


def build_pourbaix_entries(surfs):
    """Wrap the surface rows as pymatgen PourbaixEntry objects.

    The two codes agree on the free-energy form: pymatgen evaluates
    energy + npH * PREFAC * pH + nPhi * V with npH = H - 2O and
    nPhi = npH - charge, which is what dg() computes on top of ΔGh.

    ΔGh is handed over directly, undoing pymatgen's water term so its
    reconstruction lands back on it. Rebuilding it from ΔGf instead would mean
    re-deriving the concentration term, which is not a flat
    const * log10(conc): a species carrying several atoms of the element takes
    that share of it, so N2(g) gets half.

    Every row goes in as an IonEntry so the labels from label.csv survive;
    solid-phase entries would be renamed after their reduced formula, which
    also collides whenever two structures share a formula. Concentration is
    pinned at 1 M — the term is already inside ΔGh, and PourbaixDiagram would
    otherwise overwrite each entry from its conc_dict.
    """
    entries = []
    for surf in surfs:
        counts = {PMG_DUMMY_ELEMENT: 1}
        for el in ('H', 'O'):
            if surf[el]:
                counts[el] = surf[el]
        ion = Ion(Composition(counts), charge=surf['e'])
        energy = surf['ΔGh'] + MU_H2O * surf['O']
        entries.append(
            PourbaixEntry(IonEntry(ion, energy, name=surf['name']), concentration=1.0)
        )
    return entries


def plot_pymatgen_diagram(surfs, color_map, args, png_name, suffix):
    """Draw the phase diagram a second time through pymatgen.

    The grid scan in find_lowest_indices() resolves boundaries to one tick;
    pymatgen solves them analytically, so this is the same thermodynamics with
    smooth phase lines.
    """
    # Under --gc the surfaces carry A and B, making ΔG quadratic in U, which
    # PourbaixEntry cannot express. Bulk species never do, so test the rows
    # rather than the flag.
    if any(surf['A'] or surf['B'] for surf in surfs):
        print(f"Skipping the pymatgen {png_name} diagram: --gc makes ΔG quadratic "
              "in U, which PourbaixEntry cannot express.")
        return

    entries = build_pourbaix_entries(surfs)
    if not entries:
        print("Skipping the pymatgen diagram: no structures to plot.")
        return

    # A name is written inside every domain, which crowds the narrow ones. Only
    # do it when a legend was requested; otherwise leave the plot bare. The
    # labels are drawn here rather than by pymatgen, whose to_latex_string()
    # reads the '+' joining two species and the '-' in names like *OH-bridge as
    # chemical formula markup and subscripts them.
    label_domains = bool(args.legend_in or args.legend_out or args.legend_up)

    fig, ax = plt.subplots(figsize=(args.figx, args.figy))
    try:
        diagram = PourbaixDiagram(
            entries,
            conc_dict={PMG_DUMMY_ELEMENT: 1.0},
            filter_solids=False,
        )
        plotter = PourbaixPlotter(diagram)
        # Water and neutral-axis lines are left to plot_water_stability() so the
        # two diagrams carry the same --OER/--HER/--line/--fill markings.
        plotter.get_pourbaix_plot(
            limits=[[args.pHmin, args.pHmax], [args.Umin, args.Umax]],
            label_domains=False,
            show_water_lines=False,
            show_neutral_axes=False,
            ax=ax,
        )
    except Exception as exc:
        plt.close(fig)
        print(f"Skipping the pymatgen diagram: {exc}")
        return

    # pymatgen only strokes the domain outlines, so fill them here using the
    # pairing plot_2d_diagram settled on, and write the label alongside.
    label_color = (parse_colors([args.label_color]) or [(0, 0, 0, 1)])[0]
    for entry, vertices in getattr(diagram, '_stable_domain_vertices', {}).items():
        color = color_map.get(entry.name)
        if color is not None:
            ax.add_patch(Polygon(vertices, closed=True, facecolor=color,
                                 edgecolor='none', zorder=0))
        if label_domains:
            ax.annotate(entry.name, np.mean(vertices, axis=0), ha='center',
                        va='center', fontsize=args.label_fontsize, color=label_color)

    # pymatgen strokes the boundaries at lw=3; match the axes frame instead.
    frame_lw = ax.spines['left'].get_linewidth()
    for line in ax.lines:
        line.set_linewidth(frame_lw)
    ax.set_xlabel('pH', fontsize=12)
    ax.set_ylabel('Potential (V vs. SHE)', fontsize=12)

    plot_water_stability(np.arange(args.pHmin, args.pHmax + args.tick, args.tick), args)
    ax.axis([args.pHmin, args.pHmax, args.Umin, args.Umax])

    # --legend-in/-out/-up all land on the same pymatgen figure, since the names
    # go inside the domains rather than into a legend box. Collapse the three to
    # a plain _legend so they do not write three identical files.
    pmg_suffix = suffix
    for tag in ('_legend_in', '_legend_out', '_legend_up'):
        pmg_suffix = pmg_suffix.replace(tag, '_legend')
    out_name = f'{png_name}_pymatgen{pmg_suffix}.pdf'
    plt.savefig(out_name, dpi=300, bbox_inches='tight')
    print(f"Pourbaix diagram saved as {out_name}")
    if args.show_fig:
        plt.show()
    plt.close()


def main(args, png_name, suffix, element_energies):
    tick = args.tick
    pHmin, pHmax = args.pHmin, args.pHmax
    pHrange = np.arange(pHmin, pHmax + tick, tick)
    Umin, Umax = args.Umin, args.Umax
    Urange = np.arange(Umin, Umax + 0.06 * 14, tick)
    target_pH = args.pH

    json_files = sorted(glob.glob(os.path.join(args.json_dir, "*.json")), key=os.path.basename)
    file_labels, file_gibbs_corrections, file_gc_params, file_oh_counts = load_labels(
        args, json_files,
    )

    if args.show_label:
        print("\nFile labels:")
        for fname, label in file_labels.items():
            print(f"{fname}: {label}")

    if args.png:
        print("\nSaving PNG images of structures...")
        for json_file in json_files:
            png_filename = os.path.basename(json_file).replace('.json', '.png')
            try:
                write(png_filename, read(json_file), rotation=args.png_rotation, show_unit_cell=0)
            except Exception as e:
                print(f"Failed to save {png_filename}: {e}")
        print("PNG image generation completed.\n")

    surfs, sorted_elements, remaining_elements, unique_elements, min_counts = build_surfs(
        json_files, file_labels, file_gibbs_corrections, file_gc_params, args,
    )

    if args.show_count:
        print("\nMinimum count of each element across all files:")
        for el in sorted_elements:
            print(f"{el}: {min_counts[el]}")

    if args.show_element:
        print("\nsorted_elements:", sorted_elements)
        print("remaining_elements:", remaining_elements)
        print("unique_elements:", unique_elements)

    ref_energies = resolve_element_energies(element_energies, unique_elements, args)
    ref_surf_idx, ref_surf = find_ref_surface(surfs, unique_elements, args)
    thermo_data = load_thermo_data(args)
    surface_formation_corrections(surfs, ref_surf, file_oh_counts, args, unique_elements, ref_energies)

    ions, solids, gases, liquids = [], [], [], []
    for el, el_ions, el_solids, el_gases, el_liquids in load_thermo_species(
        unique_elements, remaining_elements, thermo_data, ref_energies, args,
    ):
        ions.extend(el_ions)
        solids.extend(el_solids)
        gases.extend(el_gases)
        liquids.extend(el_liquids)

        if args.hybrid and not args.no_bulk:
            el_bulks = el_ions + el_solids + el_gases + el_liquids
            bulk_colors = plot_bulk_diagram(
                el, el_bulks,
                pHrange, Urange, pHmin, pHmax, Umin, Umax, args, suffix,
            )
            if not args.no_pymatgen and bulk_colors:
                plot_pymatgen_diagram(
                    el_bulks, bulk_colors, args, f'pourbaix_bulk_{el}', suffix,
                )

    save_surfs = surfs.copy()
    nsurfs = len(surfs)
    cols = surf_df_columns(remaining_elements, args.gc)

    print()
    print(format_df_for_display(pd.DataFrame(surfs, columns=cols)))
    print(f"Surfs: {nsurfs} entries\n")

    if args.hybrid:
        species_cols = ['ΔGf', 'ΔGh', 'e'] + remaining_elements + ['conc', 'name']
        for name, data in [('Ions', ions), ('Solids', solids), ('Gases', gases), ('Liquids', liquids)]:
            if len(data) > 0:
                print(format_df_for_display(pd.DataFrame(data, columns=species_cols)))
            print(f"{name}: {len(data)} entries\n")

    if args.hybrid:
        max_counts, new_surfs = build_hybrid_combinations(surfs, nsurfs, unique_elements, ions, solids, gases, liquids)
        surfs.extend(new_surfs)
        surfs = [surf for surf in surfs if all(surf[el] == max_counts[el] for el in unique_elements)]
        nsurfs = len(surfs)

        print(format_df_for_display(pd.DataFrame(surfs, columns=cols)))
        print(f"After adding combinations: {nsurfs} surfs entries\n")

    new_ref_surf, new_ref_surf_idx = find_new_ref_surface(surfs, ref_surf, unique_elements)
    ref_surf_idx = new_ref_surf_idx
    ref_surf = new_ref_surf
    lowest_surfaces = find_lowest_indices(surfs, pHrange, Urange, surfs[ref_surf_idx])

    if args.show_min_coord:
        min_coords = {}
        for j in range(len(pHrange)):
            for i in range(len(Urange)):
                sid = int(lowest_surfaces[i, j])
                x, y = pHrange[j], Urange[i]
                if sid not in min_coords or x < min_coords[sid][0] or (x == min_coords[sid][0] and y < min_coords[sid][1]):
                    min_coords[sid] = (x, y)
        for sid in sorted(min_coords):
            x, y = min_coords[sid]
            print(f"Surface {sid}: x = {x:.2f}, y = {y:.2f}, name = {surfs[int(sid)]['name']}")

    unique_ids = []
    for i in range(len(Urange)):
        for j in range(len(pHrange)):
            surf_id = int(lowest_surfaces[i, j])
            if surf_id not in unique_ids:
                unique_ids.append(surf_id)

    color_map = plot_2d_diagram(
        surfs, save_surfs, lowest_surfaces, pHrange, Urange,
        pHmin, pHmax, Umin, Umax, ref_surf_idx, args, png_name, suffix,
    )
    plot_1d_diagram(
        surfs, nsurfs, Urange, Umin, Umax, target_pH, surfs[ref_surf_idx],
        unique_ids, args, png_name, suffix, tick,
    )
    if not args.no_pymatgen:
        plot_pymatgen_diagram(surfs, color_map, args, png_name, suffix)

if __name__ == "__main__":
    cli_args = parse_args()
    gas_energies, element_energies = load_ref_energies(cli_args)
    init_thermo_constants(gas_energies)
    output_png_name, output_suffix = build_output_names(cli_args)
    main(cli_args, output_png_name, output_suffix, element_energies)
