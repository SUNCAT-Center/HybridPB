"""Rebuild reference_energies.jsonc from the reference_energies/ structure database.

Database layout:

    reference_energies/<functional>/<formula>.json

Each file is an ASE-readable structure with an attached calculator (e.g. a VASP
`final_with_calculator.json`). The reference energy per element atom is

    mu_M = E_total/n_M - dGf/n_M - (n_O/n_M)*g_O - (n_H/n_M)*g_H

    g_O = g_H2O - g_H2 + dGf(H2O)   (O chemical potential vs. 1/2 O2)
    g_H = g_H2 / 2

where dGf is the experimental formation energy of the reduced compound taken
from thermodynamic_data.jsonc (kcal/mol), and the H2 / H2O energies come from
the same functional block of reference_energies.jsonc. A pure element has
dGf = 0 and no O/H, so it reduces to E_total/n_M.

Gibbs corrections of the compound are assumed to be zero (G_corr = 0).

Elements with no folder in the database keep their current values.
"""

import argparse
import glob
import json
import os
from collections import Counter
from math import gcd

from ase.io import read
from mendeleev import element
from pymatgen.core.ion import Ion

import HybridPourbaix as hpb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CALMOL = 23.061  # kcal/mol per eV


def atomic_number(symbol):
    """Atomic number, used to keep elements in periodic-table order."""
    return element(symbol).atomic_number


def formula_string(reduced, el):
    """Reduced composition as a formula string, element first, then O, then H."""
    order = [el] + [symbol for symbol in ('O', 'H') if symbol in reduced and symbol != el]
    return ''.join(f"{symbol}{reduced[symbol] if reduced[symbol] > 1 else ''}" for symbol in order)


def same_composition(formula, reduced):
    """True if `formula` describes the same reduced composition."""
    try:
        composition = Ion.from_formula(formula)
    except Exception:
        return False
    counts = {str(symbol): int(composition[str(symbol)]) for symbol in composition.elements}
    return reduce_counts(counts) == reduced


def reduce_counts(counts):
    """Reduce an element count dict to its smallest integer ratio."""
    divisor = 0
    for value in counts.values():
        divisor = gcd(divisor, int(value))
    divisor = divisor or 1
    return {el: int(value) // divisor for el, value in counts.items()}


def find_formation_energy(counts, el, thermo_data):
    """Look up dGf (kcal/mol per formula unit) of the reduced compound.

    Returns (dGf, formula key, metal atoms per formula unit).
    """
    reduced = reduce_counts(counts)
    solids = (thermo_data.get(el) or {}).get('solids', {})

    for formula, energy in solids.items():
        try:
            ion = Ion.from_formula(formula)
        except Exception:
            continue
        key_counts = {symbol: int(ion[symbol]) for symbol in reduced}
        if key_counts == reduced and int(ion.num_atoms) == sum(reduced.values()):
            if energy is None:
                raise LookupError(
                    f"'{formula}' is still a null placeholder in "
                    f"thermodynamic_data.jsonc['{el}']['solids']"
                )
            return float(energy), formula, reduced[el]

    pretty = ''.join(f"{sym}{n if n > 1 else ''}" for sym, n in sorted(reduced.items()))
    raise LookupError(
        f"no solid matching {pretty} in thermodynamic_data.jsonc['{el}']['solids'] "
        f"(available: {', '.join(solids) or 'none'})"
    )


def read_ueff(path, el):
    """Effective Hubbard U (U - J) applied to `el`, or None if +U was not used.

    ldau_luj usually lists every element the workflow knows about; an entry with
    L < 0 means no +U was actually applied to that element.
    """
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return None

    for entry in data.values():
        if not isinstance(entry, dict):
            continue
        params = entry.get('calculator_parameters') or {}
        if not params.get('ldau'):
            continue
        luj = (params.get('ldau_luj') or {}).get(el)
        if not luj or luj.get('L', -1) < 0:
            continue
        ueff = float(luj.get('U', 0.0)) - float(luj.get('J', 0.0))
        if ueff:
            return ueff
    return None


def structure_reference_energy(path, el, thermo_data, g_o):
    """Reference energy per `el` atom for one database structure."""
    atoms = read(path)
    counts = Counter(atoms.get_chemical_symbols())

    n_el = counts.get(el, 0)
    if n_el == 0:
        raise SystemExit(f"{path}: contains no '{el}' atoms")

    extra = set(counts) - {el, 'O', 'H'}
    if extra:
        raise SystemExit(
            f"{path}: contains {', '.join(sorted(extra))} besides {el}/O/H; "
            "only O and H can be referenced away"
        )

    n_o, n_h = counts.get('O', 0), counts.get('H', 0)
    if n_o or n_h:
        # G_corr of the compound is assumed to be zero
        dgf_kcal, formula, el_per_unit = find_formation_energy(counts, el, thermo_data)
        dgf = dgf_kcal / el_per_unit / CALMOL
    else:
        dgf, formula = 0.0, atoms.get_chemical_formula()

    mu = atoms.get_potential_energy() / n_el - dgf - (n_o / n_el) * g_o - (n_h / n_el) * hpb.gh

    cell_formula = atoms.get_chemical_formula()
    source = cell_formula if same_composition(formula, reduce_counts(counts)) \
        else f"{cell_formula} via {formula}"
    return mu, source


def scan_database(db_dir, functional, thermo_data, g_o):
    """Return ({element: {formula: (energy, source)}}, [skipped], {element: {formula: Ueff}}).

    Every *.json in the functional folder is one reference compound; the element
    it references is the single constituent that is not O or H.
    """
    func_dir = os.path.join(db_dir, functional)
    entries = {}
    skipped = []
    ueffs = {}

    for path in sorted(glob.glob(os.path.join(func_dir, '*.json'))):
        name = os.path.basename(path)
        model = os.path.splitext(name)[0]
        counts = Counter(read(path).get_chemical_symbols())

        candidates = sorted(set(counts) - {'O', 'H'})
        if not candidates:
            print(f"  SKIPPED {name}: only O/H; gas references belong in the JSONC file")
            continue
        if len(candidates) > 1:
            print(f"  SKIPPED {name}: {', '.join(candidates)} are all non-O/H, "
                  "so the referenced element is ambiguous")
            continue
        el = candidates[0]

        reduced = reduce_counts(counts)
        if not same_composition(model, reduced):
            print(f"  WARNING {name}: filename does not match its composition; "
                  f"expected {formula_string(reduced, el)}.json")

        try:
            entries.setdefault(el, {})[model] = structure_reference_energy(
                path, el, thermo_data, g_o,
            )
        except LookupError as exc:
            skipped.append(f"{model}: {exc}")
            continue

        ueff = read_ueff(path, el)
        if ueff is not None:
            ueffs.setdefault(el, {})[model] = ueff

    return {el: models for el, models in entries.items() if models}, skipped, ueffs


def merge_elements(existing, derived):
    """Merge derived values over the existing element block, keeping the rest."""
    merged = {}
    provenance = {}

    for el, models in derived.items():
        merged[el] = {model: energy for model, (energy, _) in models.items()}
        provenance[el] = ', '.join(f"{model}: {formula}" for model, (_, formula) in models.items())

    for el, value in existing.items():
        if el not in merged:
            merged[el] = value
            provenance[el] = 'not in database, kept as-is'

    order = sorted(merged, key=atomic_number)
    return {el: merged[el] for el in order}, provenance


def format_number(value):
    """Print a float without scientific notation or trailing-zero noise."""
    return f"{value:.8f}".rstrip('0').rstrip('.')


def ueff_comment(el_ueffs):
    """`// Fe2O3: Ueff = 4.30 eV` note for references calculated with +U."""
    if not el_ueffs:
        return ''
    notes = ', '.join(
        f"{model}: Ueff = {ueff:g}{'' if ueff % 1 else '.0'} eV"
        for model, ueff in el_ueffs.items()
    )
    return f'  // {notes}'


def render_elements(elements, indent, ueffs):
    """Render the elements block, one element per line."""
    lines = []
    items = list(elements.items())
    for i, (el, value) in enumerate(items):
        comma = ',' if i < len(items) - 1 else ''
        note = ueff_comment(ueffs.get(el, {}))
        if isinstance(value, dict):
            inner = ', '.join(f'"{model}": {format_number(energy)}' for model, energy in value.items())
            lines.append(f'{indent}"{el}": {{ {inner} }}{comma}{note}')
        else:
            lines.append(f'{indent}"{el}": {format_number(float(value))}{comma}{note}')
    return lines


def render_file(data, ueffs_by_functional):
    """Render the whole reference_energies.jsonc text."""
    lines = ['{']
    functionals = list(data)

    for f_i, functional in enumerate(functionals):
        block = data[functional] or {}
        gases = block.get('gases') or {}
        elements = block.get('elements') or {}
        f_comma = ',' if f_i < len(functionals) - 1 else ''

        lines.append(f'  "{functional}": {{')

        gas_items = list(gases.items())
        if gas_items:
            lines.append('    "gases": {')
            for g_i, (name, energy) in enumerate(gas_items):
                g_comma = ',' if g_i < len(gas_items) - 1 else ''
                lines.append(f'      "{name}": {format_number(float(energy))}{g_comma}')
            lines.append('    },')
        else:
            lines.append('    "gases": {},')

        if elements:
            lines.append('    "elements": {')
            lines.extend(render_elements(elements, '      ',
                                         ueffs_by_functional.get(functional, {})))
            lines.append('    }')
        else:
            lines.append('    "elements": {}')

        lines.append(f'  }}{f_comma}')

    lines.append('}')
    return '\n'.join(lines) + '\n'


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db', type=str, default=os.path.join(SCRIPT_DIR, 'reference_energies'),
                        help='Structure database directory (default: reference_energies/)')
    parser.add_argument('--out', type=str, default=os.path.join(SCRIPT_DIR, 'reference_energies.jsonc'),
                        help='JSONC file to rebuild (default: reference_energies.jsonc)')
    parser.add_argument('--thermo-data', type=str, default=os.path.join(SCRIPT_DIR, 'thermodynamic_data.jsonc'),
                        help='Formation energies used for compound references '
                             '(default: thermodynamic_data.jsonc)')
    parser.add_argument('--functional', type=str, nargs='+',
                        help='Functionals to rebuild (default: every subfolder of --db)')
    parser.add_argument('--show-source', action='store_true',
                        help='Print which structure each reference energy came from')
    parser.add_argument('--dry-run', action='store_true', help='Print the result without writing')
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.out):
        raise SystemExit(f"{args.out} not found; it supplies the H2/H2O energies")
    data = hpb.load_jsonc(args.out)

    ueffs_by_functional = {}

    if not os.path.exists(args.thermo_data):
        raise SystemExit(f"{args.thermo_data} not found; it supplies the formation energies")
    thermo_data = hpb.load_jsonc(args.thermo_data)

    if args.functional:
        functionals = args.functional
    else:
        functionals = sorted(
            os.path.basename(d.rstrip('/')) for d in glob.glob(os.path.join(args.db, '*/'))
        )
    if not functionals:
        raise SystemExit(f"No functional folders found in {args.db}")

    for functional in functionals:
        block = data.get(functional) or {}
        gases = block.get('gases') or {}
        missing = [gas for gas in ('H2', 'H2O') if gas not in gases]
        if missing:
            raise SystemExit(
                f"Functional '{functional}' in {args.out} is missing gas energies: "
                f"{', '.join(missing)}. Add them before rebuilding."
            )

        hpb.init_thermo_constants(gases)
        # hpb.go references O to water; add the water formation energy to
        # reference it to 1/2 O2, which is what the tabulated dGf assume
        g_o = hpb.go + hpb.water_formation_energy

        derived, skipped, ueffs = scan_database(args.db, functional, thermo_data, g_o)
        ueffs_by_functional[functional] = ueffs
        if not derived:
            print(f"{functional}: no structures found in {os.path.join(args.db, functional)}, skipped")
            continue

        merged, provenance = merge_elements(block.get('elements') or {}, derived)
        data[functional] = {'gases': gases, 'elements': merged}

        n_models = sum(len(models) for models in derived.values())
        print(f"{functional}: {len(derived)} elements, {n_models} structures "
              f"(g_O = {g_o:.6f} eV, g_H = {hpb.gh:.6f} eV)")
        for message in skipped:
            print(f"  SKIPPED {message}")
        if args.show_source:
            for el in sorted(derived, key=atomic_number):
                print(f"  {el}: {provenance[el]}")

    text = render_file(data, ueffs_by_functional)

    if args.dry_run:
        print()
        print(text)
        return

    with open(args.out, 'w') as f:
        f.write(text)
    print(f"\nWrote {args.out}")


if __name__ == '__main__':
    main()
