#!/usr/bin/env python3
"""Export converged cowork slabs into this example directory as label.csv's json files.

Source: /pscratch/sd/j/jiuy97/7_prediction/cowork/{5_CoRuO2, 6_CoRuO2+Ov}.  This is the
4-TRILAYER RuO2(110) tree rebuilt in September 2026, NOT the earlier 5-layer `_alloy`
model -- cowork/_alloy and cowork/_vac are the old tree and are not read.

NAMING.  label.csv uses ``{Co-cus}-{Ru-cus}`` with `x` for a bare cus, an `H-` prefix when
both bridging oxygens are protonated, and a `v-` prefix for the cus-row oxygen vacancy.
The order of the two slots is not cosmetic: the two cus metals of a doped cell
are inequivalent, one being the dopant and the other host Ru, so `O-OH` (O on Co)
and `OH-O` (O on Ru) are different surfaces.

TILT VARIANTS.  4_OH_OH/{1_cus,2_brg} and each pair inside 5_O_OH are the SAME surface
started from two guesses for where a dangling proton points.  Only the lower-energy one
belongs in a Pourbaix diagram, so they collapse to a single json here.  Site variants
(1_dopant vs 2_host) do NOT collapse -- those are different surfaces.

1_x_x/brg IS NOT A POURBAIX ENTRY.  It is exported as `brg-x-x.json` / `v-brg-x-x.json`
and deliberately left OUT of label.csv.  Here `1_x_x/brg` puts the Co DOPANT on a bridging site instead of cus.  Every adsorbate state in this
directory was computed with the dopant on cus, so letting the brg isomer win the `clean` slot
would build a diagram in which that changes as soon as anything adsorbs -- implying a
migration that is not a surface reaction.  What the pair IS good for is the site
preference itself, which the run prints as a dE line.

CONVERGENCE.  By default only runs with a DONE file are exported.  A Pourbaix diagram is
built from total energies, and an unconverged energy is not comparable with a converged
one; a half-relaxed slab would silently move a phase boundary.  Everything else is
reported as skipped so the gaps are visible rather than assumed absent.

--provisional falls back to the LAST IONIC STEP of an unconverged run (`OUTCAR@-1`) for
anything without a DONE, so a rough diagram can be drawn before the campaign finishes.
Those energies are NOT converged and the phase boundaries drawn from them will move; the
run log marks every such row with [PROVISIONAL] and prints its max force so it is obvious
which ones are placeholders.  Do not publish a diagram built this way.

Usage:
    python export_from_cowork.py                  # converged only
    python export_from_cowork.py --provisional    # + last ionic step of unfinished runs
    python export_from_cowork.py --check          # report only, write nothing
"""
import argparse
import os

import numpy as np
from ase.io import read, write

COWORK = "/pscratch/sd/j/jiuy97/7_prediction/cowork"
DEST = os.path.dirname(os.path.abspath(__file__))
STRUCT = "final_with_calculator.json"      # the run's own output, carries the energy

# target json stem -> the cowork state(s) that produce it.  A list of more than one means
# tilt variants: same surface, lowest energy wins.
STATES = {
    "x-x":       ["1_x_x"],
    "OH-x":      ["2_x_OH/1_dopant"],
    "x-OH":      ["2_x_OH/2_host"],
    "O-x":       ["3_x_O/1_dopant"],
    "x-O":       ["3_x_O/2_host"],
    "OH-OH":     ["4_OH_OH/1_cus", "4_OH_OH/2_brg"],
    "O-OH":      ["5_O_OH/1_O_OH_cus", "5_O_OH/2_O_OH_brg"],
    "OH-O":      ["5_O_OH/3_OH_O_cus", "5_O_OH/4_OH_O_brg"],
    "O-O":       ["6_O_O"],
    "H-x-x":     ["7_H_x_x"],
    "H-OH-OH":   ["8_H_OH_OH"],
    "H-O-O":     ["9_H_O_O"],
    "H2O-H2O":   ["0_H2O_H2O"],
}

# Exported for reference but kept OUT of the Pourbaix diagram -- see the module docstring.
# HybridPourbaix.py globs *.json in its --json-dir and uses label.csv only to attach
# labels, so leaving a file out of label.csv does NOT leave it out of the diagram: it
# enters unlabelled.  Keeping these in a subdirectory is what actually excludes them.
EXTRA_DIR = "brg_site"
EXTRA = {
    "brg-x-x": ["1_x_x/brg"],
}

SYSTEMS = [("", "5_CoRuO2"), ("v-", "6_CoRuO2+Ov")]

# One deliberate exception to the naming: the bare pristine surface is `clean.json`,
# not `x-x.json`.  The vacancy one stays `v-x-x.json`.
RENAME = {"x-x": "clean"}


# ---- oxygen-vacancy integrity ------------------------------------------------
# A +Ov run can relax back to a FILLED lattice: an adsorbed cus O drops into the
# vacancy and heals it.  The cell composition never changes, so nothing downstream
# notices -- the row simply stops being the surface its label claims, and it then
# competes in the diagram against the real v- states at the wrong structure.  It is
# the O-bearing states that do it (3_x_O, 6_O_O, 9_H_O_O, the O half of 5_O_OH),
# which is chemically what you would expect: the cus oxygen heals the vacancy.
#
# Every cowork structure is built on the same frozen bottom bilayer -- the fixed
# coordinates agree to 0.000000 A across the whole tree -- so the vacancy lattice
# site can be compared directly, with no alignment step.
PRISTINE = os.path.join(COWORK, "1_RuO2", "1_x_x", "start.traj")
VAC_IDX = {"": 22, "brg": 21}   # which pristine O the vacancy removes:
                                # 22 = in-plane cus-row O (z 22.32); 21 = bridging O (z 23.58)
VAC_TOL = 1.0                   # A.  Refilled sites land within 0.41 of the site; the
                                # nearest intact one is 1.39, an adsorbate leaning in from
                                # a full A higher in z.  The gap is wide and unambiguous.

_pristine = read(PRISTINE)
_cell = np.array(_pristine.get_cell())


def vacancy_filled(atoms, which):
    """True if an O has moved into the vacancy site -- i.e. this is no longer a +Ov surface."""
    site = _pristine.positions[VAC_IDX[which]]
    sym = atoms.get_chemical_symbols()
    O = np.array([atoms.positions[i] for i in range(len(atoms)) if sym[i] == "O"])
    f = np.linalg.solve(_cell.T, (O - site).T).T
    f -= np.round(f)
    return float(np.linalg.norm(f @ _cell, axis=1).min()) < VAC_TOL


def converged(d):
    """(atoms, energy) of a finished run, or None."""
    if not os.path.isfile(os.path.join(d, "DONE")):
        return None
    p = os.path.join(d, STRUCT)
    if not os.path.isfile(p):
        return None
    try:
        a = read(p)
        return a, a.get_potential_energy()
    except Exception:
        return None


def last_step(d):
    """(atoms, energy, fmax) of the last completed ionic step of an unfinished run.

    Read from OUTCAR rather than restart.json: restart.json carries the geometry a
    resubmission would START from, which after get_restart3 is the same thing but before
    it is the previous seed.  OUTCAR@-1 is unambiguous -- the furthest the run actually got.
    """
    p = os.path.join(d, "OUTCAR")
    if not os.path.isfile(p):
        return None
    try:
        a = read(p, index=-1)
        return a, a.get_potential_energy(), float(np.linalg.norm(a.get_forces(), axis=1).max())
    except Exception:
        return None


def pick(system, sources, provisional, vac=None):
    """(energy, source, atoms, tag) of the state to export, or None.

    `vac` is the vacancy kind to police ("" for the cus-row vacancy, "brg" for the
    bridging one) or None for a pristine system.  A candidate whose vacancy has been
    healed is dropped here rather than exported, and is reported by the caller.
    """
    cands = []
    for s in sources:
        got = converged(os.path.join(COWORK, system, s))
        if got is not None:
            cands.append((got[1], s, got[0], ""))
    if not cands and provisional:
        for s in sources:
            got = last_step(os.path.join(COWORK, system, s))
            if got is not None:
                cands.append((got[1], s, got[0], f"  [PROVISIONAL, fmax {got[2]:.3f}]"))
    if vac is not None:
        kept = [c for c in cands if not vacancy_filled(c[2], vac)]
        if len(kept) != len(cands):
            rejected.extend(f"{system}/{c[1]}" for c in cands if vacancy_filled(c[2], vac))
        cands = kept
    if not cands:
        return None
    best = min(cands, key=lambda c: c[0])
    if len(cands) > 1:
        best = best[:3] + (best[3] + f"  [{len(cands)}개 중 최저]",)
    return best


ap = argparse.ArgumentParser()
ap.add_argument("--check", action="store_true", help="report only, write nothing")
ap.add_argument("--provisional", action="store_true",
                help="also export OUTCAR@-1 of unfinished runs (NOT converged)")
args = ap.parse_args()

n_out = n_skip = n_prov = 0
rejected = []
healed = []
site_pref = []
for prefix, system in SYSTEMS:
    print(f"=== {system}")
    for table, subdir, tag_extra in ((STATES, "", ""),
                                     (EXTRA, EXTRA_DIR, "   ← 진단용, 다이어그램 제외")):
        for stem, sources in table.items():
            name = RENAME[stem] if prefix == "" and stem in RENAME else stem
            target = os.path.join(subdir, f"{prefix}{name}.json")
            vac = None if prefix == "" else ("brg" if table is EXTRA else "")
            got = pick(system, sources, args.provisional, vac)
            if got is None:
                # A source whose directory carries the _filled suffix relaxed its oxygen
                # vacancy away: cowork itself records the verdict, so say so rather than
                # reporting the state as merely absent.  vacancy_filled() above still
                # guards runs that have converged since the directory was last renamed.
                why = []
                for s in sources:
                    if os.path.isdir(os.path.join(COWORK, system, s + "_filled")):
                        why.append(f"{s}_filled")
                        healed.append(f"{system}/{s}")
                note = f"  빈자리 치유: {', '.join(why)}" if why else ""
                print(f"  {target:<16} -- 없음  ({', '.join(sources)}){note}")
                n_skip += 1
                continue
            e, s, atoms, tag = got
            if not args.check:
                out = os.path.join(DEST, target)
                os.makedirs(os.path.dirname(out), exist_ok=True)
                write(out, atoms)
            print(f"  {target:<16} <- {s:<22} E = {e:>12.4f} eV{tag}{tag_extra}")
            n_out += 1
            if "PROVISIONAL" in tag:
                n_prov += 1

    # site preference: the same cell with the dopant on cus vs on brg
    cus = pick(system, ["1_x_x"], args.provisional)
    brg = pick(system, ["1_x_x/brg"], args.provisional)
    if cus and brg:
        site_pref.append((system, cus[0], brg[0]))

if site_pref:
    print("\n=== 도판트 자리 선호 (1_x_x vs 1_x_x/brg)")
    for system, ecus, ebrg in site_pref:
        d = ebrg - ecus
        who = "cus" if d > 0 else "brg"
        print(f"  {system:<14} cus {ecus:>12.4f}   brg {ebrg:>12.4f}   "
              f"dE = {d:+.3f} eV  → {who} 쪽이 안정")

if rejected or healed:
    both = sorted(set(rejected) | set(healed))
    print(f"\n=== 빈자리가 다시 채워져 제외한 계산 {len(both)}개")
    for r in both:
        print(f"  {r}")
    print("  (흡착 O가 빈자리로 내려앉아 +Ov 표면이 아니게 된 것들;")
    print("   cowork 쪽 폴더에는 _filled 접미사가 붙어 있습니다)")

print(f"\n{'검사만' if args.check else '내보냄'} {n_out}개 "
      f"(그중 잠정 {n_prov}개), 없음 {n_skip}개  (총 {n_out + n_skip})")
if n_prov:
    print("⚠ [PROVISIONAL] 행은 수렴하지 않은 마지막 이온 스텝입니다 — 상경계가 움직입니다.")
if n_skip or n_prov:
    print("계산이 더 끝나면 다시 실행하면 됩니다 — 기존 파일은 덮어씁니다.")
