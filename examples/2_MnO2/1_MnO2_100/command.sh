# hybrid
mkdir -p OH/
mv OH/OHBridge_OHTop.json .
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray darkgray whitesmoke gold white lemonchiffon \
--colors-2d dodgerblue lightskyblue

# hybrid without MnO4-
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray white gold '#FFBE02' orange lemonchiffon \
--colors-2d dodgerblue lightskyblue \
--thermo-data ./thermodynamic_data.jsonc --suffix MnO4

# hybrid with OHBridge_OHTop.json
mv OHBridge_OHTop.json OH/
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray darkgray whitesmoke lightgray white lemonchiffon \
--colors-2d dodgerblue lightskyblue --suffix OH