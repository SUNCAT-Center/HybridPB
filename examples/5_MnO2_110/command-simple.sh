# hybrid
mkdir -p OH/
mv OH/OHBridge_OHTop.json .
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --OER --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray gold whitesmoke white lightgray lemonchiffon \
--colors-2d dodgerblue lightskyblue --suffix small

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray gold whitesmoke white lightgray lemonchiffon \
--colors-2d dodgerblue lightskyblue --figx 6 --figy 6 --suffix large

# hybrid without MnO4-
mv OH/OHBridge_OHTop.json .
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray gold white lemonchiffon whitesmoke \
--colors-2d dodgerblue lightskyblue \
--thermo-data ./thermodynamic_data.jsonc --figx 6 --figy 6 --suffix MnO4

# hybrid with OHBridge_OHTop.json
mv OHBridge_OHTop.json OH/
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray gray white orange whitesmoke lightgray darkgray gold \
--colors-2d dodgerblue lightskyblue --figx 6 --figy 6 --suffix OH
mv OH/OHBridge_OHTop.json .
