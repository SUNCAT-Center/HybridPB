# hybrid (bulk)
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --HER --OER \
--Umin -0.5 --Umax 2.0 --suffix K --colors-bulk white gainsboro

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --HER --OER \
--Umin -0.5 --Umax 2.0 --suffix Mn --colors-bulk white gainsboro silver gray dimgray silver gainsboro white

# hybrid
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --OER --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray darkgray whitesmoke gold white lemonchiffon \
--colors-2d dodgerblue lightskyblue --suffix small

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray darkgray whitesmoke gold white lemonchiffon \
--colors-2d dodgerblue lightskyblue --figx 6 --figy 6 --suffix large

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin 0.0 --Umax 2.0 --figx 2.5 --figy 2.5 \
--colors-bulk white whitesmoke lightgray darkgray lightgray gold white \
--colors-2d dodgerblue lightskyblue --suffix toc

# hybrid without MnO4-
mv OH/OHBridge_OHTop.json .
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray gold '#FFBE02' orange lemonchiffon \
--colors-2d dodgerblue lightskyblue \
--thermo-data ./thermodynamic_data.jsonc --figx 6 --figy 6 --suffix MnO4

# hybrid with OHBridge_OHTop.json
mv OHBridge_OHTop.json OH/
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk --Umin -0.5 --Umax 2.0 \
--colors-bulk white whitesmoke lightgray darkgray gray darkgray whitesmoke lightgray white lemonchiffon \
--colors-2d dodgerblue lightskyblue --figx 6 --figy 6 --suffix OH
mv OH/OHBridge_OHTop.json .
