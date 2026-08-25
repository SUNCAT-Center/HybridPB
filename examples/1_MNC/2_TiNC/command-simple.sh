# surface (OxHy adsorbates)
mkdir -p no3rr/
mv *.json no3rr/
mv no3rr/h*.json no3rr/oh*.json no3rr/o.json .
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --line -0.7 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--colors-2d '#5c90c2' '#f0f9db' '#f88c51' --suffix OxHy

# surface (NOxHy adsorbates)
mv no3rr/* .
mv h*.json oh*.json o.json no3rr/
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --line -0.7 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap-2d PiYG_r --cmin-2d 0.15 --cmax-2d 0.85 --cgap-2d 0.0 --suffix NOxHy

# hybrid (bulk)
mv no3rr/* .
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --OER --HER --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --hybrid --png \
--colors-bulk white darkgray gainsboro white darkgray

# hybrid
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --line -0.7 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --hybrid --no-bulk \
--cmap RdYlBu_r --cmin 0.15 --cmax 0.85 --cgap 0.0 --colors-2d '#D14895'