python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 \
--cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 \
--colors-bulk white whitesmoke white

python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --hybrid --no-bulk \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 \
--cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 \
--colors-bulk white whitesmoke white --legend-out
