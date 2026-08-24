# surface
mkdir -p vac/
mv vac*.json vac/
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 \
--colors-2d '#73ADD1' '#B9E0ED' '#F0F9DB' '#FEE79B' '#FDAD60' --Gmin -3 --Gmax 5

# surface (gc)
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 --gc \
--colors-2d '#73ADD1' '#B9E0ED' royalblue '#FEE79B' '#FDAD60' --Gmin -3 --Gmax 5

# hybrid (bulk)
mv vac/vac*.json .
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --HER --OER --Umin -1.0 --Umax 2.5 --hybrid --png \
--colors-bulk white whitesmoke gray gainsboro silver white dimgray --Gmin -3 --Gmax 5 --ref-json clean-is.json

# hybrid
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 --hybrid --no-bulk --suffix Fe2O3 \
--colors-bulk white gainsboro darkgray dimgray white gainsboro dimgray gray \
--colors-2d '#73ADD1' '#B9E0ED' '#F0F9DB' --ref-json clean-is.json --Gmin -3 --Gmax 5

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 --hybrid --no-bulk --thermo-data ./thermodynamic_data.jsonc \
--colors-bulk white gainsboro gainsboro white darkgray dimgray gray \
--colors-2d '#73ADD1' '#B9E0ED' '#F0F9DB' --ref-json clean-is.json --Gmin -3 --Gmax 5

# gc + hybrid
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 --hybrid --no-bulk --gc --suffix Fe2O3 \
--colors-bulk white gainsboro darkgray gainsboro dimgray white \
--colors-2d '#73ADD1' royalblue '#B9E0ED' '#F0F9DB' --ref-json clean-is.json --Gmin -3 --Gmax 5

python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --Umin -1.5 --Umax 2.0 --hybrid --no-bulk --gc --thermo-data ./thermodynamic_data.jsonc \
--colors-bulk white gainsboro gainsboro gray dimgray white \
--colors-2d '#73ADD1' '#B9E0ED' --ref-json clean-is.json --Gmin -3 --Gmax 5
