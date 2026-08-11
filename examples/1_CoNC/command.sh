# surface (without vac*.json files)
mkdir -p vac/
mv vac*.json vac/
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --OER --HER --png \
--colors-2d '#73ADD1' '#B9E0ED' '#FEE79B' '#FDAD60'

# hybrid (conc)
mv vac/vac*.json .
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --OER --HER --hybrid --concentration 1e-4 --suffix conc4 \
--colors-bulk white gray lightgray darkgray white \
--colors-2d '#B9E0ED' '#FEE79B' '#FDAD60'
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --OER --HER --hybrid --concentration 1e-6 --suffix conc6 \
--colors-bulk white lightgray gray lightgray darkgray white \
--colors-2d '#B9E0ED' '#FEE79B' '#FDAD60'
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-out --OER --HER --hybrid --concentration 1e-8 --suffix conc8 \
--colors-bulk white lightgray gray lightgray darkgray white \
--colors-2d '#B9E0ED' '#FEE79B' '#FDAD60'

# hybrid (pH)
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-in --OER --HER --hybrid --no-bulk --Gmin -4 --Gmax 12 --ref-json clean.json --pH 0
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-in --OER --HER --hybrid --no-bulk --Gmin -4 --Gmax 12 --ref-json clean.json --pH 7
python ../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --legend-in --OER --HER --hybrid --no-bulk --Gmin -4 --Gmax 12 --ref-json clean.json --pH 14 \
--colors-bulk white lightgray gray lightgray darkgray white \
--colors-2d '#B9E0ED' '#FEE79B' '#FDAD60'