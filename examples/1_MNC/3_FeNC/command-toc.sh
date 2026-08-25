mkdir -p vac/
mv vac*.json vac/
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --fill 0.72 0.92 --Umin -1.0 --Umax 2.0 --figx 3 --figy 3 \
--colors-2d '#73ADD1' '#B9E0ED' '#F0F9DB' '#FEE79B' '#FDAD60' --suffix toc

mv vac/vac*.json .
python ../../../HybridPourbaix.py --ref-energies ./reference_energies.jsonc --fill 0.72 0.92 --Umin -1.0 --Umax 2.0 --figx 3 --figy 3 \
--hybrid --no-bulk --gc --thermo-data ./thermodynamic_data.jsonc \
--colors-bulk white gainsboro gainsboro gray dimgray white \
--colors-2d '#73ADD1' '#B9E0ED' --suffix toc