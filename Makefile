configure:
	julia --project="." -e 'using Pkg; Pkg.instantiate()'

prepare:
	mkdir -p data/yamanishi2008/{DF,DD}
	sh src/prepare/prepare_maccs.sh data/yamanishi2008/SMILES/nr.smi
	mv data/yamanishi2008/SMILES/nr.maccs.txt data/yamanishi2008/DF/nr.maccs.txt
	julia --project="." src/prepare/TverskyMatrix.jl -i data/yamanishi2008/DF/nr.maccs.txt -o data/yamanishi2008/DD/nr.maccs_tanimoto.txt --named

test:
	mkdir -p data/example/DF
	sh src/prepare/prepare_maccs.sh data/example/SMILES/example.smi
	mv data/example/SMILES/example.maccs.txt data/example/DF/example.maccs.txt
	mkdir -p results/example/data/DD
	julia --project="." src/prepare/TverskyMatrix.jl --MD data/example/DF/example.maccs.txt --ND data/yamanishi2008/DF/nr.maccs.txt -o results/example/data/DD/example_nr.maccs_tanimoto.txt --named
	julia --project="." src/predict/SimSpread.jl --dt-train data/yamanishi2008/DT/nr_DT.txt --dd-train data/yamanishi2008/DD/nr.maccs_tanimoto.txt --dd-query results/example/data/DD/example_nr.maccs_tanimoto.txt -o results/example/preds.txt -c 0.1

clean:
	rm -rf data/yamanishi2008/{DF,DD}
	rm -rf data/example/DF
	rm -rf results

