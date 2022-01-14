
# MINING-D

MINING-D is a tool (written in Python 3.6.5) for inference of D genes using Rep-Seq data from diverse species. 
It takes as input a fasta file with consensus CDR3s from an immunosequencing dataset (see paper) and writes the inferred D genes in the output file in the fasta format. 

There are two ways of running MINING-D. 

- **From command line** - To use default MINING-D parameters, use 

    
    ``` 
    $ python MINING_D.py -i <input_CDR_file> -o <output_file> 
    ``` 
    
    To check all available options, use `python MINING-D.py --help`. 

    ```
    $ python MINING_D.py --help 
    Usage: MINING_D.py [OPTIONS]

    Options:
    -i, --input_cdr_file TEXT  Input file with consensus CDR3s  [required]
    -o, --output_file TEXT     Output file to store inferred D genes  [required]
    -k INTEGER                 Length of starting k-mers  [default: 10]
    -n, --num_k_mers INTEGER   Number of most abundant k-mers to extend [default: 300]
    -p, --p_val_th FLOAT       p-value threshold for the extension procedure [default: 4.5e-36]
    -t, --n_cores INTEGER      Number of cores  [default: 10]
    --bidir_alpha FLOAT        Alpha for filtering bidirectional extensions (see paper)  [default: 0.5]
    -g, --cliq_th INTEGER      Similarity metric threshold for generating graphs [default: 2]
    --help                     Show this message and exit.
    ```
		
- **Jupyter Notebook** - If you would like to run MINING-D in an interactive environment, there is a Jupyter notebook named "*interactive\_MINING\_D.ipynb*" in the repo. Please specify the input file, output file, and other parameters at the top of the notebook like shown in the example. 


**Dependencies**

- Biopython 
- Networkx<2.4
- Joblib
- NumPy
- SciPy
- [click](https://palletsprojects.com/p/click/)

