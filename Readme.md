#MINING-D

MINING-D is a tool (written in Python 3.6.5) for inference of D genes using Rep-Seq data from diverse species. 

Usage with default parameters

~~~ shell
$ python MINING_D.py -i input_CDR_file -o output_file
~~~

**Dependencies**

- Biopython 
- Networkx
- Joblib
- NumPy
- SciPy


To check all options for MINING-D, use
 
~~~ shell
$ python MINING_D.py --help
~~~

You should see something like this

~~~ shell
Usage: MINING_D.py [OPTIONS]

Options:
  -i, --input_cdr_file TEXT  Input file with consensus CDR3s  [required]
  -o, --output_file TEXT     Output file to store inferred D genes  [required]
  -k INTEGER                 Length of starting k-mers  [default: 10]
  -n, --num_k_mers INTEGER   Number of most abundant k-mers to extend
                             [default: 300]
  -p, --p_val_th FLOAT       p-value threshold for the extension procedure
                             [default: 4.5e-36]
  -t, --n_cores INTEGER      Number of cores  [default: 10]
  --bidir_alpha FLOAT        Alpha for filtering bidirectional extensions (see
                             paper)  [default: 0.5]
  -g, --cliq_th INTEGER      Similarity metric threshold for generating graphs
                             [default: 2]
  --help                     Show this message and exit.
~~~
