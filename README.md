# P-SOET

This program runs the Projected Site-Occupation Embedding Theory on the 1D Hubbard model:

https://arxiv.org/abs/1902.05747

It is a projected version (using Schmidt decomposition) of Site-Occupation Embedding Theory, which articles can be find at:

https://arxiv.org/abs/1409.2326

https://arxiv.org/abs/1602.02547

https://arxiv.org/abs/1710.03125

https://arxiv.org/abs/1806.06057


# Requirements

- Python 2.7 or higher
- LAPACK and BLAS libraries
- Intel fortran compiler (or make the appropriate change in the Makefile)
- For multiple impurity sites : 
Density Matrix Renormalization Group (Block-1.5 code https://sanshar.github.io/Block/overview.html)
and its own requirements.

# Installation

Clone the complete repository:
```
$ git clone https://github.com/bsenjean/P-SOET.git
```

Build the code:
```
$ cd code/
$ vi Makefile
```

Change the Makefile according to your preferences. Change compiler, add FLAGS, and change the path to the librairies.

```
$ make
```

```
$ cd ../work_dir/
$ vi script.py
```

In ```script.py``` you will see two environment variable to set, ```PSOET_DIR``` and ```DMRG_code```. Set those variables in your PATH before running the script:
```
$ export PSOET_DIR=/path/to/P-SOET/
$ export DMRG_code=/path/to/block-1.5.3/
```
In HPC Strasbourg:
```
$ export DMRG_code=/usr/local/quant/block-1.5.3/
```

# Use

```
$ python script.py
```

In ```script.py```, you have to specify the following:

- Number of sites
- Number of electrons
- On-site Coulomb repulsion U
- Nearest neighbor hopping integral t
- Number of impurity sites
- Approximation (either iBALDA, 2LBALDA or SIAMBALDA. Refer to https://arxiv.org/abs/1806.06057 for their meaning)
- The code_directory

You can also specify optionally:
- The number of renormalized states to keep in DMRG
- The number of maximum iterations
- The threshold for convergence
- The use of a global chemical potential (for now, only for a single impurity problem)

# How to cite

Cite the following paper:

- B. Senjean. Projected site-occupation embedding theory. arXiv:1902.05747 

We would be grateful if you could cite the following articles if you found this program useful:

- B. Senjean, M. Tsuchiizu, V. Robert, E. Fromager. Local density approximation in site occupation
embedding theory. Mol. Phys., 115, 48 (2017).

- B. Senjean, N. Nakatani, M. Tsuchiizu, E. Fromager. Site-Occupation Embedding Theory using
Bethe Ansatz Local Density Approximations. Phys. Rev. B, 97, 235105 (2018).

- B. Senjean, N. Nakatani, M. Tsuchiizu, E. Fromager. Multiple impurities and combined local
density approximations in Site-Occupation Embedding Theory. Theor. Chem. Acc. 137: 169 (2018)
