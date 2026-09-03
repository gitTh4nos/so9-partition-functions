# so9-partition-functions
This repository contains the Mathematica and Python code used in the project along with the necessary dependencies required for running the scripts.

> "The growth of $SO(9)$ super-representations in Type II string theory"

## Abstract
We study the growth of massive $SO(9)$ super-representations in type II and type I string theories.  We do this first directly by finding an empirical formula that specifies which representations appear at any level, and then compute the multiplicities from a refined partition function evaluated over finite fields up to level 501 for individual representations, and level 226 for all representations.   We then derive asymptotic formulae for the growth of any representation, which are constructed by integrating about the peaks of the refined partition function.  Using a Rademacher sum we can find very accurate approximations for the multiplicities of the representations.  Unlike the superstring partition function, which only receives contributions from the odd Rademacher terms, the multiplicities for any representation have contributions from both even and odd terms.  Finally, we apply the same techniques to the ordinary representations, where a simplified refined partition function allows for better computational speed and simpler expressions for the asymptotic approximations.

## Repository structure 
```text
so9-partition-functions/
├── Notebooks/
│   ├── asymp_analysis.nb             # computes asymptotic behavior for super-representations and compares it to actual multiplicities
│   ├── asymp_analysis_ordinary.nb    # similar to asymp_analysis.nb but for ordinary representations
│   ├── demo-script.nb                # Mathematica version of the python code partition_func_weyl_sparse.py (no optimizations included)
│   └── plots_sec3.nb                 # provides the data for the table and makes the 3D plot(s) shown in section 3
│
├── Python/
│   ├── 50_lowest_reps.py
│   ├── 50_lowest_reps_ord.py
│   ├── partition_func_weyl_sparse.py       # generates multiplicities for super-representations
│   ├── partition_func_weyl_sparse_ord.py   # computes multiplicities for ordinary representations
│   └── representation_finite_n_multirep.py # computes multiplicities for 1 or more super-representations
│
├── Scripts/
│   ├── coeff_gen_1.wl                # generating coefficients used in the 1st terms in the Rademacher sums
│   ├── coeff_gen_2.wl                # generating coefficients used in the 2nd terms in the Rademacher sums
│   ├── coeff_gen_3.wl                .
│   ├── coeff_gen_4.wl                .
│   ├── coeff_gen_5.wl                .
│   ├── demo-script.wl                # script version of demo-script.nb
│   └── math_Code_simplified_batch.wl # generates multiplicities for super-representations, up to level ~100
│
└── README.md
```
Data files produced by the python code and the Wolfram scripts are provided in `DataFiles/`
and the reader can skip directly to the notebooks if desired.

## Requirements 
### Mathematica 
- Wolfram Mathematica 15.0 (tested)
### Python 
- Python 3.14.4 (tested)
- No third-party packages are required.
## Usage
### Python Scripts
Here are some instructions on how to run the python scripts from the **root directory** of this repository. Note that depending on how Python is installed and configured on the operating system you might need to use the command `python` or `python3` in order to run the scripts. Also all the commands below **assume that Python is available in your system** `PATH`.
 - For the `partition_func_weyl_sparse.py` run: 
```bash 
    python Python/partition_func_weyl_sparse.py N
```
where $N$ is the level (Note that the for these purposes, the first massive level is at $N=0$, so $N=100$ means the 101st level). Also one can run: 
 ```bash
    python Python/partition_func_weyl_sparse.py --help
 ```
to see information about the variety of flags that can be used in this script.
 - For the `representation_finite_n_multirep.py` run: 
```bash
    python Python/representation_finite_n_multirep.py N L M R S
```
where $L$, $M$, $R$, and $S$ are the four Dynkin indices.  For multiple files run: 
```bash 
    python Python/representation_finite_n_multirep.py N --rep-file FILE
```
An example of `FILE` and the format is in `50_lowest_reps.py`. Again, you can run: 
```bash 
    python Python/representation_finite_n_multirep.py --help
```
to see all flags available.
 - For the `partition_func_weyl_sparse_ord.py` run:
```bash
    python Python/partition_func_weyl_sparse_ord.py N
```
where $N$ is the level. To get the multiplicities for a single representation run: 
```bash 
    python Python/partition_func_weyl_sparse_ord.py N --rep L,M,R,S
```
To get the multiplicities for multiple representations stored in a file run: 
```bash 
    python Python/partition_func_weyl_sparse_ord.py N --reps-file FILE  
```
Note that the format for the representations is different than in the super-case (there are commas between the indices).  An example file with the format is `50_lowest_reps_ord.py`. Finally run:
```bash 
    python Python/partition_func_weyl_sparse_ord.py --help
```
to see all flags.
### Wolfram Scripts
To run from terminal on a Mac, add the following to your `.zprofile` file :
```bash
    PATH=$PATH:/Applications/Wolfram.app/Contents/MacOS
``` 
To run from terminal on Windows, you have to have the directory containing `WolframKernel.exe` to your Windows `PATH`. For the `math_Code_simplified_batch.wl`, `coeff_gen_*.wl` wolfram scripts instructions on how to run are given in the beginning of the files. 
To run `demo-script.wl` from the **root directory** of the repository type: 
```bash 
    wolfram -script Scripts/demo-script.wl N Path/to/outputfile
```
where $N$ is the level and `outputfile` is the location where you want to save the output in a `.m` file.
## Large files and Git LFS
This repository uses [Git Large File Storage (Git LFS)](https://git-lfs.com/) for several large Mathematica data files.

To obtain the complete files, make sure [Git](https://git-scm.com/) and [Git LFS](https://git-lfs.com/) are installed before cloning the repository:

```bash
git lfs install
git clone https://github.com/gitTh4nos/so9-partition-functions.git
cd so9-partition-functions
git lfs pull
```

The repository's large files include:

* `DataFiles/Ordinary/level_201.m`
* `DataFiles/Super/level_200.m`
* `DataFiles/Super/level_225.m`

**Note:** Downloading the repository using GitHub's **Code → Download ZIP** option may not provide the actual Git LFS files. For the complete dataset, clone the repository with Git LFS installed or download the corresponding release/archive from Zenodo.

## Useful links
* **Paper:** arxiv (will be added later on...)
* **Code and data:** [Zenodo](https://doi.org/10.5281/zenodo.21628613)


## Citation
To be added later on...
## License 
This work is licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

See the full license here: https://creativecommons.org/licenses/by/4.0/
