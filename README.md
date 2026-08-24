# so9-partition-functions
This repository contains the Mathematica and Python code used in the project along with the necessary dependencies required for running the scripts.

> "The growth of $SO(9)$ super-representations in Type II string theory"

## Abstract
We study the growth of massive $SO(9)$ super-representations in type II and type I string theories.  We do this first directly by finding an empirical formula that specifies which representations appear at any level, and then compute the multiplicities from a refined partition function evaluated over finite fields up to level 501 for individual representations, and level 226 for all representations.   We then derive asymptotic formulae for the growth of any representation, which are constructed by integrating about the peaks of the refined partition function.  Using a Rademacher sum we can find very accurate approximations for the multiplicities of the representations.  Unlike the superstring partition function, which only receives contributions from the odd Rademacher terms, the multiplicities for any representation have contributions from both even and odd terms.  Finally, we apply the same techniques to the ordinary representations, where a simplified refined partition function allows for better computational speed and somewhat simplified expressions for the asymptotic approximations.

## Repository structure 
...
## Requirements 
### Mathematica 
- Wolfram Mathematica 15.0 (tested)
### Python 
- Python 3.14.4 (tested)
- No third-party packages are required.
## Usage
### Python
Here are some instructions on how to run the python scripts from the **root directory** of this repository:
 - For the `partition_func_weyl_sparse.py` run: 
```bash 
    python Python/partition_func_weyl_sparse.py N
```
or 
```bash
    Python/partition_func_weyl_sparse.py N
``` 
where N is the level (Note that the for these purposes, the first massive level is at N=0, so N=100 means the 101st level). Also one can run: 
 ```bash
    python partition_func_weyl_sparse.py --help
 ```
or 
```bash 
    python3 partition_func_weyl_sparse.py --help
```
to see information about the variety of flags that can be used in this script.
## Citation
...
## License 
...
