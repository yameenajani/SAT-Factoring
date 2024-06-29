**A Hybrid SAT and Lattice Reduction Approach for Integer Factorization**

This repository contains the code used in this project. The dependencies and steps to run the code are provided below.

The logic for the CAS implementation integrated with the SAT solver is available in [`maplesat/core/Solver.cc`](maplesat/core/Solver.cc) from `line 452` to `line 907`.

The scripts used to create and modify the SAT instances are available in the [`tests/scripts/`](tests/scripts/) directory. The purpose of each script is briefly outlined below -
* [`gen_nums.py`](tests/scripts/gen_nums.py) - Used to generate data (RSA modulus).
* [`generate_cubes.py`](tests/scripts/generate_cubes.py) - Selects random bits of the primes to leak (for varying percentages).
* [`modify_cnf.py`](tests/scripts/modify_cnf.py) - Modifies the base SAT encoding to include known bits and adds clauses encoding the decryption exponent `d` if the low exponent RSA encoding is being used.

The [`tests/data/`](tests/data/) directory contains some test cases.

<br>

__<u>Dependencies</u>__

The following are the dependencies required to run the code:
* [GMP](https://gmplib.org/)
* [fpLLL](https://github.com/fplll/fplll)
* [FLINT](https://flintlib.org/)
* [gcc](https://gcc.gnu.org/)
* [Glasgow Haskell Compiler](https://www.haskell.org/ghc/)
* [Python](https://www.python.org/)

<br>

__<u>Folder structure for tests</u>__

- `tests/`
    - `data/`
        - `{bitsize}/` (bitsize of primes; example: 128)
            - `data_{instance_num}.json`
    - `instances/`
        - `{bitsize}/` (bitsize of primes; example: 128)
            - `cubes/unbal` (contains json files for each instance containing the random leaked bits for different % of known bits)
                - `cubes_{instance_num}.json`
    - `solvers/`
        - `maplesat` (compiled version of the solver)

<br>

> [!IMPORTANT]
> <br>The _bitsize_ parameter corresponds to the number of bits in the primes, not the RSA modulus. If _bitsize_ is 128, then `len(p) = len(q) = 128` and `len(N) = 256` (or `255`).

__<u>How to run</u>__

To run the provided test cases, `cd` into the [`tests/`](tests/) directory and run the following:
```
./run_test.sh -b {bitsize} -p {percent} -i {instance number}
```
If you want to include the decryption exponent `d` in the SAT instance, use the `-d` flag:
```
./run_test.sh -b {bitsize} -p {percent} -i {instance number} -d
```
If you want to use the SAT+CAS approach that calls Coppersmith from within the solver, use the `-l` flag:
```
./run_test.sh -b {bitsize} -p {percent} -i {instance number} -l
```

The output files will be present in the [`tests/outputs/`](tests/outputs/) directory and will be available to view once execution of the corresponding method is complete.
