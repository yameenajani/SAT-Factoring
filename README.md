**A Hybrid SAT and Lattice Reduction Approach for Integer Factorization**

This repository contains the code used in this project. The dependencies and steps to run the code are provided below.

The logic for the CAS implementation integrated with the SAT solver is available in [`core/Solver.cc`](core/Solver.cc) from `line 450` to `line 973`.

The scripts used to create and modify the instances are available in the [`tests/scripts/`](tests/scripts/) directory. The purpose of each script is briefly outlined below -
* [`gen_nums.py`](tests/scripts/gen_nums.py) - Used to generate data (primes, instance link from CNF Generator).
* [`get_cnf.py`](tests/scripts/get_cnf.py) - Downloads the instance using the link and does basic manipulation (assigns extra MSB variables of prime to 0)
* [`update_cnf.py`](tests/scripts/update_cnf.py) - Updates the instance by assigning random known variables of the primes based on the percentage provided.
* [`add_constraints.py`](tests/scripts/add_constraints.py) - Adds additional clauses in the case where `d` is introduced.

The [`tests/data/`](tests/data/) and [`tests/instances`](tests/instances/) directories contain some test cases.

<br>

__<u>Dependencies</u>__

The following are the dependencies required to run the code -
* [GMP](https://gmplib.org/)
* [fpLLL](https://github.com/fplll/fplll)
* [FLINT](https://flintlib.org/)
* [gcc](https://gcc.gnu.org/)
* [SageMath](https://www.sagemath.org/) (Optional. Required only to run custom test cases.)
* [BeautifulSoup4](https://pypi.org/project/beautifulsoup4/) (Optional. Required only to run custom test cases.)

<br>

__<u>Folder structure for tests</u>__

- `tests/`
    - `data/`
        - `{bitsize}/` (bitsize of primes; example: 128)
            - `data.json`
    - `instances/`
        - `{bitsize}/` (bitsize of primes; example: 128)
            - `original_instances/` (contains the original instance without any modifications)
                - `instance.cnf`
            - `cubes/` (contains json files for each instance containing the random leaked bits for different % of known bits)
                - `cube.json`
            - `temp/` (contains modified instances with random known bits of p, q)
                - `{percent}/` (the percent of random known bits)
                    - `instance.cnf`
            - `temp_withd/` (contains modified instances with random known bits of p, q along with the low public exponent encoding)
                - `{percent}/` (the percent of random known bits)
                    - `instance.cnf`
    - `solvers/`
        - `maplesat_cb` (compiled version of the solver)

<br>

> [!IMPORTANT]
> <br>The _bitsize_ correspond to the number of bits in the primes and not the semi-prime. If the _bitsize_ is 128 then `len(p) = len(q) = 128` and `len(N) = 256`.

__<u>How to run</u>__

First, load all dependencies and compile the solver by running -
```
./compile_solver.sh
```
This will put the compiled version of the solver in the [`tests/solvers/`](tests/solvers/) directory.

To run the provided test cases, `cd` into the [`tests/`](tests/) directory and run the following -
```
./run_test_case.sh -b {bitsize} -p {percent} -i {instance number}
```
For example, if you wish to run the code on [`tests/instances/128/75/instance_1.cnf`](tests/instances/128/75/instance_1.cnf) then you would run -
```
./run_test_case.sh -b 128 -p 75 -i 1
```
If you want to include `d`, you can use the `-d` flag and run the following command -
```
./run_test_case.sh -b {bitsize} -p {percent} -i {instance number} -d
```

> [!NOTE]
> <br>You will not be able to run instances that don't exist. You will need to create new instances if you want to run custom test cases. The `run_test_case.sh` script will handle this for you. Just provide an instance number that does not exist to create a new test case.

<br>

<!-- To generate custom test cases, first you need to create the data file for the test case. To do this, `cd` into the [`tests/scripts/`](tests/scripts/) directory and run the following commands in sequence -
```
./gen_nums.py {bitsize} {number_of_instaces_to_generate}

./get_instance.sh -b {bitsize} -n {number_of_instaces_to_generate}
```
For example, if you wish to generate 3 data files for 256-bit N then you would run -
```
./gen_nums.py 128 3

./get_instance.sh -b 128 -n 3
```
After this you can simply run the [`run_test_case.sh`](tests/run_test_case.sh) script with correct parameters as mentioned above.

<br> -->

In both cases the instance is solved using SAT+CAS as well as standalone SAT methods. The output files will be present in the [`tests/outputs/`](tests/outputs/) directory and will be available to view once execution of the corresponding method is complete.