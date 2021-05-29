# An Improved Approximation Algorithm for the Reversal and Transposition Distance Considering Gene Order and Intergenic Sizes

This project is divided into two folders:
- Datasets: This folder contains the DS1 and DS2 datasets.
- Algorithms: This folder contains the source code in Python of two algorithms for the Sorting by Intergenic Reversals and Transpositions (SbIRT) problem.
-- 4SbIRT.py: A 4-approximation algorithm for SbIRT problem.
-- 3SbIRT.py: A 3-approximation algorithm for SbIRT problem.

Running the algorithms with the default implementation:
```sh
python3 4SbIRT.py --file <input_file>
python3 3SbIRT.py --file <input_file>
```

Running the algorithms with the greedy strategy activated:
```sh
python3 4SbIRT.py --file <input_file> --greedy T
python3 3SbIRT.py --file <input_file> --greedy T
```