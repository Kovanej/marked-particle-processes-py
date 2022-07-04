# marked-particle-processes-py
(To be) library/module for simulating and/or statistical analysis of marked particle processes realizations.
Currently only for 2D simulations.

## Usage

### Installation
#### 1) Cloning the project
* Clone the repo in branch `main`.
* Install the required packages (stored in `requirements.txt` with corresponding versions needed).
* Setup a config (see below).
* Run `main.py`.
#### 2) Install directly

## Config (WIP)
The "minimal" requirements for config values are parameters:
* `process_type` (str) - either `"ball"` or `"segment"`
* `intensity` (float) - intensity of a ground Poisson point process of germs (expected number in a unit-measure Borel subset of $R^d$).