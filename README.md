# marked-particle-processes-py
(To be) library/module for simulating and/or statistical analysis of marked particle processes realizations.
Currently only for 2D simulations.

### Installation
#### 1) Cloning the project
* Asure that you have the `python` installed in your device.
* Clone the repo in branch `main`.
* Install the required packages (stored in `requirements.txt` with corresponding versions needed).
* Setup a config (see below).
* Run `main.py`.
#### 2) Install directly
* TBD

## Config (WIP)
The "minimal" requirements for config values are parameters:
* `process_type` (str) - either `"ball"` or `"segment"`
* `intensity` (float) - intensity of a ground Poisson point process of germs (expected number in a unit-measure Borel subset of $\mathbb R^d$).
Example of "minimal config":
```
{
  "process_type": ball,
  "intensity": 100
}
```
"Maximal config options" can be found in `config_ext.json`, yet for some of the functionalities the code has to be updated (WIP).
