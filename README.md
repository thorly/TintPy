[![Language](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)

## TintPy ##

The Thorly InSAR Tools in Python (TintPy) is an package for Interferometric Synthetic Aperture Radar (InSAR) time series processing. It is mainly based on [GAMMA](https://www.gamma-rs.ch/no_cache/software.html), [MintPy](https://github.com/insarlab/MintPy), [StaMPS](https://github.com/dbekaert/StaMPS) and [GMT](https://www.generic-mapping-tools.org/).

## Setup TintPy ##

To use the package, you need to setup the environment a) by adding ${TINTPY_HOME} to your ${PYTHONPATH} to make TintPy importable in Python and b) by adding ${TINTPY_HOME} to your ${PATH} to make application scripts executable in command line, as shown bellow.

Add to your ~/.bashrc file for bash user. Source the file for the first time. It will be sourced automatically next time when you login.

```Bash
# >>> TintPy >>>
export TINTPY_HOME=/home/ly/TintPy
export PATH=${PATH}:${TINTPY_HOME}/tintpy
export PYTHONPATH=${PYTHONPATH}:${TINTPY_HOME}/tintpy
# <<< TintPy <<<
```
