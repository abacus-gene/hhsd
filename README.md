# HHSD
`hhsd` can be used to perform Hierarchical Heuristic Species Delimitation from multilocus genetic data.
---
## Installation instructions
### Prerequisites
- Users must have python 3.9+ installed. Check your current version with `python --version`. 
### - Linux 
Installation on linux can be accomplshed via the following set of commands:
```
git clone https://github.com/abacus-gene/hhsd
```
```
cd hhsd
```
```
pip install .
```
### - macOS
Installation on macOS is also relatively straightforward. Make sure that you are using the pip version specific to python3.
```
git clone https://github.com/abacus-gene/hhsd
```
```
cd hhsd
```
```
pip3 install .
```
### - Windows 
Windows installations are slightly more complicated. 
- If git is aldready installed, run the command: `git clone https://github.com/abacus-gene/hhsd`
- If not, download the code by visiting https://github.com/abacus-gene/hhsd/archive/refs/heads/main.zip. Rename the zip file to `hhsd.zip`, and then uncompress.

Then, install the program as expected:
```
cd hhsd
```
```
pip install .
```

### Checking your install
To check if hhsd was installed successfully, open your command line of choice, and type:
```
hhsd
```
This should result in a greeting message similar to the following:
```
hhsd version 0.9.9
12 cores available
specify control file for analysis with --cfile
```
---

## Detailed instructions
Details on the operation of the program, including best practices, and a thorough explanation of parameter syntax can be found in the manual.
