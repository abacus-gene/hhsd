# HHSD
`hhsd` can be used to perform Hierarchical Heuristic Species Delimitation from multilocus genetic data.
---
## Installation instructions
### Prerequisites
- Users must have python 3.9+ installed. Check your current version with `python --version`. 
- Users must **create** and **activate** a dedicated **virtual environment** for the installation of the package. See: https://docs.python.org/3/tutorial/venv.html for detailed instuctions on how this can be accomplished. 
- Create a dedicated folder (in the directory where you normally install your programs), and move your terminal to that folder.
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
hhsd version 1.1.0
12 cores available
specify control file for analysis with --cfile
```
---

## Detailed instructions
Details on the operation of the program, including best practices, and a thorough explanation of parameter syntax can be found in [https://github.com/abacus-gene/hhsd/blob/main/manual.pdf](the manual).
