# HHSD
`hhsd` can be used to perform Hierarchical Heuristic Species Delimitation from multilocus genetic data.
## Installation instructions
### Prerequisites
- Users must install BPP 4.6+ before installing `hhsd`. check your currently installed version with the `bpp` commmand, and follow the instructions at https://github.com/bpp/bpp if this requirement is not met.
- Users must have python 3.9+ installed. Check your current version with `python --version`. 
### - Linux 
Installation on linux can be accomplshed via the following set of commands:
```
git clone https://github.com/abacus-gene/hhsd
cd hhsd
pip install .
```
### - macOS
Installation on macOS is also relatively straightforward. Make sure that you are using the pip version specific to python3.
```
git clone https://github.com/abacus-gene/hhsd
cd hhsd
pip3 install .
```
### - Windows 
Windows installations are slightly more complicated. First install the program as expected, but do not close the command prompt.
```
git clone https://github.com/abacus-gene/hhsd
cd hhsd
pip install .
```
The output of pip may include a warning message about the executable file being placed in a directory not in $PATH. Copy the directory in the warning message, and navigate to the location using the file explorer. Then copy the `hhsd.exe` file to the directory where `bpp` is installed. If you set up your $PATH correctly for bpp, hhsd can inherit it. 

### Checking your install
To check if hhsd was installed successfully, open your command line of choice, and type:
```
hhsd
```
This should result in a greeting message similar to the following:
```
hhsd version 0.9.7
12 cores available
specify control file for analysis with --cfile
```
