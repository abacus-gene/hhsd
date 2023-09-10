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
- If not, download the code by visiting https://github.com/abacus-gene/hhsd/archive/refs/heads/main.zip.

Then, install the program as expected, but do not close the command prompt.
```
cd hhsd
```
```
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
## Trial run
Following the install, users can perform a qucik trial run to familiarize themselves with the program. 
### The dataset:
This simulated dataset was originally introduced in https://doi.org/10.1093/sysbio/syy051. There are five populations, with 𝐴, 𝐵, 𝐶, 𝐷 representing geographical populations of a species with a wide geographic distribution, while 𝑋 is a new species that split off from population 𝐴. The data consisted of 𝐿 = 100 simulated loci, with two diploid sequences sampled per species per locus, and 500 sites in the sequence.
### Instructions:
1) Investigate the contents of the `examples/simulated_abcdx` folder using a graphical file browser or the command line:
    - The PHYLIP formatted multiple sequence alignment (MSA) file, `MySeq.txt`
    - The Imap file specifying the assignment of individuals to populations in the guide tree, `MyImap.tzt`
    - The control file used to specify the analytic procedure, `cf_sim_merge.txt`
2) Run the example analysis by navigating to to the `examples/simulated_abcdx` folder using the command line, and running:
```
hhsd --cfile cf_sim_merge.txt
```
3) Inspect the output appearing in the command line, and the output files written to the `examples/simulated_abcdx/res_sim_merge` folder.
