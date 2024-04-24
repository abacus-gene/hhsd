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

## Further details
Further details on the operation of the program, including best practices, and a detailed explanation of parameter syntax can be found in the manual.
