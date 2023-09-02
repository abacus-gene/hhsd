'''
THIS MODULE CONTAINS THE FUNCTIONS NECESSARY TO VERIFY IF THE BPP SPECIFIC PARAMETERS THAT ARE PROVIDED
ARE COMPATIBLE WITH THE FORMATTING STANDARDS OF BPP, THE DATASET, AND THE COMPUTER
'''

import os
import sys

from module_helper import check_numeric
from module_msa_imap import alignfile_to_MSA


# check if the nloci parameter is an int
def check_nloci(
        nloci
        ):

    if nloci != None:
        if not check_numeric(nloci, "1<=x<100000", "i"):
            sys.exit("mcfile error: 'nloci' must be a positive integer.")

# check that the number of loci to check is less than or equal to the loci in the MSA
def check_nloci_msa_compat(
        input_nloci, 
        seqfile
        ):

    if input_nloci != None:
        true_nloci = len(alignfile_to_MSA(seqfile))
        user_nloci = int(input_nloci)

        if user_nloci >= true_nloci:
            sys.exit(f"mcfile error: 'nloci' ({input_nloci}) larger than number of loci in seqfile ({true_nloci}).")

# check that the phase parameter is correctly formatted
def check_phase(
        phase,
        ):

    if phase != None:
        if phase not in ["0", "1"]:
            sys.exit("mcfile error: phasing for all sequences is specifed as a single digit, with 0 (unphased) or 1 (phased).")
            
        

# check if the threads parameter passed to BPP is correctly specified
def check_threads(
        threads
        ):

    '''
    This quite complex function checks if the number, starting point, and offset of the 
    requested threads is actually compatible with the number of threads available
    to the machine. This is beacuse if the user requests more cores than available, or tries to
    request that BPP pins threads to nonexistent cores, BPP will crash. 
    '''

    if threads != None:
        # find the available number of CPU cores
        n_cpu = int(os.cpu_count())
        try:
            th = threads.split()
        except:
            sys.exit("mcfile error: 'threads' incorrectly formatted. refer to BPP manual")

        # check if threads is 3 integers, and all values are at least 1 (there is no such thing as 0 threads)
        if all(check_numeric(num, "0<x<1024", "i") for num in th):
            th = [int(x) for x in th]
            # if only one integer is given, that is interpreted as the number of threads requested
            if len(th) == 1:
                # check that the number of threads requested <= threads in the CPU
                if n_cpu < th[0]:
                    sys.exit(f"mcfile error: more 'threads' requested ({threads}) than available on computer ({n_cpu}).\ndecresase thread count.")
            elif len(th) == 2:
                # check that the requested offset and the number of threads still fits the CPU
                if n_cpu < (th[0] + (th[1]-1)):
                    sys.exit(f"mcfile error: 'threads' implies more cores ({th[0] + (th[1]-1)}) than available on computer ({n_cpu}).\ndecrease thread count and/or offset.")
            elif len(th) == 3:
                # check that all the requested threads still fit the CPU
                if n_cpu < ((th[1]-1) + (th[2]*th[0])):
                    sys.exit(f"mcfile error: 'threads' implies more cores ({(th[1]-1) + (th[2]*th[0])}) than available on computer ({n_cpu}).\ndecrease thread count and/or offset and/or interval")
        else:
            sys.exit("mcfile error: 'threads' should only contain integers. refer to BPP manual")
        

# check that the number of threads requested <= the number of loci in the MSA
def check_threads_msa_compat(
        input_threads, 
        seqfile
        ):

    if input_threads != None:
        n_threads = int(input_threads.split()[0])
        true_nloci = len(alignfile_to_MSA(seqfile))

        if n_threads > true_nloci:
            sys.exit(f"mcfile error: more 'threads' requested ({n_threads}) than the number of loci in seqfile ({true_nloci}).\ndecrease thread count.")

# check that the number of threads requested <= the number of loci specified by the user
def check_threads_nloci_compat(
        input_threads, 
        input_nloci
        ):

    if input_nloci != None:
        if input_threads != None:
            n_threads = int(input_threads.split()[0])
            if n_threads > input_nloci:
                sys.exit(f"mcfile error: more 'threads' requested ({n_threads}) than 'nloci' ({input_nloci}).\ndecrease thread count.")

# check if the locusrate parameter is correctly formatted
def check_locusrate(
        locusrate
        ):

    if locusrate != None:
        try:
            lr_par = locusrate.split()
            asd = lr_par[0]
        except:
            sys.exit("mcfile error: 'locusrate' incorrectly formatted. refer to BPP manual")

        if len(lr_par) != 1 and lr_par[0] == "0":
            sys.exit("mcfile error: if 'locusrate' begins with 0, there can be no further subparameters.")          
        
        elif len(lr_par) == 4:
            if not (lr_par[0] == "1" and all(check_numeric(value, "0<=x<150") for value in lr_par[1:3])):
                sys.exit("mcfile error: 'locusrate' with four subparameters incorrectly formatted. refer to BPP manual")
        
        elif len(lr_par) == 5:
            if not (lr_par[0] == "1" and all(check_numeric(value, "0<=x<150") for value in lr_par[1:3]) and (lr_par[4] in ["iid", "dir"])):
                sys.exit("mcfile error: 'locusrate' with five subparameters incorrectly formatted. refer to BPP manual")
        

# check that the number of burnin samples meets the minimum requirement
def check_burnin(
        burnin,
        ):

    if burnin == None:
        sys.exit("mcfile error: 'burnin' not specified.")
    elif not check_numeric(burnin, "200<=x", "i"):
        sys.exit("mcfile error: 'burnin' must be integer value >= 200")

# check that the number of samples meets the minimum requirement
def check_nsample(
        nsample,
        ):

    if nsample == None:
        sys.exit("mcfile error: 'nsample' not specified.")
    elif not check_numeric(nsample, "1000<=x", "i"):
        sys.exit("mcfile error: 'nsample' must be integer value >= 1000")

# check that the sampling frequency is in the requried range
def check_sampfreq(
        sampfreq,
        ):

    if sampfreq != None:
        if not check_numeric(sampfreq, "0<x<=100", "i"):
            sys.exit("mcfile error: 'sampfreq' must be integer value between 0 and 100")

# check if the data cleaning parameter is correctly specified
def check_cleandata(
        cleandata,
        ):

    if cleandata != None:
        if cleandata not in ["0", "1"]:
            sys.exit("mcfile error: 'cleandata' must be 0 or 1.")

# check if the seed is an int value
def check_seed(
        seed
        ):

    if seed != None:
        if not check_numeric(seed, "0<x<10000000000", "i"):
            sys.exit("mcfile error: 'seed' parameter must be integer value")

# check if the finetune parameter passed to BPP is correctly specified
def check_finetune(
        finetune
        ):

    if finetune != None:
        try:
            ft = finetune.split()
            # check if the first finetune is 0 or 1 followed by a ":", and the remainders are integers
            if ft[0] == "0" or ft[0] == "1":
                if not all(check_numeric(value, "0<=x<=10") for value in ft[1:]):
                    sys.exit("mcfile error: 'finetune' parameter incorrectly formatted. refer to BPP manual")
        except:
            sys.exit("mcfile error: 'finetune' parameter incorrectly formatted. refer to BPP manual")

# check if the tau prior parameter passed to BPP is correctly specified
def check_tauprior(
        tauprior
        ):

    if tauprior != None:
        try:
            t = tauprior.split()
            # check if the tau is two floats, and the second is less than 1
            if not (all(isinstance(float(x), float) for x in t) and len(t) == 2 and 0 < float(t[1]) < 1):
                sys.exit("mcfile error: 'tauprior' incorrectly formatted. refer to BPP manual")
        except:
            sys.exit("mcfile error: 'tauprior' incorrectly formatted. refer to BPP manual")

# check if the theta prior parameter passed to BPP is correctly specified
def check_thetaprior(
        thetaprior
        ):

    if thetaprior != None:
        try:
            t = thetaprior.split()
            # check if the theta is two floats followed by "e", the second is less than 1
            if not (all(isinstance(float(x), float) for x in t[0:1]) and len(t) == 3 and 0 < float(t[1]) < 1 and (t[2] == "e" or t[2] == "E")):
                sys.exit("mcfile error: 'thetaprior' incorrectly formatted. refer to BPP manual")
        except:
            sys.exit("mcfile error: 'thetaprior' incorrectly formatted. refer to BPP manual")

# check the migration rate prior
def check_migprior(
        migprior,
        ):

    if migprior != None:
        try:
            mp = migprior.split()
            if len(mp) != 2 or (not all(check_numeric(value, "0<=x<=50") for value in mp)):
                sys.exit(f"mcfile error: 'migprior' incorrectly formatted as '{migprior}'. please refer to BPP manual for formatting instructions")
        except:
            sys.exit(f"mcfile error: 'migprior' incorrectly formatted as '{migprior}'. please refer to BPP manual for formatting instructions")
            