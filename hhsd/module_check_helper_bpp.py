'''
FUNCTIONS NECESSARY TO VERIFY IF THE BPP SPECIFIC PARAMETERS THAT ARE PROVIDED
ARE COMPATIBLE WITH THE FORMATTING STANDARDS OF BPP, THE DATASET, AND THE COMPUTER
'''

import os
import sys

from .module_helper import check_numeric
from .module_msa_imap import alignfile_to_MSA


# check if the nloci parameter is an int
def check_nloci(
        nloci
        ):

    if nloci != None:
        if not check_numeric(nloci, "1<=x<100000", "i"):
            sys.exit("ParameterFormattingError: 'nloci' must be a positive integer.")

# check that the number of loci to check is less than or equal to the loci in the MSA
def check_nloci_msa_compat(
        input_nloci, 
        seqfile
        ):

    if input_nloci != None:
        true_nloci = len(alignfile_to_MSA(seqfile))
        user_nloci = int(input_nloci)

        if user_nloci >= true_nloci:
            sys.exit(f"ParameterIncompatibilityError: 'nloci' ({input_nloci}) larger than number of loci in seqfile ({true_nloci}).")

# check that the phase parameter is correctly formatted
def check_phase(
        phase,
        ):

    if phase != None:
        if phase not in ["0", "1"]:
            sys.exit("ParameterFormattingError: phasing for all sequences is specifed as a single digit, with 0 (unphased) or 1 (phased).")
            
        

# check if the threads parameter passed to BPP is correctly specified
def check_threads(
        threads
        ):

    n_cpu = int(os.cpu_count())

    '''
    This quite complex function checks if the number, starting point, and offset of the 
    requested threads is actually compatible with the number of threads available
    to the machine. This is beacuse if the user requests more cores than available, or tries to
    request that BPP pins threads to nonexistent cores, BPP will crash. 
    '''

    if threads == None:
        sys.exit(f"MissingParameterError: 'threads' must be specified for optimal performance.\n{n_cpu} cores are available.")

    if threads != None:
        try:
            th = threads.split()
        except:
            sys.exit("ParameterFormattingError: 'threads' incorrectly formatted. Refer to section 4.7 of the manual.")

        # check if threads is 3 integers, and all values are at least 1 (there is no such thing as 0 threads)
        if all(check_numeric(num, "0<x<1024", "i") for num in th):
            th = [int(x) for x in th]
            # if only one integer is given, that is interpreted as the number of threads requested
            if len(th) == 1:
                # check that the number of threads requested <= threads in the CPU
                if n_cpu < th[0]:
                    sys.exit(f"ResourceError: more 'threads' requested ({threads}) than available on computer ({n_cpu}).\nDecresase thread count.")
            elif len(th) == 2:
                # check that the requested offset and the number of threads still fits the CPU
                if n_cpu < (th[0] + (th[1]-1)):
                    sys.exit(f"ResourceError: 'threads' implies more cores ({th[0] + (th[1]-1)}) than available on computer ({n_cpu}).\nDecrease thread count and/or offset.")
            elif len(th) == 3:
                # check that all the requested threads still fit the CPU
                if n_cpu < ((th[1]-1) + (th[2]*th[0])):
                    sys.exit(f"ResourceError: 'threads' implies more cores ({(th[1]-1) + (th[2]*th[0])}) than available on computer ({n_cpu}).\nDecrease thread count and/or offset and/or interval")
        else:
            sys.exit("ParameterFormattingError: 'threads' should only contain integers. Refer to section 4.7 of the manual.")
        

# check that the number of threads requested <= the number of loci in the MSA
def check_threads_msa_compat(
        input_threads, 
        seqfile
        ):

    if input_threads != None:
        n_threads = int(input_threads.split()[0])
        true_nloci = len(alignfile_to_MSA(seqfile))

        if n_threads > true_nloci:
            sys.exit(f"ParameterIncompatibilityError: more 'threads' requested ({n_threads}) than the number of loci in seqfile ({true_nloci}).\ndecrease thread count.")

# check that the number of threads requested <= the number of loci specified by the user
def check_threads_nloci_compat(
        input_threads, 
        input_nloci
        ):

    if input_nloci != None:
        if input_threads != None:
            n_threads = int(input_threads.split()[0])
            if n_threads > int(input_nloci):
                sys.exit(f"ParameterIncompatibilityError: more 'threads' requested ({n_threads}) than 'nloci' ({input_nloci}).\ndecrease thread count.")

# check if the locusrate parameter is correctly formatted
def check_locusrate(
        locusrate
        ):

    if locusrate != None:
        try:
            lr_par = locusrate.split()
            asd = lr_par[0]
        except:
            sys.exit("LocusRateError: 'locusrate' incorrectly formatted. refer to BPP manual")

        if len(lr_par) != 1 and lr_par[0] == "0":
            sys.exit("LocusRateError: if 'locusrate' begins with 0, there can be no further subparameters.")          
        
        elif len(lr_par) == 4:
            if not (lr_par[0] == "1" and all(check_numeric(value, "0<=x<150") for value in lr_par[1:3])):
                sys.exit("LocusRateError: 'locusrate' with four subparameters incorrectly formatted. refer to BPP manual")
        
        elif len(lr_par) == 5:
            if not (lr_par[0] == "1" and all(check_numeric(value, "0<=x<150") for value in lr_par[1:3]) and (lr_par[4] in ["iid", "dir"])):
                sys.exit("LocusRateError: 'locusrate' with five subparameters incorrectly formatted. refer to BPP manual")
        

# check that the number of burnin samples meets the minimum requirement
def check_burnin(
        burnin,
        ):

    if burnin == None:
        sys.exit("MissingParameterError: 'burnin' not specified.")
    elif not check_numeric(burnin, "200<=x", "i"):
        sys.exit("McmcParameterError: 'burnin' must be integer value >= 200")

# check that the number of samples meets the minimum requirement
def check_nsample(
        nsample,
        ):

    if nsample == None:
        sys.exit("MissingParameterError: 'nsample' not specified.")
    elif not check_numeric(nsample, "1000<=x", "i"):
        sys.exit("McmcParameterError: 'nsample' must be integer value >= 1000")

# check that the sampling frequency is in the requried range
def check_sampfreq(
        sampfreq,
        ):

    if sampfreq != None:
        if not check_numeric(sampfreq, "0<x<=100", "i"):
            sys.exit("McmcParameterError: 'sampfreq' must be integer value between 0 and 100")

# check if the data cleaning parameter is correctly specified
def check_cleandata(
        cleandata,
        ):

    if cleandata != None:
        if cleandata not in ["0", "1"]:
            sys.exit("ParameterFormattingError: 'cleandata' must be 0 or 1.")

# check if the seed is an int value
def check_seed(
        seed
        ):

    if seed != None:
        if not check_numeric(seed, "0<x<10000000000", "i"):
            sys.exit("ParameterFormattingError: 'seed' parameter must be positive integer value")



# check if the tau prior parameter passed to BPP is correctly specified
def check_tauprior(
        tauprior
        ):

    default_error_msg = "PriorError: 'tauprior' expectes the following syntax:\n    tauprior = invgamma alpha beta\n    tauprior = gamma alpha beta"

    if tauprior != None:
        try:
            t = tauprior.split()
        except:
            sys.exit(default_error_msg)
            
        if len(t) <= 2 or len(t) == 3 and t[0] not in ["invgamma", "gamma"] or len(t) > 3:
            sys.exit(default_error_msg)

        if len(t) == 3 and t[0] == "invgamma":
            if not (check_numeric(t[1], "1<x<100") and check_numeric(t[2], "0<x<100")):
                sys.exit("PriorError: 'alpha' parameter of inverse gamma for 'tauprior' must be > 1, and 'beta' must be > 0.")
        
        elif len(t) == 3 and t[0] == "gamma":
            if not (check_numeric(t[1], "0<x<10000") and check_numeric(t[2], "0<x<10000")):
                sys.exit("PriorError: 'alpha' and 'beta' parameters of gamma for 'tauprior' must be > 0.")

        

# check if the theta prior parameter passed to BPP is correctly specified
def check_thetaprior(
        thetaprior
        ):
    
    default_error_msg = "PriorError: 'thetaprior' expectes the following syntax:\n    thetaprior = invgamma alpha beta\n    thetaprior = gamma alpha beta"

    if thetaprior != None:
        try:
            th = thetaprior.split()
        except:
            sys.exit(default_error_msg)
            
        if len(th) <= 2 or len(th) == 3 and th[0] not in ["invgamma", "gamma"] or len(th) > 3:
            sys.exit(default_error_msg)

        if len(th) == 3 and th[0] == "invgamma":
            if not (check_numeric(th[1], "2<x<100") and check_numeric(th[2], "0<x<100")):
                sys.exit("PriorError: 'alpha' parameter of inverse gamma for 'thetaprior' must be > 2, and 'beta' must be > 0.")
        
        elif len(th) == 3 and th[0] == "gamma":
            if not (check_numeric(th[1], "0<x<10000") and check_numeric(th[2], "0<x<10000")):
                sys.exit("PriorError: 'alpha' and 'beta' parameters of gamma for 'thetaprior' must be > 0.")
        

# check the migration rate prior
def check_wprior(
        wprior,
        ):

    if wprior != None:
        try:
            mp = wprior.split()
            if len(mp) != 2 or (not all(check_numeric(value, "0<=x<=50") for value in mp)):
                sys.exit(f"PriorError: 'wprior' incorrectly formatted as '{wprior}'. Refer to section 4.4 of the manual.")
        except:
            sys.exit(f"PriorError: 'wprior' incorrectly formatted as '{wprior}'. Refer to section 4.4 of the manual.")
            

# # check if the finetune parameter passed to BPP is correctly specified
# def check_finetune(
#         finetune
#         ):

#     if finetune != None:
#         try:
#             ft = finetune.split()
#             # check if the first finetune is 0 or 1 followed by a ":", and the remainders are integers
#             if ft[0] == "0" or ft[0] == "1":
#                 if not all(check_numeric(value, "0<=x<=10") for value in ft[1:]):
#                     sys.exit("cfile error: 'finetune' parameter incorrectly formatted. refer to BPP manual")
#         except:
#             sys.exit("cfile error: 'finetune' parameter incorrectly formatted. refer to BPP manual")
