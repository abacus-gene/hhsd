'''
SPECIALIZED FUNCTIONS FOR INTERPRETING, AND ANALYSING THE 
CONTENTS OF SEQUENCE ALIGNMENTS AND IMAP FILES.
'''

import io
import re
import sys
from collections import Counter
from itertools import combinations
from itertools import product
from typing import Literal, Union, Dict

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from .customtypehints import ImapIndPop, ImapPopInd, Filename, NodeName, CfileParam, NewickTree
from .module_helper import readlines, remove_empty_rows
from .data_dicts import distance_dict, avail_chars
from .module_tree import get_first_split_populations

## IO HELPER FUNCTIONS


def imapfile_read(
        imap_filename:      Filename,
        output_type:        Literal['indpop','popind',], 
        ) ->                Union[ImapIndPop, ImapPopInd]:

    '''
    Read the imap text file into various dictionary forms.
    '''

    # ingest the imap file and remove all empty rows
    imap_lines = remove_empty_rows(readlines(imap_filename))

    # create the list of individual IDs and population IDs
    imap_indiv  = ([i.split(None, 2)[0] for i in imap_lines])
    imap_pop    = ([i.split(None, 2)[1] for i in imap_lines])

    # set up output type

    # dict with key = individual ids, value = population assignments
    if   output_type == "indpop":
        out:ImapIndPop = dict(zip(imap_indiv, imap_pop))

    # dict with key = pop names, value = IDs in that pop
    elif output_type == "popind":
        out:ImapPopInd = {}
        for index, label in enumerate(imap_pop):
            if label in out:
                out[label].append(imap_indiv[index])
            else:
                out[label] = [imap_indiv[index]]

    return out

def imapfile_write(
        input_indpop_dict:  ImapIndPop,
        imap_filename:      Filename,
        ):
    
    '''
    Write an imap dict to a text Imap file.
    '''

    rows = [f"{individual}\t{input_indpop_dict[individual]}" for individual in input_indpop_dict]
    output_file = '\n'.join(rows)
    
    with open(imap_filename, 'w') as f:
        f.write(output_file)


def alignfile_to_MSA(
        align_file:         Filename
        ) ->                list[MultipleSeqAlignment]:

    '''
    Return a properly filtered BioPython MSA object when pointed to a valid alignment file
    '''

    align_raw = readlines(align_file)
    
    # remove all empty rows, and superflous whitespaces from the numerical parameter rows, as this confuses biopython
    align_noempty = remove_empty_rows(align_raw)
    align_fixednum = []
    for line in align_noempty:
        # only applies to numeric lines, ignores others
        if re.search("[0-9]+[\s]+[0-9]+", line):
            line_fix = re.sub(" +", " ", line)
            line_fix = line_fix.strip()
            align_fixednum.append(line_fix)
        else:
            align_fixednum.append(line)
    
    # translate list to stringIO object that can be parsed by biopython
    align_str = "\n".join(align_fixednum)
    buf = io.StringIO(align_str)

    # read in the data from the buffered filtered version
    alignment_list = list(AlignIO.parse(buf,"phylip-relaxed")) #wrapped in list as AlignIO.parse returns generator
    
    return alignment_list



## FUNCTIONS FOR GENERATING THE SPECIES&TREE LINES OF THE BPP CONTROL FILE

def count_seq_per_pop(
        indpop_imap:        ImapIndPop, 
        input_MSA_list:     list[MultipleSeqAlignment]
        ) ->                Dict[NodeName, int]:

    '''
    This function counts the maximum number of sequences at a single loci that are associated with a 
    population in the sequence alignment. This output is used in "autoPopParam" to automatically 
    generate the "species&tree" lines for the BPP control file. The function is also used in 
    "check_GuideTree_Imap_compat" to ensure that each population has at least two haploid sequences associated with it.
    '''

    # create empty dict to hold results
    maxcounts = dict.fromkeys({value: key for key, value in indpop_imap.items()}, 0)
    
    for curr_locus in input_MSA_list:
        # filter out a list of individual ids at a given locus
        curr_id_list = [seq.id for seq in curr_locus]
        curr_id_list = [id.split("^")[1] for id in curr_id_list]
        # replace individual ids with their population code
        curr_pop_list = [indpop_imap[id] for id in curr_id_list]
        # count the number of sequences associated with each population, and keep track of the highest value
        counts = dict(Counter(curr_pop_list))
        for key in counts:
            if counts[key] > maxcounts[key]:
                maxcounts[key] = counts[key]
    
    return maxcounts

def auto_pop_param(
        indpop_imap:        ImapIndPop, 
        seqfile:            Filename,
        phase:              Literal[0, 1], 
        ) ->                CfileParam:

    '''
    Generate the lines of the bpp control file corresponding to population numbers and sizes.
    This function reads the alignment and the IMAP. It then extracts the population labels, and uses "count_Seq_Per_Pop" 
    to count the number of  sequences associated with that population in the alignment. This data is then formatted to 
    comply with the "species&tree" row of the BPP control file.
    '''

    # invert indpop into popind type dict
    popind_dict = {value: key for key, value in indpop_imap.items()}

    # load alignment
    alignment = alignfile_to_MSA(seqfile)
    
    # row describing the number and name of populations
    n_pops = str(len(popind_dict))
    pop_names = str(popind_dict.keys())[11:-2]
    pop_names = pop_names.replace(",",""); pop_names = pop_names.replace("'","")
    
    # row describing the ,maximum number of of sequences/loci for each population 
    maxcounts = count_seq_per_pop(indpop_imap, alignment)
    n_in_pop = str([maxcounts[key] for key in maxcounts])[1:-1].replace(","," ")

    # row describing phasing
    phase_row = f"{phase} "*int(n_pops)


    # final output, formatted to comply with BPP control dict standards
    rows = {"species&tree": f"{n_pops} {pop_names}", 
            "popsizes"    : n_in_pop,
            "newick"      : None,
            "phase"       : phase_row
           }

    return rows

def auto_nloci(
        seqfile:        Filename
        ) ->            int:

    '''
    Count the number of loci in the alignment
    '''

    alignment = alignfile_to_MSA(seqfile)
    nloci = str(len(alignment))
    
    return nloci





## FUNCTIONS FOR GENERATING THE TAU AND THETA PRIOR
'''
Implements an identical method to the Minimalist BPP web app of Prof. Bruce Rannala for
the automatic generation of tau and theta prior values.
'''

# return the alignment only containing individuals from a given population
def get_population_locus_alignment(locus, indpop_dict, population):
    # add the sequences belonging to the current population to a temp aligment
    temp_aligment = MultipleSeqAlignment([])
    for sequence in locus:
        id = sequence.id
        id = id.split("^")[1]
        if indpop_dict[id] == population:
            temp_aligment.append(sequence)

    return temp_aligment

# return the alignment only containing individuals from a given populations
def get_multi_population_locus_alignment(locus, indpop_dict, population_list):
    # add the sequences belonging to the current population to a temp aligment
    temp_aligment = MultipleSeqAlignment([])
    for sequence in locus:
        id = sequence.id
        id = id.split("^")[1]
        if indpop_dict[id] in population_list:
            temp_aligment.append(sequence)

    return temp_aligment

# calculate the pairwise distance between two sequences at shared known characters
def pairwise_dist   (
        seq_1: str, 
        seq_2: str
        ) -> float:

    '''
    This custom pairwise distance function is due to the fact that the standard identity
    based paiwise distance calculation of BioPython does not ignore gaps and unknowns, 
    and also does not work for phased diploid sequences. As such, a custom function was
    necessary to have identical results to those produced by the Minimalist BPP web app 
    of Prof. Bruce Rannala.

    The pairwise distances are calculated by: 
        1) Filtering down both alignments to sites with shared non-N IUPAC codes
        2) Using a lookup table from "data_dicts" to get the distance for each pair
        3) Averaging the results
    '''

    # collect sites where both sequences have comparable IUPAC codes
    seq_1_correct = set([i for i, x in enumerate(seq_1) if x in avail_chars])
    seq_2_correct = set([i for i, x in enumerate(seq_2) if x in avail_chars])
    overlap = list(seq_1_correct.intersection(seq_2_correct))
    # check for edge case where '???' characters overlap throghout the alignment, leading to a zero length overlap of valid characters
    if len(overlap) > 0:
        alignlist = [f"{seq_1[i]}{seq_2[i]}" for i in overlap]
        
        # at each site, use the lookup table to measure the distance
        dist_persite = [distance_dict[pair] for pair in alignlist]

        return float(np.round(np.average(dist_persite), decimals = 4))

    else:
        return np.nan


# return a list of all paiwise distances in an alignment

def get_Distance_list(
        input_MSA: MultipleSeqAlignment
        ) -> list[float]:

    '''
    This function finds all the unique sequence pairs in an MSA which should have
    their distances measured. It then iterates through these pairs, using 
    "pairwise_dist" to get a pairwise distance for each.
    '''

    # isolate only the sequence strings
    seqlist = [str(sequence.seq) for sequence in input_MSA]
    
    # procude the list corresponding to which two lists will be compared in which order
        # the convoluted order is implemented to match up with the BioPython "DistanceMatrix"
    seq_com = sorted(list(combinations(list(range(len(seqlist))), 2)), key=lambda x: x[1])

    # go through the list of combinations, and measure the distance for each
    dist_list = [pairwise_dist(seqlist[seq_com[i][0]], seqlist[seq_com[i][1]]) for i, _ in enumerate(seq_com)]

    return dist_list

def get_two_pop_distance_list(
        input_MSA_l: MultipleSeqAlignment,
        input_MSA_r: MultipleSeqAlignment,
        ) -> list[float]:


    # isolate only the sequence strings
    seqlist_l = [str(sequence.seq) for sequence in input_MSA_l]
    seqlist_r = [str(sequence.seq) for sequence in input_MSA_r]
    
    seq_com = list(product(np.arange(len(seqlist_l)),np.arange(len(seqlist_r))))

    # go through the list of combinations, and measure the distance for each
    dist_list = [pairwise_dist(seqlist_l[seq_com[i][0]], seqlist_r[seq_com[i][1]]) for i, _ in enumerate(seq_com)]

    return dist_list

# measure the average within population pairwise distance in a MSA
def distance_within_pop(alignment, indpop_dict, population):
    per_locus_dist = []
    per_locus_len = []

    for locus in alignment:
        # add the sequences belonging to the current population to a temp aligment
        temp_aligment = get_population_locus_alignment(locus, indpop_dict, population)
        
        # the distance can only be calculated for more than 2 sequences
        if len(temp_aligment) >= 2:
            
            # get pairwise distances within the temp alignment
            dist_list = get_Distance_list(temp_aligment)
            
            # get average pairwise distance for the given locus
            if not np.isnan(dist_list).all():
                # append to final list
                per_locus_len.append(temp_aligment.get_alignment_length())
                per_locus_dist.append(np.nanmean(dist_list))

        ## ADD FEATURE FOR PHASED HAPLOID

    # calculate the locus length weigthed within population average
    try:
        return np.average(per_locus_dist, weights = per_locus_len)
    except:
        # this happens if only a single phased sequence is provided, then inter pop distance cannot be assessed
        sys.exit("AutoPriorError: Automatic inference of theta prior failed.\nProvide theta prior manually.")

# measure the average within population pairwise distance in a MSA
def distance_between_pop(alignment, indpop_dict, pop_l, pop_r):
    per_locus_dist = []
    per_locus_len  = []

    for locus in alignment:
        # add the sequences belonging to the current population to a temp aligment
        temp_aligment_l = get_multi_population_locus_alignment(locus, indpop_dict, pop_l)
        temp_aligment_r = get_multi_population_locus_alignment(locus, indpop_dict, pop_r)
        
        # the distance can only be calculated for more than 2 sequences
        if len(temp_aligment_l) >= 1 and len(temp_aligment_r) >= 1:
            
            dist_list = get_two_pop_distance_list(temp_aligment_l, temp_aligment_r)

            # get average pairwise distance for the given locus
            if not np.isnan(dist_list).all():
                # append to final list
                per_locus_dist.append(np.nanmean(dist_list))
                per_locus_len.append(temp_aligment_l.get_alignment_length())

        ## ADD FEATURE FOR PHASED HAPLOID

    # calculate the locus length weigthed within population average
    try:
        return np.average(per_locus_dist, weights = per_locus_len)
    except:
        # this happens if only a single phased sequence is provided, then inter pop distance cannot be assessed
        sys.exit("AutoPriorError: Automatic inference of tau prior failed.\nProvide tau prior manually.")     


# automatically generates the tau and theta prior lines of the BPP control file

def auto_prior(
        imapfile:       Filename,
        seqfile:        Filename,
        tree_newick:    NewickTree,
        tau_prior,     
        theta_prior
        ) ->            CfileParam:

    """
    'imapfile' is the path to the Imap file, which specifies the mapping of individuals to populations

    'seqfile' is the path to the MSA.

    This function generates tau and theta priors using the method implemented in Minimalist BPP by Prof. Bruce Rannala. 
    Both tau and theta priors are inverse gammas with a wide alpha parameter (3) the function calculates a suitable
    mean for these inverse gamme distributions by examining the distances in the alignment.

    Theta is calculated to give a mean which is the average of within population average pairwise distances. This
    value is expected to be a good estimate of the average effective population size. 

    Tau is estimated as the mean sequence distance between two sequences on opposide sides of the first split in the tree. 

    The final values are formatted to comply with the "tauprior" and "thetaprior" lines of the BPP control file
    """

    alignment = alignfile_to_MSA(seqfile)
    
    indpop_dict = imapfile_read(imapfile, "indpop")
    populations = list(set(indpop_dict.values()))

    autopriors = {}

    ## THETA CALCULATION
    # calculation of locus length weighted average pairwise distances within each population
    if theta_prior is None:
        dist_pop = {}
        for population in populations:
            dst = distance_within_pop(alignment, indpop_dict, population)
            if dst != None:
                dist_pop[population] = dst
            
        # final theta calculation    
        D = np.average(list(dist_pop.values()))
        theta_alpha = 3
        theta_beta = np.round(2*D, decimals = 4)
        
        # write to priors dict
        autopriors['thetaprior'] = f"{theta_alpha} {theta_beta} e"
    else:
        autopriors["thetaprior"] = None
    
    ## TAU CALCULATION
    if tau_prior is None:
        pop_l, pop_r = get_first_split_populations(tree_newick)
        M = distance_between_pop(alignment, indpop_dict, pop_l, pop_r)
        
        tau_alpha = 3
        tau_beta = np.round(2*M, decimals = 4)
        
        # write to priors dict
        autopriors['tauprior'] = f"{tau_alpha} {tau_beta}"
    else:
        autopriors["tauprior"] = None

    return autopriors