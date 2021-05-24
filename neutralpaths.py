#!/bin/env python3

import RNA
import RNAblueprint as rbp
from azoarcus.pseudoknots import pKiss
import argparse
import math
from Bio import SeqIO, Seq
from concurrent.futures import as_completed, ProcessPoolExecutor
from tqdm import tqdm


TARGET_GISSD = "...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).."
NATIVE_GISSD = "GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG"

# this tests the estimated hamming distance between 
# two random sequences compatible to a common structure
# by globally sampling using RNAblueprint and 
# computing the relative error to the estimate
#
# estimate expected random hamming from components (i.e. conditioned to compatible sequences):
# 0.75 * unpaired +
# 5/6 * base-pairs (1/3 * (3/5 * 2 + 2/5) + 2/3 * (4/5 * 2 + 1/5)) ## = 13/18 * 2 * base-pairs
# 5 of 6 base-pairs are wrong, i.e. have at least one wrong nucleotide
# in 1/3 of those cases, the correct base-pair would've been GU or UG
# in which 3 of 5 wrong cases would have 2 incorrect nucleotides and 2 just one
# similarly in 2/3 correct would have been either, GC, CG, AU or UA
# in which 4/5 incorrect had 2 wrong nucleotides and 1/5 just 1.
# 
# G-C \       / C-G
#  |   |     |   |
# G-U   > = <   U-G
#  |   |     |   |
# A-U /       \ U-A
# | denotes one changed nucleotide,  }={ mirroring base-pairs causing 2 changed nucleotides
# for example, if A-U is correct. 5/6 are wrong, of which only 1/5 hase one wrong nucleotide: G-U
#
# sequence input is only used to produce a structure
def test_estimate(sequence, samples=1000):
    _, structure = pKiss(sequence)
    dg = rbp.DependencyGraphMT([structure])
    
    components = dg.number_of_connected_components()
    bpairs = len(sequence) - components
    unpaired = components - bpairs
    
    estimate = 3/4 * unpaired + 13/18 * 2 * bpairs

    avg = 0

    print(f"length: {len(sequence)}, components: {components}, unpaired: {unpaired}, pairs: {bpairs}, samples: {samples}\n")
    #print("i; hamming; relerror")

    for i in range(0, samples):
        # we're sampling globally here
        dg.sample()
        candidate1 = dg.get_sequence()
        dg.sample()
        candidate2 = dg.get_sequence()

        hamming = RNA.hamming_distance(candidate1, candidate2)
        avg += hamming
        #print(i, hamming, abs(hamming-estimate)/estimate, sep="; ")

    avg /= samples
    print(f"mean Hamming distance of sampled pairs: {avg}\n")

    print("expected Hamming distance estimated from components:")
    print(f"d = {estimate}, rel. error = {abs(avg-estimate)/estimate}\n")
    print("expected Hamming distance via 0.75*n:")
    print(f"d = {0.75*len(sequence)}, rel. error = {abs(avg-0.75*len(sequence))/(0.75*len(sequence))}")







def neutralpath(sequence, stop_at=1000, rid=None):
    _, structure = pKiss(sequence)

    dg = rbp.DependencyGraphMT([structure])
    # ideally I'd iterate over all compatible mutations instead of uniform sampling
    #compatible_mut_count = dg.sample_clocal()
    #
    # number of connected components = sequence length - number of base pairs in structure
    # makes more sense to compare against than sequence length (= max. hamming distance)?
    components = dg.number_of_connected_components()
    dg.set_sequence(sequence)
    dg.set_history_size(2)
    
    trials = 0
    neutral_length = 0

    # this would be better, but takes longer (of course)
    # and we're sampling uniformly random anyway
    #max_stop_at = max(compatible_mut_count, stop_at)
    max_stop_at = max(components, stop_at)

    # I'm essentially doing the same as Grüner 1996 with reduced strings
    # we're moving on 1-point and base-pair mutations
    # with uniform sampling we're missing neutral neighbors sooner or later    
    while trials < max_stop_at and neutral_length < len(sequence):
        dg.sample_clocal()
        new_seq = dg.get_sequence()
        new_hamming = RNA.hamming_distance(new_seq, sequence)

        # neutral paths are required to increase distance to starting point.
        if new_hamming > neutral_length:
            #print(f"{rid}:{step}<{max_stop_at}, {neutral_length}, {components}")

            _, new_structure = pKiss(new_seq)

            if new_structure == structure:
                # this matches graphical distance "==" hamming distance
                neutral_length = new_hamming
                # alternatively, if we were just counting steps
                #neutral_length += 1
                
                # reset counter for stop condition
                trials = 0
            else:
                # increase counter for stop condition and retry
                trials += 1
                dg.revert_sequence()
    
        else:
            dg.revert_sequence()

    return neutral_length, components, structure, dg.get_sequence(), rid



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute neutral paths for input sequences with pKiss predictions")
    parser.add_argument("-S", "--sequence", dest="sequence", metavar="sequence",
                        action="store",
                        help=("RNA sequence to produce neutral path for."),
                        default=NATIVE_GISSD)
    parser.add_argument("-t", "--threshold", dest="threshold",
                        type=int,
                        action="store",
                        help=("Threshold to specify when to stop trying. Values below number of compatible point or pair mutations will be ignored."),
                        default=0
                        )
    parser.add_argument("-T", "--test", dest="test",
                        type=int,
                        action="store",
                        nargs="?",
                        help=("test estimation of expected hamming distance of compatible sequences from number of components in structure dependency graph."),
                        const=100
                        )
    parser.add_argument("-i", dest="sequence_file", 
                        metavar="sequences.fasta",
                        type=argparse.FileType('r'),
                        help=("fasta file containing multiple sequences for batch processing"),
                        default=None)
    parser.add_argument("-j", dest="number_of_jobs", metavar="2",
                        type=int,
                        help=("Number of threads to use for batch processing of a .fasta file"),
                        default=2)


    args = parser.parse_args()

    if not args.test:
        if not args.sequence_file:
            length, components, structure, lastseq, _ = neutralpath(args.sequence, args.threshold)
            print(f"neutral path length = {length}")
            print(f"number mutable points/pairs: {components}")
            print(f"structure: {structure}")
            print(f"last sequence: {lastseq}")
        else:
            indexedseqs = SeqIO.index(args.sequence_file.name, "fasta")
            totalseqs = len(indexedseqs)
            indexedseqs.close()

            pool = ProcessPoolExecutor(args.number_of_jobs)
            futs = []
            for record in SeqIO.parse(args.sequence_file, "fasta"):
                futs.append(pool.submit(neutralpath, str(record.seq), stop_at=args.threshold, rid=record))

            print("ID; plength; components; relative; lastsequence; structure")
            with tqdm(as_completed(futs), total=totalseqs) as processing:
                for completed in processing:
                    try:
                        neutral_length, components, structure, lastseq, rec = completed.result()
                        print(rec.id, 
                            neutral_length, 
                            components, 
                            neutral_length/len(lastseq),
                            lastseq, 
                            structure, 
                            sep="; ")
                    except Exception as e:
                        print(e)
    else:
        test_estimate(args.sequence, args.test)
    

