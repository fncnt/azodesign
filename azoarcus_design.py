#!/bin/env python3

import RNA
import sys
import argparse
from Bio import SeqIO, Seq
from concurrent.futures import as_completed, ProcessPoolExecutor
from numpy import argmin
from tqdm import tqdm
from uuid import uuid4 # random UUIDs

from utils import pipeline, pseudoknots

# This currently relies on the example from https://github.com/fncnt/librna-sys
try:
    from libbpdist import bp_distance_pk
except ImportError as err:
    bp_distance_pk = RNA.bp_distance
    print(f"WARNING: {err}!", file=sys.stderr)
    print("WARNING: Falling back to base pair distance without pseudoknots.", file=sys.stderr)



NATIVE_GISSD = "GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG"
TARGET_GISSD = "...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).."
SPLIT_GISSD = ["...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..))).............))))))((((...(((....)))..))))..............((.(((((....))))).....))..",
               "...(((((((..((....)).))))))).............((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...(.(((((.........((((...(((....)))..))))......))))))..((.(((((....))))).....)).."
              ]
CONSTRAINTS_COMPLETE = "GUGNCNNNNNNNNNGAAANNNNNNNNGNNANNNNNNCNAAUNCGNCNNNNNCUAAGNNNNNNNNNNNNNNUAUGNNNNNGNCGNNCCANNNNNNNNNNNNNNNNNNNNNNNNGGNGUAGAGACUANNNGNNNNNNNNCUAAGNNNNNNNNUAUGNNNNNNNCAUAGUCCNNNNNNNNNNGAAANNNNNNNNNNNNNG"
CONSTRAINTS_MIN = "GUGNNNNNNNNNNNGAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAGACUANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAGUCCNNNNNNNNNNNNNNNNNNNNNNNNNNNG"
CONSTRAINTS_COMPLETE_ALT = "GUGNCNNNNNNNNNGAAANNNNNNNNGNNANNNNNNCNAAUNCGNCNNNNNCUAAGNNNNNNNNNNNNNNUAUGNNNNNGNCGNNCCANNNNNNNNNNNNNNNNNNNNNNNNGGNGUANNNNNNNNNNGNNNNNNNNCUAAGNNNNNNNNUAUGNNNNNNNCANNNNNNNNNNNNNNNNGAAANNNNNNNNNNNNNG"
CONSTRAINTS_MIN_ALT = "GUGNNNNNNNNNNNGAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNG"

def stats(record, targets, hc, native=None, pktarget=TARGET_GISSD):
    seq = str(record.seq)

    fc = RNA.fold_compound(seq)
    if hc:
        fc.hc_add_from_db(hc)
    fc.pf()
    fc.bpp()

    structure, mfe = fc.mfe()
    
    # covers both constrained and alternative approach
    ed = max([fc.ensemble_defect(t) for t in targets])

    # ViennaRNA (at least <=2.4.18) strips pseudoknots automatically
    bpdist = RNA.bp_distance(pktarget, structure)
    eos = fc.eval_structure(pktarget)

    if native:
        hm = RNA.hamming_distance(seq, native)
    else:
        hm = 0

    pkplex = pseudoknots.RNAPKplex(seq)[0]
    bp_pkplex = bp_distance_pk(pkplex[1], pktarget)
    pkiss = pseudoknots.pKiss(seq)
    bp_pkiss = bp_distance_pk(pkiss[1], pktarget)

    _stats = (record.id, seq, structure, bpdist, mfe, eos, ed, hm, *pkplex, bp_pkplex, *pkiss, bp_pkiss)
    return _stats

def main():
    parser = argparse.ArgumentParser(description="Design Azoarcus-like sequences or compute properties of already designed sequences.")
    parser.add_argument("-i", dest="sequence_file", 
                        metavar="sequences.fasta",
                        type=argparse.FileType('r'),
                        help=("fasta file containing multiple sequences to be evaluated. "), 
                        default=None)
    parser.add_argument("-n", dest="number_of_seqs", metavar="1",
                        type=int,
                        help=("Number of sequences to design. "
                              "Ignored if \'-i\' is used."),
                        default=1)
    parser.add_argument("-j", dest="number_of_jobs", metavar="1",
                        type=int,
                        help=("Number of jobs for parallel design runs. "
                              "If \'-j\' is larger than \'-n\', the latter will be used instead."),
                        default=2)
    parser.add_argument("-r", "--relax-acceptance", dest="relax", metavar="0.001",
                        type=float,
                        help=("Relax acceptance of design candidates at each iteration; Candidates get accepted with this probability, regardless whether they minimize the objective function."),
                        default=0.001)
    parser.add_argument("-s", "--stop-at", dest="stop_at", metavar="0.1",
                        type=float,
                        help=("Threshold at which score minimization should stop."),
                        default=0.1)
    parser.add_argument("-S", "--start-seq", dest="start_seq", metavar="sequence",
                        action="store",
                        help=("Start adaptive walk at specified sequence instead of random sequence."),
                        default=None)
    parser.add_argument("-T", "--target-structure", dest="targets", metavar="target",
                        action="append",
                        help=("Target structure to design for. Can be used multiple times for multiple targets."),
                        )
    parser.add_argument("-C", "--sequence-constraint", dest="seq_constraint", metavar="constraint",
                        action="store",
                        nargs = "?",
                        help=("String containing nucleotide identities to preserve. Valid characters: A, C, G, U, N, where N denotes unconstrained positions."),
                        default="",
                        const=CONSTRAINTS_COMPLETE)
    parser.add_argument("-P", "--no-PK-constraint", dest="pk_constraint",
                        action="store_false",
                        help=("Alternative approach without structural constraints: Use this if specified sequence constraints do not include positions of a pseudoknot."))
    

    args = parser.parse_args()

    # compute properties for sequences in fasta file
    # this assumes one target with pseudoknot specified via -T
    if args.sequence_file:

        indexedseqs = SeqIO.index(args.sequence_file.name, "fasta")
        totalseqs = len(indexedseqs)
        indexedseqs.close()

        # separate "split" target structures from pseudoknotted target structure
        targets = []
        pktargets = []

        for structure in args.targets:
            _, hc = pseudoknots.hc_from_db(structure)

            if hc is None:
                targets.append(structure)
            else:
                pktargets.append(structure)

        # at least warn that a PK target is required for stats()
        # better solution would be merging split target structures
        if len(pktargets) != 1:
            print(f"WARNING: -i requires 1 target structure with pseudoknots!", file=sys.stderr)

        _, hc = pseudoknots.hc_from_db(pktargets[0])

        if not args.pk_constraint:
            hc = None

        # use -S to specify reference sequence;
        # falls back to NATIVE_GISSD
        if args.start_seq:
            nativeseq = args.start_seq
        else:
            nativeseq = NATIVE_GISSD


        stats_ppool = ProcessPoolExecutor(args.number_of_jobs)
        stats_futures = []
        for record in SeqIO.parse(args.sequence_file, "fasta"):
            stats_futures.append(stats_ppool.submit(stats, record, targets, hc, nativeseq, pktargets[0]))

        with tqdm(as_completed(stats_futures), total=totalseqs) as processing:
            print("ID; Sequence; Structure; BP (target); MFE; EOS (target); Score (ned|maxned); Hamming (native); MFE (PKplex); Structure (PKplex); BP (PKplex); MFE (pKiss); Structure (pKiss); BP (pKiss)")

            for completed in processing:
                try:
                    completed_stats = completed.result()
                except Exception as e:
                    print(e)
                else:
                    print(*completed_stats, sep="; ")

    else:
        # this starts the design pipeline

        if args.number_of_jobs > args.number_of_seqs:
            jobs = args.number_of_seqs
        else:
            jobs = args.number_of_jobs


        pk_free_target, pk_hc = pseudoknots.hc_from_db(args.targets[0])

        if not args.pk_constraint:
            pk_hc = None

        sequence_futures = pipeline.generate_sequences(jobs, args.number_of_seqs, pipeline.adaptivewalk,
            args.targets,
            args.seq_constraint,
            pk_hc,
            args.stop_at,
            args.relax,
            args.start_seq)
        
        with tqdm(as_completed(sequence_futures), total=args.number_of_seqs) as processing:

            for completed in processing:
                score, seq = completed.result()

                fc = RNA.fold_compound(seq)

                if pk_hc:
                    fc.hc_add_from_db(pk_hc)

                fc.pf()
                fc.bpp()
                structure, mfe = fc.mfe()

                print(f">{uuid4()}", seq, sep="\n")


if __name__ == "__main__":
    main()
