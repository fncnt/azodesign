from RNAblueprint import DependencyGraphMT
import RNA
import subprocess
from concurrent.futures import ProcessPoolExecutor
from random import random


def adaptivewalk(targets = [], seqconstraint = None, strconstraint = None, stop_at_score = 0.1, acceptance_prob = 0.001, startseq = None):
    dgraph = DependencyGraphMT(targets, seqconstraint)
    
    if startseq:
        dgraph.set_sequence(startseq)
        
    # initial scores
    best_dist = len(targets[0])
    best_score = 1
    
    best_seq = None

    while True:
        dgraph.sample_clocal()
        seq = dgraph.get_sequence()
        
        fc = RNA.fold_compound(seq)

        if strconstraint:
            fc.hc_add_from_db(strconstraint)

        ## required to compute ensemble defect
        fc.pf()
        fc.bpp()

        structure, mfe = fc.mfe()

        # this corresponds to the constrained approach for (len(targets) == 1 and strconstraint is not None)
        # and otherwise the alternative approach
        score = max([fc.ensemble_defect(target) for target in targets])
        
        dists = [RNA.bp_distance(structure, target) for target in targets]
        min_dist = min(dists)

        if (score < best_score or min_dist <= best_dist) or random() < acceptance_prob:
        # alternative approach corresponds to:
        #if (score < best_score) or random() < acceptance_prob:
            best_score = score
            best_seq = seq
            # only constrained approach uses base pair distance as secondary objective
            # although it would not make a huge difference
            if seqconstraint:
                best_dist = min_dist
            
        else:
            dgraph.revert_sequence()

        # Stop conditions
        if best_score <= stop_at_score or min_dist <= 2:
            break

    return best_score, best_seq


# Returns list of futures
# calling result() returns value of func
def generate_sequences(jobs=1, n_seq=1, func=None, *args, **kwargs):
    processpool = ProcessPoolExecutor(jobs)
    sequence_futures = []

    for i in range(0, n_seq):
        sequence_futures.append(processpool.submit(func, *args, **kwargs))
    
    return sequence_futures


# executing this file will read a fasta file and assign random UUIDs to its entries
# if header line is empty.
if __name__ == "__main__":
    from uuid import uuid4
    from Bio import SeqIO, Seq
    import sys

    seqfile = sys.argv[1]

    records = list(SeqIO.parse(seqfile, "fasta"))
    for record in records:
        if not record.id:
            record.id = str(uuid4())
        
    SeqIO.write(records, seqfile, "fasta-2line")