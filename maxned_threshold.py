#!/bin/env python3

import dotplot as dplot
import RNA
import numpy as np

pkstruct = "...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).."
pksplit = ["...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..))).............))))))((((...(((....)))..))))..............((.(((((....))))).....))..",
           "...(((((((..((....)).))))))).............((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...(.(((((.........((((...(((....)))..))))......))))))..((.(((((....))))).....)).."]

# AFAIK, ViennaRNA does not provide python bindings to use its ensemble defect
# implementation with a custom base pair probability matrix.
# Instead, this is a small function that does this.
# You may obtain a BPP matrix using dbn_to_bppm() from dotplot.py
def ensemble_defect(structure, bppmatrix):
    n, _ = bppmatrix.shape

    if not n == len(structure) + 1:
        return None
    else:
        ed = n - 1
        pt = RNA.ptable(structure, RNA.BRACKETS_ANY)

        for i in range(1, n):
                if pt[i] == 0:
                    # position i of target structure is unpaired.
                    # dplot.dbn_to_bppm() does not contain unpaired probabilities;
                    # compute it from row/column sum.
                    unpaired = 1 - np.sum(bppmatrix[i])
                    ed -= unpaired
                else:
                    # position i of target structure is paired
                    ed -= bppmatrix[i, pt[i]]
        
        return ed

# See Section 4.2.1, equation (10) for the definition of
# "maximum (normalized) ensemble defect" (maxNED) of two targets.
#
# CONTEXT:
# Minimizing maxNED is approximately equivalent to 
# minimizing the NED for two target structures simultaneously.
# This introduces an inherent conflict to the minimization as 
# sooner or later, minimizing for either one of the targets
# increases the NED of the other one, depending on how similar they are.
# Therefore, maxNED cannot reach zero. The minimum possible value 
# of maxNED depends on the similarity of both target structures
# (and the base pair probability matrix of the design candidate sequence).
#
# This code estimates the minimum possible maxNED value for the specific
# example of the Azoarcus target structure.
#
# A few scenarios are modelled under the assumption that the design
# candidate sequence has only a single structure in its Boltzmann ensemble.
# This essentially corresponds to NED=0 if this structure is the target.
#
# 1. BPP matrix computed using partition function algorithm with pseudoknots
#       - all non-zero entries of the matrix are equal to one.
#       - single target structure with pseudoknot: pkstruct
#       - (this is the least relevant for the design approach taken)
#
# 2. BPP matrix computed using partition function algorithm w/o pseudoknots
#   a.  - all non-zero entries are equal to one
#       - BPP matrix corresponds to either one structure in pksplit
#   b.  - all non-zero entries are equal to one except entries corresponding
#         to conflicting pairs (here: P3 and P7)
#       - the latter are equal to 0.5 
#         (for simplicity; fortunately, P3 and P7 are equal in number of pairs)  
#       - i.e. BPP matrix corresponds to superposition of pksplit   
# 
# The minimum possible value of maxNED is equal to the base pair distance
# of both target structures divided by sequence length.
# However, I stopped the design process at ~ twice the base pair distance because
#   - reducing the ensemble to a single structure seemed to be counter-productive
#     (similar to constrained approach)
#   - designs close to either one of the targets were considered valuable.
if __name__ == "__main__":
    # scenario 1
    pkbppm = dplot.dbn_to_bppm(pkstruct)
    # scenario 2a
    bppms = [dplot.dbn_to_bppm(structure) for structure in pksplit]
    # scenario 2b
    spbppm = (bppms[0] + bppms[1])/2

    # annotation
    scenarios = ["scenario 1  ", "scenario 2a ", "scenarion 2a'", "scenario 2b "]
    targets = ["(.)..(...)..", "(.)....[...]", "(.)..(.[.).]"]

    for scenario, bppm in zip(scenarios, [pkbppm, *bppms, spbppm]):
        print(scenario)
        print("target", "\tED", "NED", sep="\t")
        for target, structure in zip(targets, pksplit + [pkstruct]):
            ed = ensemble_defect(structure, bppm)
            print(target, ed, ed/len(structure), sep="\t")
