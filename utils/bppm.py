import RNA
import numpy as np


# does not handle bulges
# this is a helper for for dotplots in order to display colored stacks on the axes
def get_stackpositions(structure, brackets = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')], names=None):
    
    stacks = []

    for op, cl in brackets:

        counter_op = 0
        counter_cl = 0
        stack = []

        i = 0
        for pos in structure:
            i += 1

            if pos == op:
                counter_op += 1
            elif counter_op > 0:
                stack.append((i-counter_op, counter_op))
                counter_op = 0
            if pos == cl:
                counter_cl += 1
            elif counter_cl > 0:
                ops = stack.pop()
                stacks.append((ops, (i-counter_cl, counter_cl), "P" + names[ops[0]-1]))
                counter_cl = 0
        
    stacks.sort(key=lambda stack: stack[0][0])

    return stacks

#helper to extract pairs from pairtable
def plist(structure):
    if RNA.__version__ < "2.4.18":
        # prior to ViennaRNA 2.4.18:
        pt = RNA.ptable_from_string(structure, RNA.BRACKETS_ANY)
    else:
        # as of 2.4.18:
        pt = RNA.ptable(structure, RNA.BRACKETS_ANY)

    pairs = []

    for cl in pt:
        if pt[cl] < cl and pt[cl] != 0:
            pairs.append((pt[cl], cl))

    return pairs

# compute base pair probability matrix using ViennaRNA or NUPACK
def get_bppm(sequence, usenupack=False, rid=None):
    if usenupack:
        import nupack
        mdl = nupack.Model(material="rna")

        bppm = nupack.pairs(strands=[sequence], model=mdl).to_array()

        bppm = np.insert(bppm, 0, 0, axis=0)
        bppm = np.insert(bppm, 0, 0, axis=1)

    else:
        fc = RNA.fold_compound(sequence)
        fc.pf()
        bppm = np.asarray(fc.bpp())
    return bppm, rid

# interpret dot-bracket notation as base pair probability matrx
def dbn_to_bppm(structure):
    n = len(structure) + 1
    mat = np.zeros((n, n))

    pl = plist(structure)

    for pi, pj in pl:
        mat[pj, pi] = 1.0
        mat[pi, pj] = 1.0

    return mat