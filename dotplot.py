#!/bin/env python3

import RNA
import argparse
from pathlib import Path
import math
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib import cm
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.gridspec import GridSpec
import numpy as np
from Bio import SeqIO, Seq
from concurrent.futures import as_completed, ProcessPoolExecutor
from tqdm import tqdm

from enum import Enum

class Highlight(Enum):
    NONE = 0
    LOWER = 1
    UPPER = 2
    BOTH = 3




TARGET_STRUCTURE_PK = "........(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((((..........)))))))).)..)))...[.[[[[[...))))))((((..((((....))).)))))......]]]]]].((..(((((....))))).....)).."
TARGET_GISSD = "...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).."
REGION_GISSD = "...2222222..22....22.2222222...333333....444444...55...555555....555555..55...444444666...6.666666....666666.6..666...7.77777...3333338888...888....888..8888......777777..99.99999....99999.....99.."
NATIVE_GISSD = "GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG"

HIGHLIGHTS_GISSD = ["......................................................................................................................(.(((((......................................))))))............................",
 "...............................((((((...........................................................................................))))))..............................................................."]

P7 = [[119,125], [164, 169]]
P7a = ((119, 1), (169, 1), "P7")
P7b = ((121, 5), (164, 5), "P7")
P3 = [[32, 37], [129, 134]]

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

def dotplot(lower, upper, nupack=False, clipat=1.0, highlight=Highlight.BOTH, legend=True, colorstacks=False):
    mat = np.zeros(lower.shape)
    n, _ = mat.shape

    for i in range(1, n):
        for j in range(i, n):
            if not i == j:
                mat[j, i] = lower[i, j]
                mat[i, j] = upper[i, j]

    patches = []
    if not highlight == Highlight.NONE:
        for h in [P3, P7]:
            centerx = (h[0][1] + h[0][0])/2
            centery = (h[1][1] + h[1][0])/2
            x = h[0][1] - h[0][0]
            y = h[1][1] - h[1][0]
            l = math.sqrt(x*x + y*y) + 5
            lower_patch = Ellipse((centerx, centery), l, 6, -45)
            upper_patch = Ellipse((centery, centerx), l, 6, -45)
            
            if highlight == Highlight.BOTH or highlight == Highlight.UPPER:
                patches.append(upper_patch)
            if highlight == Highlight.BOTH or highlight == Highlight.LOWER:
                patches.append(lower_patch)

    if colorstacks:
        stacks = get_stackpositions(args.target, [('(', ')')], REGION_GISSD)
        stacks.append(P7a)
        stacks.append(P7b)
        stack_patches = []
        stack_colors = {"P2": "#d70b17",
                        "P3": "#e95f0a",#"#d13910",
                        "P4": "#f0ca0d",
                        "P5": "#f0ca0d",
                        "P6": "#1d9131",
                        "P7": "#176da6",
                        "P8": "#1a2bb3",#"#28348e",
                        "P9": "#2a107b"}#"#352e8a"}
        for stack in stacks:
            xop = stack[0][0] - 0.5
            wop = stack[0][1]
            yop = -3
            hop = 3
            xcl = stack[1][0] - 0.5
            wcl = stack[1][1]
            ycl = -3
            hcl = 3
            stack_patches.append((stack[-1], Rectangle((xop, yop), wop, hop, clip_on=False)))
            stack_patches.append((stack[-1], Rectangle((n, xop), hop, wop, clip_on=False)))
            stack_patches.append((stack[-1], Rectangle((xcl, ycl), wcl, hcl, clip_on=False)))
            stack_patches.append((stack[-1], Rectangle((n, xcl), hcl, wcl, clip_on=False)))
        
    fig = plt.figure()
    mf = fig.add_subplot(111)
    mf.set_xlim([0, n])
    mf.set_ylim([n, 0])
    #mf.set_xlabel(r"$j$", labelpad=0.01, clip_on=False, fontweight="extra bold")
    #mf.set_ylabel(r"$i$", labelpad=0.01, clip_on=False, fontweight="extra bold")

    for patch in patches:
        patch.set_facecolor("none")
        patch.set_linewidth(1)
        patch.set_edgecolor("purple")
        patch.set_alpha(0.5)
        mf.add_patch(patch)

    if colorstacks:
        xbox = Rectangle((0,-3), n, 3, clip_on=False)
        ybox = Rectangle((n,0), 3, n, clip_on=False)
        xbox.set_facecolor("#dddddd")
        ybox.set_facecolor("#dddddd")
        mf.add_patch(xbox)
        mf.add_patch(ybox)

        for stp in stack_patches:
            ptch = stp[1]
            ptch.set_facecolor(stack_colors[stp[0]])
            ptch.set_edgecolor("none")
            ptch.set_alpha(0.7)
            mf.add_patch(ptch)


    if clipat > 0.0:
        cmap = cm.get_cmap("copper_r")
    else:
        cmap = cm.get_cmap("Greys")
    #cmap.set_bad("black", clipat)

    mat = np.ma.masked_where(mat < clipat, mat)

    img = mf.imshow(mat, cmap=cmap, norm=Normalize(vmin=clipat, vmax=1.0, clip=False), interpolation="nearest")
    mf.plot([n, 1], [n, 1], linestyle="--", color="grey", linewidth=0.5)
    mf.grid(True, color="grey", linestyle="--", linewidth=0.3)
    if legend:
        cb = fig.colorbar(img, shrink=0.4, location="right", pad=0.02)
        cb.ax.set_title(r"      $P_{i,j}$", pad=10)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce dot plots using partition functions from ViennaRNA or NUPACK")
    parser.add_argument("-T", "--target-structure", dest="target", metavar="target",
                        action="store",
                        help=("Target structure to display in lower diagonal half."),
                        default=TARGET_GISSD)
    parser.add_argument("-S", "--sequence", dest="sequence", metavar="sequence",
                        action="store",
                        help=("RNA sequence to compute partition function for (plotted in upper diagonal half)."),
                        default=NATIVE_GISSD)
    parser.add_argument("-s", "--is-structure", dest="is_structure",
                        action="store_true",
                        help=("this allows providing a structure with -S for exact comparisons.")
                        )
    parser.add_argument("-N", "--nupack", dest="nupack",
                        action="store_true",
                        help=("Use NUPACK instead of ViennaRNA for partition function computation.")
                        )
    parser.add_argument("-o", "--output", dest="outname",
                        action="store",
                        help=("Save to file.")
                        )
    parser.add_argument("-c", "--clip", dest="clip",
                        type=float,
                        action="store",
                        help=("clip values at this threshold so that small probabilities are visible against the target structure."),
                        default=0.0
                        )
    parser.add_argument("-i", dest="sequence_file", 
                        metavar="sequences.fasta",
                        type=argparse.FileType('r'),
                        help=("fasta file containing multiple sequences for batch processing"),
                        default=None)
    parser.add_argument("-z", "--id-contains-sequence", dest="idtarget",
                        action="store_true",
                        help=("fasta ID has format : '>ID:sequence'. This is compatible with -s.")
                        )
    parser.add_argument("-j", dest="number_of_jobs", metavar="1",
                        type=int,
                        help=("Number of threads to use when calculating partition functions of a .fasta file"),
                        default=2)
    parser.add_argument("-H", "--highlight", dest="highlight", metavar="NONE|BOTH|LOWER|UPPER",
                        action="store",
                        nargs='?',
                        help=("Highlight (hardcoded) positions."),
                        default = Highlight.NONE.name,
                        const = Highlight.BOTH.name,
                        )
    parser.add_argument("-d", "--dpi", dest="dpi", metavar="600",
                        type=int,
                        help=("matplotlib is a bit weird when embedding bitmaps into pdf files. increase this if your pixels look strange."),
                        default=600)
    parser.add_argument("-y", "--color-stacks", dest="colorstacks",
                        action="store_true",
                        help=("Use hardcoded colors to indicate native Azoarcus stacks")
                        )
    
    

    args = parser.parse_args()


    if not args.sequence_file:
        if not args.is_structure:
            bppm, _ = get_bppm(args.sequence, args.nupack)
        else:
            bppm = dbn_to_bppm(args.sequence)

        lbppm = dbn_to_bppm(args.target)

        dotplot(lbppm, bppm, args.nupack, args.clip, Highlight[args.highlight.upper()], not args.is_structure, args.colorstacks)

        if args.outname:
            plt.savefig(args.outname, dpi=args.dpi)
        else:
            plt.show()

            
    else:
        indexedseqs = SeqIO.index(args.sequence_file.name, "fasta")
        totalseqs = len(indexedseqs)
        indexedseqs.close()

        if not args.idtarget:
            lbppm = dbn_to_bppm(args.target)

        pool = ProcessPoolExecutor(args.number_of_jobs)
        futs = []
        for record in SeqIO.parse(args.sequence_file, "fasta"):
            futs.append(pool.submit(get_bppm, str(record.seq), args.nupack, record.id))

        with tqdm(as_completed(futs), total=totalseqs) as processing:
            for completed in processing:
                try:
                    bppm, rid = completed.result()

                    if args.idtarget:
                        spl = rid.split(":")
                        rid = spl[0]
                        target = spl[1]
                        
                        if args.is_structure:
                            lbppm = dbn_to_bppm(target)
                        else:
                            lbppm, _ = get_bppm(target, args.nupack)

                    dotplot(lbppm, bppm, args.nupack, args.clip, Highlight[args.highlight.upper()], not args.is_structure, args.colorstacks)

                    fpath = Path(args.outname)
                    
                    new_stem = rid + "-" + str(fpath.stem)
                    plt.savefig(fpath.with_stem(new_stem), dpi=args.dpi)

                except Exception as e:
                    print(e)
    

