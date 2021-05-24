#!/bin/env python3

import RNA
from azoarcus import pseudoknots
import argparse
import sys

# This currently relies on the example from https://github.com/fncnt/librna-sys
try:
    from libbpdist import bp_distance_pk
except ImportError as err:
    bp_distance_pk = RNA.bp_distance
    print(f"WARNING: {err}!", file=sys.stderr)
    print("WARNING: Falling back to base pair distance without pseudoknots.", file=sys.stderr)


TARGET_GISSD = "...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).."
NATIVE_GISSD = "GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test different energy parameter sets with pKiss, RNAPKplex, RNAfold")
    parser.add_argument("-T", "--target-structure", dest="target", metavar="target",
                        action="store",
                        help=("Target structure to use."),
                        default=TARGET_GISSD)
    parser.add_argument("-S", "--sequence", dest="sequence", metavar="sequence",
                        action="store",
                        help=("Sequence to predict structure for"),
                        default=NATIVE_GISSD)
    parser.add_argument("-P", "--params", dest="params", metavar="turner2004.par",
                        nargs="+",
                        help=("Energy parameter files to use for prediction"),
                        )

    args = parser.parse_args()
    print("Parameters", "Structure", "Energy", "BP distance", "Tool", "Variant", sep="; ")
    for param in args.params:
        RNA.params_load(param)
        fc = RNA.fold_compound(args.sequence)
        structure, mfe = fc.mfe()
        print(param,
            structure,
            mfe,
            RNA.bp_distance(structure, args.target),
            "RNAfold",
            "default",
            sep="; ")

        _, hc = pseudoknots.hc_from_db(args.target)
        fc.hc_add_from_db(hc)

        structure, mfe = fc.mfe()
        print(param,
            structure,
            mfe,
            RNA.bp_distance(structure, args.target),
            "RNAfold",
            "modified",
            sep="; ")

        mfe, structure = pseudoknots.pKiss(args.sequence, penalty=9.0, params=param)
        print(param,
            structure,
            mfe,
            bp_distance_pk(structure, args.target),
            "pKiss",
            "default",
            sep="; ")

        mfe, structure = pseudoknots.pKiss(args.sequence, penalty=9.8, params=param)
        print(param,
            structure,
            mfe,
            bp_distance_pk(structure, args.target),
            "pKiss",
            "modified",
            sep="; ")

        mfestructure = pseudoknots.RNAPKplex(args.sequence, pk_penalty=8.1, params=param)[0]
        print(param,
            mfestructure[1],
            mfestructure[0],
            bp_distance_pk(mfestructure[1], args.target),
            "RNAPKplex",
            "default",
            sep="; ")

        mfestructure = pseudoknots.RNAPKplex(args.sequence, pk_penalty=0.0, params=param)[0]
        print(param,
            mfestructure[1],
            mfestructure[0],
            bp_distance_pk(mfestructure[1], args.target),
            "RNAPKplex",
            "modified",
            sep="; ")

    
    