from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import os
from tqdm import tqdm
from uuid import uuid4

# single-purpose wrapper function for pKiss
def pKiss(sequence, mode = "mfe", strategy = "P", penalty = 9.8, params=None):
    cmd = ["pKiss", "--mode", mode, "--strategy", strategy, "--Hpenalty", str(penalty)]

    if params:
        cmd.append("--param")
        cmd.append(params)

    cmd.append(sequence)

    completed = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf8")

    output = completed.stdout.split("\n")[2].split("  ")
    # return score (EOS), structure
    return float(output[0]), output[1] 

def RNAPKplex(sequence, pk_penalty = 0.0, subopt = 0.0, params = None, prefix = ""):
    indata = "\n".join([sequence, "@"])

    cmd = [f"{prefix}RNAPKplex", "-e", str(pk_penalty), "-s", str(subopt)]

    if params:
        cmd.append("-P")
        cmd.append(params)

    completed = subprocess.run(cmd, input=indata, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf8")

    output = [line.split(" (") for line in completed.stdout.split("\n")[1:-1]]
    
    return [(float(structure[1].strip("() ")), structure[0]) for structure in output]

def hc_from_db(structure):
    constraint = structure
    pk_free = structure

    for bracket in "[{<]}>":
        constraint = constraint.replace(bracket, "x")
        pk_free = pk_free.replace(bracket, ".")

    for paren in "()":
        constraint = constraint.replace(paren, ".")

    if not "x" in constraint:
        constraint = None

    return pk_free, constraint
