# Designing RNA Sequences With *Azoarcus*-like Secondary Structures

This is an implementation of a pipeline for sequence design given a target structure with simple pseudoknots, without explicitly modelling pseudoknots in the design process.

## Prerequisites

Besides [packages installable via `pip`](requirements.txt), the following is required:

- [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/) (Python bindings and `RNAPKplex`)
- [`pKiss`](https://github.com/jlab/fold-grammars)
- [`RNAblueprint`](https://github.com/ViennaRNA/RNAblueprint)
- the proof-of-concept implementation [`libbpdist`](https://github.com/fncnt/librna-sys)

## Caveats

Although the commandline options of the scripts are documented and accessible via `-h`, usage is not always ergonomic.

## Usage & Structure

### Design Pipeline

The design pipeline is implemented in [`azoarcus_design.py`](azoarcus_design.py) (and [`pipeline.py`](utils/pipeline.py)).
Per default, its output is in `FASTA` format printed to `stdout`.
To compute properties of designed sequences, use `-i`, specify reference data via `-S`, `-T`.

#### Constrained Approach

This is the default design approach.
Example usage:
```sh
python3 azoarcus_design.py -n 10 -j 4 -s 0.05 -T <target structure> -C <sequence constraints> > sequencedesigns.fasta
python3 azoarcus_design.py -j 4 -s -i sequencedesigns.fasta -T <target structure> -S <reference sequence> > sequencedesigns.csv
```
`-s 0.05` adjusts the stopping threshold for the objective function.
Structural constraints are extracted from pseudoknots in the target structure; those should be reflected in sequence constraints specified via `-C`.

Example input data:
```sh
GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG # reference sequence
...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).. # target structure
GUGNCNNNNNNNNNGAAANNNNNNNNGNNANNNNNNCNAAUNCGNCNNNNNCUAAGNNNNNNNNNNNNNNUAUGNNNNNGNCGNNCCANNNNNNNNNNNNNNNNNNNNNNNNGGNGUAGAGACUANNNGNNNNNNNNCUAAGNNNNNNNNUAUGNNNNNNNCAUAGUCCNNNNNNNNNNGAAANNNNNNNNNNNNNG # complete sequence constraints
GUGNNNNNNNNNNNGAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAGACUANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAGUCCNNNNNNNNNNNNNNNNNNNNNNNNNNNG # minimal sequence constraints
```

#### Alternative Approach

The alternative approach does not require structural constraints but two target structures obtained by *splitting* a pseudoknotted structure.
This variant of the pipeline is toggled using `-P`.
Example usage:
```sh
python3 azoarcus_design.py -P -n 10 -j 4 -s 0.15 -T <first split structure> -T <second split structure> -C <sequence constraints> > sequencedesigns.fasta
python3 azoarcus_design.py -P -j 4 -s -i sequencedesigns.fasta -T <original target structure> -T <first split structure> -T <second split structure> -S <reference sequence> > sequencedesigns.csv
```
`-s 0.15` is a higher threshold than before (the reason is illustrated in [`maxned_threshold.py`](maxned_threshold.py)).
Structural constraints are extracted from pseudoknots in the target structure; those should be reflected in sequence constraints specified via `-C`.

Example input data:
```sh
GUGCCUUGCGCCGGGAAACCACGCAAGGGAUGGUGUCAAAUUCGGCGAAACCUAAGCGCCCGCCCGGGCGUAUGGCAACGCCGAGCCAAGCUUCGGCGCCUGCGCCGAUGAAGGUGUAGAGACUAGACGGCACCCACCUAAGGCAAACGCUAUGGUGAAGGCAUAGUCCAGGGAGUGGCGAAAGUCACACAAACCGG # reference sequence
...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...[.[[[[[...))))))((((...(((....)))..))))......]]]]]]..((.(((((....))))).....)).. # target structure
...(((((((..((....)).)))))))...((((((....((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..))).............))))))((((...(((....)))..))))..............((.(((((....))))).....)).. # first split structure
...(((((((..((....)).))))))).............((((((...((...((((((....))))))..))...))))))(((...(.((((((....)))))).)..)))...(.(((((.........((((...(((....)))..))))......))))))..((.(((((....))))).....)).. # second split structure
GUGNCNNNNNNNNNGAAANNNNNNNNGNNANNNNNNCNAAUNCGNCNNNNNCUAAGNNNNNNNNNNNNNNUAUGNNNNNGNCGNNCCANNNNNNNNNNNNNNNNNNNNNNNNGGNGUANNNNNNNNNNGNNNNNNNNCUAAGNNNNNNNNUAUGNNNNNNNCANNNNNNNNNNNNNNNNGAAANNNNNNNNNNNNNG # complete-alt sequence constraints
GUGNNNNNNNNNNNGAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNG # minimal-alt sequence constraints
```

### Neutral Paths

[`neutralpaths.py`](neutralpaths.py) computes neutral path lengths using `pKiss` structure prediction.

Although not required for neutral paths, the option`-T` computes the expected Hamming distance between
two random sequences compatible to a shared structure, by both uniform sampling and an estimation using the number of paired and unpaired positions.

### Energy Parameter Sets

[`eparams.py`](eparams.py) may be used to assess different energy parameter sets (in format compatible to `ViennaRNA`) using `RNAfold`, `pKiss`and `RNAPKplex`.
Note that the parameter files are not provided here but are available in `ViennaRNA` and [here (requires conversion using `RNAparconv`)](https://www.cs.ubc.ca/labs/beta/Projects/RNA-Params/).
For this, `pKiss` was used at version `2.2.12` and [`RNAPKplex` was patched](https://github.com/ViennaRNA/ViennaRNA/compare/77861405002d93a35cec3f615e2d1a5d210964d8...7a7e84ae8f6954dff43cc31d581b0bcc63b8a1e1) (as of `ViennaRNA 2.4.18`, `RNAPKplex` seems to work better).

### Dot Plots

[`dotplot.py`](dotplot.py) was used to produce dotplots of base pair probability matrices and single secondary structures.
