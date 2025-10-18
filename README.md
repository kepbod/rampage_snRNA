# Pipeline to identify expressed snRNAs using RAMPAGE

## A schematic flow shows the pipeline

![workflow](https://github.com/kepbod/rampage_snRNA/blob/master/workflow.jpeg)

## Prerequisites

### Softwares

* [F-seq](http://fureylab.web.unc.edu/software/fseq/)

### Python libraries 

* [docopt](http://docopt.org/)
* [seqlib](https://github.com/kepbod/seqlib)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [pybedtools](https://daler.github.io/pybedtools/)
* [numpy](http://www.numpy.org/)
* [joblib](https://joblib.readthedocs.io/en/latest/)

## Usage

### Step 1:

Fetch proper read pairs and remove PCR dunplicates.

```
Usage: rm_pcr.py [options] <rampage>...

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a snRNA                       snRNA annotations (BED format).
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
```

* Inputs: BAM files of RAMPAGE (`<rampage>...`)
* Output: A output folder containing relevant files (`-o OUTPUT`)
    * `rampage_plus_5end.bed`: BED file of the 5' end of plus strand read pairs
    * `rampage_minus_5end.bed`: BED file of the 5' end of minus strand read pairs
    * `rampage_link.bed`: BED file linking the 5' and 3' ends of read pairs

Note:
1. If there are multiple RAMPAGE BAM files (different replicates) derived from the same samples, you could simply list them afterwards.
2. You could run with multiple threads using `-p THREAD`.

Example: 

```
rm_pcr.py -a snRNA_annotations.bed -o rampage_snRNA_peak rampage_rep1.bam rampage_rep2.bam
```

### Step 2:

Call peaks using 5' end of RAMPAGE read pairs

```
Usage: call_peak.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    -p PERCENT                     Retained percent of reads in resized peaks.
                                   [default: 0.7]
```

* Input: the output folder created by `rm_pcr.py`
* Output: `rampage_peaks.txt` under the input folder

Format of `rampage_peaks.txt`:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome                    |
| Start_Fseq  | Start of F-seq peak region    |
| End_Fseq    | End of F-seq peak region      |
| Name        | peak                          |
| Score       | 0                             |
| Strand      | Strand of peak                |
| Peak        | peak site                     |
| Height      | Height of peak site           |
| Peak_reads  | Reads of (peak site ± 2 bp)   |
| Total       | Total reads of peak region    |
| Start       | Start of peak region          |
| End         | End of peak region            |


Note:

1. You could set feature length for F-seq peak calling using `-l LENGTH`.
3. You could run with multiple threads using `-p THREAD`.

Example: 

```
call_peak.py rampage_snRNA_peak
```

### Step 3:

Assign multiple mapping reads using EM algorithm

```
Usage: assign_read.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    --tol TOL                      Termination for EM [default: 1e-9].
    --maxiter ITER                 Maximum iteration steps [default: 2000].
```

* Input: the output folder created by `call_peak.py`
* Output: `rampage_peak_entropy.txt` under the input folder

Format of `rampage_peak_entropy.txt`:


| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome                    |
| Start       | Start of peak region          |
| End         | End of peak region            |
| Name        | peak                          |
| Score       | 0                             |
| Strand      | Strand of peak                |
| Peak        | peak site                     |
| Start_Fseq  | Start of F-seq peak region    |
| End_Fseq    | End of F-seq peak region      |
| Total       | Total reads of peak region    |
| Entropy     | Entropy of RAMPAGE peak       |
| Positions of 3' ends | 3' ends of read pairs in the peak| 
|Read counts of each 3' end|Read counts of each 3' end  in the peak| 

Note:

1. You could set the termination cutoff for EM algorithm using `--tol TOL`.
2. You could set the maximux iteration steps of EM algorithm using `--maxiter ITER`.

Example: 

```
entropy.py rampage_snRNA_peak
```

### Step 4:

Annotate expressed snRNAs

```
Usage: annotate_snRNA.py [options] -a snRNA <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a snRNA                       snRNA annotations (BED format).
    --extend length                snRNA extended length. [default: 50]
    --span span                    span cutoff. [default: 1000]
    --coverage coverage            Coverage cutoff. [default: 0.5]
    -o out                         Output file. [default: snRNA_peak.txt]
```

* Input: the output folder created by `assign_read.py`
* Output: `snRNA_peak.txt` 

Format of `snRNA_peak.txt`:

The first thirteen columns are the same as `rampage_peak_entropy.txt`.

The additional six columns are listed below:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome of snRNA           |
| Start       | Start of snRNA                |
| End         | End of snRNA                  |
| Name        | Name of snRNA                 |
| Score       | 0                             |
| Strand      | Strand of snRNA               |

Note:

1. You could set the snRNA annotation extension length using `--extend length`.
2. You could set RAMPAGE effective length cutoff using `--span span`.
3. You could set coverage cutoff of snRNA annotation using `--coverage coverage`.

Example: 

```
annotate_snRNA.py -a snRNA_annotations.bed rampage_snRNA_peak
```

## License

Copyright (C) Xiao-Ou Zhang. See the [LICENSE](https://github.com/kepbod/rampage_te_tss/blob/master/LICENSE) file for license rights and limitations (MIT).