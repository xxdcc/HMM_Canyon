# HMM_Canyon

## User-friendly tools for exploring DNA methylation Canyon
HMM_Canyon addresses the challenge of de novo identifying large under methylated region (Canyons) through WGBS data. HMM_canyons contains useful modules to pre-process the mapped methylation profile data from BSMAP or Bismark, reporting the Canyons (or UMRs), which is useful for discovery of novel epigenetic features for DNA methylation.

For support, questions, or feature requests contact: wl1@bcm.edu

## Citation
Jeong M, Sun D, Luo M, Huang Y, Challen GA, Rodriguez B, Zhang X, Chavez L, Wang H, Hannah R, Kim SB, Yang L, Ko M, Chen R, Göttgens B, Lee JS, Gunaratne P, Godley LA, Darlington GJ, Rao A, Li W, Goodell MA. [Large conserved domains of low DNA methylation maintained by Dnmt3a](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3920905/)

## Installatioin
```C++
git clone https://github.com/xxdcc/HMM_Canyon.git
cd HMM_Canyon
g++ -o oumr hmm.cpp oumr.cpp  -Wall -I/usr/local/include -I/path_for_boost/boost/1.55.0/ -g -O3 -L/spath_for_boost/boost/1.55.0/stage/lib/ -lboost_program_options
```
## Data Preparation
Please use the following command to convert the BSMAP or Bismark mapped methylation profile to the format for HMM_Canyon.
```
python methratio2mcall.py -o output_reformated_file input_methylation_file
```
It should look like this toy sample 

chrom|start|end|ratio|totalC|methC|strand|next|Plus|totalC.1|methC.1|Minus|totalC.2|methC.2|localSeq
---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
chr1|10468|10470|0.0|3.0|0|-|G|+|0.0|0|-|3.0|0|CG
chr1|10471|10473|0.0|1.0|0|+|G|+|1.0|0|-|0.0|0|CG
chr1|10484|10486|0.332|3.0|1|B|G|+|1.0|1|-|2.0|0|CG
chr1|10489|10491|0.833|6.0|5|B|G|+|1.0|1|-|5.0|4|CG
chr1|10493|10495|0.571|7.0|4|B|G|+|1.0|0|-|6.0|4|CG
chr1|10497|10499|0.6|10.0|6|B|G|+|2.0|0|-|8.0|6|CG
chr1|10525|10527|0.846|13.0|11|B|G|+|4.0|3|-|9.0|8|CG
chr1|10542|10544|0.667|15.0|10|B|G|+|6.0|2|-|9.0|8|CG
chr1|10563|10565|0.818|11.0|9|B|G|+|2.0|1|-|9.0|8|CG

## Usage
```
installation_path/canyon -m output_reformated_file -o output_canyon_file
```

