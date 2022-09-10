# Bioinformatics-Algorithms
Implementation of multiple bio algorithms for aligning and matching RNA and DNA sequences

## 1. STAR algorithm MSA
Implementation of multiple alignments with the help of the Star-Alignment algorithm and improving the alignments in a block-based manner.

1. Compute pairwise similarities (create a distance matrix using multiple correlations between pairs)
2. Select center $s_c$ that maximizes $\sum_{i!=s}{S(s_c,s_i)}$
3. Add sequences in decreasing order of similarity to center $s_c$
4. Score this multiple alignment
5. Find the blocks that have the potential to be improved, and with the help of the alignment from the previous few steps, get the alignment and replace the block. If the alignment score improves, replace this block permanently.

Sample input:
```
4 
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

## 2. PSSM profile
Building a PSSM profile of an MSA and using that profile to find the most similar part of a long sequence.

Sample input:
```
4
HVLIP
H-MIP
HVL-P
LVLIP
LIVPHHVPIPVLVIHPVLPPHIVLHHIHVHIHLPVLHIVHHLVIHLHPIVL
```
Sample output:
```
H-L-P
```

## 3. Semi-Global Alignment
Implementation of a Semi-Global alignment that does not consider all the gaps in the first and last alignment in scoring. The scoring is calculated usning PAM50 matrix and the help of Dynamic programming.

Sample input:
```
HEAGAWGHE
PAWHEA
```

Sample output:
```
20
HEAGAWGHE-
---PAW-HEA
```
