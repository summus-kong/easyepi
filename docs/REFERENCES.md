# Primary method references

Cite the tools actually used in an analysis. These references support method choices and interpretation; they do not validate Easyepi's implementation or establish equivalence.

1. Buenrostro JD et al. Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nature Methods* (2013). https://doi.org/10.1038/nmeth.2688
2. Zhang Y et al. Model-based Analysis of ChIP-Seq (MACS). *Genome Biology* (2008). https://doi.org/10.1186/gb-2008-9-9-r137
3. Ernst J, Kellis M. ChromHMM: automating chromatin-state discovery and characterization. *Nature Methods* (2012). https://doi.org/10.1038/nmeth.1906
4. Bailey TL et al. MEME SUITE: tools for motif discovery and searching. *Nucleic Acids Research* (2009). https://doi.org/10.1093/nar/gkp335
5. Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif. *Bioinformatics* (2011). https://doi.org/10.1093/bioinformatics/btr064
6. Cock PJA et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* (2009). https://doi.org/10.1093/bioinformatics/btp163
7. Badia-i-Mompel P et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances* (2022). https://doi.org/10.1093/bioadv/vbac016
8. Schep AN et al. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nature Methods* (2017). https://doi.org/10.1038/nmeth.4401
9. Ramirez F et al. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research* (2016). https://doi.org/10.1093/nar/gkw257
10. Bentsen M et al. ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation. *Nature Communications* (2020). https://doi.org/10.1038/s41467-020-18035-1
11. Pedregosa F et al. Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research* 12, 2825-2830 (2011). https://jmlr.org/papers/v12/pedregosa11a.html
12. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Implementation references checked during development:

- [Official ChromHMM manual](https://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf): BinarizeBed positional inputs, peak mode, bin size, LearnModel seed and processors.
- decoupler 1.9.2 distributed Python source (`method_ulm.py`): `run_ulm(mat, net, source, target, weight, min_n, use_raw)` returns estimate and p-value DataFrames.
- [MACS3 callpeak documentation](https://macs3-project.github.io/MACS/docs/callpeak.html).
- [Biopython motif documentation](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html).
- [FIMO documentation](https://meme-suite.org/meme/doc/fimo.html).

The ChromHMM paper DOI is **10.1038/nmeth.1906**; references to a different DOI in earlier planning notes should not be reused.
