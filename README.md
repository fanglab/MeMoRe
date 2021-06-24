# memore

`MeMoRe` stand for Methylation Motif Refiner. This is a tool designed to validate methylation motifs detected with SMRT and ONT sequencing.

## Authors' notes
We are actively developing `MeMoRe` to facilitate usage and broaden features. All feedback is more than welcome. You can reach us on twitter ([@iamfanggang](https://twitter.com/iamfanggang) and [@AlanTourancheau](https://twitter.com/AlanTourancheau)) or directly through the [GitHub issues system](https://github.com/fanglab/MeMoRe/issues).

## Content
+ [Usage](#Usage)
+ [Tool showcase](#Tool-showcase)
  + [Principle](#Principle)
  + [Analysis of SMRT results](#Analysis-of-SMRT-results)
  + [Analysis of ONT results](#Analysis-of-ONT-results)
+ [Documentation](#Documentation)
+ [Citation](#Citation)

## Usage
To apply `MeMoRe` on your own dataset or use our testing dataset please follow: https://fanglab-tools.shinyapps.io/MeMoRe/. `MeMoRe` is designed to be user freindly and therefore is distributed as a [Shiny] application through [shinyapps.io] hosting solution. `MeMoRe` can be run from any system, without computational resources restriction.

## Tool showcase
To showcase the toolbox applications and facilitate an understanding of the methods, we build-in example datasets for SMRT and ONT analyses in [MeMoRe app](https://fanglab-tools.shinyapps.io/MeMoRe/). Typical analysis results for both SMRT and ONT datasets are presented below.

### Principle
In prokaryote and archaea, DNA modification are motif-driven, meaning that nearly all occurrences of the same sequence motif(s) will be modified. This property can be used to refine the motifs discovered from `SMRT Portal/Link Base Modification Analysis <https://www.pacb.com/support/software-downloads/>`_ or `nanodisco <https://github.com/fanglab/nanodisco>`_ pipelines.

### Analysis of SMRT results
In SMRT sequencing, DNA methylation affect the kinetics of the polymerases used for the sythesis of the SMRTBell templates. The changes of polymerase's kinetics are observed through the Inter-Pulse Duration (IPD) metric which are compared to prediction from an *in silico* model at each genomic position. The resulting metric is called the IPD ratio (IPD native/IPD *in silico*). For 6mA and 4mC DNA modification, the IPD ratio increase on top of the methylated positions while an IPD ratio of 1 means no kinetic change. However, 5mC do not typically produce detectable signal and cannot be reliably found from SMRT data.

![Output SMRT](/docs/figures/GTATAC_5_combined.png "C. perfringens's GTAT6mAC methylation motif results")

### Analysis of ONT results

![Output ONT](/docs/figures/GTATAC_5_ont.png "C. perfringens's GTAT6mAC methylation motif results")

## Documentation
For a comprehensive description of `MeMoRe` including a detailed tutorial, please consult the [complete documentation][Full Documentation].

## Citation

A preprint is in preparation. Meanwhile, please cite the GitHub repository: https://github.com/fanglab/MeMoRe.

[Shiny]: https://shiny.rstudio.com/
[shinyapps.io]: https://www.shinyapps.io/
[Full Documentation]: https://MeMoRe.readthedocs.io/en/latest/overview.html