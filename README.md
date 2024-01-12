# NEET - Nanopore Error pattern Exploration Toolkit

The Nanopore Error pattern Exploration Toolkit (`NEET`) provides a range of functionalities that provide an easily accessible and interactive analysis approach for (systematic) base-calling errors in direct RNA nanopore sequencing data. The implemented modules include options for condensing, visualizing and differentiating error features contained in direct RNA sequencing data - including mismatch, deletion and insertion rates, among others.

![Main workflow overview](https://github.com/dietvin/neet/blob/main/images/workflow_overview.png)

## Installation
It is recommended to use Conda or Mamba for installation:
```
conda install neet 
```
For more information and alternative installation approaches refer to the [Wiki](https://github.com/dietvin/neet/wiki/01-Installation).

## General usage
Once installed `NEET` can be accessed via the terminal:
```
neet --help
```
Individual modules can be accessed as follows:
```
neet [MODULE] --help
```
Available modules are: Pileup Extractor (`extractor`), Summary (`summmary`), Position-of-Interest Analyzer (`analyze_poi`), Two-Sample Extractor (`twosample`), Position Summary (`pos_summary`), Filter (`filter`) and Bedops (`bedops`). 
A detailled description of all available modules is provided in the [Wiki](https://github.com/dietvin/neet/wiki/02-Modules). The Wiki also provides [walkthroughs](https://github.com/dietvin/neet/wiki/03-Example-Workflows) for some possible use cases.
