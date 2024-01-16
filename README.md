# NEET - Nanopore Error pattern Exploration Toolkit

The Nanopore Error pattern Exploration Toolkit (`NEET`) provides a range of functionalities that provide an easily accessible and interactive analysis approach for (systematic) base-calling errors in direct RNA nanopore sequencing data. The implemented modules include options for condensing, visualizing and differentiating error features contained in direct RNA sequencing data - including mismatch, deletion and insertion rates, among others.

## Installation
It is recommended to use Conda/Mamba and pip for installation:
```
conda create -y -n neet python==3.10 && conda activate neet
pip install neet-nanopore
neet --version
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
Available modules are: Pileup Extractor (`extractor`), Summary (`summmary`), Position-of-Interest Analyzer (`analyzepoi`), Two-Sample Extractor (`twosample`), Position Summary (`possummary`), Filter (`filter`) and Bedops (`bedops`). 
The modules are interconnected to provide different workflows for different scenarios. An overview of the main workflows is given below:

![Main workflow overview](https://github.com/dietvin/neet/blob/main/images/workflow_overview.jpg)

A detailled description of all available modules is provided in the [Wiki](https://github.com/dietvin/neet/wiki/02-Modules). The Wiki also provides [walkthroughs](https://github.com/dietvin/neet/wiki/03-Example-Workflows) for some possible use cases.
