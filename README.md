# usecases_simairr

This repository contains scripts for reproducing the case studies and analyses of simAIRR manuscript. In addition, this package provides console scripts that can be run through command-line to learn custom models of realistic receptor sharing for different AIRR loci (i.e. to generate custom pgen_count_map files based on different datasets).

## Instructions

- It is recommended to install both `simAIRR` and `usecases_simairr` in a new conda environment. It is also required that [compAIRR](https://github.com/uio-bmi/compairr) and [immuneML](https://github.com/uio-bmi/immuneML) are installed. The instructions further below assumes that `simAIRR`, `compAIRR`, `immuneML` and `usecases_simairr` are installed.
- All the scripts in `usecases_simairr` are exported as console scripts to be run on command-line interface. The code chunks below show how the console scripts are to be run.

### Custom models of realistic receptor sharing for different AIRR loci

- Note that in simAIRR package, we have provided empirical models for realistic receptor sequence sharing just for TCRB loci, based on the public availability of a large experimental dataset. For other AIRR loci, custom models will be added into future versions of simAIRR when large-scale experimental data for other loci becomes available. For now users are encouraged to generate their own custom models of realistic receptor sharing for other loci (if sufficiently large data available to them). Note that, based on our analyses, a sample size of 200 (100 each of cases and controls) approximates well the realistic receptor sequence sharing. 

- Through this package, we provide console scripts that can be run through command line to learn custom models for realistic receptor sequence sharing. Each console script requires some input parameters from the users. To know which input arguments are required for each console script, we recommend using the `--help` argument of each console script. For instance, `concat_airr --help`. 

- Note that the following console scripts are expected to be run in a sequence.

| Command | Description                                                                                                                                                                         |
| --- |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `concat_airr` | Concatenate the repertoire files into two separate files that will be needed by [compAIRR](https://github.com/uio-bmi/compairr) tool to count the frequencies of each unique clones |
| `run_compairr_seqcounts` | Run [compAIRR](https://github.com/uio-bmi/compairr) tool to count the frequencies of each unique clones.                                                                            |
| `compute_pgen_public` | Compute the generation probabilities of each sequence using OLGA tool.                                                                                                              |
| `compute_pval` | Computes p-values as the probability of observing a public sequence in the same or higher number of repertoires as it was observed in for a given dataset.                          |
| `concat_pdata` | Concatenate the files generated through previous step to produce one large file that will be used in the next step                                                                  |
| `gen_pgen_count_map` | Generate pgen_count_map files in the required format desired by simAIRR                                                                                                             |


### case study-1

The ```usecases_simairr/kmers_usecase``` directory contains yaml specification files and a tsv file containing sequences that will be used in the examples below. The yaml specification files contain paths that will no longer be valid when you are trying to reproduce the code on your machines. Please change the paths to suit your machine paths.

The following code line will read-in a config file and generate multiple config files with varying simulation parameters that are compatible with simAIRR. Each simAIRR-compatible config file will be run using simAIRR to generate different simulated datasets with varying phenotype burdens (witness rates) replicated in three independent runs. (i..e three independent simulated datasets will be generated with similar simulation parameters). 

```
$ multi_simairr -i kmers_multi_simairr_config.yaml
```
The simAIRR config files that will be generated are shown in ```usecases_simairr/kmers_usecase/simulation_scripts```

The following code line will train a machine learning model on all the simulated datasets generated with above code line.

```
$ multi_ml -s /path/to/kmers_usecase -m logistic_kmer.yaml
```
### case study-2

The ```usecases_simairr/fullseq_usecase``` directory contains yaml specification files and a tsv file containing sequences that will be used in the examples below. The yaml specification files contain paths that will no longer be valid when you are trying to reproduce the code on your machines. Please change the paths to suit your machine paths.

The following code line will read-in a config file and generate multiple config files with varying simulation parameters that are compatible with simAIRR. Each simAIRR-compatible config file will be run using simAIRR to generate different simulated datasets with varying phenotype burdens (witness rates) replicated in three independent runs. (i..e three independent simulated datasets will be generated with similar simulation parameters). Since we also vary sample size in this case study, we run the following three code lines for the three different sample sizes.

```
$ multi_simairr -i samplesize_200_fullseq_multi_simairr_config.yaml
$ multi_simairr -i samplesize_400_fullseq_multi_simairr_config.yaml
$ multi_simairr -i samplesize_600_fullseq_multi_simairr_config.yaml
```

The simAIRR config files that will be generated are shown in ```usecases_simairr/fullseq_usecase/simulation_scripts```

The following code lines will train machine learning models on all the simulated datasets generated with above code line. Note that the code line is shown for all the three lines in the same line for brevity. In practice, those were three different lines for each run.

```
$ multi_ml -s /path/to/sample_size_200/{run1, run2, run3} -m ml_sample_size_200.yaml
$ multi_ml -s /path/to/sample_size_400/{run1, run2, run3} -m ml_sample_size_400.yaml
$ multi_ml -s /path/to/sample_size_600/{run1, run2, run3} -m ml_sample_size_600.yaml
```

