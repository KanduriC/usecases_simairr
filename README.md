# usecases_simairr

This repository contains scripts for reproducing the case studies and analyses of simAIRR manuscript.

## Instructions

- It is recommended to install both `simAIRR` and `usecases_simairr` in a new conda environment. It is also required that [compAIRR](https://github.com/uio-bmi/compairr) and [immuneML](https://github.com/uio-bmi/immuneML) are installed. The instructions further below assumes that `simAIRR`, `compAIRR`, `immuneML` and `usecases_simairr` are installed.
- All the scripts in `usecases_simairr` are exported as console scripts to be run on command-line interface. The code chunks below show how the console scripts are to be run.

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

