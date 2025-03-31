# Quik (Experiments)

## Name
Experiments for the publication "Fast barcode calling based on k-mer distances" by Uphoff et al. (2025)

## Description
The Java code of this repository contains the source code of the candidate set experiments and error model analysis of the paper.

## Installation
### Requirements
- Java 11
- Maven

### Compilation
To compile only the Java code of k-mer based calling algorithm, go to the root folder of the repository and run
```bash
mvn clean package
```
It should create a bin/ directory containing the executable .jar file (something like "barcode-calling-1.0.jar").



For the candidate set experiments, you need to run the algorithm by Press (2022) which is implemented in Python. For this, you have to set up the Python listener so that our tool can access Press' code. Run the "install_press.sh" script in the folder ./src/press_2022 of the repository. It will install the correct CUDA version and set up a virtual environment (venv) with all required packages. To start the listener, enter the virtual environment using
```bash
source src/press_2022/venv/bin/activate
```
and then start the communicator script with
```bash
python3 src/press_2022/java_communicator.py
```
Now the Java code should be able to call using Press' algorithm and run the experiments from the paper.

## Usage
### Simulating reads
To generate artificial erroneous reads according to the error model used by Uphoff et al., run

```bash
java -cp JAR_FILE.jar de.uni_halle.barcode_calling.experiments.data.Main <ARGUMENTS>
```

Possible arguments are:

| Argument         | Tags       | Description          |
|:-----------------|:----------:|:---------------------|
| -r --reads       | [required] |  Set output file for generated reads. |
| -l --labels      | [required] |  Set output file for the labels of the generated reads. |
| -b --barcodes    | [required] |  Set source file or directory for barcodes. |
| -n --number      | [required] |  Set the number of reads to by generated. |
| -s --pSub        | [required] |  Set the probability per base for substitution. |
| -i --pIns        | [required] |  Set the probability per base for insertion. |
| -d --pDel        | [required] |  Set the probability per base for deletion. |
| -l --log         |            |  Set output file or directory for logging. In case of directory, a default file name is chosen. If not provided, logs are only printed on the command line. |
| -h --help        |            |  Show possible arguments. |


### Experiments
To run a single experiment, go

```bash
java -cp JAR_FILE.jar de.uni_halle.barcode_calling.experiments.Main <ARGUMENTS>
```

Possible arguments are:

| Argument         | Tags       | Description          |
|:-----------------|:----------:|:---------------------|
| -o --out         | [required] |  Set output file or directory for experiment results. In case of directory, a default file name is chosen.|
| -b --barcodes    | [required] |  Set source file or directory for barcodes. In case of directory, all files in the directory are interpreted as source files. Required even for experiments without barodes.|
| -t --threads     |            |  Set the number of threads that should be used for calculations. Defaults to the number of cores + 1.|
| -l --log         |            |  Set output file or directory for logging. In case of directory, a default file name is chosen. If not provided, logs are only printed on the command line.|
| -e --experiment  | [required] |  Set the type of experiment you want to run. See the table below for available experiments and details. |
| -h --help        |            |  Show possible arguments.|

Available experiments are:
| Experiment               | Description          |
|:-------------------------|:---------------------|
| error-models                       | Compare error model used in Press (2022) to the error model by Uphoff et al. (2025). |
| pseudo-distance                    | Compare hit rate of k-mer pseudo-distances with k = 4, ..., 7 to the triage proposed by Press (2022) on specified barcode sets. ATTENTION: Requires python listener for Press' code to already be running. |
| pseudo-distance-operation-count    | Count inner loop iterations of k-mer pseudo-distances with k = 4, ..., 7 on specified barcode sets. |



## Authors
Authors: Riko Corwin Uphoff, Steffen Schüler, Ivo Große and Matthias Müller-Hannemann
