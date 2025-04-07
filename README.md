# Quik Barcode Calling


 DNA barcodes, which are short DNA strings, are regularly used as tags in pooled sequencing experiments to enable the identification of reads originating from the same sample. A crucial task in the subsequent analysis of pooled sequences is barcode calling, where one must identify the corresponding barcode for each read.
This task is computationally challenging when the probability of synthesis and sequencing errors is high, like in photolithographic microarray synthesis.
Here, we propose a simple, yet highly effective new barcode calling approach that uses a filtering technique based on precomputed k-mer lists.
This approach has a slightly higher accuracy than the state-of-the-art approach, is more than 500 times
faster than that, and thus allows barcode calling for one million barcodes and one billion reads per day on a server GPU.


## Table of Contents
- [Repository Structure](#repository-structure)
- [Installation and Dependencies](#installation)
- [Benchmark Experiment](#benchmark-experiment)
- [License](#license)
- [Citation](#citation)


## Repository Structure

```
📦 Quik
├── 📁 data/ # Barcode read file required to reproduce our experiments data 
├── 📁 src/ # Source code for analyses and experiments 
├── 📁 results/ # Output data, plots, or tables 
├── 📄 README.md # This file 
├── 📄 CMakeLists.txt # CMake file for installation
└── 📄 LICENSE # License file
```

 
## Installation

The software has been developed for Linux and has been tested on an Ubuntu 24.04 system. The following steps are required for installation:

1. Install software packages:

       sudo apt install git cmake g++ libomp-dev nvidia-cuda-toolkit

2. Checkout the project.

       git clone https://github.com/uni-halle/quick.git

3. Compile the source files. Quik comes with a CMake Script that should work for various operating systems. CMake will automatically detect whether all mandatory and optional libraries are available at your system.

       cd quik
       mkdir build
       cd build
       cmake -DCMAKE_BUILD_TYPE=Release ..
       make

The `build` directory should now contain several binaries including `benchmark_barcode_calling`, used to compare the accurary and efficiency of various of our barcode calling approaches.

## Benchmark Experiment

After building the code, the sample binary `benchmark_barcode_calling` can be used to reproduce some of the results from our article. 

The benchmark program tries to assign each read in the `<read_file>` to some barcode in the `<barcode_file>`. To assess the accuracy of the barcode assignment, a `<label_file>` is required in which the original barcode is specified for each read. 

The benchmark program will run several of our k-mer distance barcode calling approaches and presents the associated precision and recall as well as the running time in milliseconds per read. Usage:

    benchmark_barcode_calling <barcode_file> <read_file> <label_file> <distance_measure> <rejection_threshold>


| Argument               | Description   |
|:---------------------|:--------------|
| `<barcode_file>`     | Text file with one barcode per line: `barcode[0]` ... `barcode[n-1]`. |
| `<read_file>`        | Text file with one read per line: `read[0]`... `read[m-1]`.           |
| `<label_file>`       | Text file with m lines of integers: `label[0]` ... `label[m-1]`. The integer `label[i]` is associated to `read[i]` and describes the index of the barcode, from which this read originated. Thus, `read[i]` originated from `barcode[label[i]]`.      |
| `<distance_measure>` | Distance measure between reads and barcodes. Must be one of the following: `LEVENSHTEIN`, `SEQUENCE_LEVENSHTEIN`. Has an effect on the accurary of the barcode calling process.     |
| `<rejection_threshold>`  | If a read's distance to the closest barcode is larger than this integer, the read is rejected and remains unassigned. Has a large effect on the accuracy of the barcode calling process.    |

The data directory contains sample files which we used to produce the results in our journal article.

## License

This project is licensed under the GPL3. See the LICENSE file for details.

## Citation

If you use this repository or the associated article, please cite it as follows:

```
@article{FastBarcodeCallingBasedOnKMerDistances2025,
  author  = {Riko Corvin Uphoff, Steffen Schüler, Ivo Grosse, Matthias Müller-Hannemann},
  title   = {Fast barcode calling based on k-mer distances},
  journal = {to be announced},
  year    = {2025},
  volume  = {X},
  pages   = {Y-Z},
  doi     = {DOI link}
}
```
