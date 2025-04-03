# Quik Barcode Calling

This repository provides the code and data required to reproduce the results presented in the journal article

> Riko Corwin Uphoff, Steffen Schüler, Ivo Grosse and Matthias Müller-Hannemann,
> Fast barcode calling based on k-mer distances,
> submitted to PNAS Nexus, 2025.


## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation and Dependencies](#installation)
- [Usage](#usage)
- [License](#license)
- [Citation](#citation)

## Overview

The abstract of the journal article gives a short summary on our work:

> DNA barcodes, which are short DNA strings, are regularly used as tags in pooled sequencing experiments to enable the identification of reads originating from the same sample.
> A crucial task in the subsequent analysis of pooled sequences is barcode calling, where one must identify the corresponding barcode for each read.
> This task is computationally challenging when the probability of synthesis and sequencing errors is high, like in photolithographic microarray synthesis.
> Identifying the most similar barcode for each read is a theoretically attractive solution for barcode calling. However, an all-to-all exact similarity calculation is practically infeasible for applications with millions of barcodes and billions of reads.
> Hence, several computational approaches for barcode calling have been developed, but the challenge of developing an efficient and precise computational approaches remains.
>Here, we propose a simple, yet highly effective new barcode calling approach that uses a filtering technique based on precomputed $k$-mer lists.
> We find that this approach has a slightly higher accuracy than the state-of-the-art approach, is more than 500 times
>faster than that, and thus allows barcode calling for one million barcodes and one billion reads per day on a server GPU.


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

## Usage

After building the code, the sample binary `benchmark_barcode_calling` can be used to reproduce some of the results from our article. 

## License

This project is licensed under the GPL3. See the LICENSE file for details.

## Citation

If you use this repository or the associated article, please cite it as follows:

```
@article{FastBarcodeCallingBasedOnKMerDistances2025,
  author  = {Riko Corvin Uphoff, Steffen Schüler, Ivo Grosse, Matthias Müller-Hannemann},
  title   = {Fast barcode calling based on k-mer distances},
  journal = {submitted to PNAS Nexus},
  year    = {2025},
  volume  = {X},
  pages   = {Y-Z},
  doi     = {DOI link}
}
```
