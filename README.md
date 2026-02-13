# Quik Barcode Calling


DNA barcodes, which are short DNA strings, are regularly used as tags in pooled sequencing experiments to enable the identification of reads originating from the same sample. A crucial task in the subsequent analysis of pooled sequences is barcode calling, where one must identify the corresponding barcode for each read.
This task is computationally challenging when the probability of synthesis and sequencing errors is high, like in photolithographic microarray synthesis.
Here, we propose a simple, yet highly effective new barcode calling approach that uses a filtering technique based on precomputed k-mer lists.
This approach has a slightly higher accuracy than the state-of-the-art approach, is more than 500 times
faster than that, and thus allows barcode calling for one million barcodes and one billion reads per day on a server GPU.


## Installation

Our software has been developed for Linux. The following steps are required to install quik on a Ubuntu 24.04 system:

1. Install software packages:

       sudo apt install git cmake g++ libomp-dev nvidia-cuda-toolkit

2. Checkout the project.

       git clone https://github.com/uni-halle/quik.git

3. Compile the source files. Quik comes with a CMake Script that should work for various environments.

   **Important**: Quik assumes a barcode set in which all barcodes have the same length. You need to specify this length during the build process using the option `-DBARCODE_LENGTH`. The following example uses barcodes of length 34.


       cd quik
       mkdir build
       cd build
       cmake -DCMAKE_BUILD_TYPE=Release -DBARCODE_LENGTH=34 ..
       make

The `build` directory should now contain three binaries:

- `quik` is the main tool used to assign reads to barcodes.
- `benchmark_barcode_calling` is an auxiliary program to test and compare the efficiency and accuracy of quik's different filtering techniques on your system.
- `simulate_errors` is an auxiliary program used to create our synthetic test data.

All binaries explain their usage when executed without arguments.

## Usage

The `quik` command line tool is used in the following way:

```bash
  quik --barcodes <FILE> --reads <FILE> [OPTIONS]

Required arguments:
  -b, --barcodes <FILE>     Barcode file in our custom BC format.

                            Such files have two line-separated fields per barcode.
                              - Field 1 begins with a '@' character and is followed
                                by a barcode identifier (similar to FASTQ).
                              - Field 2 is the raw sequence letters.

  -r, --reads <FILE>        Read file in FASTQ format.

Optional arguments:
  -d, --distance <STRING>   Distance measure specification

                            Possible values:

                              levenshtein
                              sequence-levenshtein

                            (default: sequence-levenshtein)

  -t, --threshold-distance  Maximum allowed distance between read and barcode.
                            A read is only assigned to a barcode if its
                            distance is smaller or equal than this integer.

                            (default: 2.147.483.647 (INT32_MAX))

  -m, --method <STRING>     Barcode matching method
                            Possible values:
                            
                              4-mer-filter
                                  Highest accuracy, slowest filtering variant

                              5-mer-filter
                                  Decent accuracy, faster than 4-mer-filter

                              6-mer-filter
                                  Low accuracy, much faster than 5-mer-filter

                              7-mer-filter
                                  Lowest accuracy, fastest method

                              7-4-mer-filter
                                  First runs the fast 7-mer-filter to assign the
                                  majority of the reads, then run the 4-mer-filter
                                  to assign the remaining reads.

                            (default: 4-mer-filter)

  -g, --gpu                 Apply all calculations on the first GPU visible to quik.
                            This is usually much faster than non-gpu mode.

  -h, --help                Show this help message and exit
  -v, --verbose             Print extra information to the standard error stream

Output:
  quik writes one line per assigned read to the standard output.
  Fields are tab-separated and have the following meaning:

    read           Sequence identifier of each read (at most once)
    barcode        Sequence identifier of closest barcode
    distance       Distance between read and barcode

  Reads with distance > threshold-distance do not occur in the output.


Examples:
  # Basic usage on GPU
  quik --barcodes barcodes.bc --reads reads.fq --gpu

  # Use 5-mer filtering and classical Levenshtein distance
  quik -b barcodes.bc -r reads.fq -m 5-mer-filter -d levenshtein
```

## Data

The data directory contains sample files that can be used for testing and which we used to produce some of the results in our journal article.
The reads and ground truth assignments in this directory have been created from the auxiliary tool `simulate_errors` using `prob=0.2`.



## License

This project is licensed under the GPL3. See the LICENSE file for details.

## Citation

If you use this repository or the associated article, please cite:

Riko Corwin Uphoff, Steffen Schüler, Ivo Grosse, Matthias Müller-Hannemann, Fast barcode calling based on k-mer distances, PNAS Nexus, 2026;, pgag001, https://doi.org/10.1093/pnasnexus/pgag001


```
@article{10.1093/pnasnexus/pgag001,
    author = {Uphoff, Riko Corwin and Schüler, Steffen and Grosse, Ivo and Müller-Hannemann, Matthias},
    title = {Fast barcode calling based on k-mer distances},
    journal = {PNAS Nexus},
    pages = {pgag001},
    year = {2026},
    month = {01},
    issn = {2752-6542},
    doi = {10.1093/pnasnexus/pgag001},
    url = {https://doi.org/10.1093/pnasnexus/pgag001},
    eprint = {https://academic.oup.com/pnasnexus/advance-article-pdf/doi/10.1093/pnasnexus/pgag001/66276546/pgag001.pdf},
}
```