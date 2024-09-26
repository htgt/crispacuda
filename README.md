# crispacuda
# 
This application searches for CRISPR off-targets on an NVIDIA GPU. crispacuda will use a dedicated GPU to to load known CRISPR sequences into memory and then perform searches on these sequences.

## Index Files

In order for crispacuda to run it needs an index file containing an array of all known CRISPRs in a genome assembly.
The files are in binary format and are created using the [CRISPR-Analyser](https://github.com/htgt/CRISPR-Analyser) utility.
Follow the [Create Index](https://github.com/htgt/CRISPR-Analyser?tab=readme-ov-file#create-index) instructions on how to create a binary index file to use in crispacuda.

## Runtime Options

crispacuda is run from the command line with the following options:

- `[-i <FILE>]` - specify the index file to search (required)
- `[-s <INT>]` - the index to start the searches from
- `[-n <INT>]` - how many off-targets to calculate from the start
- `[-m <INT>]` - the maximum number of mismatches to records results for. Defaults to 4
- `[-q]` - do not report a list of off-targets
- `[-d <INT>]` - the GPU device to use
- `[-z]` - specify CRISPRs as strings rather than IDs
- `[-h]` - print a help message and GPU information


