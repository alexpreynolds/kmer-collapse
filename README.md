# kmer-collapse
Collapse a set of redundant kmers to use IUPAC degenerate bases

## Overview

Given an input set of kmers, find the smallest set of kmers that encapsulates all diversity in the input set using IUPAC degenerate bases. This aims to solve the problem described here: https://www.biostars.org/p/9498272/

## Usage

Install the [marisa-trie](https://pypi.org/project/marisa-trie/) library, if necessary.

Modify the script's `input` variable to specify desired sequences, and then run the script:

```
$ python kmer-collapse.py
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "ACAAAAAAAA",
        "AGAAAAAAAA"
    ],
    "encoded_output": [
        "WAAAAAAAAA",
        "ASAAAAAAAA"
    ]
}
```

## Notes

This has not been tested with any kmer sets but those examples provided. However, it aims to be scalable by pruning combinations of sub-kmers along the way, which would otherwise yield incorrect encodings. This also uses a trie for faster prefix testing. If futher performance is needed, some easy wins would be to cache sub-kmer tests, since most of these test outcomes would be redundant.

Additionally, no error checking is done on the input kmer alphabet or on the consistency of kmer lengths. It may be useful to validate input before using this script.

## Examples

These examples are available from the script by uncommenting the relevant `input`.

### A

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA"
    ],
    "encoded_output": [
        "WAAAAAAAAA"
    ]
}
```

### B

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "GCGAAAAAAA"
    ],
    "encoded_output": [
        "GCGAAAAAAA",
        "WAAAAAAAAA"
    ]
}
```

### C 

```
{
    "input": [
        "AAAAAAAAAA"
    ],
    "encoded_output": [
        "AAAAAAAAAA"
    ]
}
```

### D

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "CAAAAAAAAA",
        "GAAAAAAAAA"
    ],
    "encoded_output": [
        "NAAAAAAAAA"
    ]
}
```

### E

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "TTAAAAAAAA",
        "ATAAAAAAAA"
    ],
    "encoded_output": [
        "WWAAAAAAAA"
    ]
}
```

### F

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "CAAAAAAAAA",
        "GAAAAAAAAA",
        "TACAGATACA",
        "AACAGAAAAA"
    ],
    "encoded_output": [
        "NAAAAAAAAA",
        "TACAGATACA",
        "AACAGAAAAA"
    ]
}
```

### G

```
{
    "input": [
        "AAAAAAAAAA",
        "TAAAAAAAAA",
        "ACAAAAAAAA",
        "AGAAAAAAAA"
    ],
    "encoded_output": [
        "ASAAAAAAAA",
        "WAAAAAAAAA"
    ]
}
```
