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

In this usage example, `WAAAAAAAAA` can expand to `AAAAAAAAAA` and `TAAAAAAAAA`, while `ASAAAAAAAA` can expand to `ACAAAAAAAA` and `AGAAAAAAAA`, covering the input set. This encoding covers the original input, when using IUPAC mapping, and is one of two solutions that is also the smallest such set.

## Notes

This aims for the smallest set of encoded strings. However, this does not report all possible such solutions in case there are more than one, as shown in the example above.

This has not been tested with any kmer sets but those examples provided. However, it aims to be scalable by pruning combinations of sub-kmers along the way that would yield incorrect encodings. This also uses a trie for fast, space-efficient prefix testing. If futher performance is needed, some easy wins would be to cache sub-kmer prefix tests, as most of these tests will be redundant.

Additionally, no error checking is done on the input kmer alphabet, the consistency of kmer lengths, or the uniqueness of kmers. It may be useful to validate input before using this script.

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
