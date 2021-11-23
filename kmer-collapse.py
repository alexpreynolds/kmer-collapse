#!/usr/bin/env python

'''
kmer-collapse.py

Starting with a set of kmers using a standard DNA alphabet (ACTG), collapse
it to an encoded, optimized set using IUPAC-derived, degenerate base alphabet.
'''

import sys
import json
import itertools as it
import marisa_trie as mt

"""
Generate all subsets of specified set. Empty set is fine.
"""
def powerset(seq):
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item

'''
We define a couple substitution objects to map between degenerate and
unambiguous bases; ref.: http://www.bioinformatics.org/sms/iupac.html
'''

def forward_substitutions():
    '''
    Return an object of tuple key to degenerate base value pairings. 
    Different tuples can point to the same degenerate base.
    '''
    raw_forward_substitutions = {
        'AG'   : 'R',
        'CT'   : 'Y',
        'GC'   : 'S',
        'AT'   : 'W',
        'GT'   : 'K',
        'AC'   : 'M',
        'CGT'  : 'B',
        'AGT'  : 'D',
        'ACT'  : 'H',
        'ACG'  : 'V',
        'ACGT' : 'N',
    }
    permuted_forward_substitutions = {}
    for key, value in raw_forward_substitutions.items():
        new_keys = list(it.permutations(list(key)))
        for new_key in new_keys:
            permuted_forward_substitutions[new_key] = value
    return permuted_forward_substitutions

def reverse_substitutions():
    '''
    Return an object of degenerate base key to unambiguous base value pairings.
    '''
    return {
        'A' : ['A'],
        'T' : ['T'],
        'C' : ['C'],
        'G' : ['G'],
        'R' : ['A', 'G'],
        'Y' : ['C', 'T'],
        'S' : ['G', 'C'],
        'W' : ['A', 'T'],
        'K' : ['G', 'T'],
        'M' : ['A', 'C'],
        'B' : ['C', 'G', 'T'],
        'D' : ['A', 'G', 'T'],
        'H' : ['A', 'C', 'T'],
        'V' : ['A', 'C', 'G'],
        'N' : ['A', 'C', 'G', 'T'],
    }

forward_subs = forward_substitutions()
reverse_subs = reverse_substitutions()

def expand_encoded_candidate(encoded_candidate):
    '''
    Given an encoded mer, expand it to all ordered combinations of unambiguous
    bases. For example, 'WW' expands to ['AA', 'AT', 'TA', 'TT']
    '''
    n_candidate = len(encoded_candidate)
    subs = []
    for c in list([x for x in encoded_candidate]):
        rrs = reverse_subs[c]
        subs.append(rrs)
    expanded_candidates = []
    i = 1
    while True:
        sub = subs.pop(0)
        if len(subs) == 0: break
        new_sub = [f'{a}{b}' for a in sub for b in subs[0]]
        i += 1
        if i == n_candidate: 
            expanded_candidates = new_sub
        else:
            subs[0] = new_sub
    expanded_candidates = list(set(expanded_candidates))
    return expanded_candidates

def test_encoded_candidate(encoded_candidate):
    '''
    Validate encoded candidate, if all expanded candicates are valid prefixes
    in global input trie.
    '''
    expanded_candidates = expand_encoded_candidate(encoded_candidate)
    for ec in expanded_candidates:
        if len(input_trie.keys(ec)) == 0:
            return False
    return True

def main():
    '''
    Set up input data with kmers to encode.
    '''
    # input = ['AAAAAAAAAA', 'TAAAAAAAAA']
    # input = ['AAAAAAAAAA','TAAAAAAAAA','GCGAAAAAAA']
    # input = ['AAAAAAAAAA']
    # input = ['AAAAAAAAAA', 'TAAAAAAAAA', 'CAAAAAAAAA', 'GAAAAAAAAA']
    # input = ['AAAAAAAAAA', 'TAAAAAAAAA', 'TTAAAAAAAA', 'ATAAAAAAAA']
    # input = ['AAAAAAAAAA', 'TAAAAAAAAA', 'CAAAAAAAAA', 'GAAAAAAAAA', 'TACAGATACA', 'AACAGAAAAA']
    input = ['AAAAAAAAAA', 'TAAAAAAAAA', 'ACAAAAAAAA', 'AGAAAAAAAA']
    k = len(input[0])
    global input_trie
    input_trie = mt.Trie(input)
    
    '''
    At each column, encode that column's bases, taken from all kmers.
    '''
    per_sequence_encoded_candidates = [{x:[]} for x in input]
    for col_idx in range(k):
        bases_to_permute = list(set([x[col_idx] for x in input]))
        bases_powerset = list(powerset(bases_to_permute))
        per_column_encodings = {}
        for bpse in bases_powerset:
            tk = tuple(bpse)
            if tk in forward_subs:
                for tkb in tk:
                    if tkb not in per_column_encodings:
                        per_column_encodings[tkb] = []
                    if forward_subs[tk] not in per_column_encodings[tkb]:
                        per_column_encodings[tkb].append(forward_subs[tk])
        for b in bases_to_permute:
            if b not in per_column_encodings:
                per_column_encodings[b] = []
            per_column_encodings[b].append(b)

        '''
        Initialize candidates per sequence, or extend after testing.
        '''
        for psec in per_sequence_encoded_candidates:
            psck = list(psec.keys())[0]
            ek = psck[col_idx]
            if col_idx == 0:
                for ekv in per_column_encodings[ek]:
                    psec[psck].append('{}'.format(ekv))
            else:
                old_candidates = psec[psck]
                new_candidates = []
                for oc in old_candidates:
                    for pcec in per_column_encodings[psck[col_idx]]:
                        new_candidate = '{}{}'.format(oc, pcec)
                        '''
                        Expand the new candidate. Test it. If the new candidate
                        passes testing when expanded, add it to the new list.
                        Continually pruning the tree of bad candidates should 
                        help limit the number of combinations to explore for 
                        longer kmers.
                        '''
                        if test_encoded_candidate(new_candidate):
                            new_candidates.append(new_candidate)
                psec[psck] = new_candidates

    '''
    Filter candidates by kmer "row" to resolve conflicts.
    '''
    kmer_available = {}
    for psec in per_sequence_encoded_candidates:
        for pseck in psec.keys():
            kmer_available[pseck] = True
    encodings = set()
    for psec in per_sequence_encoded_candidates:
        for pseck, psecv in psec.items():
            encoded_candidates = psecv
            for ec in encoded_candidates:
                eecs = expand_encoded_candidate(ec)
                '''
                Prime the pump.
                '''
                if len(encodings) == 0:
                    for eec in eecs:
                        encodings.add(ec)
                        kmer_available[eec] = False
                    break
                '''
                If we pick an encoding that would generate a kmer
                that had previously been generated, this would give 
                an incorrect result. So we walk through the expanded
                set of unambiguous kmers and add an encoding if all 
                its expanded kmers are available for inclusion in the 
                result set.
                '''
                expansion_test_passed = True
                for eec in eecs:
                    if not kmer_available[eec]:
                        expansion_test_passed = False
                if expansion_test_passed:
                    for eec in eecs:
                        encodings.add(ec)
                        kmer_available[eec] = False
    encodings = list(encodings)

    '''
    Write encoded output.
    '''
    result = {
        "input" : input,
        "encoded_output" : encodings
    }
    sys.stdout.write('{}\n'.format(json.dumps(result, indent=4)))

if __name__ == '__main__':
    main()
