#!/usr/bin/env python

import sys
import subprocess
import re
import os.path
from collections import namedtuple


TestCase = namedtuple('TestCase', 'input alignment score')

test_cases = [
    TestCase(
        input = { 'a' : 'ACGATA', 'b' : 'AGTGGTA'
                , 'l' : 1, 'm' : 2, 't' : 0.1
                , 'A' : 0.25, 'C' : 0.25
                , 'G' : 0.25, 'T' : 0.25 },
        alignment = ('A-CGATA', 'AGTGGTA'),
        score = 4.341e-12),
    TestCase(
        input = { 'a' : 'ACATA', 'b' : 'CAATT'
                , 'l' : 1, 'm' : 2, 't' : 0.1
                , 'A' : 0.25, 'C' : 0.25
                , 'G' : 0.25, 'T' : 0.25 },
        alignment = ('-ACATA', 'CA-ATT'),
        score = 3.81782e-10),
    TestCase(
        input = { 'a' : 'ACATA', 'b' : 'CAATT'
                , 'l' : 1, 'm' : 2, 't' : 0.1
                , 'A' : 0.27, 'C' : 0.24
                , 'G' : 0.26, 'T' : 0.23 },
        alignment = ('AC-ATA', '-CAATT'),
        score = 4.21821e-10),
    TestCase(
        input = { 'a' : 'ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA'
                , 'b' : 'AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA'
                , 'l' : 1, 'm' : 2, 't' : 0.1
                , 'A' : 0.25, 'C' : 0.25
                , 'G' : 0.25, 'T' : 0.25 },
        alignment =
            ('ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA'
           , 'A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A'),
        score = 8.54034e-78),
    TestCase(
        input = { 'a' : 'ATGGGTGATGTTGAGAAAGGCAAGAAGATTTTTATTATGAAGTGTTCCCAGTGCCACACCGTTGAAAAGGGAGGCAAGCACAAGACTGGGCCAAATCTCCATGGTCTCTTTGGGCGGAAGACAGGTCAGGCCCCTGGATACTCTTACACAGCCGCCAATAAGAACAAAGGCATCATCTGGGGAGAGGATACACTGATGGAGTATTTGGAGAATCCCAAGAAGTACATCCCTGGAACAAAAATGATCTTTGTCGGCATTAAGAAGAAGGAAGAAAGGGCAGACTTAATAGCTTATCTCAAAAAAGCTACTAATGAGTAA' # Human
                , 'b' : 'ATGGGTGATGTTGAAAAAGGCAAGAAGATTTTTGTTCAGAAGTGTGCCCAGTGCCACACTGTGGAAAAGGGAGGCAAGCATAAGACTGGACCAAATCTCCACGGTCTGTTCGGGCGGAAGACAGGCCAGGCTGCTGGATTCTCTTACACAGATGCCAACAAGAACAAAGGCATCACCTGGGGAGAGGATACCCTGATGGAGTATTTGGAGAATCCCAAAAAGTACATCCCTGGAACAAAAATGATCTTCGCTGGAATTAAGAAGAAGGGAGAAAGGGCAGACCTAATAGCTTATCTTAAAAAGGCTACTAATGAGTAA' # Mouse
                , 'l' : 1, 'm' : 2, 't' : 0.1
                , 'A' : 0.27, 'C' : 0.24
                , 'G' : 0.26, 'T' : 0.23 },
        alignment = ('ATGGGTGATGTTGAGAAAGGCAAGAAGATTTTTATT-ATGAAGTGTTCCCAGTGCCACACCGTTGAAAAGGGAGGCAAGCACAAGACTGGGCCAAATCTCCATGGTCTCTTTGGGCGGAAGACAGGTCAGGCCCCTGGATACTCTTACACAGCCGCCAATAAGAACAAAGGCATCATCTGGGGAGAGGATACACTGATGGAGTATTTGGAGAATCCCAAGAAGTACATCCCTGGAACAAAAATGATCTTTG-TCGGCATTAAGAAGAAGGAAGAAAGGGCAGACTTAATAGCTTATCTCAAAAAAGCTACTAATGAGTAA', 'ATGGGTGATGTTGAAAAAGGCAAGAAGATTTTTGTTCA-GAAGTGTGCCCAGTGCCACACTGTGGAAAAGGGAGGCAAGCATAAGACTGGACCAAATCTCCACGGTCTGTTCGGGCGGAAGACAGGCCAGGCTGCTGGATTCTCTTACACAGATGCCAACAAGAACAAAGGCATCACCTGGGGAGAGGATACCCTGATGGAGTATTTGGAGAATCCCAAAAAGTACATCCCTGGAACAAAAATGATCTTCGCT-GGAATTAAGAAGAAGGGAGAAAGGGCAGACCTAATAGCTTATCTTAAAAAGGCTACTAATGAGTAA'),
        score = 0)
    ]

def input_args(test_case):
    return [x for k, v in test_case.input.items() for x in ['-' + k, str(v)]]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('Please specify the relative or absolute path to mlalign as the first and only parameter!')

    mlalign = os.path.abspath(sys.argv[1])

    print('Testing executable at ' + mlalign)
    print('')

    regexp = re.compile(r'length (\d+) with score ([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?):\s*([aAcCgGtT-]+)\s+([aAcCgGtT-]+)')
    for i, test_case in enumerate(test_cases):
        ea, eb = test_case.alignment
        assert len(ea) == len(eb), 'The expected alignment has differing length in test case ' + str(test_case)
        output = subprocess.check_output([mlalign] + input_args(test_case))
        match = regexp.search(output.decode('utf-8'))
        if not match:
            sys.exit('Output format could not be matched against for test case ' + str(test_case))

        n, score, aa, ab = match.group(1, 2, 4, 5) # group 3 is part of the float matching regex, so we just ignore it

        assert len(aa) == len(ab), 'The actual alignment has differing length. {0} vs {1}'.format(len(aa), len(ab))
        assert int(n) == len(aa), 'The length of the actual alignment {0} and the announced length {1} mismatch.'.format(len(aa), int(n))
        assert int(n) == len(ea), 'The length of the expected alignment {0} and the announced length {1} mismatch.'.format(len(ea), int(n))

        a_equal = [x == y for x, y in zip(ea, aa)]
        if not all(a_equal):
            assert False, 'The topmost actual alignment and the expected alignment differ at index {0}'.format(a_equal.index(False))

        b_equal = [x == y for x, y in zip(eb, ab)]
        if not all(b_equal):
            assert False, 'The bottommost actual alignment and the expected alignment differ at index {0}'.format(b_equal.index(False))

        assert test_case.score == float(score), 'The actual score {0} did not match the expected score of {1}'.format(float(score), test_case.score)

        print('  Passed test case {0}.'.format(i))

    print('')
    print('Passed all {0} test cases!'.format(len(test_cases)))
