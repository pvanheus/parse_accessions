#!/usr/bin/env python

"""Functions for parsing accessions and returning their source database and data type."""
# Copyright (c) 2019, Peter van Heusden <pvh@sanbi.ac.za>
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Peter van Heusden BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import print_function
import json
import re

class RuleMatcher(object):
    """Superclass for matcher objects that have a match() method to return whether they match an accession."""

    def matches(self, accession):
        """matches(self, accession)

        Returns true of accession matches this rule."""
        pass

class RangeMatcher(RuleMatcher):
    """RuleMatcher subclass that can match an accession whose prefix is in a range of letter prefixes."""
    
    def __init__(self, prefix, suffix_length):
        # this is not a greedy match so we need to specify that it is followed by a number
        self.re = re.compile(f'{prefix}[A-Z]{suffix_length}' + r'\d')


    def matches(self, accession):
        return self.re.match(accession) is not None

def make_range_matcher(spec):
    """make_range_matcher(spec)

    spec - string (e.g. 'BAAA-BZZZ')

    Turns a range specification into a object that can match an accession that falls
    within that range. Returns a RangeMatcher() object"""
    if '-' in spec:
        (start, end) = spec.split('-')
    elif '_' in spec:
        (start, end) = spec.split('_')
    else:
        raise ValueError('require specification with - or _ in it, got: ' + spec)
    assert end.endswith('Z'), f'Range specifier has unknown format: {spec}'
    for i, _ in enumerate(start):
        if start[i] != end[i]:
            break
    prefix_length = i
    assert prefix_length < len(start), f'Prefix not found in range specifier {spec}'
    suffix_length = len(start) - prefix_length
    return (prefix_length + suffix_length, RangeMatcher(start[:prefix_length], suffix_length))


REFSEQ_PREFIX_RE = re.compile('[A-Z]{2}_')
def build_accession_parser(rules_file):
    """build_accession_parser(rules_file)

    rules_file - file object open for reading

    Builds a rule parse usable by match_accession() in this module from the accession rules
    described in rules_file. These rules can be downloaded by the scrape_accession_rules
    script."""

    rules_data = json.load(rules_file)
    rules_by_prefix_len = {}
    for prefix_list, database, molecule_type, type_description in rules_data:
        for prefix in prefix_list:
            prefix_length = len(prefix)
            if REFSEQ_PREFIX_RE.match(prefix) is not None:
                # RefSeq whose accessions start with XX_ has its own rules
                if 'RefSeq' not in rules_by_prefix_len:
                    rules_by_prefix_len['RefSeq'] = []
                rules_by_prefix_len['RefSeq'].append((prefix, database, molecule_type, type_description))
            elif '-' in prefix or '_' in prefix:
                (prefix_length, matcher) = make_range_matcher(prefix)
                if prefix_length not in rules_by_prefix_len:
                    rules_by_prefix_len[prefix_length] = []
                rules_by_prefix_len[prefix_length].append((matcher, database, molecule_type, type_description))
            else:
                if prefix_length not in rules_by_prefix_len:
                    rules_by_prefix_len[prefix_length] = []
                rules_by_prefix_len[prefix_length].append((prefix, database, molecule_type, type_description))
    return rules_by_prefix_len

# letter_re is a greedy match so we do not need to specify that the letter prefix is followed by a number
LETTER_RE = re.compile('[A-Z]+')
def match_accession(accession, rules_by_prefix_len):
    """match_accession(accession, rules_by_prefix_len)
    
    Returns the tuple (database, accession_type, type_description) for accession using the rules specified by
    rules_by_prefix_len. The database is one of 'GenBank', 'NCBI', 'GenBank and DDBJ' , 'DDBJ', and 'EMBL',
    accession_type is one of 'nucleotide', 'protein', 'WGS' and 'MGA' and the type_description is the description
    of the type of data associated with that accession."""

    letter_match = LETTER_RE.match(accession)
    if letter_match is None:
        raise ValueError('an accession number must start with at least one capital letter, this does not: ' + accession)
    letter_prefix = letter_match.group(0)
    letter_match_length = len(letter_prefix)
    accession_type = ''
    if letter_match_length == 0:
        # this should never happen
        raise ValueError('an accession number must start with at least one capital letter, this does not: ' + accession)
    if letter_match_length < 3:
        accession_type = 'nucleotide'
    elif letter_match_length < 4:
        if letter_prefix in ('SRA', 'SRP', 'SRX', 'SRR', 'SRS', 'SRZ',
                             'ERA', 'ERP', 'ERX', 'ERR', 'ERS', 'ERZ',
                             'DRA', 'DRP', 'DRX', 'DRR', 'DRS', 'DRZ'):
            accession_type = 'SRA'
        else:
            accession_type = 'protein'
    elif letter_match_length in (4, 6):
        accession_type = 'WGS'
    elif letter_match_length == 5:
        accession_type = 'MGA'
    else:
        raise ValueError('an accession number must start with less than 7 capital letters, this does not: ' + accession_type)

    if accession[2] == '_':
        # details from https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
        # thanks to Torsten Seemann
        for prefix, database, molecule_type, type_description in rules_by_prefix_len['RefSeq']:
            accession_type = molecule_type
            if accession[:3] == prefix:
                return (database, accession_type, 'RefSeq: ' + type_description)

    rules = rules_by_prefix_len[letter_match_length]
    for rule in rules:
        (matcher, database, _, type_description) = rule
        if (isinstance(matcher, RuleMatcher) and matcher.matches(accession)) or letter_prefix == matcher:
            return (database, accession_type, type_description)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('rules_data_file', type=argparse.FileType())
    parser.add_argument('accession')
    args = parser.parse_args()
    rules_by_prefix_len = build_accession_parser(args.rules_data_file)
    print(match_accession(args.accession, rules_by_prefix_len))
