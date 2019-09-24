#!/usr/bin/env python

import json
import re

class RuleMatcher(object):
    def matches(this, accession):
        pass

class RangeMatcher(RuleMatcher):
    
    def __init__(self, prefix, suffix_length):
        self.re = re.compile(r'{prefix}[A-Z]{suffix_length}')
    
    def matches(this, accession):
        return self.re.match(accession) is not None

def make_range_matcher(spec):
    if '-' in spec:
        (start, end) = spec.split('-')
    elif '_' in spec:
        (start, end) = spec.split('_')
    else:
        raise ValueError('require specification with - or _ in it, got: ' + spec)
    assert end.endswith('Z'), f'Range specifier has unknown format: {spec}'
    for i in range(len(start)):
        if start[i] != end[i]:
            break
    prefix_length = i
    assert prefix_length < len(start), f'Prefix not found in range specifier {spec}'
    suffix_length = len(start) - prefix_length
    return (prefix_length + suffix_length, RangeMatcher(start[:prefix_length], suffix_length))


def build_accession_parser(rules_file):
    rules_data = json.load(rules_file)
    rules_by_prefix_len = {}
    for prefix_list, database, type_description in rules_data:
        for prefix in prefix_list:
            prefix_length = len(prefix)
            if '-' in prefix or '_' in prefix:
                (prefix_length, matcher) = make_range_matcher(prefix)
                if prefix_length not in rules_by_prefix_len:
                    rules_by_prefix_len[prefix_length] = []
                rules_by_prefix_len[prefix_length].append((matcher, database, type_description))
            else:
                if prefix_length not in rules_by_prefix_len:
                    rules_by_prefix_len[prefix_length] = []
                rules_by_prefix_len[prefix_length].append((prefix, database, type_description))
    return rules_by_prefix_len

letter_re = re.compile('[A-Z]+')
def match_accession(accession, rules_by_prefix_len):
    if accession[2] == '_':
        return ('NCBI', 'RefSeq')
    letter_match = letter_re.match(accession)
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
        accession_type = 'protein'
    elif letter_match_length == 4 or letter_match_length == 6:
        accession_type = 'WGS'
    elif letter_match_length == 5:
        accession_type = 'MGA'
    else:
        raise ValueError('an accession number must start with less than 7 capital letters, this does not: ' + accession_type)
    
    rules = rules_by_prefix_len[letter_match_length]
    for rule in rules:
        (matcher, database, type_description) = rule
        if (isinstance(matcher, RuleMatcher) and matcher.matches(accession)) or letter_prefix == matcher:
            return (database, type_description)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('rules_data_file', type=argparse.FileType())
    parser.add_argument('accession')
    args = parser.parse_args()
    rules_by_prefix_len = build_accession_parser(args.rules_data_file)
    print(match_accession(args.accession, rules_by_prefix_len))