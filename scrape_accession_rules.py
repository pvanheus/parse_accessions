#!/usr/bin/env python
"""Script for parsing accession description page and saving to JSON.

Usage: scrape_accession_rules.py <output file>

e.g.

scrape_accession_rules.py accession_rules.json
"""

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

import argparse
import json
from typing import List, Tuple, TextIO
RulesData = List[Tuple[List[str], str, str, str]]
from bs4 import BeautifulSoup
import requests


def parse_rules(text: str) -> RulesData:
    """parse_rules(text: str) -> RulesData:
       
       parse the rules from the NCBI webpage, returns a list of lists, which each inner list having the
       elements:

       prefix - a list of accession prefixes that this rule applies to
       database - name of the database or databases associated with this rule
       type_description - description of the type of data associated with this rule
       """
    soup = BeautifulSoup(text, 'html.parser')
    data = []
    for table in soup.find_all('table', cellpadding='3'):
        for row in table.find_all('tr'):
            (prefix, database, type_description, _) = row.text.split('\n')
            if prefix.strip() == 'Prefix':
                continue
            if ',' in prefix:
                prefix = [p.strip() for p in prefix.split(',')]
            else:
                prefix = [prefix]
            data.append((prefix, database, 'unknown', type_description))
    return data


def parse_refseq_rules(text: str) -> RulesData:
    """parse_refseq_rules(text: str) -> RulesData
    
    Parse the rules for NCBI RefSeq (thanks to Torsten Seemann for pointing them out)"""
    database = 'NCBI'
    soup = BeautifulSoup(text, 'html.parser')
    table = soup.find('div', id='ch18.T.refseq_accession_numbers_and_mole').table
    data = []
    for row in table.tbody.find_all('tr'):
        prefix = row.td.text
        molecule_type = row.td.next_sibling.text
        type_description = row.td.next_sibling.next_sibling.text
        data.append(([prefix], database, molecule_type, type_description))
    return data


def fetch(url: str = 'https://www.ncbi.nlm.nih.gov/Sequin/acc.html') -> str:
    """fetch(url:str) -> str

    Fetches the accession type description page (by default from https://www.ncbi.nlm.nih.gov/Sequin/acc.html)"""
    response = requests.get(url)
    if response.status_code != requests.codes['ok']:
        raise requests.exceptions.RequestException(f'Failed to download {url}', response=response)
    return response.text


def save_data(data: RulesData, output: TextIO) -> None:
    """save_data(data: RulesData, str, str], output: TextIO)

    Saves the data from parsing the accession description page to a JSON format file. output must be a
    file open for writing."""
    json.dump(data, output, indent=2)
    output.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Save rules from NCBI website')
    parser.add_argument('output_file', type=argparse.FileType('w'))
    args = parser.parse_args()
    data = parse_rules(fetch())
    refseq_data = parse_refseq_rules(fetch(url='https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly'))
    data.extend(refseq_data)
    save_data(data, args.output_file)
