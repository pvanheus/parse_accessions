#!/usr/bin/env python

import argparse
import json
from typing import List, TextIO

from bs4 import BeautifulSoup
import requests


# def iterate_id_range(spec: str) -> List[str]:
#     (start, end) = spec.split('-')
#     assert end.endswith('Z'), f'Range specifier has unknown format: {spec}'
#     for i in range(len(start)):
#         if start[i] != end[i]:
#             break
#     prefix_length = i
#     assert prefix_length < len(start), f'Prefix not found in range specifier {spec}'
#     suffix_length = len(start) - prefix_length

def parse_rules(text: str) -> List[List[str]]:
    """parse_rules(text: str):
       parse the rules from the NCBI webpage"""
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
            data.append([prefix, database, type_description])
    return data


def fetch(url='https://www.ncbi.nlm.nih.gov/Sequin/acc.html') -> str:
    response = requests.get(url)
    if response.status_code == requests.codes['ok']:
        return response.text
    else:
        raise requests.exceptions.RequestException(f'Failed to download {url}', response=response)


def save_data(data: List[List[str]], output: TextIO) -> None:
    json.dump(data, output, indent=2)
    output.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Save rules from NCBI website')
    parser.add_argument('output_file', type=argparse.FileType('w'))
    args = parser.parse_args()
    save_data(parse_rules(fetch()), args.output_file)
