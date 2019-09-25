import sys,os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import pytest

from parse_accession import build_accession_parser, match_accession

ACCESSION_RULES_FILENAME = os.path.realpath(os.path.dirname(__file__)+"/../accession_rules.json")


@pytest.fixture()
def matching_rules():
    rules_by_prefix_len = build_accession_parser(open(ACCESSION_RULES_FILENAME))
    return rules_by_prefix_len


def test_match_accession(matching_rules):
    (database, accession_type, description) = match_accession('AE014297', matching_rules)
    assert database == 'GenBank'
    assert accession_type == 'nucleotide'
    assert description == 'Genome project data'

def test_match_refseq(matching_rules):
    (database, accession_type, description) = match_accession('NC_000962', matching_rules)
    assert database == 'NCBI'
    assert accession_type == 'Genomic'
    assert description == 'RefSeq: Complete genomic molecule, usually reference assembly'

def test_match_sra(matching_rules):
    (database, accession_type, description) = match_accession('ERR5413122', matching_rules)
    assert database == 'EMBL'
    assert accession_type == 'SRA'
    assert description == 'ERA run accession: A Run is an object that contains actual sequencing data for a particular sequencing experiment. Experiments may contain many Runs depending on the number of sequencing instrument runs that were needed.'    