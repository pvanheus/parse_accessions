This project contains two sets of code:

Firstly, parser for the accession format described here: https://www.ncbi.nlm.nih.gov/Sequin/acc.html called 
`scrape_accession_rules.py`. This can be run as 

```
scrape_accession_rules.py accession_rules.json
```

to scrape the rules at that website and save them in JSON format. The file generated by this parser is
used in the second script, `parse_ncbi_accession.py` which can be used as a Python module but also
as a command line tool e.g. `parse_ncbi_accession.py accession_rules.json AE014297`.

The `parse_ncbi_accession.py` script aims to work with all Python versions. The `scrape_accession_rules.py`
script requires at least Python 3.6 and the modules specified in `conda_requirements.txt`.