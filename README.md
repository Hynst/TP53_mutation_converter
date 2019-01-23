# TP53 mutation converter

This is wrapper script of functions from python package hgvs (https://www.ncbi.nlm.nih.gov/pubmed/30129167). hgvs is used to parse, format, validate, normalize, and map biological sequence variants according to recommendations of the Human Genome Variation Society.

To use this script follow instructions.

## Clone git repository
Clone this repository to your local directory (no need for installation)

## Instal python packages
Install following packages to your Python2.7 enviroment (python2+ is required fro hgvs package)

1) hgvs (https://hgvs.readthedocs.io/en/stable/installation.html)
2) pandas (http://pandas.pydata.org/pandas-docs/stable/install.html)

## Run script
Run script as:
```
python2.7 /path/to/script_directory/Convert_HGVS.py \
  input_table.txt \
  conversion_direction_tag \ # string "toHGVS" or "fromHGVS"
  output_name_prefix # in example "example_output"
```
Input table should have same header names as examples inputs in git repository
