# diamondpy
This is python binding for diamond--a sequence alignment BLAST software implements in C++.

## install

- from source code
```bash
git clone https://github.com/GCS-ZHN/diamondpy.git
cd diamondpy
conda create -n diamondpy python=3.8
# zlib is required by diamond.
conda install -c anaconda zlib
pip install -v .
```

- from pypi
```bahs
pip install diamondpy
```

## Usage
```python
from diamondpy import Diamond

# create a object
diamond = Diamond(
    database="database.dmnd",
    n_threads=4
)

# make db if you don't create it or just download one from websites
diamond.makedb("database.fasta")
print(diamond.version)

# print database statistic info
diamond.dbinfo()

# run blast for proteins
diamond.blastp(
    query="test_proteins.fasta",
    out="test_output"
)
```