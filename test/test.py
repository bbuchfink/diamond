from diamond4py import Diamond

diamond = Diamond(
    database="test/cafa4.dmnd",
    n_threads=4,
)
print(diamond.version)
diamond.dbinfo()
diamond.blastp(
    query="test/test_proteins.fasta",
    out="test/test_blastp_output"
)
# diamond.makedb("/home/hongning/AnnoPRO/test_proteins.fasta")
print("done")