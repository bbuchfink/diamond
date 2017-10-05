@test "ARDB: Database creation" {
    ../diamond makedb --in inputs/ARDB.faa -d ARDB
}

@test "ARDB_small: Alignment serial" {
     ../diamond blastx -d ARDB -q inputs/small.fastq.gz -o ARDB_small_s.m8 -p 1
     echo "b4783ca2281c76907f5645ee2b535fc6 ARDB_small_s.m8" |  md5sum -c -
     rm ARDB_small_s.m8
}


@test "ARDB_small: Alignment parallel" {
     ../diamond blastx -d ARDB -q inputs/small.fastq.gz -o ARDB_small_p.m8
     echo "b4783ca2281c76907f5645ee2b535fc6 ARDB_small_p.m8" |  md5sum -c -
     rm ARDB_small_p.m8
}

@test "ARDB: Clean up" {
	rm -f ARDB.dmnd
}
