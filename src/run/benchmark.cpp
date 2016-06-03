#include "../basic/sequence.h"
#include "../data/sequence_set.h"
#include "../dp/floating_sw.h"
#include "../util/Timer.h"

void benchmark_sw()
{
	static const unsigned n = 100000;

	/*
	> d2va1a_ c.73.1.0 (A:) automated matches {Ureaplasma parvum [TaxId: 
134821]}
Length=234

 Score = 26.2 bits (56),  Expect = 1.1
 Identities = 18/66 (27%), Positives = 28/66 (42%), Gaps = 4/66 (6%)

Query  24  QADATVATFFNGIDMPNQTNKTAA--FLCAALGGPNAWTGRNLKE--VHANMGVSNAQFT  79
           Q D+++  F    D+  Q  K +    +   LGG N W G   KE  +  N+  +     
Sbjct  16  QNDSSIIDFIKINDLAEQIEKISKKYIVSIVLGGGNIWRGSIAKELDMDRNLADNMGMMA  75

Query  80  TVIGHL  85
           T+I  L
Sbjct  76  TIINGL  81	*/

	vector<Letter> s1 = sequence::from_string("SLFEQLGGQAAVQAVTAQFYANIQADATVATFFNGIDMPNQTNKTAAFLCAALGGPNAWTGRNLKEVHANMGVSNAQFTTVIGHLRSALTGAGVAAALVEQTVAVAETVRGDVVTV");
	vector<Letter> s2 = sequence::from_string("RKQRIVIKISGACLKQNDSSIIDFIKINDLAEQIEKISKKYIVSIVLGGGNIWRGSIAKELDMDRNLADNMGMMATIINGLALENALNHLNVNTIVLSAIKCDKLVHESSANNIKKAIEKEQVMIFVAGTGFPYFTTDSCAAIRAAETESSIILMGKNGVDGVYDSDPKINPNAQFYEHITFNMALTQNLKVMDATALALCQENNINLLVFNIDKPNAIVDVLEKKNKYTIVSK");
	Sequence_set ss;
	ss.push_back(s1);
	ss.push_back(s2);
	ss.finish_reserve();

	unsigned long cell_updates = 0;
	local_match hsp(0, 0, &ss[1][16]);
	vector<char> transcript_buf;
	Timer t;
	t.start();

	for (unsigned i = 0; i < n; ++i) {

		floating_sw(&ss[0][24],
			hsp,
			32,
			config.xdrop,
			config.gap_open + config.gap_extend,
			config.gap_extend,
			transcript_buf,
			cell_updates,
			Score_only());

	}
	t.stop();

	cout << hsp.score_ << ' ' << cell_updates << endl;
	cout << (double)cell_updates / 1e9 / t.getElapsedTimeInSec();

}