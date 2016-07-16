/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "../basic/sequence.h"
#include "../data/sequence_set.h"
#include "../dp/floating_sw.h"
#include "../util/Timer.h"
#include "../dp/dp.h"
#include "../align/align.h"

void benchmark_sw()
{
	static const unsigned n = 1;

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

	local_match hsp2(0, 0, &ss[1][0]);
	//greedy_align(ss[0], hsp2);
	
	Timer t;
	t.start();

	for (unsigned i = 0; i < n; ++i) {

		greedy_align(ss[0], ss[1]);

	}
	t.stop();

	cout << " n/sec=" << (double)n / t.getElapsedTimeInSec() << endl;
	return;

	uint64_t cell_updates = 0;
	local_match hsp(0, 0, &ss[1][16]);
	
	{
		Timer t;
		t.start();

		for (unsigned i = 0; i < n; ++i) {

			floating_sw(&ss[0][24],
				hsp.subject_,
				hsp,
				32,
				config.xdrop,
				config.gap_open + config.gap_extend,
				config.gap_extend,
				cell_updates,
				hsp.query_anchor_,
				hsp.subject_anchor,
				Score_only());

		}
		t.stop();

		cout << hsp.score << ' ' << cell_updates << endl;
		cout << "gcups=" << (double)cell_updates / 1e9 / t.getElapsedTimeInSec() << " n/sec=" << (double)n / t.getElapsedTimeInSec() << endl;
	}

	/*Query= d1mpxa2 c.69.1.21 (A:24-404) Alpha-amino acid ester hydrolase
	{Xanthomonas citri [TaxId: 346]}

	> sp|Q9L9D7|COCE_RHOSM Cocaine esterase OS=Rhodococcus sp. (strain
	MB1 Bresler) GN=cocE PE=1 SV=1
	Length=574

	Score = 94.0 bits (232),  Expect = 1e-019
	Identities = 103/380 (27%), Positives = 157/380 (41%), Gaps = 51/380 (13%)

	Query  20   NDYIKREVMIPMRDGVKLHTVIVLPKGAKNAPIVLTRTPYDASGRTERLA-SPHMKDLLS  78
	            N  +   VM+PMRDGV+L   +  P      P++L R PYD   + +  A S    + L
	Sbjct  5    NYSVASNVMVPMRDGVRLAVDLYRPDADGPVPVLLVRNPYD---KFDVFAWSTQSTNWLE  61

	Query  79   AGDDVFVEGGYIRVFQDVRGKYGSEGDYVMTRPLRGPLNPSEVDHATDAWDTIDWLVKNV  138
	                 FV  GY  V QD RG + SEG++V             VD   DA DT+ W+++
	Sbjct  62   -----FVRDGYAVVIQDTRGLFASEGEFV-----------PHVDDEADAEDTLSWILEQ-  104

	Query  139  SESNGKVGMIGSSYEGFTVVMALTNPHPALKVAVPESPMIDGWMGDDWFNYGAFRQVNFD  198
	            +  +G VGM G SY G T   A  +    LK   P     D +    W  YG    ++ +
	Sbjct  105  AWCDGNVGMFGVSYLGVTQWQAAVSGVGGLKAIAPSMASADLYRA-PW--YGPGGALSVE  161

	Query  199  YFTGQLSKRGKGAGIARQG--HDDYSNFLQ-AGSAGDFAKAAGLEQL----------PW-  244
	               G  +  G G   +R     +D ++F+Q A    D A AA +  L          PW
	Sbjct  162  ALLGWSALIGTGLITSRSDARPEDAADFVQLAAILNDVAGAASVTPLAEQPLLGRLIPWV  221

	Query  245  WHKLTEHAAYDAFWQEQALDKVMA--RTPLKVPTMWLQGLWDQEDMWGAIHSYAAMEPRD  302
	              ++ +H   D  WQ  +L + +    TP  +   W  G   +     ++ ++ A+
	Sbjct  222  IDQVVDHPDNDESWQSISLFERLGGLATPALITAGWYDGFVGE-----SLRTFVAV----  272

	Query  303  KRNTLNYLVMGPWRHSQVNYDGSALGALNFEGDTARQFRHDVLRPFFDQYL-VDGAPKAD  361
	            K N    LV+GPW HS +    +A            Q    + + FFD++L  +    A
	Sbjct  273  KDNADARLVVGPWSHSNLT-GRNADRKFGIAATYPIQEATTMHKAFFDRHLRGETDALAG  331

	Query  362  TPPVFIYNTGENHWDRLKAW  381
	             P V ++  G + W     W
	Sbjct  332  VPKVRLFVMGIDEWRDETDW  351

	*/

}