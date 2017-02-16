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
#include "../align/extend_ungapped.h"
#include "../search/sse_dist.h"
#include "../dp/score_profile.h"
#include "../output/output_format.h"

void benchmark_cmp()
{
#ifdef __SSE2__
	const size_t n = 1000000000llu;
	__m128i r1 = _mm_set_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
	const __m128i r2 = _mm_set_epi8(0, 2, 3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 13, 14, 0, 16);
	Timer t;
	t.start();
	unsigned x = 0;
	for (size_t i = 0; i < n; ++i) {
		r1 = _mm_set_epi32(x, x, x, x);
		//x += popcount32(_mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2)));
		x += _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
	}
	cout << "x=" << x << " t=" << t.getElapsedTimeInMicroSec() * 1000 / n << endl;
#endif
}

int xdrop_ungapped2(const Letter *query, const Letter *subject)
{
	int score(0), st(0);

	const Letter *q(query), *s(subject);

	st = score;
	while (score - st < config.raw_ungapped_xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;

		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;

		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
	}
	return score;
}

int xdrop_window(const Letter *query, const Letter *subject)
{
	static const int window = 40;
	int score(0), st(0), n = 0;

	const Letter *q(query), *s(subject);

	st = score;
	while (n < window
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		//score = std::max(score, st);
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		//score = std::max(score, st);
		++q;
		++s;
		++n;

		st += score_matrix(*q, *s);
		//score = std::max(score, st);
		++q;
		++s;
		++n;
	}
	return st;
}

void benchmark_ungapped(const Sequence_set &ss, unsigned qa, unsigned sa)
{
	static const size_t n = 100000000llu;
	Timer t;
	t.start();

	const Letter *q = &ss[0][qa], *s = &ss[1][sa];
	int score=0;

	uint64_t mask = 0;
	for (unsigned i = 0; i < 64; ++i) {
		if (q[i] == '\xff' || s[i] == '\xff')
			break;
		if (q[i] == s[i])
			mask |= 1llu << i;
	}

	for (size_t i = 0; i < n; ++i) {

		//score += xdrop_window(q, s);
		//score += binary_ungapped(mask);

	}
	t.stop();

	cout << score << endl;
	cout << "t=" << t.getElapsedTimeInMicroSec() * 1000 / n << " ns" << endl;
}

void benchmark_greedy(const Sequence_set &ss, unsigned qa, unsigned sa)
{
	static const unsigned n = 100000;
	vector<Diagonal_segment> d;
	d.push_back(ungapped_extension(sa, qa, ss[0], ss[1]));
	Long_score_profile qp(ss[0]);
	//greedy_align(ss[0], qp, ss[1], d[0], true);
	//greedy_align(ss[0], qp, ss[1], qa, sa, true);
	Hsp_data hsp;
	greedy_align2(ss[0], qp, ss[1], d, true, hsp);
	Text_buffer buf;
	Pairwise_format().print_match(Hsp_context(hsp, 0, ss[0], ss[0], "", 0, 0, "", 0, 0, 0), buf);
	buf << '\0';
	cout << buf.get_begin();

	Timer t;
	t.start();

	for (unsigned i = 0; i < n; ++i) {

		//greedy_align(ss[0], qp, ss[1], d[0], false);
		//greedy_align(ss[0], qp, ss[1], qa, sa, false);
		greedy_align2(ss[0], qp, ss[1], d, false, hsp);

	}
	t.stop();

	cout << " usec=" << t.getElapsedTimeInSec() / (double)n * 1000000.0 << endl;
	cout << "t=" << t.getElapsedTimeInMicroSec() << endl;
}

void benchmark_floating(const Sequence_set &ss, unsigned qa, unsigned sa)
{
	static const unsigned n = 10000;
	uint64_t cell_updates = 0;
	local_match hsp(0, 0, &ss[1][sa]);

	{
		Timer t;
		t.start();

		for (unsigned i = 0; i < n; ++i) {

			floating_sw(&ss[0][qa],
				hsp.subject_,
				hsp,
				32,
				score_matrix.rawscore(config.gapped_xdrop),
				config.gap_open + config.gap_extend,
				config.gap_extend,
				cell_updates,
				hsp.query_anchor_,
				hsp.subject_anchor,
				No_score_correction(),
				Score_only());

		}
		t.stop();

		cout << hsp.score << ' ' << cell_updates << endl;
		cout << "gcups=" << (double)cell_updates / 1e9 / t.getElapsedTimeInSec() << " n/sec=" << (double)n / t.getElapsedTimeInSec() << endl;
	}
}

void benchmark_sw()
{
	Sequence_set ss;
	vector<Letter> s1, s2;
	unsigned qa = 0, sa = 0;
	goto aln1;	
	
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


	s1 = sequence::from_string("SLFEQLGGQAAVQAVTAQFYANIQADATVATFFNGIDMPNQTNKTAAFLCAALGGPNAWTGRNLKEVHANMGVSNAQFTTVIGHLRSALTGAGVAAALVEQTVAVAETVRGDVVTV");
	s2 = sequence::from_string("RKQRIVIKISGACLKQNDSSIIDFIKINDLAEQIEKISKKYIVSIVLGGGNIWRGSIAKELDMDRNLADNMGMMATIINGLALENALNHLNVNTIVLSAIKCDKLVHESSANNIKKAIEKEQVMIFVAGTGFPYFTTDSCAAIRAAETESSIILMGKNGVDGVYDSDPKINPNAQFYEHITFNMALTQNLKVMDATALALCQENNINLLVFNIDKPNAIVDVLEKKNKYTIVSK");
	qa = 23;
	sa = 15;
	goto ende;	

	/*

	Query= d1g2na_ a.123.1.1 (A:) Ultraspiracle protein, usp {Tobacco budworm
	(Heliothis virescens) [TaxId: 7102]}

	> sp|Q6DHP9|RXRGB_DANRE Retinoic acid receptor RXR-gamma-B OS=Danio
	rerio GN=rxrgb PE=2 SV=1
	Length=452

	Score = 189 bits (479),  Expect = 6e-055
	Identities = 101/249 (41%), Positives = 153/249 (61%), Gaps = 24/249 (10%)

	Query  4    QELSIERLLEMESLVADPSEEFQFLRVGPDSNVPPKFRAPVSSLCQIGNKQIAALVVWAR  63
	            +++ ++++L+ E  V   +E +       +S+       PV+++C   +KQ+  LV WA+
	Sbjct  221  EDMPVDKILDAELSVEPKTETYT------ESSPSNSTNDPVTNICHAADKQLFTLVEWAK  274

	Query  64   DIPHFSQLEMEDQILLIKGSWNELLLFAIAWRSMEFLTEERDGVDGTGNRTTSPPQLMCL  123
	             IPHFS L ++DQ++L++  WNELL+ + + RS+      +DG+               L
	Sbjct  275  RIPHFSDLPLDDQVILLRAGWNELLIASFSHRSITV----KDGI--------------LL  316

	Query  124  MPGMTLHRNSALQAGVGQIFDRVLSELSLKMRTLRVDQAEYVALKAIILLNPDVKGLKNR  183
	              G+ +HR+SA  AGVG IF+RVL+EL  KM+ +++D+ E   L+AI+L NPD KGL N
	Sbjct  317  GTGLHVHRSSAHSAGVGSIFNRVLTELVSKMKDMQMDKTELGCLRAIVLFNPDAKGLSNS  376

	Query  184  QEVEVLREKMFLCLDEYCRRSRSSEEGRFAALLLRLPALRSISLKSFEHLFFFHLVADTS  243
	             EVE LREK++  L+ Y ++    + GRFA LLLRLPALRSI LK  EHLFFF L+ DT
	Sbjct  377  LEVEALREKVYASLETYTKQKYPDQPGRFAKLLLRLPALRSIGLKCLEHLFFFKLIGDTP  436

	Query  244  IAGYIRDAL  252
	            I  ++ + L
	Sbjct  437  IDTFLMEML  445
	*/


	s1 = sequence::from_string("aavqelsierllemeslvadpseefqflrvgpdsnvppkfrapvsslcqignkqiaalvv\
wardiphfsqlemedqillikgswnelllfaiawrsmeflteerdgvdgtgnrttsppql\
mclmpgmtlhrnsalqagvgqifdrvlselslkmrtlrvdqaeyvalkaiillnpdvkgl\
knrqevevlrekmflcldeycrrsrsseegrfaalllrlpalrsislksfehlfffhlva\
dtsiagyirdalrnha");

	s2 = sequence::from_string("MDTHDTYLHLHSSPLNSSPSQPPVMSSMVGHPSVISSSRPLPSPMSTLGSSMNGLPSPYS\
VITPSLSSPSISLPSTPSMGFNTLNSPQMNSLSMNGNEDIKPPPGLAPLGNMSSYQCTSP\
GSLSKHICAICGDRSSGKHYGVYSCEGCKGFFKRTIRKDLTYTCRDIKECLIDKRQRNRC\
QYCRYQKCLAMGMKREAVQEERQRGKEKSDTEVETTSRFNEDMPVDKILDAELSVEPKTE\
TYTESSPSNSTNDPVTNICHAADKQLFTLVEWAKRIPHFSDLPLDDQVILLRAGWNELLI\
ASFSHRSITVKDGILLGTGLHVHRSSAHSAGVGSIFNRVLTELVSKMKDMQMDKTELGCL\
RAIVLFNPDAKGLSNSLEVEALREKVYASLETYTKQKYPDQPGRFAKLLLRLPALRSIGL\
KCLEHLFFFKLIGDTPIDTFLMEMLEAPHQIT");

	qa = 3;
	sa = 220;

	goto ende;

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

aln1:

	s1 = sequence::from_string("tspmtpditgkpfvaadasndyikrevmipmrdgvklhtvivlpkgaknapivltrtpyd\
asgrterlasphmkdllsagddvfveggyirvfqdvrgkygsegdyvmtrplrgplnpse\
vdhatdawdtidwlvknvsesngkvgmigssyegftvvmaltnphpalkvavpespmidg\
wmgddwfnygafrqvnfdyftgqlskrgkgagiarqghddysnflqagsagdfakaagle\
qlpwwhkltehaaydafwqeqaldkvmartplkvptmwlqglwdqedmwgaihsyaamep\
rdkrntlnylvmgpwrhsqvnydgsalgalnfegdtarqfrhdvlrpffdqylvdgapka\
dtppvfiyntgenhwdrlkaw");

	s2 = sequence::from_string("MVDGNYSVASNVMVPMRDGVRLAVDLYRPDADGPVPVLLVRNPYDKFDVFAWSTQSTNWL\
EFVRDGYAVVIQDTRGLFASEGEFVPHVDDEADAEDTLSWILEQAWCDGNVGMFGVSYLG\
VTQWQAAVSGVGGLKAIAPSMASADLYRAPWYGPGGALSVEALLGWSALIGTGLITSRSD\
ARPEDAADFVQLAAILNDVAGAASVTPLAEQPLLGRLIPWVIDQVVDHPDNDESWQSISL\
FERLGGLATPALITAGWYDGFVGESLRTFVAVKDNADARLVVGPWSHSNLTGRNADRKFG\
IAATYPIQEATTMHKAFFDRHLRGETDALAGVPKVRLFVMGIDEWRDETDWPLPDTAYTP\
FYLGGSGAANTSTGGGTLSTSISGTESADTYLYDPADPVPSLGGTLLFHNGDNGPADQRP\
IHDRDDVLCYSTEVLTDPVEVTGTVSARLFVSSSAVDTDFTAKLVDVFPDGRAIALCDGI\
VRMRYRETLVNPTLIEAGEIYEVAIDMLATSNVFLPGHRIMVQVSSSNFPKYDRNSNTGG\
VIAREQLEEMCTAVNRIHRGPEHPSHIVLPIIKR");

	qa = 19;
	sa = 4;

	ende:
	ss.push_back(s1);
	ss.push_back(s2);
	ss.finish_reserve();

	//benchmark_floating(ss, qa, sa);
	benchmark_greedy(ss, qa, sa);
	//benchmark_cmp();
	//benchmark_ungapped(ss, qa, sa);

}