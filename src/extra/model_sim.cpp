/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <time.h>
#include <fstream>
#include "../basic/sequence.h"
#include "../basic/shape_config.h"
#include "extra.h"

struct MSS_model
{
	bool hit(const Letter *x, const Letter *y, size_t l) const
	{
		for (unsigned s = 0; s < shapes.count(); ++s) {
			if (l < shapes[s].length_) continue;
			for (unsigned i = 0; i <= l - shapes[s].length_; ++i) {
				uint64_t a, b;
				shapes[s].set_seed(a, &x[i]);
				shapes[s].set_seed(b, &y[i]);
				if (a == b)
					return true;
			}
		}
		return false;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length());
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i))
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i))
				return true;
		}
		return false;
	}
};

struct MSS_twohit_model
{
	unsigned hit(const Letter *x, const Letter *y, size_t l) const
	{
		unsigned n = 0;
		for (unsigned s = 0; s < shapes.count(); ++s) {
			if (l < shapes[s].length_) continue;
			for (unsigned i = 0; i <= l - shapes[s].length_; ++i) {
				uint64_t a, b;
				shapes[s].set_seed(a, &x[i]);
				shapes[s].set_seed(b, &y[i]);
				if (a == b)
					++n;
			}
		}
		return n;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length()) >= 2;
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		unsigned n = 0;
		for (unsigned i = 0; i < l; ++i) {
			n += hit(&x[i], y.data(), l - i);
			if (i > 0)
				n += hit(x.data(), &y[i], l - i);
		}
		return n >= 2;
	}
};

struct MSS_twohit_diagonal_model
{
	unsigned hit(const Letter *x, const Letter *y, size_t l) const
	{
		unsigned n = 0;
		for (unsigned s = 0; s < shapes.count(); ++s) {
			if (l < shapes[s].length_) continue;
			for (unsigned i = 0; i <= l - shapes[s].length_; ++i) {
				uint64_t a, b;
				shapes[s].set_seed(a, &x[i]);
				shapes[s].set_seed(b, &y[i]);
				if (a == b)
					++n;
			}
		}
		return n;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length()) >= 2;
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i) >= 2)
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i) >= 2)
				return true;
		}
		return false;
	}
};

template<int window>
struct MSS_twohit_window_model
{
	bool hit(const Letter *x, const Letter *y, size_t l) const
	{
		int last = -2 * window;
		for (int i = 0; i <= l; ++i) {
			for (unsigned s = 0; s < shapes.count(); ++s) {
				if (l - i < shapes[s].length_) continue;
				uint64_t a, b;
				shapes[s].set_seed(a, &x[i]);
				shapes[s].set_seed(b, &y[i]);
				if (a == b) {
					if (i != last && i - last < window)
						return true;
					last = i;
				}
			}
		}
		return false;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length());
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i))
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i))
				return true;
		}
		return false;
	}
};

template<int window, int k, int mismatch>
struct Mismatch_twohit_window_model
{
	bool hit(const Letter *x, const Letter *y, size_t l) const
	{
		int last = -2 * window;
		if (l < k)
			return false;
		for (int i = 0; i <= l - k; ++i) {
			unsigned n = 0;
			for (int j = 0; j < k; ++j)
				if (x[i + j] != y[i + j])
					++n;
			if (n <= mismatch) {
				if (i != last && i - last < window)
					return true;
				last = i;
			}
		}
		return false;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length());
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i))
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i))
				return true;
		}
		return false;
	}
};

template<int window>
struct MSS_twohit_window_neighborhood_model
{
	bool hit(const Letter *x, const Letter *y, size_t l) const
	{
		int last = -2 * window;
		for (int i = 0; i <= (int)l; ++i) {
			for (unsigned s = 0; s < shapes.count(); ++s) {
				if (l - i < shapes[s].length_) continue;
				if(shapes[s].score(&x[i], &y[i]) >= config.neighborhood_score) {
					if (i != last && i - last < window)
						return true;
					last = i;
				}
			}
		}
		return false;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length());
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i))
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i))
				return true;
		}
		return false;
	}
};

struct MSS_neighborhood_model
{
	bool hit(const Letter *x, const Letter *y, size_t l) const
	{
		for (size_t i = 0; i <= l; ++i) {
			for (unsigned s = 0; s < shapes.count(); ++s) {
				if (l - i < shapes[s].length_) continue;
				if (shapes[s].score(&x[i], &y[i]) >= config.neighborhood_score)
					return true;
			}
		}
		return false;
	}
	bool hit(const sequence &x, const sequence &y) const
	{
		return hit(x.data(), y.data(), x.length());
	}
	bool hit_shifted(const sequence &x, const sequence &y) const
	{
		const size_t l = x.length();
		for (unsigned i = 0; i < l; ++i) {
			if (hit(&x[i], y.data(), l - i))
				return true;
			if (i > 0 && hit(x.data(), &y[i], l - i))
				return true;
		}
		return false;
	}
};

void get_random_seq(vector<Letter> &seq)
{
	for (vector<Letter>::iterator i = seq.begin(); i != seq.end(); ++i)
		*i = get_distribution<20>(background_freq);
}

const double subst_freq[20][20] = {
	{ 0, 0.051599, 0.0319706, 0.0472072, 0.0259005, 0.042533, 0.0794012, 0.104447, 0.0190464, 0.0529913, 0.0806505, 0.0561607, 0.023709, 0.0243562, 0.0447981, 0.119439, 0.0710634, 0.00618423, 0.019181, 0.099362, },
	{ 0.0836991,0,0.0557221,0.0526752,0.00570959,0.0834397,0.0992451,0.0474531,0.0360983,0.0290852,0.0591073,0.18726,0.0187483,0.0177066,0.03043,0.062954,0.0596622,0.00896391,0.0252332,0.036807, },
	{ 0.0688047,0.0737048,0,0.134642,0.0104452,0.0556231,0.0811654,0.0970697,0.0503946,0.0214787,0.0370703,0.069479,0.0150052,0.0171298,0.0255345,0.114536,0.0729503,0.00573275,0.021037,0.0281966, },
	{ 0.0734774,0.0538258,0.10774,0,0.00541203,0.0638292,0.220103,0.0852837,0.0313663,0.0145824,0.0242017,0.0687326,0.00866069,0.0125848,0.0451219,0.0885509,0.0562425,0.0041692,0.0139378,0.0221784, },
	{ 0.186323,0.0324277,0.0277663,0.0177638,0,0.0153757,0.0222231,0.052226,0.0127764,0.0743205,0.11608,0.021819,0.0292451,0.0436977,0.0221956,0.0899075,0.0771679,0.00626418,0.0267835,0.125637, },
	{ 0.0815328,0.112847,0.0568952,0.074416,0.00530242,0,0.150972,0.0456165,0.0364876,0.0241864,0.0523059,0.119369,0.0224538,0.0139128,0.0272128,0.0645326,0.0511373,0.00542559,0.0193477,0.0360457, },
	{ 0.0918357,0.0766689,0.0556527,0.174685,0.00514321,0.105233,0,0.0497034,0.0284655,0.0215927,0.0364467,0.102734,0.0130324,0.0127301,0.0374955,0.0725947,0.0566717,0.00542622,0.0173958,0.0364932, },
	{ 0.20364,0.0508512,0.0818626,0.0839228,0.0149571,0.0404233,0.0622734,0,0.0241239,0.0237986,0.038645,0.0567032,0.015483,0.0207614,0.0397716,0.122147,0.0536346,0.00808141,0.0181211,0.0407985, },
	{ 0.0699983,0.0907505,0.0813351,0.0777822,0.0081228,0.0669371,0.0813018,0.0560113,0,0.0264461,0.0551068,0.0740711,0.0180528,0.03949,0.025928,0.0666684,0.0459882,0.0124538,0.0684967,0.0350589, },
	{ 0.0758071,0.0218901,0.0124197,0.0102327,0.0139658,0.0158577,0.0178954,0.0154911,0.00998824,0,0.261344,0.0185933,0.0553554,0.0532969,0.0164271,0.0248502,0.0451354,0.00741016,0.0241487,0.299892, },
	{ 0.0827941,0.0356979,0.0164472,0.0130886,0.0158172,0.0256513,0.0247873,0.0195319,0.0124665,0.218909,0,0.0302946,0.0895811,0.0919587,0.0199965,0.030285,0.046212,0.0136618,0.0364274,0.176391, },
	{ 0.0781029,0.208329,0.0570106,0.0577245,0.00610149,0.0855974,0.10978,0.0459911,0.0280067,0.0234476,0.0526726,0,0.0186556,0.014372,0.0345505,0.0674287,0.0530917,0.00620819,0.0164983,0.0364304, },
	{ 0.0820435,0.0371522,0.0193431,0.0143829,0.0129935,0.039314,0.0273425,0.02189,0.0123121,0.135883,0.247222,0.0313609,0,0.0639739,0.0174525,0.0360643,0.0483175,0.0103003,0.0273692,0.115283, },
	{ 0.0624752,0.023874,0.0168958,0.0146007,0.0120188,0.0162837,0.0221759,0.0282569,0.025181,0.107736,0.228084,0.018739,0.0498004,0,0.0143663,0.0337512,0.0352976,0.040237,0.152406,0.09782, },
	{ 0.137934,0.0566305,0.0383101,0.0798596,0.00858529,0.0448655,0.091169,0.0649277,0.0214561,0.0353504,0.0589727,0.0747337,0.0158675,0.0192857,0,0.0958676,0.0712531,0.00607841,0.0187647,0.0600885, },
	{ 0.161678,0.0503686,0.0675213,0.0725064,0.0170918,0.0439576,0.0688274,0.0920208,0.023245,0.0261139,0.039281,0.0556799,0.0163333,0.0207061,0.0417712,0,0.140562,0.00529629,0.0168799,0.0401595, },
	{ 0.116544,0.0543551,0.0502298,0.051343,0.0162852,0.0415458,0.0696321,0.0441794,0.0200499,0.0539524,0.072698,0.0528105,0.0237411,0.0226826,0.0364311,0.149828,0,0.00676344,0.0208007,0.0961279, },
	{ 0.0616437,0.051267,0.0204352,0.0156793,0.0128913,0.0281531,0.0336645,0.0438026,0.0318158,0.0565597,0.135468,0.0329638,0.0367357,0.153548,0.018681,0.0425701,0.0372178,0,0.126995,0.0599093, },
	{ 0.0629104,0.0481246,0.0279043,0.0248727,0.0132862,0.0319342,0.0336743,0.0255073,0.0636279,0.0693312,0.123129,0.03326,0.0336581,0.1923,0.0183125,0.0405417,0.039996,0.0424192,0,0.0752103, },
	{ 0.122866,0.0265131,0.0140538,0.0140485,0.020662,0.0200589,0.0307184,0.022685,0.010223,0.253759,0.190633,0.0272914,0.0436694,0.0420276,0.0241076,0.034372,0.0705984,0.00696795,0.0247446,0, } };

void get_related_seq(const sequence &seq, vector<Letter> &out, double id)
{
	for (unsigned i = 0; i < seq.length(); ++i)
		if ((double)rand() / RAND_MAX < id)
			out[i] = seq[i];
		else
			out[i] = get_distribution<20>(subst_freq[(size_t)seq[i]]);
}

template<typename _model>
void model_sim()
{
	static const unsigned seq_len = 300;
	static const double id = 0.25;
	srand((unsigned)time(0));

	vector<Letter> query (seq_len), subject (seq_len);
	const sequence q(query), s(subject);
	
	size_t n = 0, hit_bg = 0, hit_fg = 0;

	while (true) {
		if (n % 100 == 0)
			get_random_seq(query);

		get_random_seq(subject);
		/*cout << q << endl;
		cout << s << endl << endl;*/
		if (_model().hit_shifted(q, s))
			++hit_bg;

		get_related_seq(q, subject, id);
		if (_model().hit(q, s))
			++hit_fg;

		if (n % 100 == 0) {
			cout << "n=" << n << " background hit%=" << percentage<double,size_t>(hit_bg, n) <<  " foreground hit%=" << percentage<double, size_t>(hit_fg, n) << endl;
		}
		++n;
	}
}

void model_sim()
{
	//Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST"); // murphy.10
	Reduction::reduction = Reduction("A K R E D N Q C G H I L V M F Y W P S T");
	model_sim<MSS_model>();
}

void model_seqs()
{
	static const unsigned seq_len = 300;
	static const double id = 0.25;
	vector<Letter> seq1(seq_len), seq2(seq_len);
	sequence s1(seq1), s2(seq2);
	const size_t count = atoi(config.seq_no[0].c_str());
	std::ofstream out(config.output_file.c_str());
	for (unsigned i = 0; i < count; ++i) {
		get_random_seq(seq1);
		out << '>' << i << endl;
		out << s1 << endl;
		get_related_seq(s1, seq2, id);
		cout << '>' << i << endl;
		cout << s2 << endl;
	}
}