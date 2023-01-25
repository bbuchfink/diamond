/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "value.h"
#include "shape_config.h"
#include "../util/sequence/translate.h"
#include "statistics.h"
#include "sequence.h"
#include "../masking/masking.h"
#include "../util/util.h"
#include "../stats/standard_matrix.h"

const char* Const::version_string = "2.1.0";
using std::string;
using std::vector;
using std::count;


const char* Const::program_name = "diamond";

AlignMode::AlignMode(unsigned mode) :
	mode(mode)
{
	sequence_type = SequenceType::amino_acid;
	switch (mode) {
	case blastx:
		input_sequence_type = SequenceType::nucleotide;
		query_contexts = 6;
		query_translated = true;
		query_len_factor = 3;
		break;
    case blastn:
        input_sequence_type = SequenceType::nucleotide;
        query_translated = false;
        query_contexts = 2;
        query_len_factor = 1;
        sequence_type = SequenceType::nucleotide;
        break;
	default:
		input_sequence_type = SequenceType::amino_acid;
		query_contexts = 1;
		query_translated = false;
		query_len_factor = 1;
	}
}

unsigned AlignMode::from_command(unsigned command)
{
	switch (command) {
	case Config::blastx:
		return blastx;
    case Config::blastn:
        return blastn;
	default:
		return blastp;
	}
}

AlignMode align_mode(AlignMode::blastp);

Statistics statistics;
ShapeConfig shapes;
unsigned shape_from, shape_to;

const Letter Translator::reverseLetter[5] = { 3, 2, 1, 0, 4 };

Letter Translator::lookup[5][5][5];
Letter Translator::lookupReverse[5][5][5];

const char* Translator::codes[] = {
	0,
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 1
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", // 2
	"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 3
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 4
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", // 5
	"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 6
	0,
	0,
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 9
	"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 10
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 11
	"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 12
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", // 13
	"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 14
	0,
	"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 16
	0,
	0,
	0,
	0,
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", // 21
	"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 22
	"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 23
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", // 24
	"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 25
	"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" // 26
};

void Translator::init(unsigned id)
{
	static const unsigned idx[] = { 2, 1, 3, 0 };
	if (id >= sizeof(codes) / sizeof(codes[0]) || codes[id] == 0)
		throw std::runtime_error("Invalid genetic code id.");
	for (unsigned i = 0; i < 5; ++i)
		for (unsigned j = 0; j < 5; ++j)
			for (unsigned k = 0; k < 5; ++k)
				if (i == 4 || j == 4 || k == 4) {
					lookup[i][j][k] = value_traits.mask_char;
					lookupReverse[i][j][k] = value_traits.mask_char;
				}
				else {
					lookup[i][j][k] = value_traits.from_char(codes[id][(int)idx[i] * 16 + (int)idx[j] * 4 + (int)idx[k]]);
					lookupReverse[i][j][k] = value_traits.from_char(codes[id][idx[(int)reverseLetter[i]] * 16 + idx[(int)reverseLetter[j]] * 4 + idx[(int)reverseLetter[k]]]);
				}
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 4; ++j) {
			if (count(lookup[i][j], lookup[i][j] + 4, lookup[i][j][0]) == 4)
				lookup[i][j][4] = lookup[i][j][0];
			if (count(lookupReverse[i][j], lookupReverse[i][j] + 4, lookupReverse[i][j][0]) == 4)
				lookupReverse[i][j][4] = lookupReverse[i][j][0];
		}
}

vector<Letter> Sequence::from_string(const char* str, const ValueTraits&vt)
{
	vector<Letter> seq;
	while (*str)
		seq.push_back(vt.from_char(*(str++)));
	return seq;
}

void Seed::enum_neighborhood(unsigned pos, int treshold, vector<Seed>& out, int score)
{
	const size_t l = data_[pos];
	score -= score_matrix(l, l);
	for (size_t i = 0; i < 20; ++i) {
		int new_score = score + score_matrix(l, i);
		data_[pos] = (Letter)i;
		if (new_score >= treshold) {
			if (pos < config.seed_weight - 1)
				enum_neighborhood(pos + 1, treshold, out, new_score);
			else
				out.push_back(*this);
		}
	}
	data_[pos] = (Letter)l;
}

void Seed::enum_neighborhood(int treshold, vector<Seed>& out)
{
	out.clear();
	enum_neighborhood(0, treshold, out, score(*this));
}

std::vector<Letter> Sequence::reverse() const {
	std::vector<Letter> v;
	v.reserve(len_);
	std::reverse_copy(data(), end(), std::back_inserter(v));
	return v;
}

void Statistics::print() const
{
	using std::endl;
	//log_stream << "Used ref size = " << data_[REF_SIZE] << endl;
	//log_stream << "Traceback errors = " << data_[BIAS_ERRORS] << endl;
	//log_stream << "Low complexity seeds  = " << data_[LOW_COMPLEXITY_SEEDS] << endl;
	log_stream << "Hits (filter stage 0) = " << data_[SEED_HITS] << endl;
	log_stream << "Hits (filter stage 1) = " << data_[TENTATIVE_MATCHES1] << " (" << data_[TENTATIVE_MATCHES1] * 100.0 / data_[SEED_HITS] << " %)" << endl;
	log_stream << "Hits (filter stage 2) = " << data_[TENTATIVE_MATCHES2] << " (" << data_[TENTATIVE_MATCHES2] * 100.0 / data_[TENTATIVE_MATCHES1] << " %)" << endl;
	log_stream << "Hits (filter stage 3) = " << data_[TENTATIVE_MATCHES3] << " (" << data_[TENTATIVE_MATCHES3] * 100.0 / data_[TENTATIVE_MATCHES2] << " %)" << endl;
	//log_stream << "Hits (filter stage 4) = " << data_[TENTATIVE_MATCHES4] << " (" << data_[TENTATIVE_MATCHES4] * 100.0 / data_[TENTATIVE_MATCHES3] << " %)" << endl;
	log_stream << "Target hits (stage 0) = " << data_[TARGET_HITS0] << endl;
	log_stream << "Target hits (stage 1) = " << data_[TARGET_HITS1] << endl;
	log_stream << "Target hits (stage 2) = " << data_[TARGET_HITS2] << endl;
	log_stream << "Target hits (stage 3) = " << data_[TARGET_HITS3] << " (" << data_[TARGET_HITS3_CBS] << " (" << (double)data_[TARGET_HITS3_CBS] * 100.0 / data_[TARGET_HITS3] << "%) with CBS)" << endl;
	log_stream << "Target hits (stage 4) = " << data_[TARGET_HITS4] << endl;
	log_stream << "Target hits (stage 5) = " << data_[TARGET_HITS5] << endl;
	log_stream << "Target hits (stage 6) = " << data_[TARGET_HITS6] << endl;
	log_stream << "Swipe realignments    = " << data_[SWIPE_REALIGN] << endl;
	if (data_[MASKED_LAZY])
		log_stream << "Lazy maskings         = " << data_[MASKED_LAZY] << endl;
	log_stream << "Matrix adjusts        = " << data_[MATRIX_ADJUST_COUNT] << endl;
	log_stream << "Extensions (8 bit)    = " << data_[EXT8] << endl;
	log_stream << "Extensions (16 bit)   = " << data_[EXT16] << endl;
	log_stream << "Extensions (32 bit)   = " << data_[EXT32] << endl;
	log_stream << "Overflows (8 bit)     = " << data_[EXT_OVERFLOW_8] << endl;
	log_stream << "Wasted (16 bit)       = " << data_[EXT_WASTED_16] << endl;
	log_stream << "Effort (Extension)    = " << 2 * data_[EXT16] + data_[EXT8] << endl;
	log_stream << "Effort (Cells)        = " << 2 * data_[DP_CELLS_16] + data_[DP_CELLS_8] << endl;
	log_stream << "Cells (8 bit)         = " << data_[DP_CELLS_8] << endl;
	log_stream << "Cells (16 bit)        = " << data_[DP_CELLS_16] << endl;
	log_stream << "SWIPE tasks           = " << data_[SWIPE_TASKS_TOTAL] << endl;
	log_stream << "SWIPE tasks (async)   = " << data_[SWIPE_TASKS_ASYNC] << endl;
	log_stream << "Trivial aln           = " << data_[TRIVIAL_ALN] << endl;
	log_stream << "Hard queries          = " << data_[HARD_QUERIES] << endl;
#ifdef DP_STAT
	log_stream << "Gross DP Cells        = " << data_[GROSS_DP_CELLS] << endl;
	log_stream << "Net DP Cells          = " << data_[NET_DP_CELLS] << " (" << data_[NET_DP_CELLS] * 100.0 / data_[GROSS_DP_CELLS] << " %)" << endl;
#endif
	log_stream << "Gapped filter (targets) = " << data_[GAPPED_FILTER_TARGETS] << endl;
	log_stream << "Gapped filter (hits) stage 1 = " << data_[GAPPED_FILTER_HITS1] << endl;
	log_stream << "Gapped filter (hits) stage 2 = " << data_[GAPPED_FILTER_HITS2] << endl;
	log_stream << "Time (Load seed hit targets) = " << (double)data_[TIME_LOAD_HIT_TARGETS] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Sort targets by score) = " << (double)data_[TIME_SORT_TARGETS_BY_SCORE] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Gapped filter)         = " << (double)data_[TIME_GAPPED_FILTER] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Matrix adjust)         = " << (double)data_[TIME_MATRIX_ADJUST] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Chaining)              = " << (double)data_[TIME_CHAINING] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (DP target sorting)     = " << (double)data_[TIME_TARGET_SORT] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Query profiles)        = " << (double)data_[TIME_PROFILE] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Smith Waterman)        = " << (double)data_[TIME_SW] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Anchored SWIPE Alloc)  = " << (double)data_[TIME_ANCHORED_SWIPE_ALLOC] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Anchored SWIPE Sort)   = " << (double)data_[TIME_ANCHORED_SWIPE_SORT] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Anchored SWIPE Add)    = " << (double)data_[TIME_ANCHORED_SWIPE_ADD] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Anchored SWIPE Output) = " << (double)data_[TIME_ANCHORED_SWIPE_OUTPUT] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Anchored SWIPE)        = " << (double)data_[TIME_ANCHORED_SWIPE] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Smith Waterman TB)     = " << (double)data_[TIME_TRACEBACK_SW] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Smith Waterman-32)     = " << (double)data_[TIME_EXT_32] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Traceback)             = " << (double)data_[TIME_TRACEBACK] / 1e6 << "s (CPU)" << endl;
	log_stream << "Time (Target parallel)       = " << (double)data_[TIME_TARGET_PARALLEL] / 1e6 << "s (wall)" << endl;
	log_stream << "Time (Load seed hits)        = " << (double)data_[TIME_LOAD_SEED_HITS] / 1e6 << "s (wall)" << endl;
	log_stream << "Time (Sort seed hits)        = " << (double)data_[TIME_SORT_SEED_HITS] / 1e6 << "s (wall)" << endl;
	log_stream << "Time (Extension)             = " << (double)data_[TIME_EXT] / 1e6 << "s (wall)" << endl;
	//log_stream << "Time (greedy extension)      = " << data_[TIME_GREEDY_EXT]/1e9 << "s" << endl;
	//log_stream << "Gapped hits = " << data_[GAPPED_HITS] << endl;
	//log_stream << "Overlap hits = " << data_[DUPLICATES] << endl;
	//log_stream << "Secondary hits = " << data_[SECONDARY_HITS] << endl;
	//log_stream << "Erased hits = " << data_[ERASED_HITS] << endl;
	//log_stream << "High similarity hits = " << data_[HIGH_SIM] << endl;
	//log_stream << "Net hits = " << data_[OUT_HITS] << endl;
	//log_stream << "Matches = " << data_[OUT_MATCHES] << endl;
	//log_stream << "Total score = " << data_[SCORE_TOTAL] << endl;
	//log_stream << "Aligned query len = " << data_[ALIGNED_QLEN] << endl;
	//log_stream << "Gapped matches = " << data_[GAPPED] << endl;
	//log_stream << "MSE = " << (double)data_[SQUARED_ERROR] / (double)data_[OUT_HITS] << endl;
	//log_stream << "Cells = " << data_[CELLS] << endl;
	verbose_stream << "Temporary disk space used (search): " << (double)data_[SEARCH_TEMP_SPACE] / (1 << 30) << " GB" << endl;
	message_stream << "Reported " << data_[PAIRWISE] << " pairwise alignments, " << data_[MATCHES] << " HSPs." << endl;
	message_stream << data_[ALIGNED] << " queries aligned." << endl;
}

Reduction::Reduction(const char* definition_string)
{
	memset(map_, 0, sizeof(map_));
	memset(map8_, 0, sizeof(map8_));
	memset(map8b_, 0, sizeof(map8b_));
	map_[(long)MASK_LETTER] = MASK_LETTER;
	map_[(long)STOP_LETTER] = MASK_LETTER;
	const vector<string> tokens(tokenize(definition_string, " "));
	size_ = (unsigned)tokens.size();
	bit_size_exact_ = log(size_) / log(2);
	bit_size_ = (int)ceil(bit_size_exact_);
	freq_.fill(0.0);
	for (unsigned i = 0; i < size_; ++i)
		for (unsigned j = 0; j < tokens[i].length(); ++j) {
			const char ch = tokens[i][j];
			const long letter = (long)value_traits.from_char(ch);
			map_[letter] = i;
			map8_[letter] = i;
			map8b_[letter] = i;
			freq_[i] += Stats::blosum62.background_freqs[letter];
		}
	for (double& f : freq_)
		f = std::log(f);
	map8_[(long)MASK_LETTER] = (Letter)size_;
	map8_[(long)STOP_LETTER] = (Letter)size_;
	map8_[(long)DELIMITER_LETTER] = (Letter)size_;
	map8b_[(long)MASK_LETTER] = (Letter)(size_ + 1);
	map8b_[(long)STOP_LETTER] = (Letter)(size_ + 1);
	map8b_[(long)DELIMITER_LETTER] = (Letter)(size_ + 1);
}

std::string Reduction::decode_seed(const uint64_t seed, const size_t len) const {
	string s(len, '-');
	uint64_t c = seed;
	for (size_t i = 0; i < len; ++i) {
		s[len - i - 1] = amino_acid_traits.alphabet[c % size_];
		c /= size_;
	}
	return s;
}