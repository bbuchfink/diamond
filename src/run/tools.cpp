/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <set>
#include <iostream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>
#include "tools.h"
#include "basic/config.h"
#include "data/sequence_set.h"
#include "data/sequence_file.h"
#include "masking/masking.h"
#include "basic/packed_transcript.h"
#include "data/dmnd/dmnd.h"
#include "util/sequence/sequence.h"
#include "basic/match.h"

using std::thread;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;
using std::runtime_error;
using std::vector;
using std::string;

void get_seq()
{
	config.database.require();
	SequenceFile* db_file = SequenceFile::auto_create({ config.database });
	db_file->get_seq();
	delete db_file;
}

void random_seqs()
{
	DatabaseFile db_file(config.database);
	Block* ref_seqs = db_file.load_seqs(std::numeric_limits<size_t>::max());
	const auto& r = ref_seqs->seqs();
	cout << "Sequences = " << r.size() << endl;
	std::set<unsigned> n;
	const size_t count = atoi(config.seq_no[0].c_str());
	while (n.size() < count)
		n.insert((rand()*RAND_MAX + rand()) % r.size());
	OutputFile out(config.output_file);
	unsigned j = 0;
	
	std::string s;
	for (std::set<unsigned>::const_iterator i = n.begin(); i != n.end(); ++i) {
		std::stringstream ss;
		ss << '>' << j++ << endl;
		if (config.reverse)
			; // ref_seqs::get()[*i].print(ss, value_traits, sequence::Reversed());
		else
			ss << r[*i];
		ss << endl;
		s = ss.str();
		out.write(s.data(), s.length());
	}
	out.close();
	delete ref_seqs;
}

void run_masker()
{
	TextInputFile f(config.single_query_file());
	vector<Letter> seq, seq2;
	string id;
	//const FASTA_format format;
	size_t letters = 0, seqs = 0, seqs_total = 0;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//while (format.get_seq(id, seq, f, value_traits)) {
	while(true) {
		cout << '>' << id << endl;
		seq2 = seq;
		Masking::get()(seq2.data(), seq2.size(), MaskingAlgo::TANTAN, 0);
		/*for (size_t i = 0; i < seq.size(); ++i) {
			char c = value_traits.alphabet[(long)seq[i]];
			if (seq2[i] == value_traits.mask_char)
				c = tolower(c);
			cout << c;
		}
		cout << endl;*/
		cout << Sequence(seq2.data(), seq2.size()) << endl;
		size_t n = 0;
		for (size_t i = 0; i < seq2.size(); ++i)
			if (seq2[i] == value_traits.mask_char) {
				++n;
			}
		letters += n;
		if (n > 0)
			++seqs;
		++seqs_total;
	}
	cerr << "#Sequences: " << seqs << "/" << seqs_total << ", #Letters: " << letters << ", t=" << (double)duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - t1).count()/1000.0 << endl;
}

void fastq2fasta()
{
	unique_ptr<TextInputFile> f(new TextInputFile(config.single_query_file()));
	vector<Letter> seq;
	string id;
//	const FASTQ_format format;
	input_value_traits = value_traits = nucleotide_traits;
	size_t n = 0, max = atoi(config.seq_no[0].c_str());
	//while (n < max && format.get_seq(id, seq, *f, value_traits)) {
	while(true) {
		cout << '>' << id << endl;
		cout << Sequence(seq.data(), seq.size()) << endl;
		++n;
	}
}

void info()
{
	vector<std::string> arch_flags;
#ifdef __SSE2__
	arch_flags.push_back("sse2");
#endif
#ifdef __SSE3__
	arch_flags.push_back("sse3");
#endif
#ifdef __SSSE3__
	arch_flags.push_back("ssse3");
#endif
#ifdef __POPCNT__
	arch_flags.push_back("popcnt");
#endif
#ifdef __ARM_NEON
	arch_flags.push_back("neon");
#endif

	cout << "Architecture flags: ";
	for (vector<std::string>::const_iterator i = arch_flags.begin(); i != arch_flags.end(); ++i)
		cout << *i << ' ';
	cout << endl;
}

void pairwise_worker(TextInputFile *in, std::mutex *input_lock, std::mutex *output_lock) {
	//FASTA_format format;
	std::string id_r, id_q;
	vector<Letter> ref, query;

	while(true) {
		input_lock->lock();
		/*if (!format.get_seq(id_r, ref, *in, value_traits)) {
			input_lock->unlock();
			return;
		}
		if (!format.get_seq(id_q, query, *in, value_traits)) {
			input_lock->unlock();
			return;
		}*/
		input_lock->unlock();
		const std::string ir = Util::Seq::seqid(id_r.c_str()), iq = Util::Seq::seqid(id_q.c_str());
		Hsp hsp(true);
		//smith_waterman(Sequence(query), Sequence(ref), hsp);
		HspContext context(hsp, 0, 0, TranslatedSequence(query), "", 0, 0, nullptr, 0, 0, Sequence());
		HspContext::Iterator it = context.begin();
		std::stringstream ss;
		while (it.good()) {
			if (it.op() == op_substitution)
				ss << ir << '\t' << iq << '\t' << it.subject_pos << '\t' << it.query_pos.translated << '\t' << it.query_char() << endl;
			else if (it.op() == op_deletion)
				ss << ir << '\t' << iq << '\t' << it.subject_pos << '\t' << "-1" << '\t' << '-' << endl;
			++it;
		}
		output_lock->lock();
		cout << ss.str();
		output_lock->unlock();
	}
}

void pairwise()
{
	input_value_traits = nucleotide_traits;
	value_traits = nucleotide_traits;
	score_matrix = ScoreMatrix("DNA", 5, 2, 0, 1);

	TextInputFile in(config.single_query_file());
	std::mutex input_lock, output_lock;
	vector<thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(pairwise_worker, &in, &input_lock, &output_lock);
	for (auto &t : threads)
		t.join();
}

void reverse() {
	input_value_traits = amino_acid_traits;
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	TextBuffer buf;
	//while (FASTA_format().get_seq(id, seq, in, value_traits)) {
	while(true) {
		buf << '>';
		buf << '\\';
		buf.write_raw(id.data(), id.size());
		buf << '\n';
		Sequence(seq).print(buf, amino_acid_traits, Sequence::Reversed());
		buf << '\n';
		buf << '\0';
		cout << buf.data();
		buf.clear();
	}
	in.close();
}