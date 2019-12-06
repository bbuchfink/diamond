#include <iostream>
#include <string>
#include "../util/io/temp_file.h"
#include "../util/io/text_input_file.h"
#include "test.h"
#include "../util/sequence/sequence.h"
#include "../util/log_stream.h"

using std::endl;
using std::string;

namespace Test {

void run() {
	const size_t prot_count = 1000, prot_len = 300, homolog_count = 19;
	std::minstd_rand0 rand_engine;
	std::uniform_real_distribution<double> id_dist(0.2, 1.0);

	task_timer timer("Generating test dataset");
	TempFile proteins;
	for (size_t i = 0; i < prot_count; ++i) {
		auto v = generate_random_seq(prot_len, rand_engine);
		sequence seq(v);
		const string id = std::to_string(i);
		Util::Sequence::format(seq, id.c_str(), nullptr, proteins, "fasta", amino_acid_traits);
		for (size_t j = 0; j < homolog_count; ++j) {
			const string id = std::to_string(i) + '_' + std::to_string(j);
			Util::Sequence::format(sequence(simulate_homolog(seq, id_dist(rand_engine), rand_engine)), id.c_str(), nullptr, proteins, "fasta", amino_acid_traits);
		}
	}
	InputFile query_file(proteins);
	query_file.close_and_delete();
}

}