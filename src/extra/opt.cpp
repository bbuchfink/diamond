/****
Copyright (c) 2017, Benjamin Buchfink
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

#include <string.h>
#include <time.h>
#include "../util/util.h"
#include "../basic/shape_config.h"
#include "extra.h"
#include "../util/log_stream.h"
#include "../util/util.h"
#include "../util/thread.h"

double Trail::background_p() const
{
	double p = 0;
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
			if (bucket[i] == bucket[j])
				p += background_freq[i] * background_freq[j];
	return pow(p, shapes[0].weight_);
}

std::ostream& operator<<(std::ostream &s, const Trail &t)
{
	const int buckets = t.buckets();
	for (int i = 0; i < buckets; ++i) {
		for (int j = 0; j < 20; ++j)
			if (t.bucket[j] == i)
				s << value_traits.alphabet[j];
		s << ' ';
	}
	return s;
}

struct Trails
{

	Trails()
	{
		for (int i = 0; i < 20; ++i)
			for (int j = 0; j < 20; ++j)
				pheromone[i][j] = 1.0;
	}

	Trail get() const
	{
		Trail t;
		int next = 0;
		double p[20];
		while (true) {

			double sum = pheromone[next][next];
			for (int i = next + 1; i < 20; ++i)
				if (t.bucket[i] == -1)
					sum += pheromone[next][i] * std::max(score_matrix(next, i) + 1, 1);

			memset(p, 0, sizeof(p));
			p[next] = pheromone[next][next] / sum;
			for (int i = next + 1; i < 20; ++i)
				if (t.bucket[i] == -1)
					p[i] = (pheromone[next][i] * std::max(score_matrix(next, i) + 1, 1)) / sum;

			int i = get_distribution<20>(p);

			if (i == next) {
				const int bucket = t.bucket[next];
				next = t.next();
				if (next != -1)
					t.bucket[next] = bucket + 1;
				else
					break;
			}
			else {
				t.bucket[i] = t.bucket[next];
				next = i;
			}

		}
		return t;
	}

	void update(const Trail &t, double sens)
	{
		const int buckets = t.buckets();
		vector<vector<int> > v(buckets);
		for (int i = 0; i < 20; ++i)
			v[t.bucket[i]].push_back(i);
		for (int i = 0; i < buckets; ++i) {
			int j;
			for (j = 0; j < v[i].size() - 1; ++j)
				pheromone[v[i][j]][v[i][j + 1]] += sens;
			pheromone[v[i][j]][v[i][j]] += sens;
		}
	}

	double pheromone[20][20];
};

const size_t max_ants = 100000;
Trail ants[max_ants];
double sens[max_ants];
Trails trails;

void get_sens_worker(vector<char>::const_iterator query, vector<char>::const_iterator query_end, vector<char>::const_iterator subject, vector<double> *sens)
{
	vector<bool> hit;

	for (; query < query_end; query += 70, subject += 70) {
		hit.insert(hit.begin(), config.n_ants, false);
		for (size_t j = 0; j <= 70 - shapes[0].length_; ++j) {
			for (size_t k = 0; k < config.n_ants; ++k)
				if (shapes[0].hit(&query[j], &subject[j], ants[k]) && !hit[k]) {
					++((*sens)[k]);
					hit[k] = true;
				}
		}
		hit.clear();
	}
}

void get_sens(const vector<char> &query, const vector<char> &subject)
{
	const size_t n_seqs = query.size() / 70;
	partition<size_t> p(n_seqs, config.threads_);
	vector<vector<double> > v(config.threads_);
	Thread_pool threads;
	for (unsigned i = 0; i < config.threads_; ++i) {
		v[i].resize(config.n_ants);
		threads.push_back(launch_thread(get_sens_worker, query.begin() + p.getMin(i) * 70, query.begin() + p.getMax(i) * 70, subject.begin() + p.getMin(i) * 70, &v[i]));
	}
	threads.join_all();
	memset(sens, 0, sizeof(sens));
	for (unsigned i = 0; i < config.threads_; ++i)
		for (unsigned j = 0; j < config.n_ants; ++j)
			sens[j] += v[i][j];

	for (size_t k = 0; k < config.n_ants; ++k)
		sens[k] /= n_seqs;
}

void opt()
{
	static const size_t region = 70;
	static const size_t count = (size_t)1e6;
	static const double id = 0.25;
	
	srand((unsigned)time(0));
	
	task_timer timer("Init");
	vector<char> query(count*region);
	get_random_seq(query);

	vector<char> subject(count*region);
	get_related_seq(sequence(query), subject, id);

	timer.go("Calculating sensitivity");
	ants[0] = Trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	get_sens(query, subject);
	timer.finish();

	double p_bg = ants[0].background_p();
	cout << "Sensitivity = " << sens[0] << endl;
	cout << "P(background) = " << p_bg << endl;	

	while (true) {
		timer.go("Setting ants");
		for (size_t i = 0; i < config.n_ants; ++i)
			ants[i] = trails.get();
		
		timer.go("Getting sensitivity");
		get_sens(query, subject);
		
		double max_sens_eff = 0;
		size_t max_ant;
		for (size_t i = 0; i < config.n_ants; ++i) {
			const double e = sens[i] * std::min(p_bg / ants[i].background_p(), 1.0);
			if (e > max_sens_eff) {
				max_sens_eff = e;
				max_ant = i;
			}
		}
		timer.finish();
		cout << "Effective sensitivity = " << max_sens_eff << endl;
		cout << "Sensitivity = " << sens[max_ant] << endl;
		cout << "P(background) = " << ants[max_ant].background_p() << endl;
		cout << ants[max_ant] << endl << endl;;

		trails.update(ants[max_ant], max_sens_eff);
	}

}