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

#include <string.h>
#include <time.h>
#include "../util/util.h"
#include "../basic/shape_config.h"
#include "extra.h"
#include "../util/log_stream.h"
#include "../util/util.h"
#include "../util/thread.h"

double Letter_trail::background_p() const
{
	double p = 0;
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
			if (bucket[i] == bucket[j])
				p += background_freq[i] * background_freq[j];
	return p;
}

double Letter_trail::foreground_p(double id) const
{
	double p = 0;
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
			if (bucket[i] == bucket[j] && i != j)
				p += background_freq[i] * subst_freq[i][j];
	return id + (1 - id)*p;
}

double background_p(const Trail& t)
{
	double p = 1.0;
	for (int pos = 0; pos < OPT_W; ++pos)
		p *= t[pos].background_p();
	return p;
}

std::ostream& operator<<(std::ostream &s, const Letter_trail &t)
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

double tau_min, tau_max;

void set_limits(double &x)
{
	x = std::min(std::max(x, tau_min), tau_max);
	//x = std::max(x, tau_min);
}

const size_t max_ants = 10000;
Trail ants[max_ants];
double sens[max_ants];

struct Trails
{

	Trails()
	{
		for (int i = 0; i < 20; ++i)
			for (int j = 0; j < 20; ++j)
				for (int p = 0; p < OPT_W; ++p)
					pheromone[p][i][j] = 100.0;
	}

	double delta_tau(int pos, int i, int j) const
	{
		//return pheromone[pos][i][j] * pow(subst_freq[i][j], config.d_exp);
		return pheromone[pos][i][j] * pow(std::max(score_matrix(i, j) + 1, 0), config.d_exp);
	}

	double delta_tau0(int pos, int i) const
	{
		return pheromone[pos][i][i] * pow(config.d_new, config.d_exp);
	}

	Letter_trail get(int pos) const
	{
		Letter_trail t;
		int next = 0;
		double p[20];
		while (true) {

			double sum = delta_tau0(pos, next);
			for (int i = next + 1; i < 20; ++i)
				if (t.bucket[i] == -1)
					sum += delta_tau(pos, next, i);

			memset(p, 0, sizeof(p));
			p[next] = delta_tau0(pos, next) / sum;
			for (int i = next + 1; i < 20; ++i)
				if (t.bucket[i] == -1)
					p[i] = delta_tau(pos, next, i) / sum;

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

	void get(Trail &out) const
	{
		for (int pos = 0; pos < OPT_W; ++pos)
			out[pos] = get(pos);
	}

	void update(const Letter_trail &t, int pos, double sens)
	{
		const int buckets = t.buckets();
		vector<vector<int> > v(buckets);
		for (int i = 0; i < 20; ++i)
			v[t.bucket[i]].push_back(i);
		for (int i = 0; i < buckets; ++i) {
			int j;
			for (j = 0; j < (int)v[i].size() - 1; ++j) {
				pheromone[pos][v[i][j]][v[i][j + 1]] += sens;
				set_limits(pheromone[pos][v[i][j]][v[i][j + 1]]);
			}
			pheromone[pos][v[i][j]][v[i][j]] += sens;
			set_limits(pheromone[pos][v[i][j]][v[i][j]]);
		}
	}

	void evaporate()
	{
		for (int pos = 0; pos < OPT_W; ++pos)
			for (int i = 0; i < 20; ++i)
				for (int j = 0; j < 20; ++j) {
					pheromone[pos][i][j] *= config.rho;
					set_limits(pheromone[pos][i][j]);
				}
	}

	void update(const Trail &t, double sens)
	{		
		for (int pos = 0; pos < OPT_W; ++pos)
			update(t[pos], pos, sens);
	}

	void update_all()
	{
		for (size_t i = 0; i < config.n_ants; ++i)
			update(ants[i], sens[i]);
	}

	double pheromone[OPT_W][20][20];
};

Trails trails;

void get_sens_worker(vector<char>::const_iterator query, vector<char>::const_iterator query_end, vector<char>::const_iterator subject, vector<double> *sens)
{
	vector<bool> hit;

	for (; query < query_end; query += 70, subject += 70) {
		hit.insert(hit.begin(), config.n_ants, false);
		for (size_t j = 0; j <= 70 - shapes[0].length_; ++j) {
			for (size_t k = 0; k < config.n_ants; ++k)
				if (shapes[1].hit(&query[j], &subject[j], ants[k]) && !hit[k]) {
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

bool hit(const sequence &seq, const vector<Letter> &v, const Trail &t)
{
	for (size_t j = 0; j <= 70 - shapes[0].length_; ++j)
		if (shapes[0].hit(&seq[j], &v[j], t))
			return true;
	return false;
}

void get_related_seq(const sequence &seq, vector<Letter> &out, double id, size_t region, const Trail &previous)
{
	vector<Letter> r(region);
	out.clear();
	for (size_t i = 0; i < seq.length(); i += region) {
		const sequence query(&seq[i], region);
		do {
			get_related_seq(query, r, id);
		} while (hit(query, r, previous));
		out.insert(out.end(), r.begin(), r.end());
	}
}

void opt()
{
	static const size_t region = 70;
	static const size_t count = (size_t)1e6;
	static const double id = 0.25;

	srand((unsigned)time(0));

	Trail previous;
	previous[0] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[1] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[2] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[3] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[4] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[5] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	previous[6] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));

	previous[0] = Letter_trail(Reduction("A K R E D N Q C G H I L V M F Y W P S T"));
	cout << previous[0].foreground_p(id) / previous[0].background_p() << endl;
	
	return;
	
	task_timer timer("Init");
	vector<char> query(count*region);
	get_random_seq(query);

	vector<char> subject(count*region);
	get_related_seq(sequence(query), subject, id, region, previous);

	timer.go("Calculating sensitivity");
	for (int pos = 0; pos < OPT_W; ++pos)
		ants[0][pos] = Letter_trail(Reduction("A KR EDNQ C G H ILVM FYW P ST"));
	get_sens(query, subject);
	timer.finish();
	
	double p_bg = background_p(ants[0]);
	cout << "Sensitivity = " << sens[0] << endl;
	cout << "P(background) = " << p_bg << endl;	
	double global_best = 0;

	while (true) {
		timer.go("Setting ants");
		for (size_t i = 0; i < config.n_ants; ++i)
			trails.get(ants[i]);
		
		timer.go("Getting sensitivity");
		get_sens(query, subject);
		
		double max_sens_eff = 0;
		size_t max_ant;
		for (size_t i = 0; i < config.n_ants; ++i) {
			const double e = sens[i] * std::min(p_bg / background_p(ants[i]), 1.0);
			sens[i] = e;
			if (e > max_sens_eff) {
				max_sens_eff = e;
				max_ant = i;
			}
		}
		timer.finish();
		global_best = std::max(global_best, max_sens_eff);
		tau_max = 1 / (1 - config.rho)*global_best;
		tau_min = tau_max*(1 - pow(config.p_best, 0.05)) / 9 / pow(config.p_best, 0.05);
		cout << "Effective sensitivity = " << max_sens_eff << ", global = " << global_best << endl;
		cout << "Sensitivity = " << sens[max_ant] << endl;
		cout << "P(background) = " << background_p(ants[max_ant]) << endl;
		cout << "tau_max = " << tau_max << " tau_min = " << tau_min << endl;
		for (int pos = 0; pos < OPT_W; ++pos)
			cout << ants[max_ant][pos] << endl;
		cout << endl;

		trails.evaporate();
		trails.update(ants[max_ant], max_sens_eff);
		//trails.update_all();
	}

}