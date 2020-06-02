/****
DIAMOND protein aligner
Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)
Code developed by Patrick Ettenhuber <patrick.ettenhuber@qiagen.com>

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

#include <map>
#include <string>
#include "../basic/score_matrix.h"
#include "../basic/match.h"

class Variable{
public:
	virtual std::string get_name() = 0;
	virtual double get(const Hsp_context& r) = 0;
};

class QueryLength: public Variable{
public:
	std::string get_name(){
		return "qlen";
	}
	double get(const Hsp_context& r){
		return r.query.source().length();
	}
};
class SubjectLength: public Variable{
public:
	std::string get_name(){
		return "slen";
	}
	double get(const Hsp_context& r){
		return r.subject_len;
	}
};
class QueryStart: public Variable{
public:
	std::string get_name(){
		return "qstart";
	}
	double get(const Hsp_context& r){
		return r.oriented_query_range().begin_ + 1;
	}
};
class QueryEnd: public Variable{
public:
	std::string get_name(){
		return "qend";
	}
	double get(const Hsp_context& r){
		return r.oriented_query_range().end_ + 1;
	}
};
class SubjectStart: public Variable{
public:
	std::string get_name(){
		return "sstart";
	}
	double get(const Hsp_context& r){
		return r.subject_range().begin_ + 1;
	}
};
class SubjectEnd: public Variable{
public:
	std::string get_name(){
		return "send";
	}
	double get(const Hsp_context& r){
		return r.subject_range().end_;
	}
};
class EValue: public Variable{
public:
	std::string get_name(){
		return "evalue";
	}
	double get(const Hsp_context& r){
		return r.evalue();
	}
};
class BitScore: public Variable{
public:
	std::string get_name(){
		return "bitscore";
	}
	double get(const Hsp_context& r){
		return r.bit_score();
	}
};
class Score: public Variable{
public:
	std::string get_name(){
		return "score";
	}
	double get(const Hsp_context& r){
		return r.score();
	}
};
class Length: public Variable{
public:
	std::string get_name(){
		return "length";
	}
	double get(const Hsp_context& r){
		return r.length();
	}
};
class PercentIdenticalMatches: public Variable{
public:
	std::string get_name(){
		return "pident";
	}
	double get(const Hsp_context& r){
		return (double)r.identities() * 100 / r.length();
	}
};
class NumberIdenticalMatches: public Variable{
public:
	std::string get_name(){
		return "nident";
	}
	double get(const Hsp_context& r){
		return r.identities();
	}
};
class NumberMismatches: public Variable{
public:
	std::string get_name(){
		return "mismatch";
	}
	double get(const Hsp_context& r){
		return r.mismatches();
	}
};
class NumberPositiveMatches: public Variable{
public:
	std::string get_name(){
		return "positive";
	}
	double get(const Hsp_context& r){
		return r.positives();
	}
};
class NumberGapOpenings: public Variable{
public:
	std::string get_name(){
		return "gapopen";
	}
	double get(const Hsp_context& r){
		return r.gap_openings();
	}
};
class NumberGaps: public Variable{
public:
	std::string get_name(){
		return "gaps";
	}
	double get(const Hsp_context& r){
		return r.gaps();
	}
};
class PercentagePositiveMatches: public Variable{
public:
	std::string get_name(){
		return "ppos";
	}
	double get(const Hsp_context& r){
		return (double)r.positives() * 100.0 / r.length();
	}
};
class QueryFrame: public Variable{
public:
	std::string get_name(){
		return "qframe";
	}
	double get(const Hsp_context& r){
		return r.blast_query_frame();
	}
};
class QueryCoveragePerHsp: public Variable{
public:
	std::string get_name(){
		return "qcovhsp";
	}
	double get(const Hsp_context& r){
		return (double)r.query_source_range().length()*100.0 / r.query.source().length();
	}
};
class SwDiff: public Variable{
public:
	std::string get_name(){
		return "swdiff";
	}
	double get(const Hsp_context& r){
		return r.sw_score() - r.bit_score();
	}
};
class Time: public Variable{
public:
	std::string get_name(){
		return "time";
	}
	double get(const Hsp_context& r){
		return r.time();
	}
};
class SubjectCoveragePerHsp: public Variable{
public:
	std::string get_name(){
		return "scovhsp";
	}
	double get(const Hsp_context& r){
		return (double)r.subject_range().length() * 100.0 / r.subject_len;
	}
};
class UngappedScore: public Variable{
public:
	std::string get_name(){
		return "ungapped_score";
	}
	double get(const Hsp_context& r){
		return score_matrix.bitscore(r.ungapped_score);
	}
};

class VariableRegistry{
private:
	VariableRegistry(){};
	// To include new variables add the instantiation here and in the clustering_variable.cpp file. Then add it to the StaticConstructor below
	static QueryLength queryLength;
	static SubjectLength subjectLength;
	static QueryStart queryStart;
	static QueryEnd queryEnd;
	static SubjectStart subjectStart;
	static SubjectEnd subjectEnd;
	static EValue eValue;
	static BitScore bitScore;
	static Score score;
	static Length length;
	static PercentIdenticalMatches percentIdenticalMatches;
	static NumberIdenticalMatches numberIdenticalMatches;
	static NumberMismatches numberMismatches;
	static NumberPositiveMatches numberPositiveMatches;
	static NumberGapOpenings numberGapOpenings;
	static NumberGaps numberGaps;
	static PercentagePositiveMatches percentagePositiveMatches;
	static QueryFrame queryFrame;
	static QueryCoveragePerHsp queryCoveragePerHsp;
	static SwDiff swdiff;
	static Time time;
	static SubjectCoveragePerHsp subjectCoveragePerHsp;
	static UngappedScore ungappedScore;
public:
	static std::map<std::string, Variable*> regMap;
	static Variable* get(std::string key){
		std::map<std::string, Variable*>::iterator ca = VariableRegistry::regMap.find(key);
		if(ca == VariableRegistry::regMap.end()){
			throw std::runtime_error(std::string("Unknown variable: ")+ key);
		}
		return ca->second;
	}
	static bool has(std::string key){
		return VariableRegistry::regMap.find(key) != VariableRegistry::regMap.end();
	}
	static vector<std::string> getKeys(){
		auto it = regMap.begin();
		vector<std::string> keys;
		while(it != regMap.end()){
			keys.push_back(it->first);
			it++;
		}
		return keys;
	}
	static struct StaticConstructor {
		StaticConstructor() {
			regMap.emplace(queryLength.get_name(), &queryLength);
			regMap.emplace(subjectLength.get_name(), &subjectLength);
			regMap.emplace(queryStart.get_name(), &queryStart);
			regMap.emplace(queryEnd.get_name(), &queryEnd);
			regMap.emplace(subjectStart.get_name(), &subjectStart);
			regMap.emplace(subjectEnd.get_name(), &subjectEnd);
			regMap.emplace(eValue.get_name(), &eValue);
			regMap.emplace(bitScore.get_name(), &bitScore);
			regMap.emplace(score.get_name(), &score);
			regMap.emplace(length.get_name(), &length);
			regMap.emplace(percentIdenticalMatches.get_name(), &percentIdenticalMatches);
			regMap.emplace(numberIdenticalMatches.get_name(), &numberIdenticalMatches);
			regMap.emplace(numberMismatches.get_name(), &numberMismatches);
			regMap.emplace(numberPositiveMatches.get_name(), &numberPositiveMatches);
			regMap.emplace(numberGapOpenings.get_name(), &numberGapOpenings);
			regMap.emplace(numberGaps.get_name(), &numberGaps);
			regMap.emplace(percentagePositiveMatches.get_name(), &percentagePositiveMatches);
			regMap.emplace(queryFrame.get_name(), &queryFrame);
			regMap.emplace(queryCoveragePerHsp.get_name(), &queryCoveragePerHsp);
			regMap.emplace(swdiff.get_name(), &swdiff);
			regMap.emplace(time.get_name(), &time);
			regMap.emplace(subjectCoveragePerHsp.get_name(), &subjectCoveragePerHsp);
			regMap.emplace(ungappedScore.get_name(), &ungappedScore);
		}
	} _staticConstructor;
};
