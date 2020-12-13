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
#include "../stats/score_matrix.h"
#include "../basic/match.h"

class Variable{
public:
	virtual double get(const Hsp_context& r) = 0;
	virtual ~Variable(){};
};

class QueryLength: public Variable{
public:
	static const std::string get_name(){
		return "qlen";
	}
	double get(const Hsp_context& r){
		return r.query.source().length();
	}
};
class SubjectLength: public Variable{
public:
	static const std::string get_name(){
		return "slen";
	}
	double get(const Hsp_context& r){
		return r.subject_len;
	}
};
class QueryStart: public Variable{
public:
	static const std::string get_name(){
		return "qstart";
	}
	double get(const Hsp_context& r){
		return r.oriented_query_range().begin_ + 1;
	}
};
class QueryEnd: public Variable{
public:
	static const std::string get_name(){
		return "qend";
	}
	double get(const Hsp_context& r){
		return r.oriented_query_range().end_ + 1;
	}
};
class SubjectStart: public Variable{
public:
	static const std::string get_name(){
		return "sstart";
	}
	double get(const Hsp_context& r){
		return r.subject_range().begin_ + 1;
	}
};
class SubjectEnd: public Variable{
public:
	static const std::string get_name(){
		return "send";
	}
	double get(const Hsp_context& r){
		return r.subject_range().end_;
	}
};
class EValue: public Variable{
public:
	static const std::string get_name(){
		return "evalue";
	}
	double get(const Hsp_context& r){
		return r.evalue();
	}
};
class BitScore: public Variable{
public:
	static const std::string get_name(){
		return "bitscore";
	}
	double get(const Hsp_context& r){
		return r.bit_score();
	}
};
class Score: public Variable{
public:
	static const std::string get_name(){
		return "score";
	}
	double get(const Hsp_context& r){
		return r.score();
	}
};
class Length: public Variable{
public:
	static const std::string get_name(){
		return "length";
	}
	double get(const Hsp_context& r){
		return r.length();
	}
};
class PercentIdenticalMatches: public Variable{
public:
	static const std::string get_name(){
		return "pident";
	}
	double get(const Hsp_context& r){
		return (double)r.identities() * 100 / r.length();
	}
};
class NumberIdenticalMatches: public Variable{
public:
	static const std::string get_name(){
		return "nident";
	}
	double get(const Hsp_context& r){
		return r.identities();
	}
};
class NumberMismatches: public Variable{
public:
	static const std::string get_name(){
		return "mismatch";
	}
	double get(const Hsp_context& r){
		return r.mismatches();
	}
};
class NumberPositiveMatches: public Variable{
public:
	static const std::string get_name(){
		return "positive";
	}
	double get(const Hsp_context& r){
		return r.positives();
	}
};
class NumberGapOpenings: public Variable{
public:
	static const std::string get_name(){
		return "gapopen";
	}
	double get(const Hsp_context& r){
		return r.gap_openings();
	}
};
class NumberGaps: public Variable{
public:
	static const std::string get_name(){
		return "gaps";
	}
	double get(const Hsp_context& r){
		return r.gaps();
	}
};
class PercentagePositiveMatches: public Variable{
public:
	static const std::string get_name(){
		return "ppos";
	}
	double get(const Hsp_context& r){
		return (double)r.positives() * 100.0 / r.length();
	}
};
class QueryFrame: public Variable{
public:
	static const std::string get_name(){
		return "qframe";
	}
	double get(const Hsp_context& r){
		return r.blast_query_frame();
	}
};
class QueryCoveragePerHsp: public Variable{
public:
	static const std::string get_name(){
		return "qcovhsp";
	}
	double get(const Hsp_context& r){
		return (double)r.query_source_range().length()*100.0 / r.query.source().length();
	}
};
class SubjectCoveragePerHsp: public Variable{
public:
	static const std::string get_name(){
		return "scovhsp";
	}
	double get(const Hsp_context& r){
		return (double)r.subject_range().length() * 100.0 / r.subject_len;
	}
};
class UngappedScore: public Variable{
public:
	static const std::string get_name(){
		return "ungapped_score";
	}
	double get(const Hsp_context& r){
		return score_matrix.bitscore(r.ungapped_score);
	}
};

class StaticVariableRegistry{
	std::map<std::string, Variable*> regMap;
public:
	StaticVariableRegistry(){
		// To include new variables add the instantiation here and in the clustering_variable.cpp file. Then add it to the StaticConstructor below
		regMap[QueryLength::get_name()] = new QueryLength();
		regMap[SubjectLength::get_name()] = new SubjectLength();
		regMap[QueryStart::get_name()] = new QueryStart();
		regMap[QueryEnd::get_name()] = new QueryEnd();
		regMap[SubjectStart::get_name()] = new SubjectStart();
		regMap[SubjectEnd::get_name()] = new SubjectEnd();
		regMap[EValue::get_name()] = new EValue();
		regMap[BitScore::get_name()] = new BitScore();
		regMap[Score::get_name()] = new Score();
		regMap[Length::get_name()] = new Length();
		regMap[PercentIdenticalMatches::get_name()] = new PercentIdenticalMatches();
		regMap[NumberIdenticalMatches::get_name()] = new NumberIdenticalMatches();
		regMap[NumberMismatches::get_name()] = new NumberMismatches();
		regMap[NumberPositiveMatches::get_name()] = new NumberPositiveMatches();
		regMap[NumberGapOpenings::get_name()] = new NumberGapOpenings();
		regMap[NumberGaps::get_name()] = new NumberGaps();
		regMap[PercentagePositiveMatches::get_name()] = new PercentagePositiveMatches();
		regMap[QueryFrame::get_name()] = new QueryFrame();
		regMap[QueryCoveragePerHsp::get_name()] = new QueryCoveragePerHsp();
		regMap[SubjectCoveragePerHsp::get_name()] = new SubjectCoveragePerHsp();
		regMap[UngappedScore::get_name()] = new UngappedScore();
	}
	~StaticVariableRegistry(){
		for(auto it = regMap.begin(); it != regMap.end(); it++){
			delete it->second;
			it->second = nullptr;
		}
	}
	Variable* get(std::string key) const {
		auto ca = regMap.find(key);
		if(ca == regMap.end()){
			throw std::runtime_error(std::string("Unknown variable: ")+ key);
		}
		return ca->second;
	}
	bool has(std::string key) const {
		return regMap.find(key) != regMap.end();
	}
	vector<std::string> getKeys() const {
		auto it = regMap.begin();
		vector<std::string> keys;
		while(it != regMap.end()){
			keys.push_back(it->first);
			it++;
		}
		return keys;
	}
};

class VariableRegistry{
	static StaticVariableRegistry vr;
public:
	static Variable* get(std::string key){
		return vr.get(key);
	}
	static bool has(std::string key){
		return vr.has(key);
	}
	static vector<std::string> getKeys(){
		return vr.getKeys();
	}
};
