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
#include "def.h"
#include "../dp/flags.h"
#include "def.h"

class Variable{
public:
	Variable(const HspValues hsp_values = HspValues::NONE, const Output::Flags flags = Output::Flags::NONE):
		hsp_values(hsp_values),
		flags(flags)
	{}
	virtual double get(const HspContext& r) const = 0;
	virtual ~Variable(){};
	const HspValues hsp_values;
	const Output::Flags flags;
};

class QueryLength: public Variable{
public:
	static const std::string get_name(){
		return "qlen";
	}
	virtual double get(const HspContext& r) const override {
		return r.query.source().length();
	}
};
class SubjectLength: public Variable{
public:
	static const std::string get_name(){
		return "slen";
	}
	virtual double get(const HspContext& r) const override {
		return r.subject_len;
	}
};
class QueryStart: public Variable{
public:
	QueryStart():
		Variable(HspValues::QUERY_START)
	{}
	static const std::string get_name(){
		return "qstart";
	}
	virtual double get(const HspContext& r) const override {
		return r.oriented_query_range().begin_ + 1;
	}
};
class QueryEnd: public Variable{
public:
	QueryEnd() :
		Variable(HspValues::QUERY_END)
	{}
	static const std::string get_name(){
		return "qend";
	}
	virtual double get(const HspContext& r) const override {
		return r.oriented_query_range().end_ + 1;
	}
};
class SubjectStart: public Variable{
public:
	SubjectStart() :
		Variable(HspValues::TARGET_START)
	{}
	static const std::string get_name(){
		return "sstart";
	}
	virtual double get(const HspContext& r) const override {
		return r.subject_range().begin_ + 1;
	}
};
class SubjectEnd: public Variable{
public:
	SubjectEnd() :
		Variable(HspValues::TARGET_END)
	{}
	static const std::string get_name(){
		return "send";
	}
	virtual double get(const HspContext& r) const override {
		return r.subject_range().end_;
	}
};
class EValue: public Variable{
public:
	static const std::string get_name(){
		return "evalue";
	}
	virtual double get(const HspContext& r) const override {
		return r.evalue();
	}
};
class BitScore: public Variable{
public:
	static const std::string get_name(){
		return "bitscore";
	}
	virtual double get(const HspContext& r) const override {
		return r.bit_score();
	}
};
class RawScore: public Variable{
public:
	static const std::string get_name(){
		return "score";
	}
	virtual double get(const HspContext& r) const override {
		return r.score();
	}
};
class Length: public Variable{
public:
	Length():
		Variable(HspValues::LENGTH)
	{}
	static const std::string get_name(){
		return "length";
	}
	virtual double get(const HspContext& r) const override {
		return r.length();
	}
};
class PercentIdenticalMatches: public Variable{
public:
	PercentIdenticalMatches():
		Variable(HspValues::LENGTH | HspValues::IDENT)
	{}
	static const std::string get_name(){
		return "pident";
	}
	virtual double get(const HspContext& r) const override {
		return (double)r.identities() * 100 / r.length();
	}
};
class NumberIdenticalMatches: public Variable{
public:
	NumberIdenticalMatches() :
		Variable(HspValues::IDENT)
	{}
	static const std::string get_name(){
		return "nident";
	}
	virtual double get(const HspContext& r) const override {
		return r.identities();
	}
};
class NumberMismatches: public Variable{
public:
	NumberMismatches() :
		Variable(HspValues::MISMATCHES)
	{}
	static const std::string get_name(){
		return "mismatch";
	}
	virtual double get(const HspContext& r) const override {
		return r.mismatches();
	}
};
class NumberPositiveMatches: public Variable{
public:
	NumberPositiveMatches() :
		Variable(HspValues::TRANSCRIPT)
	{}
	static const std::string get_name(){
		return "positive";
	}
	virtual double get(const HspContext& r) const override {
		return r.positives();
	}
};
class NumberGapOpenings: public Variable{
public:
	NumberGapOpenings() :
		Variable(HspValues::GAP_OPENINGS)
	{}
	static const std::string get_name(){
		return "gapopen";
	}
	virtual double get(const HspContext& r) const override {
		return r.gap_openings();
	}
};
class NumberGaps: public Variable{
public:
	NumberGaps() :
		Variable(HspValues::GAPS)
	{}
	static const std::string get_name(){
		return "gaps";
	}
	virtual double get(const HspContext& r) const override {
		return r.gaps();
	}
};
class PercentagePositiveMatches: public Variable{
public:
	PercentagePositiveMatches() :
		Variable(HspValues::TRANSCRIPT)
	{}
	static const std::string get_name(){
		return "ppos";
	}
	virtual double get(const HspContext& r) const override {
		return (double)r.positives() * 100.0 / r.length();
	}
};
class QueryFrame: public Variable{
public:
	static const std::string get_name(){
		return "qframe";
	}
	virtual double get(const HspContext& r) const override {
		return r.blast_query_frame();
	}
};
class QueryCoveragePerHsp: public Variable{
public:
	QueryCoveragePerHsp() :
		Variable(HspValues::QUERY_COORDS)
	{}
	static const std::string get_name(){
		return "qcovhsp";
	}
	virtual double get(const HspContext& r) const override {
		return (double)r.query_source_range().length()*100.0 / r.query.source().length();
	}
};
class SubjectCoveragePerHsp: public Variable{
public:
	SubjectCoveragePerHsp() :
		Variable(HspValues::TARGET_COORDS)
	{}
	static const std::string get_name(){
		return "scovhsp";
	}
	virtual double get(const HspContext& r) const override {
		return (double)r.subject_range().length() * 100.0 / r.subject_len;
	}
};

struct NormalizedBitScoreGlobal : public Variable {
	NormalizedBitScoreGlobal():
		Variable(HspValues::NONE, Output::Flags::SELF_ALN_SCORES)
	{}
	static const std::string get_name() {
		return "normalized_bitscore_global";
	}
	virtual double get(const HspContext& r) const override {
		return r.bit_score() / std::max(r.query_self_aln_score, r.target_self_aln_score) * 100;
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
		regMap[RawScore::get_name()] = new RawScore();
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
		regMap[NormalizedBitScoreGlobal::get_name()] = new NormalizedBitScoreGlobal();
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
	std::vector<std::string> getKeys() const {
		auto it = regMap.begin();
		std::vector<std::string> keys;
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
	static std::vector<std::string> getKeys(){
		return vr.getKeys();
	}
};
