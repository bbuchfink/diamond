/****
DIAMOND protein aligner
Copyright (C) 2020 Patrick Ettenhuber <pettenhuber@gmail.com>

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

#include "output_format.h"
#include "clustering_variables.h"


class RecursiveParser{
	const Hsp_context& r; 
	const char * expression_to_parse;
	char peek()
	{
		return *expression_to_parse;
	}
	char peek(uint32_t ahead)
	{
		return *(expression_to_parse+ahead);
	}
	char get()
	{
		return *expression_to_parse++;
	}
	void advance(uint32_t ahead)
	{
		expression_to_parse+=ahead;
	}
	std::string variable()
	{
		vector<char> v;
		while ( (peek() >= 'a' && peek() <= 'z') ||  (peek() >= 'A' && peek() <= 'Z') || peek() == '_')
		{
			v.push_back(get());
		}
		return std::string(v.begin(), v.end());
	}
	uint32_t integer()
	{
		uint32_t result = get() - '0';
		while (peek() >= '0' && peek() <= '9')
		{
			result = 10*result + get() - '0';
		}
		return result;
	}
	double number()
	{
		double result = integer();
		if(peek() == '.'){
			advance(1); // '.'
			const char *pos = expression_to_parse;
			double decimals = integer() * 1.0 / (expression_to_parse - pos);
			result += decimals;
		}
		return result;
	}

	double factor()
	{
		if (peek() >= '0' && peek() <= '9') {
			return number();
		}
		else if (peek() == '(') {
			advance(1); // '('
			double result = expression();
			advance(1); // ')'
			return result;
		}
		else if (peek() == '-') {
			advance(1);
			return -factor();
		}
		else if (peek() == 'm' && peek(1)=='a' && peek(2)=='x' && peek(3)=='('){
			advance(4); // 'max('
			double result1 = expression();
			advance(1); // ','
			double result2 = expression();
			advance(1); // ')'
			return std::max(result1, result2);
		}
		else if (peek() == 'm' && peek(1)=='i' && peek(2)=='n' && peek(3)=='('){
			advance(4); // 'min('
			double result1 = expression();
			advance(1); // ','
			double result2 = expression();
			advance(1); // ')'
			return std::min(result1, result2);
		}
		else if (peek() == 'e' && peek(1)=='x' && peek(2)=='p' && peek(3)=='('){
			advance(4); // 'exp('
			double result1 = expression();
			advance(1); // ')'
			return std::exp(result1);
		}
		else {
			return VariableRegistry::get(variable())->get(r);
		}
	}

	double term()
	{
		double result = factor();
		while (peek() == '*' || peek() == '/'){
			if (get() == '*'){
				result *= factor();
			}
			else{
				result /= factor();
			}
		}
		return result;
	}

	double expression() 
	{
		double result = term();
		while (peek() == '+' || peek() == '-'){
			if (get() == '+'){
				result += term();
			}
			else{
				result -= term();
			}
		}
		return result;
	}

public:
	static const vector<std::string> variables;
	RecursiveParser(const Hsp_context& r, const char * c): r(r), expression_to_parse(c){}
	double evaluate(){
		return expression();
	}
};


void Clustering_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const 
{
	out.write((uint32_t)query_num);
}

void Clustering_format::print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) 
{
	if (r.query_id <= r.subject_id) {
		out.write((uint32_t)r.subject_id);
		RecursiveParser rp(r, format->c_str());
		out.write(rp.evaluate());
	}
}
