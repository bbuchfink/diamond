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
#pragma once
#include "clustering_variables.h"
#include <cmath>
#include <cstdlib>
#include <functional>
#include <stdexcept>

class RecursiveParser {
	const HspContext* r; 
	const char * expression_to_parse;
	std::vector<Variable*> vars_;
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
	std::function<bool(double a, double b)> relation(){
		if(peek() == '>' && peek(1) == '='){
			advance(2);
			return [](double a, double b){return a>=b;};
		}
		else if(peek() == '<' && peek(1) == '='){
			advance(2);
			return [](double a, double b){return a<=b;};
		}
		else if(peek() == '=' && peek(1) == '='){
			advance(2);
			return [](double a, double b){return a==b;};
		}
		else if(peek() == '='){
			advance(1);
			return [](double a, double b){return a==b;};
		}
		else if(peek() == '>'){
			advance(1);
			return [](double a, double b){return a>b;};
		}
		else if(peek() == '<'){
			advance(1);
			return [](double a, double b){return a<b;};
		}
		std::runtime_error(std::string("Error while evaluating the expression ")+std::string(expression_to_parse)+std::string(" Unknown relation: ")+peek()+peek(1));
		return [](double a, double b){ return false; };
	}
	std::string variable()
	{
		std::vector<char> v;
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
		else if (peek() == 'l' && peek(1)=='o' && peek(2)=='g' && peek(3)=='('){
			advance(4); // 'log('
			double result1 = expression();
			advance(1); // ')'
			return std::log(result1);
		}
		else if (peek() == 'I' && peek(1)=='('){
			advance(2); // 'I('
			double result1 = expression();
			std::function<bool(double, double)> rel = relation();
			double result2 = expression();
			advance(1); // ')'
			return rel(result1, result2);
		}
		else {
			auto v = VariableRegistry::get(variable());
			if(!r) {
				vars_.push_back(v);
				return 4;
			}
			else {
				return v->get(*r);
			}
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
	RecursiveParser(const HspContext* r, const char * c): r(r), expression_to_parse(c)
	{}
	double evaluate() {
		return expression();
	}
	std::vector<Variable*> variables() const {
		return vars_;
	}
	const std::string static clean_expression(const std::string* const expression){
		std::string cleanedExpression = std::string(*expression);
		cleanedExpression.erase(std::remove_if(cleanedExpression.begin(), cleanedExpression.end(), [](unsigned char c){return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\v' || c == '\f';}), cleanedExpression.end());
		return cleanedExpression;
	}
};

