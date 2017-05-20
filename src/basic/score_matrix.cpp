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

#include <string>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include "score_matrix.h"
#include "../blast/raw_scoremat.h"

using std::string;

Score_matrix score_matrix;

/** Number of statistical parameters in each row of the precomputed tables. */
#define BLAST_NUM_STAT_VALUES 11  /**< originally 8, now 11 to support Spouge's FSC. see notes below */

/** Holds values (gap-opening, extension, etc.) for a matrix. */
typedef double array_of_8[BLAST_NUM_STAT_VALUES];

const double INT2_MAX = std::numeric_limits<double>::max();

#define BLOSUM45_VALUES_MAX 14 /**< Number of different combinations supported for BLOSUM45. */
static array_of_8 blosum45_values[BLOSUM45_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7, 0.641318, 9.611060, 9.611060 },
	{ 13, 3, (double)INT2_MAX, 0.207, 0.049, 0.14, 1.5, -22, 0.671128, 35.855900, 35.963900 },
	{ 12, 3, (double)INT2_MAX, 0.199, 0.039, 0.11, 1.8, -34, 0.691530, 45.693600, 45.851700 },
	{ 11, 3, (double)INT2_MAX, 0.190, 0.031, 0.095, 2.0, -38, 0.691181, 62.874100, 63.103700 },
	{ 10, 3, (double)INT2_MAX, 0.179, 0.023, 0.075, 2.4, -51, 0.710529, 88.286800, 88.639100 },
	{ 16, 2, (double)INT2_MAX, 0.210, 0.051, 0.14, 1.5, -24, 0.666680, 36.279800, 36.452400 },
	{ 15, 2, (double)INT2_MAX, 0.203, 0.041, 0.12, 1.7, -31, 0.673871, 44.825700, 45.060400 },
	{ 14, 2, (double)INT2_MAX, 0.195, 0.032, 0.10, 1.9, -36, 0.685753, 60.736200, 61.102300 },
	{ 13, 2, (double)INT2_MAX, 0.185, 0.024, 0.084, 2.2, -45, 0.698480, 85.148100, 85.689400 },
	{ 12, 2, (double)INT2_MAX, 0.171, 0.016, 0.061, 2.8, -65, 0.713429, 127.758000, 128.582000 },
	{ 19, 1, (double)INT2_MAX, 0.205, 0.040, 0.11, 1.9, -43, 0.672302, 53.071400, 53.828200 },
	{ 18, 1, (double)INT2_MAX, 0.198, 0.032, 0.10, 2.0, -43, 0.682580, 72.342400, 73.403900 },
	{ 17, 1, (double)INT2_MAX, 0.189, 0.024, 0.079, 2.4, -57, 0.695035, 103.055000, 104.721000 },
	{ 16, 1, (double)INT2_MAX, 0.176, 0.016, 0.063, 2.8, -67, 0.712966, 170.100000, 173.003000 },
};  /**< Supported values (gap-existence, extension, etc.) for BLOSUM45. */

#define BLOSUM50_VALUES_MAX 16 /**< Number of different combinations supported for BLOSUM50. */
static array_of_8 blosum50_values[BLOSUM50_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0, 0.609639, 5.388310, 5.388310 },
	{ 13, 3, (double)INT2_MAX, 0.212, 0.063, 0.19, 1.1, -16, 0.639287, 18.113800, 18.202800 },
	{ 12, 3, (double)INT2_MAX, 0.206, 0.055, 0.17, 1.2, -18, 0.644715, 22.654600, 22.777700 },
	{ 11, 3, (double)INT2_MAX, 0.197, 0.042, 0.14, 1.4, -25, 0.656327, 29.861100, 30.045700 },
	{ 10, 3, (double)INT2_MAX, 0.186, 0.031, 0.11, 1.7, -34, 0.671150, 42.393800, 42.674000 },
	{ 9, 3, (double)INT2_MAX, 0.172, 0.022, 0.082, 2.1, -48, 0.694326, 66.069600, 66.516400 },
	{ 16, 2, (double)INT2_MAX, 0.215, 0.066, 0.20, 1.05, -15, 0.633899, 17.951800, 18.092100 },
	{ 15, 2, (double)INT2_MAX, 0.210, 0.058, 0.17, 1.2, -20, 0.641985, 21.940100, 22.141800 },
	{ 14, 2, (double)INT2_MAX, 0.202, 0.045, 0.14, 1.4, -27, 0.650682, 28.681200, 28.961900 },
	{ 13, 2, (double)INT2_MAX, 0.193, 0.035, 0.12, 1.6, -32, 0.660984, 42.059500, 42.471600 },
	{ 12, 2, (double)INT2_MAX, 0.181, 0.025, 0.095, 1.9, -41, 0.678090, 63.747600, 64.397300 },
	{ 19, 1, (double)INT2_MAX, 0.212, 0.057, 0.18, 1.2, -21, 0.635714, 26.311200, 26.923300 },
	{ 18, 1, (double)INT2_MAX, 0.207, 0.050, 0.15, 1.4, -28, 0.643523, 34.903700, 35.734800 },
	{ 17, 1, (double)INT2_MAX, 0.198, 0.037, 0.12, 1.6, -33, 0.654504, 48.895800, 50.148600 },
	{ 16, 1, (double)INT2_MAX, 0.186, 0.025, 0.10, 1.9, -42, 0.667750, 76.469100, 78.443000 },
	{ 15, 1, (double)INT2_MAX, 0.171, 0.015, 0.063, 2.7, -76, 0.694575, 140.053000, 144.160000 },
};  /**< Supported values (gap-existence, extension, etc.) for BLOSUM50. */

#define BLOSUM62_VALUES_MAX 12 /**< Number of different combinations supported for BLOSUM62. */
static array_of_8 blosum62_values[BLOSUM62_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660 },
	{ 11, 2, (double)INT2_MAX, 0.297, 0.082, 0.27, 1.1, -10, 0.641766, 12.673800, 12.757600 },
	{ 10, 2, (double)INT2_MAX, 0.291, 0.075, 0.23, 1.3, -15, 0.649362, 16.474000, 16.602600 },
	{ 9, 2, (double)INT2_MAX, 0.279, 0.058, 0.19, 1.5, -19, 0.659245, 22.751900, 22.950000 },
	{ 8, 2, (double)INT2_MAX, 0.264, 0.045, 0.15, 1.8, -26, 0.672692, 35.483800, 35.821300 },
	{ 7, 2, (double)INT2_MAX, 0.239, 0.027, 0.10, 2.5, -46, 0.702056, 61.238300, 61.886000 },
	{ 6, 2, (double)INT2_MAX, 0.201, 0.012, 0.061, 3.3, -58, 0.740802, 140.417000, 141.882000 },
	{ 13, 1, (double)INT2_MAX, 0.292, 0.071, 0.23, 1.2, -11, 0.647715, 19.506300, 19.893100 },
	{ 12, 1, (double)INT2_MAX, 0.283, 0.059, 0.19, 1.5, -19, 0.656391, 27.856200, 28.469900 },
	{ 11, 1, (double)INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200 },
	{ 10, 1, (double)INT2_MAX, 0.243, 0.024, 0.10, 2.5, -44, 0.693267, 83.178700, 85.065600 },
	{ 9, 1, (double)INT2_MAX, 0.206, 0.010, 0.052, 4.0, -87, 0.731887, 210.333000, 214.842000 },
}; /**< Supported values (gap-existence, extension, etc.) for BLOSUM62. */

#define BLOSUM80_VALUES_MAX 10 /**< Number of different combinations supported for BLOSUM80. */
static array_of_8 blosum80_values[BLOSUM80_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6, 0.564057, 1.918130, 1.918130 },
	{ 25, 2, (double)INT2_MAX, 0.342, 0.17, 0.66, 0.52, -1.6, 0.563956, 1.731000, 1.731300 },
	{ 13, 2, (double)INT2_MAX, 0.336, 0.15, 0.57, 0.59, -3, 0.570979, 2.673470, 2.692300 },
	{ 9, 2, (double)INT2_MAX, 0.319, 0.11, 0.42, 0.76, -6, 0.587837, 5.576090, 5.667860 },
	{ 8, 2, (double)INT2_MAX, 0.308, 0.090, 0.35, 0.89, -9, 0.597556, 7.536950, 7.686230 },
	{ 7, 2, (double)INT2_MAX, 0.293, 0.070, 0.27, 1.1, -14, 0.615254, 11.586600, 11.840400 },
	{ 6, 2, (double)INT2_MAX, 0.268, 0.045, 0.19, 1.4, -19, 0.644054, 19.958100, 20.441200 },
	{ 11, 1, (double)INT2_MAX, 0.314, 0.095, 0.35, 0.90, -9, 0.590702, 8.808610, 9.223320 },
	{ 10, 1, (double)INT2_MAX, 0.299, 0.071, 0.27, 1.1, -14, 0.609620, 13.833800, 14.533400 },
	{ 9, 1, (double)INT2_MAX, 0.279, 0.048, 0.20, 1.4, -19, 0.623800, 24.252000, 25.490400 },
}; /**< Supported values (gap-existence, extension, etc.) for BLOSUM80. */

#define BLOSUM90_VALUES_MAX 8 /**< Number of different combinations supported for BLOSUM90. */
static array_of_8 blosum90_values[BLOSUM90_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4 , 0.544178, 1.377760, 1.377760 },
	{ 9, 2, (double)INT2_MAX, 0.310, 0.12, 0.46, 0.67, -6 , 0.570267, 4.232290, 4.334170 },
	{ 8, 2, (double)INT2_MAX, 0.300, 0.099, 0.39, 0.76, -7, 0.581580, 5.797020, 5.961420 },
	{ 7, 2, (double)INT2_MAX, 0.283, 0.072, 0.30, 0.93, -11, 0.600024, 9.040880, 9.321600 },
	{ 6, 2, (double)INT2_MAX, 0.259, 0.048, 0.22, 1.2, -16, 0.629344, 16.024400, 16.531600 },
	{ 11, 1, (double)INT2_MAX, 0.302, 0.093, 0.39, 0.78, -8, 0.576919, 7.143250, 7.619190 },
	{ 10, 1, (double)INT2_MAX, 0.290, 0.075, 0.28, 1.04, -15, 0.591366, 11.483900, 12.269800 },
	{ 9, 1, (double)INT2_MAX, 0.265, 0.044, 0.20, 1.3, -19, 0.613013, 21.408300, 22.840900 },
};  /**< Supported values (gap-existence, extension, etc.) for BLOSUM90. */

#define PAM250_VALUES_MAX 16 /**< Number of different combinations supported for PAM250. */
static array_of_8 pam250_values[PAM250_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0, 0.660059, 11.754300, 11.754300 },
	{ 15, 3, (double)INT2_MAX, 0.205, 0.049, 0.13, 1.6, -23, 0.687656, 34.578400, 34.928000 },
	{ 14, 3, (double)INT2_MAX, 0.200, 0.043, 0.12, 1.7, -26, 0.689768, 43.353000, 43.443800 },
	{ 13, 3, (double)INT2_MAX, 0.194, 0.036, 0.10, 1.9, -31, 0.697431, 50.948500, 51.081700 },
	{ 12, 3, (double)INT2_MAX, 0.186, 0.029, 0.085, 2.2, -41, 0.704565, 69.606500, 69.793600 },
	{ 11, 3, (double)INT2_MAX, 0.174, 0.020, 0.070, 2.5, -48, 0.722438, 98.653500, 98.927100 },
	{ 17, 2, (double)INT2_MAX, 0.204, 0.047, 0.12, 1.7, -28, 0.684799, 41.583800, 41.735800 },
	{ 16, 2, (double)INT2_MAX, 0.198, 0.038, 0.11, 1.8, -29, 0.691098, 51.635200, 51.843900 },
	{ 15, 2, (double)INT2_MAX, 0.191, 0.031, 0.087, 2.2, -44, 0.699051, 67.256700, 67.558500 },
	{ 14, 2, (double)INT2_MAX, 0.182, 0.024, 0.073, 2.5, -53, 0.714103, 96.315100, 96.756800 },
	{ 13, 2, (double)INT2_MAX, 0.171, 0.017, 0.059, 2.9, -64, 0.728738, 135.653000, 136.339000 },
	{ 21, 1, (double)INT2_MAX, 0.205, 0.045, 0.11, 1.8, -34, 0.683265, 48.728200, 49.218800 },
	{ 20, 1, (double)INT2_MAX, 0.199, 0.037, 0.10, 1.9, -35, 0.689380, 60.832000, 61.514100 },
	{ 19, 1, (double)INT2_MAX, 0.192, 0.029, 0.083, 2.3, -52, 0.696344, 84.019700, 84.985600 },
	{ 18, 1, (double)INT2_MAX, 0.183, 0.021, 0.070, 2.6, -60, 0.710525, 113.829000, 115.184000 },
	{ 17, 1, (double)INT2_MAX, 0.171, 0.014, 0.052, 3.3, -86, 0.727000, 175.071000, 177.196000 },
}; /**< Supported values (gap-existence, extension, etc.) for PAM250. */

#define PAM30_VALUES_MAX 7 /**< Number of different combinations supported for PAM30. */
static array_of_8 pam30_values[PAM30_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818 },
	{ 7, 2, (double)INT2_MAX, 0.305, 0.15, 0.87, 0.35, -3, 0.479087, 1.014010, 1.162730 },
	{ 6, 2, (double)INT2_MAX, 0.287, 0.11, 0.68, 0.42, -4, 0.499980, 1.688060, 1.951430 },
	{ 5, 2, (double)INT2_MAX, 0.264, 0.079, 0.45, 0.59, -7, 0.533009, 3.377010, 3.871950 },
	{ 10, 1, (double)INT2_MAX, 0.309, 0.15, 0.88, 0.35, -3, 0.474741, 1.372050, 1.788770 },
	{ 9, 1, (double)INT2_MAX, 0.294, 0.11, 0.61, 0.48, -6, 0.492716, 2.463920, 3.186150 },
	{ 8, 1, (double)INT2_MAX, 0.270, 0.072, 0.40, 0.68, -10, 0.521286, 5.368130, 6.763480 },
}; /**< Supported values (gap-existence, extension, etc.) for PAM30. */


#define PAM70_VALUES_MAX 7 /**< Number of different combinations supported for PAM70. */
static array_of_8 pam70_values[PAM70_VALUES_MAX] = {
	{ (double)INT2_MAX, (double)INT2_MAX, (double)INT2_MAX, 0.3345, 0.229, 1.029, 0.3250,   -0.7, 0.511296, 0.633439, 0.633439 },
	{ 8, 2, (double)INT2_MAX, 0.301, 0.12, 0.54, 0.56, -5, 0.549019, 2.881650, 3.025710 },
	{ 7, 2, (double)INT2_MAX, 0.286, 0.093, 0.43, 0.67, -7, 0.565659, 4.534540, 4.785780 },
	{ 6, 2, (double)INT2_MAX, 0.264, 0.064, 0.29, 0.90, -12, 0.596330, 7.942630, 8.402720 },
	{ 11, 1, (double)INT2_MAX, 0.305, 0.12, 0.52, 0.59, -6, 0.543514, 3.681400, 4.108020 },
	{ 10, 1, (double)INT2_MAX, 0.291, 0.091, 0.41, 0.71, -9, 0.560723, 6.002970, 6.716570 },
	{ 9, 1, (double)INT2_MAX, 0.270, 0.060, 0.28, 0.97, -14, 0.585186, 11.360800, 12.636700 },
}; /**< Supported values (gap-existence, extension, etc.) for PAM70. */

struct Matrix_info
{
	const char *name;
	const array_of_8 *constants;
	const char *scores;
	const unsigned count;
	const int default_gap_open, default_gap_extend;

	static const Matrix_info& get(const string &name)
	{
		for (unsigned i = 0; i < sizeof(matrices) / sizeof(matrices[0]); ++i)
			if (name == matrices[i].name)
				return matrices[i];
		throw std::runtime_error("Invalid scoring matrix: " + name);
		return matrices[0];
	}

	const double* get_constants(int gap_open, int gap_extend) const
	{
		for (unsigned i = 0; i < count; ++i)
			if (constants[i][0] == gap_open && constants[i][1] == gap_extend)
				return constants[i];
		throw std::runtime_error("Invalid gap open and/or gap extend scores.");
		return 0;
	}

	static const Matrix_info matrices[8];
};

const Matrix_info Matrix_info::matrices[] = {
	{ "BLOSUM45", blosum45_values, (const char*)NCBISM_Blosum45.scores, BLOSUM45_VALUES_MAX, 14, 2 },
	{ "BLOSUM50", blosum50_values, (const char*)NCBISM_Blosum50.scores, BLOSUM50_VALUES_MAX, 13, 2 },
	{ "BLOSUM62", blosum62_values, (const char*)NCBISM_Blosum62.scores, BLOSUM62_VALUES_MAX, 11, 1 },
	{ "BLOSUM80", blosum80_values, (const char*)NCBISM_Blosum80.scores, BLOSUM80_VALUES_MAX, 10, 1 },
	{ "BLOSUM90", blosum90_values, (const char*)NCBISM_Blosum90.scores, BLOSUM90_VALUES_MAX, 10, 1 },
	{ "PAM70", pam70_values, (const char*)NCBISM_Pam70.scores, PAM70_VALUES_MAX, 10, 1 },
	{ "PAM250", pam250_values, (const char*)NCBISM_Pam250.scores, PAM250_VALUES_MAX, 14, 2 },
	{ "PAM30", pam30_values, (const char*)NCBISM_Pam30.scores, PAM30_VALUES_MAX, 9, 1 }
};

Score_matrix::Score_matrix(const string & matrix, int gap_open, int gap_extend, int reward, int penalty):
	gap_open_ (gap_open == -1 ? Matrix_info::get(matrix).default_gap_open : gap_open),
	gap_extend_ (gap_extend == -1 ? Matrix_info::get(matrix).default_gap_extend : gap_extend),
	constants_ (Matrix_info::get(matrix).get_constants(gap_open_, gap_extend_)),
	name_(matrix),
	matrix8_(Matrix_info::get(matrix).scores),
	bias_((char)(-low_score())),
	matrix8u_(Matrix_info::get(matrix).scores, bias_),
	matrix16_(Matrix_info::get(matrix).scores)
{ }

char Score_matrix::low_score() const
{
	char low = std::numeric_limits<char>::max();
	for (Letter i = 0; i < (char)value_traits.alphabet_size; ++i)
		for (Letter j = i + 1; j < (char)value_traits.alphabet_size; ++j)
			low = std::min(low, (char)this->operator()(i, j));
	return low;
}

double Score_matrix::avg_id_score() const
{
	double s = 0;
	for (int i = 0; i < 20; ++i)
		s += this->operator()(i, i);
	return s / 20;
}

std::ostream& operator<<(std::ostream& s, const Score_matrix &m)
{
	s << "(Matrix=" << m.name_
		<< " Lambda=" << m.lambda()
		<< " K=" << m.k()
		<< " Penalties=" << m.gap_open_
		<< '/' << m.gap_extend_ << ')';
	return s;
}

const char* custom_scores(const string &matrix_file)
{
	static char scores[25 * 25];
	string l, s;
	std::stringstream ss;
	vector<Letter> pos;
	unsigned n = 0;
	memset(scores, 0xff, sizeof(scores));
	if (matrix_file == "")
		return scores;
	std::ifstream f(matrix_file.c_str());
	while (!f.eof()) {
		std::getline(f, l);
		if (l[0] == '#')
			continue;
		if (pos.size() == 0) {
			for (string::const_iterator i = l.begin(); i != l.end(); ++i)
				if (*i == ' ' || *i == '\t')
					continue;
				else
					pos.push_back(value_traits.from_char(*i));
		}
		else {
			if (n >= pos.size())
				break;
			ss << l;
			if (value_traits.from_char(ss.get()) != pos[n])
				throw std::runtime_error("Invalid custom scoring matrix file format.");
			for (unsigned i = 0; i < pos.size(); ++i) {
				int score;
				ss >> score;
				scores[(int)pos[n] * 25 + (int)pos[i]] = score;
			}
			ss.clear();
			++n;
		}
	}
	return scores;
}

Score_matrix::Score_matrix(const string &matrix_file, double lambda, double K, int gap_open, int gap_extend):
	gap_open_(gap_open),
	gap_extend_(gap_extend),
	name_("custom"),
	matrix8_(custom_scores(matrix_file)),
	bias_((char)(-low_score())),
	matrix8u_(custom_scores(matrix_file), bias_),
	matrix16_(custom_scores(matrix_file))
{
	static double constants[5];
	constants[3] = lambda;
	constants[4] = K;
	constants_ = constants;
}