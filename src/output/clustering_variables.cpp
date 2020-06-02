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

#include "clustering_variables.h"

std::map<std::string, Variable*> VariableRegistry::regMap;
QueryLength VariableRegistry::queryLength;
SubjectLength VariableRegistry::subjectLength;
QueryStart VariableRegistry::queryStart;
QueryEnd VariableRegistry::queryEnd;
SubjectStart VariableRegistry::subjectStart;
SubjectEnd VariableRegistry::subjectEnd;
EValue VariableRegistry::eValue;
BitScore VariableRegistry::bitScore;
Score VariableRegistry::score;
Length VariableRegistry::length;
PercentIdenticalMatches VariableRegistry::percentIdenticalMatches;
NumberIdenticalMatches VariableRegistry::numberIdenticalMatches;
NumberMismatches VariableRegistry::numberMismatches;
NumberPositiveMatches VariableRegistry::numberPositiveMatches;
NumberGapOpenings VariableRegistry::numberGapOpenings;
NumberGaps VariableRegistry::numberGaps;
PercentagePositiveMatches VariableRegistry::percentagePositiveMatches;
QueryFrame VariableRegistry::queryFrame;
QueryCoveragePerHsp VariableRegistry::queryCoveragePerHsp;
SwDiff VariableRegistry::swdiff;
Time VariableRegistry::time;
SubjectCoveragePerHsp VariableRegistry::subjectCoveragePerHsp;
UngappedScore VariableRegistry::ungappedScore;
VariableRegistry::StaticConstructor VariableRegistry::_staticConstructor;

