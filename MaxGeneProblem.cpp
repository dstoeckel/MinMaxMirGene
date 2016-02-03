/*
 * MinMaxMirnaGene - A program for computing optimal miRNA-gene covers.
 * Copyright (C) 2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *
 * MinMaxMirnaGene is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MinMaxMirnaGene is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MinMaxMirnaGene. If not, see <http://www.gnu.org/licenses/>.
 */
#include "MaxGeneProblem.h"

#include <array>
#include <algorithm>

MaxGeneProblem::MaxGeneProblem(const TargetMappings& mappings,
                               std::size_t num_mirnas)
    : ILPProblem(mappings), num_mirnas_(num_mirnas)
{
	createProblem_();
}

MaxGeneProblem::MaxGeneProblem(TargetMappings&& mappings,
                               std::size_t num_mirnas)
    : ILPProblem(std::move(mappings)), num_mirnas_(num_mirnas)
{
	createProblem_();
}

void MaxGeneProblem::setNumMirna(size_t k)
{
	const int idx = 0;
	const double value = k;
	handleCPLEXError_(CPXchgrhs(env_, lp_, 1, &idx, &value));
}

void MaxGeneProblem::createObjectiveFunction_()
{
	const std::size_t nvar = mappings_.numGenes() + mappings_.numMirnas();

	int status = CPXchgobjsen(env_, lp_, CPX_MAX);
	handleCPLEXError_(status);

	std::vector<char> ctype(nvar, 'B');
	std::vector<double> row(nvar, 1.0);

	std::fill_n(row.begin(), mappings_.numMirnas(), 0.0);

	status = CPXnewcols(env_, lp_, nvar, &row[0], nullptr, nullptr, &ctype[0],
	                    nullptr);
	handleCPLEXError_(status);
}

void MaxGeneProblem::createConstraints_()
{
	createNumMirnaConstraint_();
	createMappingConstraints_();
}

void MaxGeneProblem::createNumMirnaConstraint_()
{
	const std::size_t nvar = mappings_.numMirnas();

	const char sense = 'E';
	const double rhs = num_mirnas_;
	std::vector<int> indices(nvar);
	std::vector<double> row(nvar);
	std::array<int, 2> rmatbeg = {0, static_cast<int>(nvar)};

	std::iota(indices.begin(), indices.end(), 0);
	std::fill(row.begin(), row.end(), 1.0);

	int status = CPXaddrows(env_, lp_, 0, 1, nvar, &rhs, &sense, &rmatbeg[0],
	                        &indices[0], &row[0], nullptr, nullptr);
	handleCPLEXError_(status);
}
