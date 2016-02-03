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
#include "MinMaxProblem.h"
#include "CPLEXException.h"

#include <algorithm>

MinMaxProblem::MinMaxProblem(const TargetMappings& mappings,
                             double mirna_weight, double gene_weight)
    : ILPProblem(mappings),
      mirna_weight_(mirna_weight),
      gene_weight_(gene_weight)
{
	createProblem_();
}

MinMaxProblem::MinMaxProblem(TargetMappings&& mappings, double mirna_weight,
                             double gene_weight)
    : ILPProblem(std::move(mappings)),
      mirna_weight_(mirna_weight),
      gene_weight_(gene_weight)
{
	createProblem_();
}

void MinMaxProblem::createObjectiveFunction_()
{
	const std::size_t nvar = mappings_.numGenes() + mappings_.numMirnas();

	int status = CPXchgobjsen(env_, lp_, CPX_MAX);
	handleCPLEXError_(status);

	std::vector<char> ctype(nvar, 'B');
	std::vector<double> row(nvar, gene_weight_);

	std::fill_n(row.begin(), mappings_.numMirnas(), -mirna_weight_);

	status = CPXnewcols(env_, lp_, nvar, &row[0], nullptr, nullptr, &ctype[0], nullptr);
	handleCPLEXError_(status);
}

void MinMaxProblem::createConstraints_() { createMappingConstraints_(); }
