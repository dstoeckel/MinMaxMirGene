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

	CPXchgobjsen(env_, lp_, CPX_MAX);
	std::vector<int> indices(nvar);
	std::vector<char> ctype(nvar, 'B');
	std::vector<double> row(nvar, gene_weight_);
	std::vector<double> lb(nvar, 0.0);
	std::vector<double> ub(nvar, 1.0);

	std::iota(indices.begin(), indices.end(), 0);
	std::fill_n(row.begin(), mappings_.numMirnas(), -mirna_weight_);

	int status = CPXnewcols(env_, lp_, nvar, &row[0], &lb[0], &ub[0], 0, 0);
	handleCPLEXError_(status);

	status = CPXchgctype(env_, lp_, indices.size(), &indices[0], &ctype[0]);
	handleCPLEXError_(status);
}

void MinMaxProblem::createConstraints_() { createMappingConstraints_(); }
