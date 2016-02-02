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
#include "ILPProblem.h"

#include "CPLEXException.h"

void ILPProblem::handleCPLEXError_(int status)
{
	if(status) {
		char buffer[CPXMESSAGEBUFSIZE];
		CPXCCHARptr errstr;
		errstr = CPXgeterrorstring(env_, status, buffer);

		if(errstr == NULL) {
			sprintf(buffer, "CPLEX Error %5d:  Unknown error code.n", status);
		}

		throw CPLEXException(buffer);
	}
}

ILPProblem::ILPProblem(const TargetMappings& mappings) : mappings_(mappings) {}

ILPProblem::ILPProblem(TargetMappings&& mappings)
    : mappings_(std::move(mappings))
{
}

size_t ILPProblem::numVariables() const { return CPXgetnumcols(env_, lp_); }

size_t ILPProblem::numConstraints() const { return CPXgetnumrows(env_, lp_); }

size_t ILPProblem::numNonZero() const { return CPXgetnumnz(env_, lp_); }

void ILPProblem::createProblem_()
{
	int status = 0;
	env_ = CPXopenCPLEX(&status);
	handleCPLEXError_(status);

	lp_ = CPXcreateprob(env_, &status, "MinMax");
	handleCPLEXError_(status);

	mappings_.finalize();
	createObjectiveFunction_();
	createConstraints_();
}

void ILPProblem::createMappingConstraints_()
{
	const size_t num_constr = mappings_.numGenes();
	const size_t num_indices = mappings_.numMappings() + mappings_.numGenes();
	std::vector<int> indices(num_indices);
	std::vector<double> row(num_indices);

	std::vector<int> rmatbeg(num_constr + 1, 0);
	std::vector<double> rhs(num_constr, 0.0);
	std::vector<char> sense(num_constr, 'G');

	size_t cur_gene = 0;
	size_t cur_constr = 0;
	size_t cur_index = 1;

	indices[0] = mappings_.numMirnas();
	row[0] = -1.0;
	for(const auto& mapping : mappings_) {
		if(mapping.gene() != cur_gene) {
			cur_gene = mapping.gene();
			indices[cur_index] = cur_gene + mappings_.numMirnas();
			row[cur_index] = -1.0;
			rmatbeg[++cur_constr] = cur_index++;
		}

		indices[cur_index] = mapping.mirna();
		row[cur_index] = 1.0;
		++cur_index;
	}

	rmatbeg[++cur_constr] = cur_index;

#ifndef NDEBUG
	std::cout << "indices: ";
	std::copy(indices.begin(), indices.end(),
	          std::ostream_iterator<int>(std::cout, ", "));
	std::cout << '\n';
	std::cout << "weights: ";
	std::copy(row.begin(), row.end(),
	          std::ostream_iterator<double>(std::cout, ", "));
	std::cout << '\n';
	std::cout << "constr ptr: ";
	std::copy(rmatbeg.begin(), rmatbeg.end(),
	          std::ostream_iterator<double>(std::cout, ", "));
	std::cout << '\n';
#endif

	int status = CPXaddrows(env_, lp_, 0, num_constr, num_indices, &rhs[0],
	                        &sense[0], &rmatbeg[0], &indices[0], &row[0], 0, 0);

	handleCPLEXError_(status);
}

ILPProblem::Result ILPProblem::solve()
{
	int status = CPXmipopt(env_, lp_);
	handleCPLEXError_(status);

	const size_t num_variables = mappings_.numMirnas() + mappings_.numGenes();
	std::vector<double> row(num_variables);
	status = CPXgetx(env_, lp_, &row[0], 0, num_variables - 1);
	handleCPLEXError_(status);

	std::vector<std::string> mirnas;
	std::vector<std::string> genes;

	for(size_t i = 0; i < mappings_.numMirnas(); ++i) {
		if(row[i] > 0.5) {
			mirnas.push_back(mappings_.mirna(i));
		}
	}

	for(size_t i = 0; i < mappings_.numGenes(); ++i) {
		if(row[i + mappings_.numMirnas()] > 0.5) {
			genes.push_back(mappings_.gene(i));
		}
	}

	return {std::move(mirnas), std::move(genes)};
}
