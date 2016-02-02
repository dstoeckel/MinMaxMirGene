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
#ifndef ILPPROBLEM_H
#define ILPPROBLEM_H

#include "TargetMappings.h"

#include <ilcplex/cplex.h>

#include <vector>
#include <utility>

class ILPProblem
{
  public:
	using Result =
	    std::pair<std::vector<std::string>, std::vector<std::string>>;
	explicit ILPProblem(const TargetMappings& mappings);
	explicit ILPProblem(TargetMappings&& mappings);

	std::size_t numVariables() const;
	std::size_t numConstraints() const;
	std::size_t numNonZero() const;

	virtual Result solve();

  protected:
	TargetMappings mappings_;
	CPXENVptr env_;
	CPXLPptr lp_;

	virtual void createProblem_();
	virtual void createObjectiveFunction_() = 0;
	virtual void createConstraints_() = 0;

	void createMappingConstraints_();
	void handleCPLEXError_(int status);
};

#endif // ILPPROBLEM_H
