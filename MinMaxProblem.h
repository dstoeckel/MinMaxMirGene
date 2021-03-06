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
#ifndef MINMAXPROBLEM_H
#define MINMAXPROBLEM_H

#include "ILPProblem.h"

#include <string>
#include <vector>

class MinMaxProblem : public ILPProblem
{
  public:
	MinMaxProblem(const TargetMappings& mappings, double mirna_weight,
	              double gene_weight);
	MinMaxProblem(TargetMappings&& mappings, double mirna_weight,
	              double gene_weight);

  private:
	virtual void createObjectiveFunction_();
	virtual void createConstraints_();

	double mirna_weight_;
	double gene_weight_;
};

#endif // MINMAXPROBLEM_H
