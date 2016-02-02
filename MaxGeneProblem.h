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
#ifndef MAXGENEPROBLEM_H
#define MAXGENEPROBLEM_H

#include "ILPProblem.h"

class MaxGeneProblem : public ILPProblem
{
  public:
	explicit MaxGeneProblem(const TargetMappings& mappings,
	                        std::size_t num_mirnas);
	explicit MaxGeneProblem(TargetMappings&& mappings, std::size_t num_mirnas);

	void setNumMirna(std::size_t k);

  protected:
	virtual void createObjectiveFunction_();
	virtual void createConstraints_();
	void createNumMirnaConstraint_();

  private:
	size_t num_mirnas_;
};

#endif // MAXGENEPROBLEM_H
