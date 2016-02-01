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
#include "TargetMappings.h"

#include <fstream>
#include <iostream>

TargetMappings readMappings(const std::string& path)
{
	TargetMappings mappings;

	std::ifstream input(path);

	if(!input) {
		std::cerr << "Could not open file '" << path << "' for reading.\n";
		exit(-1);
	}

	std::string mirna;
	std::string gene;
	input >> mirna >> gene;
	while(input) {
		mappings.add(mirna, gene);
		input >> mirna >> gene;
	}

	return mappings;
}

void writeList(const std::string& path, const std::vector<std::string>& items)
{
	std::ofstream output(path);

	if(!output) {
		std::cerr << "Could not open file '" << path << "' for writing.\n";
		exit(-3);
	}

	for(const auto& item : items) {
		output << item << '\n';
	}
}

int main(int argc, char* argv[])
{
	if(argc <= 5) {
		std::cerr << "Not enough arguments supplied. Usage:\n\t" << argv[0]
		          << " input.txt mirna_weight gene_weight mirnas.out genes.out";
		return -4;
	}

	double mirnaWeight = 0.0;
	double geneWeight = 0.0;

	try {
		mirnaWeight = std::stof(argv[2]);
		geneWeight = std::stof(argv[3]);
	} catch(const std::exception& e) {
		std::cerr << "Error converting argument to a number\n";
		return -6;
	}

	auto mappings = readMappings(argv[1]);
	std::cout << "Read " << mappings.numMappings() << " mappings between "
	          << mappings.numMirnas() << " miRNAs and " << mappings.numGenes()
	          << " genes.\n";

	try {
		MinMaxProblem problem(mappings, mirnaWeight, geneWeight);
		std::cout << "Created ILP formulation with " << problem.numVariables()
		          << " variables, " << problem.numConstraints()
		          << " constraints, and " << problem.numNonZero()
		          << " non-zero entries.\n";

		std::vector<std::string> mirnas;
		std::vector<std::string> genes;
		std::tie(mirnas, genes) = problem.solve();

		std::cout << "Solution contains " << mirnas.size() << " miRNAs and "
		          << genes.size() << " genes.\n";

		writeList(argv[4], mirnas);
		writeList(argv[5], genes);
	} catch(const std::exception& e) {
		std::cerr << "Error during ILP computation: " << e.what() << std::endl;
		return -7;
	}

	return 0;
}

