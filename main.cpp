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
#include "CPLEXException.h"
#include "MaxGeneProblem.h"
#include "MinMaxProblem.h"
#include "TargetMappings.h"

#include <cstring>
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
		exit(-8);
	}

	for(const auto& item : items) {
		output << item << '\n';
	}
}

const char* commandList()
{
	return "\tminmax\n\tmaxgene\n\tmaxgene-curve\n\tminmirna\n\tminmirna-curve";
}

void printILPStatistics(const ILPProblem& problem)
{
	std::cout << "Created ILP formulation with " << problem.numVariables()
	          << " variables, " << problem.numConstraints()
	          << " constraints, and " << problem.numNonZero()
	          << " non-zero entries.\n";
}

void solveProblem(ILPProblem& problem, const std::string& mpath,
                  const std::string& gpath)
{
	printILPStatistics(problem);

	std::vector<std::string> mirnas;
	std::vector<std::string> genes;
	std::tie(mirnas, genes) = problem.solve();

	std::cout << "Solution contains " << mirnas.size() << " miRNAs and "
	          << genes.size() << " genes.\n";

	writeList(mpath, mirnas);
	writeList(gpath, genes);
}

int minMax(int argc, char* argv[], TargetMappings& mappings)
{
	if(argc <= 6) {
		std::cerr << "Not enough arguments supplied. Usage:\n\t" << argv[0]
		          << " minmax mappings.txt mirna_weight gene_weight mirnas.out "
		             "genes.out\n";
		return -3;
	}

	double mirnaWeight = 0.0;
	double geneWeight = 0.0;

	try {
		mirnaWeight = std::stof(argv[3]);
		geneWeight = std::stof(argv[4]);
	} catch(const std::exception& e) {
		std::cerr << "Error converting argument to a number\n";
		return -4;
	}

	MinMaxProblem problem(mappings, mirnaWeight, geneWeight);
	solveProblem(problem, argv[5], argv[6]);

	return 0;
}

int maxGene(int argc, char* argv[], TargetMappings& mappings)
{
	if(argc <= 5) {
		std::cerr << "Not enough arguments supplied. Usage:\n\t" << argv[0]
		          << " minmax mappings.txt num_mirna mirnas.out genes.out\n";
		return -3;
	}

	int num_mirnas = 0;
	try {
		num_mirnas = std::stoi(argv[3]);
	} catch(const std::exception& e) {
		std::cerr << "Could not convert " << argv[3] << " to a number\n";
		return -6;
	}

	MaxGeneProblem problem(mappings, num_mirnas);
	solveProblem(problem, argv[4], argv[5]);

	return 0;
}

int maxGeneCurve(int argc, char* argv[], TargetMappings& mappings)
{
	std::ofstream curve(argv[3]);

	if(!curve) {
		std::cerr << "Could not open file '" << argv[3] << "' for writing!\n";
		return -9;
	}

	MaxGeneProblem problem(mappings, 0);

	std::size_t i = 1;
	for(; i < mappings.numMirnas(); ++i) {
		problem.setNumMirna(i);
		auto result = problem.solve();

		if(result.second.size() == mappings.numGenes()) {
			break;
		}
		curve << i << '\t' << result.second.size() << '\n';
	}

	for(; i < mappings.numMirnas(); ++i) {
		curve << i << '\t' << mappings.numGenes() << '\n';
	}

	return 0;
}

int minMirna(int argc, char* argv[], TargetMappings& mappings)
{
	std::cerr << "Not yet implemented!\n";
	return 0;
}

int minMirnaCurve(int argc, char* argv[], TargetMappings& mappings)
{
	std::cerr << "Not yet implemented!\n";
	return 0;
}

int dispatchCLIArguments(int argc, char* argv[], TargetMappings& mappings)
{
	if(strcmp(argv[1], "minmax") == 0) {
		return minMax(argc, argv, mappings);
	} else if(strcmp(argv[1], "maxgene") == 0) {
		return maxGene(argc, argv, mappings);
	} else if(strcmp(argv[1], "maxgene-curve") == 0) {
		return maxGeneCurve(argc, argv, mappings);
	} else if(strcmp(argv[1], "minmirna") == 0) {
		return minMirna(argc, argv, mappings);
	} else if(strcmp(argv[1], "minmirna-curve") == 0) {
		return minMirnaCurve(argc, argv, mappings);
	} else {
		std::cerr << "Unknown command '" << argv[1]
		          << "'.\n\nAvailable commands:\n" << commandList() << '\n';
		return -2;
	}
}

void printMappingStatistics(const TargetMappings& mappings)
{
	std::cout << "Read " << mappings.numMappings() << " mappings between "
	          << mappings.numMirnas() << " miRNAs and " << mappings.numGenes()
	          << " genes.\n";
}

int main(int argc, char* argv[])
{
	if(argc <= 2) {
		std::cerr << "Not enough arguments supplied. Usage:\n\n\t" << argv[0]
		          << " [command] mappings.txt [...]\n\nAvailable commands:\n"
		          << commandList() << '\n';
		return -1;
	}

	auto mappings = readMappings(argv[2]);
	printMappingStatistics(mappings);

	try {
		return dispatchCLIArguments(argc, argv, mappings);
	} catch(const CPLEXException& e) {
		std::cerr << "Error during ILP computation: " << e.what() << std::endl;
		return -7;
	}

	return 0;
}
