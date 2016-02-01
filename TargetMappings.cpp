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
#include "TargetMappings.h"

#include <algorithm>

void TargetMappings::add(const std::string& mirna, const std::string& gene)
{
	size_t m = findOrCreate_(mirna, mirna_to_id_, mirna_names_);
	size_t g = findOrCreate_(gene, gene_to_id_, gene_names_);

	mappings_.emplace_back(m, g);
}

size_t TargetMappings::numMirnas() const { return mirna_names_.size(); }

size_t TargetMappings::numGenes() const { return gene_names_.size(); }

size_t TargetMappings::numMappings() const { return mappings_.size(); }

void TargetMappings::finalize()
{
	auto comp = [](const TargetMapping& a, const TargetMapping& b) {
		if(a.gene() == b.gene()) {
			return a.mirna() < b.mirna();
		}

		return a.gene() < b.gene();
	};

	std::sort(mappings_.begin(), mappings_.end(), comp);
	auto new_end = std::unique(mappings_.begin(), mappings_.end());

	mappings_.erase(new_end, mappings_.end());
}

size_t TargetMappings::findOrCreate_(const std::string& s,
                                     std::unordered_map<std::string, int>& m,
                                     std::vector<std::string>& v) const
{
	auto res = m.find(s);

	if(res == m.end()) {
		res = m.emplace(s, v.size()).first;
		v.push_back(s);
	}

	return res->second;
}
