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
#ifndef MINMAX_TARGET_MAPPINGS_H
#define MINMAX_TARGET_MAPPINGS_H

#include <string>
#include <unordered_map>
#include <vector>

class TargetMapping
{
  public:
	TargetMapping(size_t m, size_t g) : mirna_(m), gene_(g) {}

	size_t mirna() const { return mirna_; }

	size_t gene() const { return gene_; }

	bool operator==(const TargetMapping& t)
	{
		return mirna_ == t.mirna_ && gene_ == t.gene_;
	}

  private:
	size_t mirna_;
	size_t gene_;
};

class TargetMappings
{
  public:
	using iterator = std::vector<TargetMapping>::iterator;
	using const_iterator = std::vector<TargetMapping>::const_iterator;

	void add(const std::string& mirna, const std::string& gene);

	iterator begin() { return mappings_.begin(); }
	iterator end() { return mappings_.end(); }

	const_iterator begin() const { return mappings_.begin(); }
	const_iterator end() const { return mappings_.end(); }

	const_iterator cbegin() const { return mappings_.begin(); }
	const_iterator cend() { return mappings_.end(); }

	size_t numGenes() const;
	size_t numMirnas() const;
	size_t numMappings() const;

	const std::string& mirna(size_t i) { return mirna_names_[i]; }
	const std::string& gene(size_t i) { return gene_names_[i]; }

	void finalize();

  private:
	std::vector<std::string> gene_names_;
	std::vector<std::string> mirna_names_;

	std::unordered_map<std::string, int> gene_to_id_;
	std::unordered_map<std::string, int> mirna_to_id_;

	std::vector<TargetMapping> mappings_;

	size_t findOrCreate_(const std::string& s,
	                     std::unordered_map<std::string, int>& m,
	                     std::vector<std::string>& v) const;
};

#endif // MINMAX_TARGET_MAPPINGS_H
