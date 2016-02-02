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
#ifndef MINMAX_CPLEXEXCEPTION_H
#define MINMAX_CPLEXEXCEPTION_H

#include <stdexcept>

class CPLEXException : public std::exception
{
  public:
	CPLEXException(const char* msg) throw() : what_(msg) {}
	~CPLEXException() throw(){};
	const char* what() const throw() { return what_.c_str(); }

  private:
	std::string what_;
};

#endif // MINMAX_CPLEXEXCEPTION_H
