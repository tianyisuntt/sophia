/*
 * Copyright 2016 TU Graz
 * Contributed by: Hannah Schreiber
 *
 * This file is part of Sophia.
 *
 * Sophia is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Sophia is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Sophia. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

//#define PHAT
//#define GUDHI

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <ios>
#include <iomanip>
#include <ctime>

#ifdef GUDHI
#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#endif

#include "hashcomplex.h"
#include "persistence.h"

#define LOG(msg) if (!silent) std::cout << msg << std::endl;

class TowerReader
{
public:
	TowerReader(bool silent);
	~TowerReader();

	enum type : int {INCLUSION, CONTRACTION, INACTIVE, COMMENT};

	void convert(std::string inputFileName);
	void convertWithVerticesOutput(std::string inputFileName, std::string *&outputStream);
	void convertWithFacesOutput(std::string inputFileName, std::string *&outputStream);
	void convertWithActivityOutput(std::string inputFileName, std::string *&outputStream);
	void convertWithVerticesOutput(std::string inputFileName, std::ofstream *outputStream);
	void convertWithFacesOutput(std::string inputFileName, std::ofstream *outputStream);
	void convertWithActivityOutput(std::string inputFileName, std::ofstream *outputStream);
	void convertAndComputePersistence(std::string inputFileName, std::string outputFileName, double chunkSize);
	void computePersistence(std::string inputFileName, std::string outputFileName, double chunkSize);
#ifdef PHAT
	void computePersistenceWithPhat(std::string inputFileName, std::string outputFileName);
#endif
#ifdef GUDHI
	void computePersistenceWithGudhi(std::string inputFileName, std::string outputFileName);
#endif

private:
	ComplexStructure *complex;
	Persistence *pers;
	bool silent;

	TowerReader::type readOperation(std::string *line, std::vector<double> *vertices, double *timestamp);
	bool readOperationsToConvert(std::string inputFileName);
	void printComplexData();
};

#endif // OPERATIONS_H
