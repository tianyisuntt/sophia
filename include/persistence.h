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

#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include <vector>
#include <list>
#include <unordered_map>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "hashcomplex.h"
#include "listcolumn.h"
#include "heapcolumn.h"

#define LIST

class Persistence
{
public:
	Persistence(double reductionInterval, std::string persistencePairsFileName);
	~Persistence();

	class BoundaryMatrix
	{
	public:
		BoundaryMatrix(std::string persistencePairsFileName);
		~BoundaryMatrix();

		void insertColumn(double insertionNumber, std::vector<double> *boundary, double timestamp);
		void reduce(double start);
		void clearOut();
		void markInactive(std::vector<double> *insertionNumbers);
		void markInactive(double insertionNumber);
		void insertVertex(double insertionNumber, double timestamp);

		double getLastInsertNumber() const;
		int getMaxDim() const;

	private:
#ifdef LIST
		typedef std::list<double> ContainerType;
		typedef ListColumn ColumnType;
#else
		typedef std::vector<double> ContainerType;
		typedef HeapColumn ColumnType;
#endif

		std::unordered_map<double, ColumnType*> *columns;
		std::unordered_map<double, double> *latest;
		std::unordered_map<double, std::pair<bool, bool>*> *isActivePositive;
		std::unordered_map<double, double> *timestamps;
		double lastInsertNumber;
		int maxDim;
		std::ofstream *persistencePairsFile;

		void clearColumn(double columnIndex);
		void printPersistencePair(int dim, double birth, double death);
	};

	bool conv_insertSimplex(std::vector<double> *simplex, double timestamp);
	void conv_contractVertices(double u, double v, double timestamp);
	void insertSimplex(double insertionNum, std::vector<double> *boundary, double timestamp);
	void finalizeReduction();
	void markInactive(double insertionNumber);

	void printComplexData();

private:
	ComplexStructure *complex;
	BoundaryMatrix *matrix;
	double reductionInterval;
	double lastReduction;

	void computePartialPersistence();
};

#endif // PERSISTENCE_H
