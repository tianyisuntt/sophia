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

#ifndef COMPLEX_H
#define COMPLEX_H

#include <unordered_map>
#include <vector>
#include <list>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <functional>

#include "complexstructure.h"

class HashComplex : public ComplexStructure
{
public:
	HashComplex();
	HashComplex(std::string *&outputString, streamType type = VERTICES);
	HashComplex(std::ofstream *outputStream, streamType type = VERTICES);
	~HashComplex();

	class Simplex
	{
	public:
		Simplex(index num, simplex_base *vertices);
		~Simplex();

		index getInsertionNum() const;
		void setInsertionNum(index value);

		void addCofacet(Simplex *coface, vertex v);
		std::unordered_map<vertex, Simplex*>* getCofacets();

		simplex_base *getVertices() const;

	private:
		index insertionNum;
		std::unordered_map<vertex, Simplex*> *cofacets;
		simplex_base* vertices;
	};

	struct KeyHasher {
		std::size_t operator()(const std::pair<simplex_base*, int> *k) const
		{
			std::size_t seed;
			if (k->second < 0) seed = k->first->size();
			else seed = k->first->size() - 1;

			for (int i = 0; i < (int)k->first->size(); i++) {
				if (i != k->second) seed ^= (std::size_t)(k->first->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

	struct SimplicesEquals : std::binary_function<const std::pair<simplex_base*, int>*, const std::pair<simplex_base*, int>*, bool> {
		bool operator()(const std::pair<simplex_base*, int> *s1, const std::pair<simplex_base*, int> *s2) const
		{
			const std::pair<simplex_base*, int> *key;
			const std::pair<simplex_base*, int> *inMap;
			simplex_base::size_type size;

			if (s1->second > -1) {
				key = s1;
				inMap = s2;
			} else {
				key = s2;
				inMap = s1;
			}

			if (key->second < 0) size = key->first->size();
			else size = key->first->size() - 1;
			if (size != inMap->first->size()) return false;
			int j = 0;
			for (simplex_base::size_type i = 0; i < size; i++){
				if (j == key->second) j++;
				if (key->first->at(j) != inMap->first->at(i)) {
					return false;
				}
				j++;
			}
			return true;
		}
	};

	double getNumberOfVertices() const;
	bool insertSimplex(simplex_base *numVertices, double timestamp) override;
	index contractVertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers) override;
	index getBoundary(simplex_base *simplex, std::vector<index> *boundary) override;
	bool isSaturated();
	bool exists(simplex_base *simplex) override;
	double getSize() const override;
	double getMaxSize() const override;
	double getTowerWidth() const override;
	double getFiltrationSize() const override;
	vertex getVertexNumber(double num) const override;

	void print(std::string outputFileName) override;

private:
	std::unordered_map<double, vertex> *vertices;
	double maxIndex;
	double maxSize;
	double towerWidth;
	int maxDim;
	std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
			KeyHasher,
			SimplicesEquals
		> *simplices;
	std::string *outputString;
	std::ofstream *outputStream;
	streamType type;

	index contractVertexTo(std::vector<simplex_base*> *acs, simplex_base *vertexToRemoveAsSimplex, double timestamp,
						  std::vector<std::vector<index>*> *boundaries, std::vector<index> *&insertionNumbers);
	bool internInsertSimplex(simplex_base *vs, double timestamp);
	void deleteSimplex(simplex_base *simplex, std::vector<index> *&insertionNumbers);
	int getVertexIndex(simplex_base *simplex, vertex v);
	simplex_base* getExtendedSimplex(simplex_base *simplex, vertex v);
	vertex getSmallestActiveClosedStar(simplex_base *v, simplex_base *u,
							   std::vector<simplex_base *> *acsActive, std::vector<simplex_base *> *acsInactive);
	vertex getSmallestActiveStar(simplex_base *v, simplex_base *u, std::queue<simplex_base*> *qv, std::queue<simplex_base*> *qu);
	void preprocessActiveClosedStar(vertex v, std::vector<simplex_base *> *&acsInactive,
									std::vector<simplex_base*> *&acsActive, std::vector<simplex_base*> *&acs);
	index internGetBoundary(simplex_base *simplex, std::vector<index> *boundary);
	void streamSimplex(Simplex *vs, double timestamp);
	void streamInactivity(index insertionNumber);
};

#endif // COMPLEX_H
