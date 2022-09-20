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

#ifndef COMPLEXSTRUCTURE
#define COMPLEXSTRUCTURE

#include <vector>
#include <string>

class ComplexStructure
{
public:
	enum streamType : int {VERTICES, FACES, ACTIVITY};

	virtual ~ComplexStructure() {}
	using vertex = double;
	using index = double;
	using simplex_base = std::vector<vertex>;

	virtual bool insertSimplex(simplex_base *numVertices, double timestamp) = 0;
	virtual index contractVertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers) = 0;
	virtual index getBoundary(simplex_base *simplex, std::vector<index> *boundary) = 0;
	virtual bool exists(simplex_base *simplex) = 0;
	virtual void print(std::string outputFileName) = 0;
	virtual double getSize() const = 0;
	virtual double getMaxSize() const = 0;
	virtual double getTowerWidth() const = 0;
	virtual double getFiltrationSize() const = 0;
	virtual vertex getVertexNumber(double num) const = 0;
};

#endif // COMPLEXSTRUCTURE

