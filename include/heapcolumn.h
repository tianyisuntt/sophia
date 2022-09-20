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

#ifndef HEAPCOLUMN_H
#define HEAPCOLUMN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

class HeapColumn
{
public:
	HeapColumn(std::vector<double> *column, int dim);
	~HeapColumn();

public:
	void add(HeapColumn *columnToAdd);
	int getDim() const;
	double getPivot();
	double popPivot();
	double getSize();
	double at(double index);
	void set(std::vector<double> *newColumn);
	void insert(double value);

	void print();

private:
	int dim;
	std::vector<double> *column;
	std::vector<double> *temp_column;
	double insertsSinceLastPrune;

	void prune();
};

#endif // HEAPCOLUMN_H
