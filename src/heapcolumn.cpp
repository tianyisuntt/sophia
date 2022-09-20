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

#include "heapcolumn.h"

HeapColumn::HeapColumn(std::vector<double> *column, int dim) : dim(dim), column(column), insertsSinceLastPrune(0)
{
	std::make_heap(this->column->begin(), this->column->end());
	temp_column = new std::vector<double>();
}

HeapColumn::~HeapColumn()
{
	delete column;
	delete temp_column;
}

void HeapColumn::add(HeapColumn *columnToAdd)
{
	double size = columnToAdd->getSize();
	for (double i = 0; i < size; i++) {
		column->push_back(columnToAdd->at(i));
		std::push_heap(column->begin(), column->end());
	}
	insertsSinceLastPrune += size;

	if (2 * insertsSinceLastPrune > (double)column->size()) prune();
}

int HeapColumn::getDim() const
{
	return dim;
}

double HeapColumn::getPivot()
{
	double pivot = popPivot();
	if (pivot != -1){
		column->push_back(pivot);
		std::push_heap(column->begin(), column->end());
	}
	return pivot;
}

double HeapColumn::getSize()
{
	prune();
	return column->size();
}

double HeapColumn::at(double index)
{
	return column->at(index);
}

void HeapColumn::set(std::vector<double> *newColumn)
{
	column->clear();
	column->insert(column->end(), newColumn->begin(), newColumn->end());
	std::make_heap(column->begin(), column->end());
}

void HeapColumn::insert(double value)
{
	column->push_back(value);
	std::push_heap(column->begin(), column->end());
}

void HeapColumn::prune()
{
	if (insertsSinceLastPrune == 0) return;

	std::vector<double> *tempCol = temp_column;
	tempCol->clear();
	double pivot = popPivot();
	while (pivot != -1) {
		tempCol->push_back(pivot);
		pivot = popPivot();
	}
	temp_column = column;
	column = tempCol;
	std::reverse(column->begin(), column->end());
	std::make_heap(column->begin(), column->end());

	insertsSinceLastPrune = 0;
}

double HeapColumn::popPivot()
{
	if (column->empty()) {
		return -1;
	} else {
		double pivot = column->front();
		std::pop_heap(column->begin(), column->end());
		column->pop_back();
		while (!column->empty() && column->front() == pivot) {
			std::pop_heap(column->begin(), column->end());
			column->pop_back();
			if (column->empty()) {
				return -1;
			} else {
				pivot = column->front();
				std::pop_heap(column->begin(), column->end());
				column->pop_back();
			}
		}
		return pivot;
	}
}

void HeapColumn::print()
{
	std::cout << "dim: " << dim << "\n";
	std::cout << "column: ";
	for (std::vector<double>::size_type i = 0; i < column->size(); i++) std::cout << column->at(i) << " ";
	std::cout << "\n";
}
