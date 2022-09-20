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

#include "listcolumn.h"

ListColumn::ListColumn(std::list<double> *column, int dim) : dim(dim), column(column)
{}

ListColumn::~ListColumn()
{
	delete column;
}

void ListColumn::add(ListColumn *columnToAdd)
{
	std::list<double>::iterator itToAdd = columnToAdd->getBeginIterator(), itTarget = column->begin();
	while (itToAdd != columnToAdd->getEndIterator() && itTarget != column->end()){
		if (*itToAdd == *itTarget){
			column->erase(itTarget++);
			itToAdd++;
		} else if (*itToAdd < *itTarget){
			column->insert(itTarget, *itToAdd);
			itToAdd++;
		} else {
			itTarget++;
		}
	}
	while (itToAdd != columnToAdd->getEndIterator()){
		column->push_back(*itToAdd);
		itToAdd++;
	}
}

int ListColumn::getDim() const
{
	return dim;
}

std::list<double>::iterator ListColumn::getBeginIterator()
{
	return column->begin();
}

std::list<double>::reverse_iterator ListColumn::getReverseBeginIterator()
{
	return column->rbegin();
}

std::list<double>::iterator ListColumn::getEndIterator()
{
	return column->end();
}

std::list<double>::reverse_iterator ListColumn::getReverseEndIterator()
{
	return column->rend();
}

void ListColumn::erase(std::list<double>::iterator *pos)
{
	column->erase(*pos);
}

double ListColumn::getPivot()
{
	if (column->empty()) return -1;
	return column->back();
}

double ListColumn::getSize()
{
	return column->size();
}

