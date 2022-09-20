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

#ifndef LISTCOLUMN_H
#define LISTCOLUMN_H

#include <iostream>
#include <list>

class ListColumn
{
public:
	ListColumn(std::list<double> *column, int dim);
	~ListColumn();

public:
	void add(ListColumn *columnToAdd);
	int getDim() const;
	std::list<double>::iterator getBeginIterator();
	std::list<double>::reverse_iterator getReverseBeginIterator();
	std::list<double>::iterator getEndIterator();
	std::list<double>::reverse_iterator getReverseEndIterator();
	void erase(std::list<double>::iterator *pos);
	double getPivot();
	double getSize();

private:
	int dim;
	std::list<double> *column;
};

#endif // LISTCOLUMN_H
