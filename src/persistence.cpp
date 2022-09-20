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

#include "persistence.h"

Persistence::Persistence(double reductionInterval, std::string persistencePairsFileName) : reductionInterval(reductionInterval), lastReduction(-1)
{
	complex = new HashComplex();
	matrix = new BoundaryMatrix(persistencePairsFileName);
	if (this->reductionInterval < 1) this->reductionInterval = 1;
}

Persistence::~Persistence()
{
	delete complex;
	delete matrix;
}

bool Persistence::conv_insertSimplex(std::vector<double> *simplex, double timestamp)
{
	std::vector<double> boundary;
	double insertionNum;

	if (!complex->insertSimplex(simplex, timestamp)) return false;
	insertionNum = complex->getBoundary(simplex, &boundary);
	if (simplex->size() == 1) {
		matrix->insertVertex(insertionNum, timestamp);
		return true;
	}
	matrix->insertColumn(insertionNum, &boundary, timestamp);

	if (fmod(insertionNum, reductionInterval) == 0) {
		computePartialPersistence();
	}

	return true;
}

void Persistence::conv_contractVertices(double u, double v, double timestamp)
{
	std::vector<std::vector<double>*> boundaries;
	std::vector<double> *insertionNumbers = new std::vector<double>();
	bool reduce = false;

	double first = complex->contractVertices(u, v, timestamp, &boundaries, insertionNumbers);

	for (std::vector<std::vector<double>*>::size_type i = 0; i < boundaries.size(); i++){
		matrix->insertColumn(first + i, boundaries.at(i), timestamp);
		delete boundaries.at(i);
		if (fmod((first + i), reductionInterval) == 0) reduce = true;
	}
	matrix->markInactive(insertionNumbers);
	delete insertionNumbers;

	if (reduce) computePartialPersistence();
}

void Persistence::insertSimplex(double insertionNum, std::vector<double> *boundary, double timestamp)
{
	if (boundary->size() == 0) {
		matrix->insertVertex(insertionNum, timestamp);
		if (insertionNum != 0 && fmod(insertionNum, reductionInterval) == 0) computePartialPersistence();
		return;
	}
	matrix->insertColumn(insertionNum, boundary, timestamp);

	if (insertionNum != 0 && fmod(insertionNum, reductionInterval) == 0) computePartialPersistence();
}

void Persistence::finalizeReduction()
{
	if (lastReduction != matrix->getLastInsertNumber()) matrix->reduce(lastReduction + 1);
	lastReduction = matrix->getLastInsertNumber();
}

void Persistence::markInactive(double insertionNumber)
{
	matrix->markInactive(insertionNumber);
}

void Persistence::printComplexData()
{
	std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "Filtration Size: " << complex->getFiltrationSize() << std::endl;
	std::cout << "Max Size: " << complex->getMaxSize() << std::endl;
	std::cout << "Tower Width: " << complex->getTowerWidth() << std::endl;
}

void Persistence::computePartialPersistence()
{
	matrix->reduce(lastReduction + 1);
	matrix->clearOut();
	lastReduction = matrix->getLastInsertNumber();
}

Persistence::BoundaryMatrix::BoundaryMatrix(std::string persistencePairsFileName) : lastInsertNumber(-1), maxDim(-1)
{
	columns = new std::unordered_map<double, ColumnType*>();
	latest = new std::unordered_map<double, double>();
	isActivePositive = new std::unordered_map<double, std::pair<bool, bool>*>();
	timestamps = new std::unordered_map<double, double>();
	persistencePairsFile = new std::ofstream(persistencePairsFileName);
	if (!persistencePairsFile->is_open()){
		std::cout << "Persistence Pairs File could not be open\n";
		exit(0);
	}
}

Persistence::BoundaryMatrix::~BoundaryMatrix()
{
	for (auto it = columns->begin(); it != columns->end(); it++){
		delete it->second;
	}
	delete columns;
	delete latest;
	for (auto it = isActivePositive->begin(); it != isActivePositive->end(); it++){
		delete it->second;
	}
	delete isActivePositive;
	persistencePairsFile->close();
	delete persistencePairsFile;
	delete timestamps;
}

void Persistence::BoundaryMatrix::insertColumn(double insertionNumber, std::vector<double> *boundary, double timestamp)
{
	ContainerType *boundaryCells = new ContainerType();
	isActivePositive->emplace(insertionNumber, new std::pair<bool, bool>(true, true));

	for (int i = 0; i < (int)boundary->size(); i++){
		boundaryCells->push_back(boundary->at(i));
	}

	columns->emplace(insertionNumber, new ColumnType(boundaryCells, boundary->size() - 1));

	lastInsertNumber = insertionNumber;
	if (maxDim < (int)boundary->size() - 1) maxDim = boundary->size() - 1;
	timestamps->emplace(insertionNumber, timestamp);
}

void Persistence::BoundaryMatrix::reduce(double start)
{
	for (int d = maxDim; d > 0; d--){
		for (double i = start; i <= lastInsertNumber; i++){
			if (columns->find(i) != columns->end() && columns->at(i)->getDim() == d){
				ColumnType *curr = columns->at(i);
				double pivot = curr->getPivot();

				while (pivot != -1 && latest->find(pivot) != latest->end()){
					curr->add(columns->at(latest->at(pivot)));
					pivot = curr->getPivot();
				}

				if (pivot != -1){
					isActivePositive->at(i)->second = false;
					latest->emplace(pivot, i);
					clearColumn(pivot);
					printPersistencePair(d - 1, pivot, i);
				} else {
					clearColumn(i);
				}
			}
		}
	}
}

void Persistence::BoundaryMatrix::clearOut()
{
	double r;
#ifndef LIST
	ContainerType tmp;
#endif

	double c;
	for (auto it = columns->begin(); it != columns->end(); it++){
		c = it->first;
		ColumnType *column = columns->at(c);
		r = column->getPivot();
		if (isActivePositive->at(r)->first){
#ifdef LIST
			ContainerType::reverse_iterator it;
			ContainerType::iterator it2;
			it = column->getReverseBeginIterator();
			it++;
			while (it != column->getReverseEndIterator()){
				if (latest->find(*it) != latest->end() && !isActivePositive->at(*it)->first){
					column->add(columns->at(latest->at(*(it--))));
				} else if (!isActivePositive->at(*it)->second && !isActivePositive->at(*it)->first) {
					it2 = (++it).base();
					it--; it--;
					column->erase(&it2);
				}
				it++;
			}
#else
			tmp.clear();
			tmp.push_back(column->popPivot());
			double max = column->popPivot();
			while (max != -1){
				if (latest->find(max) != latest->end() && !isActivePositive->at(max)->first){
					column->insert(max);
					column->add(columns->at(latest->at(max)));
				} else if (isActivePositive->at(max)->second || isActivePositive->at(max)->first) {
					tmp.push_back(max);
				}
				max = column->popPivot();
			}
			std::reverse(tmp.begin(), tmp.end());
			column->set(&tmp);
#endif
		}
	}

	for (auto it = columns->begin(), next_it = columns->begin(); it != columns->end(); it = next_it)
	{
		next_it = it; ++next_it;
		c = it->first;
		r = columns->at(c)->getPivot();
		if (!isActivePositive->at(r)->first)
		{
			latest->erase(r);
			delete isActivePositive->at(r);
			isActivePositive->erase(r);
			clearColumn(c);
			if (!isActivePositive->at(c)->first) {
				delete isActivePositive->at(c);
				isActivePositive->erase(c);
			}
		}
	}

	for (auto it = isActivePositive->begin(), next_it = isActivePositive->begin(); it != isActivePositive->end(); it = next_it)
	{
		next_it = it; ++next_it;
		c = it->first;
		if (columns->find(c) == columns->end() && !it->second->second && !it->second->first)
		{
			delete isActivePositive->at(c);
			isActivePositive->erase(c);
		}
	}
}

void Persistence::BoundaryMatrix::markInactive(std::vector<double> *insertionNumbers)
{
	for (std::vector<double>::size_type i = 0; i < insertionNumbers->size(); i++){
		isActivePositive->at(insertionNumbers->at(i))->first = false;
	}
}

void Persistence::BoundaryMatrix::markInactive(double insertionNumber)
{
	isActivePositive->at(insertionNumber)->first = false;
}

void Persistence::BoundaryMatrix::insertVertex(double insertionNumber, double timestamp)
{
	isActivePositive->emplace(insertionNumber, new std::pair<bool, bool>(true, true));
	timestamps->emplace(insertionNumber, timestamp);
}

void Persistence::BoundaryMatrix::printPersistencePair(int dim, double birth, double death)
{
	if (timestamps->at(birth) != timestamps->at(death)) *persistencePairsFile << std::setprecision(std::numeric_limits<double>::digits10 + 1)
																			  << dim << " " << timestamps->at(birth) << " " << timestamps->at(death) << std::endl;
	timestamps->erase(birth);
	timestamps->erase(death);
}

void Persistence::BoundaryMatrix::clearColumn(double columnIndex)
{
	if (columns->find(columnIndex) == columns->end()) return;
	delete columns->at(columnIndex);
	columns->erase(columnIndex);
}

int Persistence::BoundaryMatrix::getMaxDim() const
{
	return maxDim;
}

double Persistence::BoundaryMatrix::getLastInsertNumber() const
{
	return lastInsertNumber;
}

