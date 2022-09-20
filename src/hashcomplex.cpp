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

#include "hashcomplex.h"

HashComplex::HashComplex() : maxIndex(-1), maxSize(0), towerWidth(0), maxDim(0), outputString(NULL), outputStream(NULL), type(VERTICES)
{
	vertices = new std::unordered_map<double, vertex>();

	simplices = new std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
			KeyHasher,
			SimplicesEquals
		>();
}

HashComplex::HashComplex(std::string *&outputString, streamType type) : maxIndex(-1), maxSize(0), towerWidth(0), maxDim(0), outputString(outputString), outputStream(NULL), type(type)
{
	vertices = new std::unordered_map<double, vertex>();

	simplices = new std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
			KeyHasher,
			SimplicesEquals
		>();
}

HashComplex::HashComplex(std::ofstream *outputStream, streamType type) : maxIndex(-1), maxSize(0), towerWidth(0), maxDim(0), outputString(NULL), outputStream(outputStream), type(type)
{
	vertices = new std::unordered_map<double, vertex>();

	simplices = new std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
			KeyHasher,
			SimplicesEquals
		>();
}

HashComplex::~HashComplex()
{
	for (auto it = simplices->begin(); it != simplices->end(); ++it){
		delete it->first;
		delete it->second;
	}
	delete simplices;
	delete vertices;
}

bool HashComplex::insertSimplex(simplex_base *numVertices, double timestamp)
{
	simplex_base *vs = new simplex_base();

	if (numVertices->size() == 1){
		vertices->emplace(numVertices->at(0), numVertices->at(0));
		vs->push_back(numVertices->at(0));
	} else {
		for (simplex_base::size_type i = 0; i < numVertices->size(); i++){
			vs->push_back(vertices->at(numVertices->at(i)));
		}
		std::sort(vs->begin(), vs->end());
	}

	if (internInsertSimplex(vs, timestamp)){
		if (towerWidth < (double)simplices->size()) towerWidth = simplices->size();
		return true;
	}

	return false;
}

bool HashComplex::internInsertSimplex(simplex_base *vs, double timestamp)
{
	std::pair<simplex_base*,int> *p = new std::pair<simplex_base*,int>(vs, -1);
	Simplex *splx = new Simplex(maxIndex + 1, vs);

	if (simplices->emplace(p, splx).second == false) {
		p->first = NULL;
		delete splx;
		delete p;
		return false;
	}

	streamSimplex(splx, timestamp);
	maxIndex++;
	if ((int)vs->size() - 1 > maxDim) maxDim = vs->size() - 1;
	if (maxSize < (double)simplices->size()) maxSize = simplices->size();

	if (vs->size() <= 1) return true;

	for (simplex_base::size_type i = 0; i < vs->size(); i++){
		p->second = i;
		simplices->at(p)->addCofacet(splx, vs->at(i));
	}
	p->second = -1;

	return true;
}

ComplexStructure::index HashComplex::contractVertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index> *> *boundaries, std::vector<index> *&inactiveInsertionNumbers)
{
	vertex nu = vertices->at(u);
	vertex nv = vertices->at(v);

	vertices->erase(v);
	if (nu == nv) return -1;

	simplex_base vv(1, nv);
	simplex_base vu(1, nu);
	std::vector<simplex_base*> *acsActive = new std::vector<simplex_base*>();
	std::vector<simplex_base*> *acsInactive = new std::vector<simplex_base*>();

	vertex num = getSmallestActiveClosedStar(&vv, &vu, acsActive, acsInactive);
	std::vector<simplex_base*> *acsExt = new std::vector<simplex_base*>();

	// ------------------------- acsInactive and acsActive are deleted in preprocessActiveClosedStar
	// ------------------------- acsExt is deleted in contractVertexTo
	if (num == nu){
		vertices->at(u) = nv;
		preprocessActiveClosedStar(nv, acsInactive, acsActive, acsExt);
		return contractVertexTo(acsExt, &vu, timestamp, boundaries, inactiveInsertionNumbers);
	} else {
		preprocessActiveClosedStar(nu, acsInactive, acsActive, acsExt);
		return contractVertexTo(acsExt, &vv, timestamp, boundaries, inactiveInsertionNumbers);
	}
}

ComplexStructure::index HashComplex::contractVertexTo(std::vector<simplex_base*> *acs, simplex_base *vertexToRemoveAsSimplex, double timestamp,
								   std::vector<std::vector<index> *> *boundaries, std::vector<index> *&insertionNumbers)
{
	index first = -1;

	for (std::vector<simplex_base*>::iterator it = acs->begin(); it != acs->end(); it++){
		if (internInsertSimplex(*it, timestamp)) {
			if (boundaries != NULL){
				std::vector<index> *boundary = new std::vector<index>();
				if (first == -1) first = maxIndex;
				internGetBoundary(*it, boundary);
				boundaries->push_back(boundary);
			}
		}
	}

	deleteSimplex(vertexToRemoveAsSimplex, insertionNumbers);

	if (towerWidth < (double)simplices->size()) towerWidth = simplices->size();

	delete acs;
	return first;
}

ComplexStructure::index HashComplex::getBoundary(simplex_base *simplex, std::vector<index> *boundary)
{
	simplex_base convertedSimplex;

	for (simplex_base::size_type i = 0; i < simplex->size(); i++){
		convertedSimplex.push_back(vertices->at(simplex->at(i)));
	}
	std::sort(convertedSimplex.begin(), convertedSimplex.end());

	return internGetBoundary(&convertedSimplex, boundary);
}

ComplexStructure::index HashComplex::internGetBoundary(HashComplex::simplex_base *simplex, std::vector<index> *boundary)
{
	std::pair<simplex_base*,int> p(simplex, -1);

	if (simplex->size() == 1) return simplices->at(&p)->getInsertionNum();

	for (simplex_base::size_type i = 0; i < simplex->size(); i++){
		p.second = i;
		boundary->push_back(simplices->at(&p)->getInsertionNum());
	}
	std::sort(boundary->begin(), boundary->end());

	p.second = -1;
	return simplices->at(&p)->getInsertionNum();
}

bool HashComplex::isSaturated()
{
	return (simplices->size() == pow(2, vertices->size()) - 1);
}

bool HashComplex::exists(simplex_base *simplex)
{
	simplex_base vs;
	std::pair<simplex_base*,int> p(&vs, -1);
	for (simplex_base::size_type i = 0; i < simplex->size(); i++){
		vs.push_back(vertices->at(simplex->at(i)));
	}
	std::sort(vs.begin(), vs.end());

	return (simplices->find(&p) != simplices->end());
}

void HashComplex::print(std::string outputFileName)
{
	if (outputString == NULL) {
		std::cout << "The print string was not initialized (NULL).\n";
		return;
	}
	std::ofstream file(outputFileName);
	if (file.is_open()){
		file << (*outputString);
		file.close();
	} else std::cout << "Unable to open file " << outputFileName << "\n";
}

double HashComplex::getSize() const
{
	return simplices->size();
}

double HashComplex::getMaxSize() const
{
	return maxSize;
}

double HashComplex::getTowerWidth() const
{
	return towerWidth;
}

double HashComplex::getFiltrationSize() const
{
	return maxIndex + 1;
}

ComplexStructure::vertex HashComplex::getVertexNumber(double num) const
{
	if (vertices->find(num) == vertices->end()) return -1;
	return vertices->at(num);
}

double HashComplex::getNumberOfVertices() const
{
	return vertices->size();
}

void HashComplex::deleteSimplex(simplex_base *simplex, std::vector<index> *&insertionNumbers)
{
	std::pair<simplex_base*,int> p(simplex, -1);
	Simplex *splx = simplices->at(&p);
	std::unordered_map<vertex, Simplex*> *cofacets = splx->getCofacets();

	if (simplex->size() > 1){
		for (int i = 0; i < (int)simplex->size(); i++){
			p.second = i;
			simplices->at(&p)->getCofacets()->erase(simplex->at(i));
		}
	}
	p.second = -1;

	while (!cofacets->empty()) deleteSimplex(cofacets->begin()->second->getVertices(), insertionNumbers);

	std::pair<simplex_base*,int> *tmp = simplices->find(&p)->first;
	if (insertionNumbers != NULL) insertionNumbers->push_back(splx->getInsertionNum());
	streamInactivity(splx->getInsertionNum());

	simplices->erase(&p);
	delete tmp;
	delete splx;
}

int HashComplex::getVertexIndex(simplex_base *simplex, vertex v)
{
	int i = 0;
	while (i < (int)simplex->size() && simplex->at(i) != v){
		i++;
	}
	if (i == (int)simplex->size()) return -1;
	return i;
}

ComplexStructure::simplex_base *HashComplex::getExtendedSimplex(simplex_base *simplex, vertex v)
{
	simplex_base *eSimp = new simplex_base();
	int i = 0;
	while (i < (int)simplex->size() && simplex->at(i) < v){
		eSimp->push_back(simplex->at(i));
		i++;
	}
	if ((i < (int)simplex->size() && simplex->at(i) != v) || i == (int)simplex->size()){
		eSimp->push_back(v);
	}
	while (i < (int)simplex->size()){
		eSimp->push_back(simplex->at(i));
		i++;
	}
	return eSimp;
}

ComplexStructure::vertex HashComplex::getSmallestActiveClosedStar(simplex_base *v, simplex_base *u,
													 std::vector<simplex_base*> *acsActive, std::vector<simplex_base*> *acsInactive)
{
	std::queue<simplex_base*> qv;
	std::queue<simplex_base*> qu;
	simplex_base *s;
	std::pair<simplex_base*,int> p;

	if (getSmallestActiveStar(v, u, &qv, &qu) == v->at(0)){
		while (!qv.empty()){
			s = qv.front();
			qv.pop();
			p.first = s;
			p.second = getVertexIndex(s, v->at(0));
			if (s->size() > 1){
				acsActive->push_back(simplices->at(&p)->getVertices());
			}
			p.second = -1;
			acsInactive->push_back(s);
		}

		return v->at(0);
	} else {
		while (!qu.empty()){
			s = qu.front();
			qu.pop();
			p.first = s;
			p.second = getVertexIndex(s, u->at(0));
			if (s->size() > 1){
				acsActive->push_back(simplices->at(&p)->getVertices());
			}
			p.second = -1;
			acsInactive->push_back(s);
		}

		return u->at(0);
	}
}

ComplexStructure::vertex HashComplex::getSmallestActiveStar(simplex_base *v, simplex_base *u, std::queue<simplex_base*> *qv, std::queue<simplex_base*> *qu)
{
	auto hash = [](simplex_base* const& k) {
		std::size_t seed = k->size();
		for (simplex_base::size_type i = 0; i < k->size(); i++) {
			seed ^= (std::size_t)(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	};
	auto comp = [](const simplex_base* s1, const simplex_base* s2) {
		simplex_base::size_type size = s1->size();
		if (size != s2->size()) return false;
		for (simplex_base::size_type i = 0; i < size; i++){
			if (s1->at(i) != s2->at(i)) return false;
		}
		return true;
	};

	std::unordered_map<simplex_base*, bool, decltype(hash), decltype(comp)> visitedV(100, hash, comp);
	std::unordered_map<simplex_base*, bool, decltype(hash), decltype(comp)> visitedU(100, hash, comp);
	std::pair<simplex_base*,int> p;
	std::unordered_map<vertex, Simplex*> *cofacetsV;
	std::unordered_map<vertex, Simplex*> *cofacetsU;
	std::unordered_map<vertex, Simplex*>::iterator vit;
	std::unordered_map<vertex, Simplex*>::iterator uit;
	std::queue<simplex_base*> tv;
	std::queue<simplex_base*> tu;

	qv->push(v);
	qu->push(u);

	p.first = v;
	p.second = -1;
	cofacetsV = simplices->at(&p)->getCofacets();
	vit = cofacetsV->begin();

	p.first = u;
	cofacetsU = simplices->at(&p)->getCofacets();
	uit = cofacetsU->begin();

	if (vit != cofacetsV->end() && vit->first == u->at(0)) vit++;
	if (uit != cofacetsU->end() && uit->first == v->at(0)) uit++;

	while (vit != cofacetsV->end() && uit != cofacetsU->end()) {
		p.first = vit->second->getVertices();
		visitedV.emplace(p.first, true);
		qv->push(p.first);
		tv.push(p.first);

		p.first = uit->second->getVertices();
		visitedU.emplace(p.first, true);
		qu->push(p.first);
		tu.push(p.first);

		vit++;
		while ((vit == cofacetsV->end() && !tv.empty()) || (vit != cofacetsV->end() && (vit->first == u->at(0) || visitedV.find(vit->second->getVertices()) != visitedV.end()))){
			if (vit == cofacetsV->end() && !tv.empty()){
				p.first = tv.front();
				tv.pop();
				cofacetsV = simplices->at(&p)->getCofacets();
				vit = cofacetsV->begin();
			}
			if (vit != cofacetsV->end() && (vit->first == u->at(0) || visitedV.find(vit->second->getVertices()) != visitedV.end())) vit++;
		}

		uit++;
		while ((uit == cofacetsU->end() && !tu.empty()) || (uit != cofacetsU->end() && (uit->first == v->at(0) || visitedU.find(uit->second->getVertices()) != visitedU.end()))){
			if (uit == cofacetsU->end() && !tu.empty()){
				p.first = tu.front();
				tu.pop();
				cofacetsU = simplices->at(&p)->getCofacets();
				uit = cofacetsU->begin();
			}
			if (uit != cofacetsU->end() && (uit->first == v->at(0) || visitedU.find(uit->second->getVertices()) != visitedU.end())) uit++;
		}
	}

	if (uit == cofacetsU->end()) return u->at(0);
	else return v->at(0);
}

//#define PRUNNING	/* deletes collapse-pairs in active closed star (with greedy alg.), because they are not needed in the filtration */

void HashComplex::preprocessActiveClosedStar(vertex v, std::vector<simplex_base*> *&acsInactive,
											 std::vector<simplex_base*> *&acsActive, std::vector<simplex_base*> *&acs)
{
#ifdef PRUNNING
	int dim = maxDim + 1;
	std::unordered_map<std::pair<simplex_base*, int>*, simplex_base*, KeyHasher, SimplicesEquals> buckets[dim + 1];
	std::unordered_map<std::pair<simplex_base*, int>*, int, KeyHasher, SimplicesEquals> cofacets;
	std::pair<simplex_base*, int> p(NULL, -1);
	std::pair<simplex_base*, int> p2(NULL, -1);
	std::pair<simplex_base*, int> *tmp = NULL;
	std::unordered_map<std::pair<simplex_base*, int>*, simplex_base*, KeyHasher, SimplicesEquals>::iterator it;
	std::unordered_map<std::pair<simplex_base*, int>*, simplex_base*, KeyHasher, SimplicesEquals>::iterator pre;
	std::unordered_map<std::pair<simplex_base*, int>*, simplex_base*, KeyHasher, SimplicesEquals>::iterator pos;
	std::unordered_map<std::pair<simplex_base*, int>*, int, KeyHasher, SimplicesEquals>::iterator find;
	bool found = false;

	for (std::vector<simplex_base*>::size_type j = 0; j < acsInactive->size(); j++){
		p.first = getExtendedSimplex(acsInactive->at(j), v);
		if (simplices->find(&p) == simplices->end()){
			tmp = new std::pair<simplex_base*, int>(p.first, -1);
			buckets[p.first->size() - 1].emplace(tmp, p.first);
			cofacets.emplace(tmp, 0);
			for (int i = 0; i < (int)tmp->first->size(); i++){
				tmp->second = i;
				find = cofacets.find(tmp);
				if (find != cofacets.end()) find->second++;
			}
		} else {
			delete p.first;
		}
	}
	delete acsInactive;

	while (dim > 0){
		it = buckets[dim].begin();
		pre  = buckets[dim].begin();
		while (it != buckets[dim].end()){
			found = false;
			p.first = it->second;
			p.second = 0;
			while (p.second < (int)p.first->size() && !found){
				pos = buckets[dim - 1].find(&p);
				if (pos != buckets[dim - 1].end()){
					p2.first = pos->second;
					if (cofacets.at(&p2) == 1) found = true;
				}
				p.second++;
			}
			if (found){
				if ((int)p2.first->size() > 1){
					for (int i = 0; i < (int)p2.first->size(); i++) {
						p2.second = i;
						find = cofacets.find(&p2);
						if (find != cofacets.end()) find->second--;
					}
				}
				p2.second = -1;
				if ((int)p.first->size() > 1){
					for (int i = 0; i < (int)p.first->size(); i++) {
						p.second = i;
						find = cofacets.find(&p);
						if (find != cofacets.end()) find->second--;
					}
				}
				p.second = -1;
				cofacets.erase(&p);
				cofacets.erase(&p2);

				delete pos->first;
				delete pos->second;
				delete it->first;
				delete it->second;

				buckets[dim - 1].erase(pos);
				if (it != buckets[dim].begin()){
					buckets[dim].erase(it);
					it = pre;
					it++;
				} else {
					buckets[dim].erase(it);
					it = buckets[dim].begin();
					pre = buckets[dim].begin();
				}
			} else {
				if (it != buckets[dim].begin()) pre++;
				it++;
			}
		}
		dim--;
	}

	p.second = -1;
	while (acsActive->size() > 0){
		p.first = getExtendedSimplex(acsActive->back(), v);
		if (simplices->find(&p) == simplices->end()){
			buckets[p.first->size() - 1].emplace(new std::pair<simplex_base*, int>(p.first, -1), p.first);
		} else {
			delete p.first;
		}
		acsActive->pop_back();
	}
	delete acsActive;

	for (int i = 0; i < maxDim + 2; i++){
		for (it = buckets[i].begin(); it != buckets[i].end(); it++){
			delete it->first;
			acs->push_back(it->second);
		}
	}
#else
	std::vector<simplex_base*>::iterator itI = acsInactive->begin();
	std::vector<simplex_base*>::iterator itA = acsActive->begin();

	while (itI != acsInactive->end() && itA != acsActive->end()){
		while (itI != acsInactive->end() && (*itI)->size() <= (*itA)->size()){
			acs->push_back(getExtendedSimplex(*itI, v));
			itI++;
		}
		while (itI != acsInactive->end() && itA != acsActive->end() && (*itI)->size() > (*itA)->size()){
			acs->push_back(getExtendedSimplex(*itA, v));
			itA++;
		}
	}
	while (itI != acsInactive->end()){
		acs->push_back(getExtendedSimplex(*itI, v));
		itI++;
	}
	while (itA != acsActive->end()){
		acs->push_back(getExtendedSimplex(*itA, v));
		itA++;
	}

	delete acsInactive;
	delete acsActive;
#endif
}

void HashComplex::streamSimplex(Simplex *vs, double timestamp)
{
	int size = (int)vs->getVertices()->size();

	if (outputString != NULL){
		std::stringstream convert;

		if (type == ACTIVITY){
			convert.str(std::string());
			convert.clear();
			convert << std::setprecision(std::numeric_limits<double>::digits10 + 1) << timestamp;
			*outputString += convert.str();
			*outputString += " i ";
			if (size > 1){
				std::vector<index> boundary(size);
				std::pair<simplex_base*,int> p(vs->getVertices(), -1);
				for (int i = 0; i < size; i++){
					p.second = i;
					boundary.at(i) = simplices->at(&p)->getInsertionNum();
				}
				p.second = -1;
				std::sort(boundary.begin(), boundary.end());
				for (int i = 0; i < size; i++){
					convert.str(std::string());
					convert.clear();
					convert << boundary.at(i);
					*outputString += convert.str();
					*outputString += " ";
				}
			}
			convert.str(std::string());
			convert.clear();
			convert << vs->getInsertionNum();
			*outputString += convert.str();
			*outputString += "\n";
		} else {
			convert << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (size - 1);
			*outputString += convert.str();
			*outputString += " ";

			if (type == FACES){
				if (size > 1){
					std::vector<index> boundary(size);
					std::pair<simplex_base*,int> p(vs->getVertices(), -1);
					for (int i = 0; i < size; i++){
						p.second = i;
						boundary.at(i) = simplices->at(&p)->getInsertionNum();
					}
					p.second = -1;
					std::sort(boundary.begin(), boundary.end());
					for (int i = 0; i < size; i++){
						convert.str(std::string());
						convert.clear();
						convert << boundary.at(i);
						*outputString += convert.str();
						*outputString += " ";
					}
				}
			} else {
				for (int i = 0; i < size; i++){
					convert.str(std::string());
					convert.clear();
					convert << vs->getVertices()->at(i);
					*outputString += convert.str();
					*outputString += " ";

				}
			}

			convert.str(std::string());
			convert.clear();
			convert << timestamp;
			*outputString += convert.str();
			*outputString += "\n";
		}
	}
	if (outputStream != NULL){
		if (type == ACTIVITY){
			*outputStream << std::setprecision(std::numeric_limits<double>::digits10 + 1) << timestamp << " i ";
			if (size > 1){
				std::vector<index> boundary(size);
				std::pair<simplex_base*,int> p(vs->getVertices(), -1);
				for (int i = 0; i < size; i++){
					p.second = i;
					boundary.at(i) = simplices->at(&p)->getInsertionNum();
				}
				p.second = -1;
				std::sort(boundary.begin(), boundary.end());
				for (int i = 0; i < size; i++){
					*outputStream << boundary.at(i) << " ";
				}
			}
			*outputStream << vs->getInsertionNum() << "\n";
		} else {
			*outputStream << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (size - 1) << " ";
			if (type == FACES){
				if (size > 1){
					std::vector<index> boundary(size);
					std::pair<simplex_base*,int> p(vs->getVertices(), -1);
					for (int i = 0; i < size; i++){
						p.second = i;
						boundary.at(i) = simplices->at(&p)->getInsertionNum();
					}
					p.second = -1;
					std::sort(boundary.begin(), boundary.end());
					for (int i = 0; i < size; i++){
						*outputStream << boundary.at(i) << " ";
					}
				}
			} else {
				for (int i = 0; i < size; i++){
					*outputStream << vs->getVertices()->at(i) << " ";
				}
			}
			*outputStream << timestamp << "\n";
		}
	}
}

void HashComplex::streamInactivity(index insertionNumber)
{
	if (type != ACTIVITY) return;

	if (outputString != NULL){
		std::stringstream convert;
		convert << insertionNumber;
		*outputString += "a ";
		*outputString += convert.str();
		*outputString += "\n";
	}
	if (outputStream != NULL){
		*outputStream << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "a " << insertionNumber << "\n";
	}
}

HashComplex::Simplex::Simplex(index num, simplex_base* vertices)
{
	insertionNum = num;
	cofacets = new std::unordered_map<vertex, Simplex*>();
	this->vertices = vertices;
}

HashComplex::Simplex::~Simplex()
{
	delete cofacets;
	delete vertices;
}

ComplexStructure::index HashComplex::Simplex::getInsertionNum() const
{
	return insertionNum;
}

void HashComplex::Simplex::setInsertionNum(index value)
{
	insertionNum = value;
}

void HashComplex::Simplex::addCofacet(Simplex *coface, vertex v)
{
	cofacets->emplace(v, coface);
}

std::unordered_map<HashComplex::vertex, HashComplex::Simplex *> *HashComplex::Simplex::getCofacets()
{
	return cofacets;
}

HashComplex::simplex_base *HashComplex::Simplex::getVertices() const
{
	return vertices;
}


