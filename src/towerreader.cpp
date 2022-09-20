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

#include "towerreader.h"

#ifdef PHAT
#include <phat/compute_persistence_pairs.h>
/*#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_heap.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/full_pivot_column.h>*/
#include <phat/representations/bit_tree_pivot_column.h>
#include <phat/algorithms/twist_reduction.h>
/*#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>*/
#include <phat/helpers/dualize.h>
#endif

TowerReader::TowerReader(bool silent) : complex(NULL), pers(NULL), silent(silent)
{ }

TowerReader::~TowerReader()
{
	if (complex != NULL) delete complex;
	if (pers != NULL) delete pers;
}

void TowerReader::convert(std::string inputFileName)
{
	complex = new HashComplex();
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithVerticesOutput(std::string inputFileName, std::string *&outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::VERTICES);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithFacesOutput(std::string inputFileName, std::string *&outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::FACES);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithActivityOutput(std::string inputFileName, std::string *&outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::ACTIVITY);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithVerticesOutput(std::string inputFileName, std::ofstream *outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::VERTICES);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithFacesOutput(std::string inputFileName, std::ofstream *outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::FACES);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertWithActivityOutput(std::string inputFileName, std::ofstream *outputStream)
{
	complex = new HashComplex(outputStream, HashComplex::streamType::ACTIVITY);
	if (readOperationsToConvert(inputFileName)) printComplexData();
}

void TowerReader::convertAndComputePersistence(std::string inputFileName, std::string outputFileName, double chunkSize)
{
	std::string line;
	std::ifstream file(inputFileName);
	pers = new Persistence(chunkSize, outputFileName);

	if (file.is_open()){
		std::vector<double> vertices;
		double timestamp = -1;
		double defaultTimestamp = 0;
		while (getline(file,line)){
			TowerReader::type type = readOperation(&line, &vertices, &timestamp);
			if (timestamp != -1) defaultTimestamp = timestamp;

			if (type == INCLUSION){
				pers->conv_insertSimplex(&vertices, defaultTimestamp);
			} else if (type == CONTRACTION) {
				pers->conv_contractVertices(vertices.at(0), vertices.at(1), defaultTimestamp);
			}

			defaultTimestamp++;
			timestamp = -1;
		}
		file.close();
	} else {
		std::cout << "Unable to open file " << inputFileName << "\n";
		return;
	}
	pers->finalizeReduction();
	if (!silent) pers->printComplexData();
}

void TowerReader::computePersistence(std::string inputFileName, std::string outputFileName, double chunkSize)
{
	std::string line;
	std::ifstream file(inputFileName);
	pers = new Persistence(chunkSize, outputFileName);

	if (file.is_open()){
		std::vector<double> boundary;
		double timestamp = -1;
		double defaultTimestamp = 0;
		double insertionNum;
		while (getline(file,line)){
			TowerReader::type type = readOperation(&line, &boundary, &timestamp);
			if (timestamp != -1) defaultTimestamp = timestamp;

			if (type == INCLUSION){
				insertionNum = boundary.back();
				boundary.pop_back();
				pers->insertSimplex(insertionNum, &boundary, defaultTimestamp);
			} else if (type == INACTIVE) {
				pers->markInactive(boundary.at(0));
			}

			defaultTimestamp++;
			timestamp = -1;
		}
		file.close();
	} else {
		std::cout << "Unable to open file " << inputFileName << "\n";
		return;
	}
	pers->finalizeReduction();
}

#ifdef PHAT
void TowerReader::computePersistenceWithPhat(std::string inputFileName, std::string outputFileName)
{
	std::string *stream = new std::string("");
	std::vector<double> timestamps;
	convertWithFacesOutput(inputFileName, stream);

	if (complex == NULL) return;

	phat::boundary_matrix<phat::bit_tree_pivot_column> matrix;

	double read_timer = omp_get_wtime();

	LOG("Reading input stream...");

	std::string cur_line;
	std::stringstream input_stream(*stream);
	delete stream;

	matrix.set_num_cols((phat::index)complex->getFiltrationSize());
	phat::column temp_col;
	phat::index cur_col = -1;
	while (getline(input_stream, cur_line)){
		cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
		if (cur_line != "" && cur_line[ 0 ] != '#') {
			cur_col++;
			std::stringstream ss(cur_line);

			int64_t temp_dim;
			ss >> temp_dim;
			matrix.set_dim(cur_col, (phat::dimension) temp_dim);

			int64_t temp_index;
			temp_col.clear();
			while (ss.good()) {
				ss >> temp_index;
				temp_col.push_back((phat::index)temp_index);
			}
			timestamps.push_back(temp_col.back());
			temp_col.pop_back();
			std::sort(temp_col.begin(), temp_col.end());
			matrix.set_col(cur_col, temp_col);
		}
	}
	input_stream.str(std::string());

	double read_time = omp_get_wtime() - read_timer;
	double read_time_rounded = floor(read_time * 10.0 + 0.5) / 10.0;
	LOG("Reading input stream took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << read_time_rounded << " s");

	double pairs_timer = omp_get_wtime();
	phat::persistence_pairs pairs;
	LOG("Computing persistence pairs ...");
	phat::compute_persistence_pairs<phat::twist_reduction>(pairs, matrix);
	double pairs_time = omp_get_wtime() - pairs_timer;
	double pairs_time_rounded = floor(pairs_time * 10.0 + 0.5) / 10.0;
	LOG("Computing persistence pairs took " << std::setiosflags( std::ios::fixed ) << std::setiosflags( std::ios::showpoint ) << std::setprecision( 1 ) << pairs_time_rounded << " s");

	double write_timer = omp_get_wtime();

	LOG("Writing output file " << outputFileName << " in ascii mode ...");

	std::ofstream output_stream(outputFileName);
	if(output_stream.fail()) {
		std::cout << "Unable to open file " << outputFileName << "\n";
		return;
	}

	pairs.sort();
	for(phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ ) {
		phat::index birth = pairs.get_pair(idx).first;
		phat::index death = pairs.get_pair(idx).second;
		if (timestamps.at(birth) != timestamps.at(death)) output_stream << timestamps.at(birth) << " " << timestamps.at(death) << std::endl;
	}

	output_stream.close();

	double write_time = omp_get_wtime() - write_timer;
	double write_time_rounded = floor(write_time * 10.0 + 0.5) / 10.0;
	LOG("Writing output file took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << write_time_rounded << " s");
}
#endif

#ifdef GUDHI
void TowerReader::computePersistenceWithGudhi(std::string inputFileName, std::string outputFileName)
{
	int p = 11;
	double min_persistence = 0;

	std::string *inputStream = new std::string("");
	convertWithVerticesOutput(inputFileName, inputStream);

	std::clock_t start = std::clock();

	LOG("Simplex_tree from stream - output_file=" << outputFileName.c_str() << " - p=" << p << " - min_persistence=" << min_persistence);

	Gudhi::Simplex_tree<> simplex_tree;
	std::stringstream simplex_tree_stream(*inputStream);
	delete inputStream;
	simplex_tree_stream >> simplex_tree;

	LOG("The complex contains " << simplex_tree.num_simplices() << " simplices - dimension " << simplex_tree.dimension());

	simplex_tree.initialize_filtration();

	Gudhi::persistent_cohomology::Persistent_cohomology<Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp> pcoh(simplex_tree);
	pcoh.init_coefficients(p);

	pcoh.compute_persistent_cohomology(min_persistence);

	LOG("Computing time: " << (std::clock() - start) / (double)CLOCKS_PER_SEC << " s");

	if (outputFileName.empty()) {
		pcoh.output_diagram();
	} else {
		std::ofstream out(outputFileName);
		pcoh.output_diagram(out);
		out.close();
	}
}
#endif

TowerReader::type TowerReader::readOperation(std::string *line, std::vector<double> *vertices, double *timestamp)
{
	TowerReader::type type;
	vertices->clear();
	double num;

	size_t next = line->find_first_not_of(' ', 0);
	size_t current = next;
	next = line->find_first_of(' ', current);
	if (next == std::string::npos) return COMMENT;
	if (line->substr(current, next - current) == "i") type = INCLUSION;
	else if (line->substr(current, next - current) == "c") type = CONTRACTION;
	else if (line->substr(current, next - current) == "a") type = INACTIVE;
	else if (line->substr(current, next - current) == "#") return COMMENT;
	else {
		*timestamp = stod(line->substr(current, next - current));
		next = line->find_first_not_of(' ', next + 1);
		current = next;
		next = line->find_first_of(' ', current);
		if (next == std::string::npos) {
			std::cout << "Operation syntaxe error in file.\n";
			exit(0);
		}
		if (line->substr(current, next - current) == "i") type = INCLUSION;
		else if (line->substr(current, next - current) == "c") type = CONTRACTION;
		else if (line->substr(current, next - current) == "a") type = INACTIVE;
		else if (line->substr(current, next - current) == "#") return COMMENT;
		else {
			std::cout << "Operation syntaxe error in file.\n";
			exit(0);
		}
	}

	next = line->find_first_not_of(' ', next + 1);
	while (next != std::string::npos){
		current = next;
		next = line->find_first_of(' ', current);
		num = stod(line->substr(current, next - current));
		vertices->push_back(num);
		if (next != std::string::npos) next = line->find_first_not_of(' ', next + 1);
	}

	return type;
}

bool TowerReader::readOperationsToConvert(std::string inputFileName)
{
	std::string line;
	std::ifstream file(inputFileName);

	if (file.is_open()){
		std::vector<double> vertices;
		double timestamp = -1;
		double defaultTimestamp = 0;
		while (getline(file,line)){
			TowerReader::type type = readOperation(&line, &vertices, &timestamp);
			if (timestamp != -1) defaultTimestamp = timestamp;

			if (type == INCLUSION){
				if (complex->insertSimplex(&vertices, defaultTimestamp)) defaultTimestamp++;
			} else if (type == CONTRACTION) {
				std::vector<double> *dummy = NULL;
				complex->contractVertices(vertices.at(0), vertices.at(1), defaultTimestamp, NULL, dummy);
				defaultTimestamp++;
			}

			timestamp = -1;
		}
		file.close();
	} else {
		std::cout << "Unable to open file " << inputFileName << "\n";
		return false;
	}

	return true;
}

void TowerReader::printComplexData()
{
	if (silent) return;
	LOG(std::setprecision(std::numeric_limits<double>::digits10 + 1) << "Filtration Size: " << complex->getFiltrationSize());
	LOG("Max Size: " << complex->getMaxSize());
	LOG("Tower Width: " << complex->getTowerWidth());
}

