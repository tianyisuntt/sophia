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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "towerreader.h"

using namespace std;

void displayUsage()
{
	cout << "Usage:\n";
	cout << "  -cv input_file_name output_file_name \t\t\t #Builds filtration in VERTICES Format\n";
	cout << "  -cf input_file_name output_file_name \t\t\t #Builds filtration in FACES Format\n";
#ifdef PHAT
	cout << "  -cphat input_file_name output_file_name \t\t #Builds filtration and computes persistence pairs with PHAT library (default parameters)\n";
#endif
#ifdef GUDHI
	cout << "  -cgudhi input_file_name output_file_name \t\t #Builds filtration and computes persistence pairs with GUDHI library (default parameters)\n";
#endif
	cout << "  -cp input_file_name output_file_name chunk_size \t #Builds filtration and computes persistence pairs\n";
	cout << "  -p input_file_name output_file_name chunk_size \t #Computes persistence pairs\n";
	cout << "Adding the argument --silent avoid output in the terminal.\n";
	cout << "\nfor more details, see Readme file\n\n";
}

int main(int argc, char *argv[])
{
	if (argc < 2 || argc > 6){
		displayUsage();
		return 0;
	}

	std::vector<char*> args;
	bool silent = false;
	for (int i = 1; i < argc; i++){
		if (strcmp(argv[i], "--silent") == 0) silent = true;
		else args.push_back(argv[i]);
	}

	TowerReader tr(silent);
	clock_t time;
	if (strcmp(args.at(0), "-cv") == 0) {
		if (args.size() != 3){
			displayUsage();
			return 0;
		}
		std::ofstream outputStream(args.at(2));
		time = clock();
		tr.convertWithVerticesOutput(args.at(1), &outputStream);
		time = clock() - time;
		outputStream.close();
	} else if (strcmp(args.at(0), "-cf") == 0) {
		if (args.size() != 3){
			displayUsage();
			return 0;
		}
		std::ofstream outputStream(args.at(2));
		time = clock();
		tr.convertWithActivityOutput(args.at(1), &outputStream);
		time = clock() - time;
		outputStream.close();
	} else
#ifdef PHAT
	if (strcmp(args.at(0), "-cphat") == 0) {
		if (args.size() != 3){
			displayUsage();
			return 0;
		}
		time = clock();
		tr.computePersistenceWithPhat(args.at(1), args.at(2));
		time = clock() - time;
	} else
#endif
#ifdef GUDHI
	if (strcmp(args.at(0), "-cgudhi") == 0) {
		if (args.size() != 3){
			displayUsage();
			return 0;
		}
		time = clock();
		tr.computePersistenceWithGudhi(args.at(1), args.at(2));
		time = clock() - time;
	} else
#endif
		if (strcmp(args.at(0), "-cp") == 0) {
		if (args.size() != 4){
			displayUsage();
			return 0;
		}
		time = clock();
		tr.convertAndComputePersistence(args.at(1), args.at(2), atof(args.at(3)));
		time = clock() - time;
	} else if (strcmp(args.at(0), "-p") == 0) {
		if (args.size() != 4){
			displayUsage();
			return 0;
		}
		time = clock();
		tr.computePersistence(args.at(1), args.at(2), atof(args.at(3)));
		time = clock() - time;
	} else {
		displayUsage();
		return 0;
	}

	if (!silent) cout << "Total elapsed time: " << std::setprecision(0) << ((float)time)/CLOCKS_PER_SEC << " seconds.\n";
	return 0;
}

