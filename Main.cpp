/*
 * Life.cpp
 *
 *  Created on: Feb 12, 2020
 *      Author: Nicolas
 */

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <math.h>
#include <fstream>
#include <vector>
#include <ctime>
#include <chrono>
#include <algorithm>
#include "Life.h"
#include "Life.cpp"

using namespace std;

int id, p;

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 1000);

	// initial matrix to start
	int rows = 10;
	int cols = 10;
	int pad_cols = cols + 2;
	int pad_rows = rows + 2;
	int n = cols * rows;
	int pad_n = pad_rows * pad_cols;
	bool periodic;
	int p_rows;
	int p_cols;
	int sub_rows;
	int sub_cols;
	int id_col;
	int id_row;
	int steps = 100;
	double prob = 0.5;

	periodic = true;
	fstream outfile;
	outfile.open("output_"+to_string(id)+".txt", fstream::out); // clears file to be written in
	outfile.close();

	// initialize board starting state
	auto* values = new int[rows * cols];
	for (int i = 0; i < n; i++) values[i] = 0;
	auto* board = new Life(rows, cols, values);
	
	//set processors their own grids.
	board->find_dimensions(p, p_rows, p_cols, sub_rows, sub_cols, id);
	board->id_to_index(id, id_row, id_col, p_cols);
	board->id_from_index(id_row, id_col, p_rows, p_cols, periodic);

	// initialize the core sub boards
	auto* s_values = new int[sub_rows * sub_cols];
	for (int i = 0; i < sub_rows * sub_cols; i++) // initialize array
		s_values[i] = 0;

	auto* sub_board = new Life(sub_rows, sub_cols, s_values);
	sub_board->randomize(prob);
	auto* sub_pad = sub_board->pad();

	// play the game
	auto t_start = chrono::high_resolution_clock::now();
	for (int i = 0; i < steps; i++){
		sub_pad->make_send_datatypes(p, p_rows, p_cols, id_row, id_col, periodic);
		sub_pad->life(*sub_pad->unpad(), *sub_pad);

		outfile.open("output_" + to_string(id) + ".txt", fstream::out | fstream::app);
		sub_pad->unpad()->display(id);
		outfile.close();
		outfile << endl;
		outfile << endl;
	}
	auto t_end = chrono::high_resolution_clock::now();

	// calculate times
	chrono::duration<double> fs = t_end - t_start;
	if (id == 0){
	    cout << "Processor 0 time: " << fs.count() << "s\n" << endl;
	// make settings file
			board->write_settings(p, p_rows, p_cols, steps);
		}
	board->write_sub_settings(sub_rows, sub_cols, id);

	delete sub_pad;
	delete[] s_values;
	delete sub_board;
	delete[] values;
	delete board;
	MPI_Finalize();
}

