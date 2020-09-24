/*
 * Life.cpp
 *
 *  Created on: Feb 16, 2020
 *      Author: Nicolas
 */

#include "Life.h"
#include <iostream>
#include <fstream>

using namespace std;
fstream outfile;

/**
 * Life class containing methods for grid set up, iterating the game and communication between processors.
 */

Life::Life(int rows, int cols, int* values)
:rows(rows), cols(cols), values(values), pad_rows(rows + 2),
 pad_cols(cols + 2), n(rows * cols), pad_n(pad_rows * pad_cols)
{}

/**
 * Writes settings of the game in a txt file for future referencing.
 *
 * @params p: number of processsors
 *		   p_rows: rows of processors
 *		   p_cols: cols of processors
 *		   steps: number of game iterations
 */
void Life::write_settings(int p, int p_rows, int p_cols, int steps){

	fstream settings;
	settings.open("settings.txt", fstream::out);
	settings << this->rows << endl;
	settings << this->cols << endl;
	settings << p << endl;
	settings << p_rows << endl;
	settings << p_cols << endl;
	settings << steps << endl;
	settings.close();
}

/**
 * Writes settings of the sub rows in a txt file for future referencing. 
 *
 * @params sub_rows: rows that processor has
 *	       sub_cols: cols that processor has
 *	       id: current processor
 */
void Life::write_sub_settings(int sub_rows, int sub_cols, int id){
	fstream sub_set;
	sub_set.open("sub_set_"+to_string(id)+".txt", fstream::out);
	sub_set << sub_rows << endl;
	sub_set << sub_cols << endl;
	sub_set.close();
}

/**
 * Writes game iteration to txt file for each processor
 *
 * @params id: current processor
 */
void Life::display(int id){ // write to file
	outfile.open("output_"+to_string(id)+".txt", fstream::out | fstream::app);
	for (int i = 0; i < this->rows; i++){
		for (int j = 0; j < this->cols; j++)
		    outfile << this->values[i*this->cols + j] << " ";
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}

/**
 * Writes game iteration to screen
 */
void Life::screen(){ // write to screen
	for (int i = 0; i < this->rows; i++){
		for (int j = 0; j < this->cols; j++)
		    cout << this->values[i*this->cols+j] << " ";
		cout << endl;
	}
}

/**
 * Pads the Life matrix
 *
 * @returns padded: padded Life object
 */
Life* Life::pad(){ // add layer of zeros around the input object
	int* pad_values = new int[this->pad_n];
	for (int i = 0; i < this->pad_n; i++) pad_values[i] = 0;
	auto* padded = new Life(this->pad_rows, this->pad_cols, pad_values);

	for (int i = 0; i < this->rows; i++){
		for (int j = 0; j < this->cols; j++)
			padded->values[(i + 1)*(this->cols + 2) + j + 1] = this->values[i*this->cols + j];
	}
	return padded;
	delete[] pad_values;
}

/**
 * Unpads the Life matrix
 *
 * @returns padded: unpadded Life object
 */
Life* Life::unpad(){ // remove of outer layer
	int* unpad_vals = new int[(this->rows - 2) * (this->cols - 2)];
	for (int i = 0; i < (this->rows - 2) * (this->cols - 2); i++) unpad_vals[i] = 0;
	auto* unpadded = new Life(this->rows - 2, this->cols - 2, unpad_vals);

	for (int i = 0; i < unpadded->rows; i++){
		for (int j = 0; j < unpadded->cols; j++)
			unpadded->values[i*unpadded->cols + j] = this->values[(i + 1)*(unpadded->cols + 2) + j + 1];
	}
	return unpadded;
	delete[] unpad_vals;
}

/**
 * Randomizes the game grid with alive and dead cells according to prob.
 *
 * @params prob: percentage of values that should be alive.
 */
void Life::randomize(double prob){
	for (int i = 0; i < this->rows * this->cols; i++)
		if (rand() % 100 > 100*prob)
		    this->values[i] = 0;
		else
			this->values[i] = 1;
}

/**
 * Rules for the next iteration of the Game of Life
 *
 * @params uni: reference to Life object
 		   s: current cell
 */
int Life::rule(Life& uni, int s){
	int neighbors = 0;
	for (int i = (0-1); i < 2; i++)
		neighbors += uni.values[s - (uni.cols) + i] + uni.values[s + (uni.cols) + i];
	neighbors += uni.values[s - 1] + uni.values[s + 1];

	// next step fertility settings
	if ((uni.values[s] == 1) && (2 <= neighbors && neighbors <= 3))
		return 1;
	else if ((uni.values[s] == 0) && (neighbors == 3))
		return 1;
	else
		return 0;
}

/**
 * Plays one iteration of the Game of Life
 *
 * @params matrix: reference to Life object
 		   init: current cell
 */
void Life::life(Life& matrix, Life& init){ // matrix is not updated, only init is
	auto* values = new int[init.n];
	for (int i = 0; i < init.n; i++) values[i] = 0;
	auto* tmp = new Life(init.rows, init.cols, values); // pass in values of whatever i want to play the game instead

	for (int i = 0; i < init.n; i++) tmp->values[i] = 0;

	for (int i = 0; i < matrix.rows; i++){
		for (int j = 0; j < matrix.cols; j++){
			tmp->values[(i+1)*init.cols + (j+1)] = rule(init, (i+1)*init.cols + (j+1));
		}
	}

	for (int i = 0; i < init.n; i++) init.values[i] = tmp->values[i]; // update init

	auto* disp = init.unpad();
//	cout << endl;
	delete disp;
	delete[] values;
	delete tmp;
}

/**
 * Finds dimensions of each processor grid and grid of processors
 *
 * @params p: number of processors
 		   p_rows: number of row of processors
 		   p_cols: number of columns of processors
 		   sub_rows: rows controlled by processor
 		   sub_cols: columns controlled by processor
 		   id: current processor
 */
void Life::find_dimensions(int p, int& p_rows, int& p_cols, int& sub_rows, int& sub_cols, int id) {
	int min_gap = max(this->rows, this->cols); // gap can never be bigger than rows or columns

	for (int i = 1; i <= p; i++){
		if (p % i == 0){
			int gap = abs(this->cols / (p / i) - this->rows / i); // most square subgrid
			if (gap < min_gap){
				min_gap = gap;
				sub_rows = this->rows / i;
				sub_cols = this->cols / (p / i);
				p_rows = i;
				p_cols = p / i;
			}
		}
	}
	int rows_add = this->rows % p_rows; // how many p_rows needs a row added

	// rows
	if (id < rows_add * p_cols)
		sub_rows++;

	// cols
	int cols_add = this->cols % p_cols; // how many p_cols needs a col added
	for (int i = 0; i < cols_add; i++)
		if (id % p_cols == i)
			sub_cols++;
}

/**
 * Obtains index of processors within the grid of processors
 *
 * @params id: current processor
 		   id_row: row that current processor is located at
 		   id_col: column that current processor is located at
 		   p_cols: number of columns of processors
 */
void Life::id_to_index(int id, int& id_row, int& id_col, int p_cols){
	id_col = id % p_cols;
	id_row = id / p_cols;
}


/**
 * Obtains ID of processors from its location in the grid of processors
 *
 * @params id_row: row that current processor is located at
 		   id_col: column that current processor is located at
 		   p_rows: number of rows of processors
 		   p_cols: number of columns of processors
 		   periodic: wrap-around or static boundaries
 */
int Life::id_from_index(int id_row, int id_col, int p_rows, int p_cols, bool periodic){
	if (periodic)
		return ((id_row + p_rows) % p_rows) * p_cols + ((id_col + p_cols) % p_cols);

	else {
		if (id_row >= p_rows || id_row < 0)
			return -1;
		if (id_col >= p_cols || id_col < 0)
			return -1;

		return id_row * p_cols + id_col;
	}
}

/**
 * Creates MPI datatypes to communicate rows, columns and corners to other processors. Then communicates and recieves them.
 *
 * @params p: number of processors
 		   p_rows: number of rows of processors
 		   p_cols: number of columns of processors
 		   id_row: row that current processor is located at
 		   id_col: column that current processor is located at
 		   periodic: wrap-around or static boundaries
 */
void Life::make_send_datatypes(int p, int p_rows, int p_cols, int id_row, int id_col, bool periodic) {
	vector<vector<int>> block_length(16); // 16 data types
	vector<vector<MPI_Aint>> addresses(16);
	vector<vector<MPI_Datatype>> typelist(16);
	vector<MPI_Datatype> types(16);
	vector<vector<int>> recv_data(16);

	// type 0, send top left corner
	block_length[0].push_back(1);
	MPI_Aint tmp0;
	MPI_Get_address(&this->values[this->cols + 1], &tmp0);
	addresses[0].push_back(tmp0);
	typelist[0].push_back(MPI_INT);

	// type 1, send top
	block_length[1].push_back(pad_cols - 2);
	MPI_Aint tmp1;
	MPI_Get_address(&this->values[this->cols + 1], &tmp1);
	addresses[1].push_back(tmp1);
	typelist[1].push_back(MPI_INT);

	// type 2, send top right corner
	block_length[2].push_back(1);
	MPI_Aint tmp2;
	MPI_Get_address(&this->values[this->cols + this->cols-2], &tmp2);
	addresses[2].push_back(tmp2);
	typelist[2].push_back(MPI_INT);

	// type 3, send left side
	for (int i = 1; i < this->rows - 1; i++){
		block_length[3].push_back(1);
		MPI_Aint tmp3;
		MPI_Get_address(&this->values[i * this->cols + 1], &tmp3);
		addresses[3].push_back(tmp3);
		typelist[3].push_back(MPI_INT);
	}

	// type 4, send right side
	for (int i = 1; i < this->rows - 1; i++){
		block_length[4].push_back(1);
		MPI_Aint tmp4;
		MPI_Get_address(&this->values[i * this->cols + (this->cols -2)], &tmp4);
		addresses[4].push_back(tmp4);
		typelist[4].push_back(MPI_INT);
	}

	// type 5, send bottom left corner
	block_length[5].push_back(1);
	MPI_Aint tmp5;
	MPI_Get_address(&this->values[((this->rows - 2) * this->cols) + 1], &tmp5);
	addresses[5].push_back(tmp5);
	typelist[5].push_back(MPI_INT);

	// type 6, send bottom row
	block_length[6].push_back(this->cols - 2);
	MPI_Aint tmp6;
	MPI_Get_address(&this->values[((this->rows - 2) * this->cols) + 1], &tmp6);
	addresses[6].push_back(tmp6);
	typelist[6].push_back(MPI_INT);

	// type 7, send bottom right corner
	block_length[7].push_back(1);
	MPI_Aint tmp7;
	MPI_Get_address(&this->values[((this->rows - 2) * this->cols) + (this->cols - 2)], &tmp7);
	addresses[7].push_back(tmp7);
	typelist[7].push_back(MPI_INT);

	// type 8, receive top left corner
	block_length[8].push_back(1);
	MPI_Aint tmp8;
	MPI_Get_address(&this->values[0], &tmp8);
	addresses[8].push_back(tmp8);
	typelist[8].push_back(MPI_INT);

	// type 9, receive top row
	block_length[9].push_back(pad_cols - 2);
	MPI_Aint tmp9;
	MPI_Get_address(&this->values[1], &tmp9);
	addresses[9].push_back(tmp9);
	typelist[9].push_back(MPI_INT);

	// type 10, receive top right corner
	block_length[10].push_back(1);
	MPI_Aint tmp10;
	MPI_Get_address(&this->values[this->cols - 1], &tmp10);
	addresses[10].push_back(tmp10);
	typelist[10].push_back(MPI_INT);

	// type 11, receive left side
	for (int i = 1; i < this->rows - 1; i++){
		block_length[11].push_back(1);
		MPI_Aint tmp11;
		MPI_Get_address(&this->values[i * this->cols], &tmp11);
		addresses[11].push_back(tmp11);
		typelist[11].push_back(MPI_INT);
	}

	// type 12, receive right side
	for (int i = 1; i < this->rows - 1; i++){
		block_length[12].push_back(1);
		MPI_Aint tmp12;
		MPI_Get_address(&this->values[i * this->cols + (this->cols - 1)], &tmp12);
		addresses[12].push_back(tmp12);
		typelist[12].push_back(MPI_INT);
	}

	// type 13, receive bottom left side
	block_length[13].push_back(1);
	MPI_Aint tmp13;
	MPI_Get_address(&this->values[(this->rows - 1) * this->cols], &tmp13);
	addresses[13].push_back(tmp13);
	typelist[13].push_back(MPI_INT);

	// type 14, receive bottom row
	block_length[14].push_back(this->cols - 2);
	MPI_Aint tmp14;
	MPI_Get_address(&this->values[((this->rows - 1) * this->cols) + 1], &tmp14);
	addresses[14].push_back(tmp14);
	typelist[14].push_back(MPI_INT);

	// type 15, receive bottom right
	block_length[15].push_back(1);
	MPI_Aint tmp15;
	MPI_Get_address(&this->values[(this->rows*this->cols)-1], &tmp15);
	addresses[15].push_back(tmp15);
	typelist[15].push_back(MPI_INT);

	// create types
	for (int i = 0; i < 16; i++){
		MPI_Type_create_struct(block_length[i].size(), &block_length[i][0], &addresses[i][0], &typelist[i][0], &types[i]);
		MPI_Type_commit(&types[i]);
	}

	// send and receive
	int cnt_type_send = 0;
	MPI_Request* request = new MPI_Request[8 * 2];
	int cnt = 0;
	int tag_num = 0;
	int iter = 0;

	for (int i = -1; i <= 1; i++){
		for (int j = -1; j <= 1; j++){
			int com_i = id_row + i;
			int com_j = id_col + j;
			int com_id = id_from_index(com_i, com_j, p_rows, p_cols, periodic);
			if (iter != 4){ // does not send or receive from itself
				if (com_id >= 0 && com_id < p){
					MPI_Isend(MPI_BOTTOM, 1, types[cnt_type_send], com_id, tag_num, MPI_COMM_WORLD, &request[cnt * 2]); // send with rank ascending
					MPI_Irecv(MPI_BOTTOM, 1, types[cnt_type_send + 8], com_id, (8 - tag_num), MPI_COMM_WORLD, &request[cnt * 2 + 1]); // recv but with rank descending from 8
					cnt++;
				}
				cnt_type_send++;
			}
			tag_num++;
			iter++;
		}
	}
	MPI_Waitall(2 * cnt, request, MPI_STATUSES_IGNORE) == MPI_SUCCESS;  // wait for all communications before continuing
	for (int i = 0; i < 16; i++)
		MPI_Type_free(&types[i]);
	delete[] request;
}
























