/*
 * Life.h
 *
 *  Created on: Feb 16, 2020
 *      Author: Nicolas
 */

#ifndef LIFE_H_
#define LIFE_H_
#pragma once

class Life {
public:
	Life(int rows, int cols, int* values);

	void write_settings(int p, int p_rows, int p_cols, int steps);
	void write_sub_settings(int sub_rows, int sub_cols, int id);

	void display(int id);
	void screen();

	Life* pad();
	Life* unpad();

	void randomize(double prob);

	int rule(Life& uni, int s);
	void life(Life& matrix, Life& init);

	void find_dimensions(int p, int& p_rows, int& p_cols, int& sub_rows, int& sub_cols, int id);
	void id_to_index(int id, int& id_row, int& id_col, int p_cols);
	int id_from_index(int id_row, int id_col, int p_rows, int p_cols, bool periodic);

	void make_send_datatypes(int p, int p_rows, int p_cols, int id_row, int id_col, bool periodic);

	int rows = -1;
	int cols = -1;
	int* values = nullptr;
	int pad_rows = rows + 2;
	int pad_cols = cols + 2;
	int n = rows * cols;
	int pad_n = pad_rows * pad_cols;

private:
};

#endif /* LIFE_H_ */
