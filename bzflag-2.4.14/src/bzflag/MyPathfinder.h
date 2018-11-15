/*
	Created by:		Brendan McCabe
	Date:			10/7/2018
	Description:	This is a class used to find a shortest path on a graph using the a* algorithm
					You must first set the dimensions of the graph and then build it.  You can then
					set the current location, target location, and the obstacles. Once everything is set
					you can run the find_path function.  To get the path iterate through the parent cells
					starting from target and ending at your position.
*/

#pragma once
#include <vector>
class MyPathfinder
{
public:

	struct Cell {
		int g_cost;
		int heuristic;
		int f_cost;
		bool obstacle;
		bool is_target;
		bool discovered;
		int position[2];
		Cell* parent;
	};

	bool cell_comp(Cell a, Cell b);

	MyPathfinder();
	~MyPathfinder();

	void set_dimensions(int rows, int columns);
	void set_target(int x, int y);
	void set_location(int x, int y);
	void set_obstacle(int x, int y);
	void set_distances(int adj, int diag);
	void find_path();
	bool valid_neighbor(int* neighbor_pos);
	void clear_graph();
	void build_graph();
	void calc_gcost(Cell* cell, Cell* current);
	void calc_fcost(Cell* cell);
	

	std::vector<std::vector<Cell*>*> graph;
	int* dimensions;
	Cell* target;
	Cell* location;
	std::vector<int*> path;
	int diagonal_dist;
	int adjacent_dist;
};

