/*
	Created by:		Brendan McCabe
	Date:			10/7/2018
	Description:	This is a class used to find a shortest path on a graph using the a* algorithm
				    You must first set the dimensions of the graph and then build it.  You can then
					set the current location, target location, and the obstacles. Once everything is set
					you can run the find_path function.  To get the path iterate through the parent cells
					starting from target and ending at your position.
*/

#include "MyPathfinder.h"
#include <iostream>
#include <algorithm>


MyPathfinder::MyPathfinder()
{
	dimensions = new int(2);
	target = NULL;
	location = NULL;
	diagonal_dist = 14;
	adjacent_dist = 10;
}


MyPathfinder::~MyPathfinder()
{
	while (!graph.empty()) {
		std::vector<Cell*>* next = graph.at(0);
		while (!next->empty()) {
			Cell* erasing = (next->at(0));
			next->erase(next->begin());
			delete erasing;
		}
		graph.erase(graph.begin());
	}
}

/*
Function name:	set_dimensions
Parameters:		rows(int): number of rows on the graph
				columns(int): number of columns on the graph
Return:			none
Description:	sets the dimensions of the graph
*/
void MyPathfinder::set_dimensions(int rows, int columns) {
	dimensions[0] = rows;
	dimensions[1] = columns;
}

/*
Function name:	set_target
Parameters:		x(int): column that the target is on
				y(int): row that the target is on
Return:			none
Description:	sets location of target
*/
void MyPathfinder::set_target(int x, int y) {
	if (x < dimensions[0] && y < dimensions[1] && x >= 0 && y >= 0) {
		target = graph.at(x)->at(y);
	}
	else {
		//CHANGE
		target = graph.at(1)->at(1);
		//CHANGE
	}
}

/*
Function name:	set_location
Parameters:		x(int): column of the current position
				y(int): row of the current position
Return:			none
Description:	sets the location of the current position
*/
void MyPathfinder::set_location(int x, int y) {
	if (x < dimensions[0] && y < dimensions[1]) {
		location = graph.at(x)->at(y);
	}
	else {
		location = NULL;
	}
}

/*
Function name:	set_obstacle
Parameters:		x(int): column that the obstacle is on
				y(int): row that the obstacle is on
Return:			none
Description:	places an obstacle on the graph
*/
void MyPathfinder::set_obstacle(int x, int y) {
	if (x < dimensions[0] && y < dimensions[1]) {
		(*(*(graph.at(x))).at(y)).obstacle = true;
	}
}

/*
Function name:	set_distances
Parameters:		adj(int): horizontal and vertical distances
				diag(int): diagonal distances
Return:			none
Description:	sets the distance from one cell to its neighbors
*/
void MyPathfinder::set_distances(int adj, int diag) {
	diagonal_dist = diag;
	adjacent_dist = adj;
}

/*
Function name:	find_path
Parameters:		none
Return:			none
Description:	finds the shortest path on the current graph. Afterwards this path
				can be traced backwards by iterating through the parents starting from
				the target cell
*/
void MyPathfinder::find_path() {
	clear_graph();
	target->is_target = true;
	int neighbors[][2] = { {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1} };
	if (!graph.empty()) {
		std::vector<Cell*> discovered;
		std::vector<Cell*> visited;
		Cell *current = graph.at(location->position[0])->at(location->position[1]);
		while (!current->is_target) {
			std::cout << current->position[0] << ", " << current->position[1] << "\n";
			for (int i = 0; i < 8; i++) {
				int neighbor_pos[2];
				neighbor_pos[0] = current->position[0] + neighbors[i][0];
				neighbor_pos[1] = current->position[1] + neighbors[i][1];
				if (valid_neighbor(neighbor_pos)) {
					Cell* next_neighbor = graph.at(neighbor_pos[0])->at(neighbor_pos[1]);
					if (!next_neighbor->obstacle && std::find(visited.begin(), visited.end(), next_neighbor) == visited.end()) {
						calc_gcost(next_neighbor, current);
						calc_fcost(next_neighbor);
						next_neighbor->heuristic = next_neighbor->f_cost + next_neighbor->g_cost;
						if(std::find(visited.begin(), visited.end(), next_neighbor) == visited.end() && std::find(discovered.begin(), discovered.end(), next_neighbor) == discovered.end())
							discovered.push_back(next_neighbor);
					}
				}
			}
			//std::cout << "END\n\n";
			visited.push_back(current);
			if (std::find(discovered.begin(), discovered.end(), current) != discovered.end()) {
				discovered.erase(std::find(discovered.begin(), discovered.end(), current));
			}
			auto smallest = std::min_element(discovered.begin(), discovered.end(),[](Cell* a, Cell* b) {return a->heuristic < b->heuristic;});
			if (smallest != discovered.end())
				current = *smallest;
			else
				break;
		}
	}
}

/*NOT USED
*/
bool MyPathfinder::cell_comp(Cell a, Cell b){
	return a.heuristic < b.heuristic; 
}

/*
Function name:	valid_neighbor
Parameters:		neighbor_pos(int): the coordinates being checked 
Return:			none
Description:	checks if coordinates are within the dimensions of the graph
*/
bool MyPathfinder::valid_neighbor(int* neighbor_pos) {
	if (	neighbor_pos[0] >= 0
			&& neighbor_pos[0] < dimensions[0]
			&& neighbor_pos[1] >= 0
			&& neighbor_pos[1] < dimensions[1])
		return true;
	return false;
}

/*
Function name:	calc_gcost
Parameters:		cell(int): neighbor cell with the gcost being updated
				current(int): the currently visited cell 
Return:			none
Description:	updates the gcost of a cell 
*/
void MyPathfinder::calc_gcost(Cell* cell, Cell* current) {
	int temp_g = 0;
	temp_g += current->g_cost;
	if (current->position[0] == cell->position[0] || current->position[1] == cell->position[1]) {
		temp_g += adjacent_dist;
	}
	else {
		temp_g += diagonal_dist;
	}
	if (cell->discovered) {
		if (cell->g_cost > temp_g) {
			cell->parent = current;
		}
		cell->g_cost = std::min(cell->g_cost, temp_g);
	}
	else {
		cell->g_cost = temp_g;
		cell->discovered = true;
		cell->parent = current;
	}
}

/*
Function name:	calc_fcost
Parameters:		cell(int): neighbor cell with the fcost being updated
Return:			none
Description:	updates the fcost of a cell
*/
void MyPathfinder::calc_fcost(Cell* cell) {
	cell->f_cost = 0;
	int xdist = abs(cell->position[0] - (*target).position[0]);
	int ydist = abs(cell->position[1] - (*target).position[1]);
	cell->f_cost = 0;
	cell->f_cost += diagonal_dist * std::min(xdist, ydist);
	cell->f_cost += adjacent_dist * abs(xdist - ydist);
}


/*
Function name:	clear_graph
Parameters:		none
Return:			none
Description:	resets all the values in the cells of the graph but does not remove obstacles
*/
void MyPathfinder::clear_graph() {
	path.clear();
	if (!graph.empty()) {
		std::vector<std::vector<Cell*>*>::iterator row;
		std::vector<Cell*>::iterator column;
		for (row = graph.begin(); row != graph.end(); row++) {
			for (column = (*(*row)).begin(); column != (*(*row)).end(); column++) {
				(*(*column)).f_cost = 0;
				(*(*column)).g_cost = 0;
				(*(*column)).heuristic = 0;
				(*(*column)).is_target = false;
				(*(*column)).discovered = false;
				(*(*column)).parent = NULL;
			}
		}
	}
}


/*
Function name:	build_graph
Parameters:		none
Return:			none
Description:	allocates memory for graph and initializes all cell values
*/
void MyPathfinder::build_graph() {
	if (dimensions[0] > 0 && dimensions[1] > 0) {
		for (int i = 0; i < dimensions[0]; i++) {
			std::vector<Cell*>* next_row = new std::vector<Cell*>;
			for (int j = 0; j < dimensions[1]; j++) {
				Cell* next_col = new Cell;
				(*next_col).f_cost = 0;
				(*next_col).g_cost = 0;
				(*next_col).heuristic = 0;
				(*next_col).is_target = false;
				(*next_col).obstacle = false;
				(*next_col).parent = NULL;
				(*next_col).position[0] = i;
				(*next_col).position[1] = j;
				(*next_col).discovered = false;
				(*next_row).push_back(next_col);
			}
			graph.push_back(next_row);
		}
	}
}

