#include "Celllist.hpp"
#include <cassert>
#include <iostream>
#include <cmath>

/*class Celllist{
	public:
	Celllist(){};
	void generateCelllist(double cutoff, double dx);
	int cellNeighbour(int c, int n);
	
	std::vector<std::vector<int>> celllist_;
	
	private:
	int ncell_;
	
};
*/
Celllist::Celllist(){};

void Celllist::generateCelllist(double cutoff, double dx){
	ncell_ = std::ceil(1./cutoff);
/*	while(ncell_*L < 1.){
		L = L+dx;
		ncell_ = 1./L;
		std::cout << "L = " << L << " ncell = " << ncell_ << std::endl;
	}
	*/
	
	//const int meshptpercell = M/ncell;
	celllist_ = std::vector<std::vector<int>> (ncell_*ncell_);

	//save index of mesh points of each cell
	for(int i = 0; i < ncell_; ++i){
		for(int j = 0; j < ncell_; ++j){
			for(int k1 = 0; k1 < cutoff/dx; ++k1){
				for(int k2 = 0; k2<cutoff/dx; ++k2){
					celllist_[i*ncell_+j].push_back(i*cutoff/dx*ncell_+j*cutoff/dx + k1*cutoff/dx + k2);
				}
			}
		}
	}
	
 }

 int Celllist::cellNeighbour(int c, int n) {
	assert(n <= 8);
	
	int col = c % ncell_;
    int row = c / ncell_;

    switch( n % 3 )
    {
		case 0:     break;
        case 1: col = (col+ncell_-1) % ncell_; break;
        case 2: col = (col       +1) % ncell_; break;
    }
    switch( n / 3 )
    {
    case 0:         break;
    case 1: row = (row+ncell_-1) % ncell_; break;
    case 2: row = (row       +1) % ncell_; break;
    }

    return ncell_*row + col;
}

 
 int Celllist::cellNeighbour_openBC(int c, int n) {
	assert(n <= 8);
	
	int col = c % ncell_;
    int row = c / ncell_;

    switch( n % 3 )
    {
		case 0:     break;
        case 1: col = (col - 1); break;
        case 2: col = (col + 1); break;
    }
    switch( n / 3 )
    {
    case 0:         break;
    case 1: row = (row - 1); break;
    case 2: row = (row + 1); break;
    }

	//if neighbour is outside of the mesh: return -1 to indicate empty ghost cell.
	//else, return index of neighbour cell.
	if(col >= ncell_ || col < 0 || row >= ncell_ || row < 0)
		return -1;
	else
		return ncell_*row + col;
}

int Celllist::getNumberOfCells(){
	return ncell_;
}
