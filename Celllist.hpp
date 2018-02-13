#ifndef CELLLIST
#define CELLLIST

#include <vector>

class Celllist{
	public:
	Celllist();
	void generateCelllist(double cutoff, double dx);
	int cellNeighbour_openBC(int c, int n);
	int cellNeighbour(int c, int n);
	std::vector<std::vector<int>> celllist_;
	int getNumberOfCells();
	
	private:
	int ncell_;

};

#endif
