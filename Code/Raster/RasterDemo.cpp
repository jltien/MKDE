/* ***
*/

// include
#include "basicFunctions.h"
#include "file_io.h"
#include "raster.h"

using namespace std;

int main(int argc, char * argv[]) {
    // 1. Make a raster
    // vars
    long nrow = 20, ncol=30;
    double xmin = -15.0, ymin = -10.0, cs = 1.0;
    gridFloat rst(1, nrow, ncol, xmin, ymin, cs); // contains to data
    // 2. Set values in raster
    long r2chng = 2, c2chng = 2;
    rst.setAllCellsToZero(true); // set all cells to 0.0, change cells with no data
    cout << "rst[" << r2chng << "][" << c2chng << "] = " << rst.getGridValue(r2chng, c2chng) << endl;
    rst.setGridValue(r2chng, c2chng, 2.0); // rw, cl, val; sets the cell value
    cout << "rst[" << r2chng << "][" << c2chng << "] = " << rst.getGridValue(r2chng, c2chng) << endl;
    rst.addValueToGridCell(r2chng, c2chng, 1.0); // rw, cl, val; adds val to current value in cell
    cout << "rst[" << r2chng << "][" << c2chng << "] = " << rst.getGridValue(r2chng, c2chng) << endl;
    // 3. Write raster to file
    rst.summary();
    rst.printESRIascii("raster_demo.asc");
    //
    return 0;
}