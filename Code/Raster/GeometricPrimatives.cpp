/******************************************************************************
* GeometricPrimatives.cpp                                                     *
*                                                                             *
* Definitions for basic geometric data types.                                 *
*                                                                             *
* Author: Jeff A. Tracey                                                      *
* Created: 29 November 2002.                                                  *
* Copyright Jeff A. Tracey 2002. All Rights Reserved.                         *
******************************************************************************/

#include "/home/jeff/projects/myCPPlib/GeometricPrimatives.h"

#ifdef _ANSI_STD_H // using ANSI/ISO standard header file names
using namespace std;
#endif


/*-----------------------------------------------------------------------------
* Definitions of Member functions for box class
-----------------------------------------------------------------------------*/

// test to see if an edge (defined by 2 point objects) intersects a box
// clearly, there is a lot that can be done with this function in the way
// of optimization
// THIS FUNCTION HAS A PROBLEM WHICH CAUSES IT TO RETURN TRUE
// IN SOME CASES WHEN IT SHOULD RETURN FALSE (1/12/07)
bool box::edgeIsInBox(point & p1, point & p2) {
    if ((pointIsInBox(p1))||(pointIsInBox(p2))) { // p1, p2, or p1 and p2 in box
        return true;
    }
    else if ((p1.get_x() < ll_corner.get_x())&&(p2.get_x() < ll_corner.get_x())) {
        return false; // to the left
    }
    else if ((p1.get_x() > ur_corner.get_x())&&(p2.get_x() > ur_corner.get_x())) {
        return false; // to the right
    }
    else if ((p1.get_y() < ll_corner.get_y())&&(p2.get_y() < ll_corner.get_y())) {
        return false; // below
    }
    else if ((p1.get_y() > ur_corner.get_y())&&(p2.get_y() > ur_corner.get_y())) {
        return false; // above
    }
    else { // test for intersections with edge
        double bboxXmin = ll_corner.get_x();
        double bboxXmax = ur_corner.get_x();
        double bboxYmin = ll_corner.get_y();
        double bboxYmax = ur_corner.get_y();
        double boxXvert[] = {bboxXmin, bboxXmax, bboxXmax, bboxXmin, bboxXmin};
        double boxYvert[] = {bboxYmin, bboxYmin, bboxYmax, bboxYmax, bboxYmin};
        edge boxEdge;
        edge testEdge(p1, p2);
        int intersectCount = 0;
        polygonPoint testP;
        for (int i = 0; i < 4; i++) {
            boxEdge.reset(boxXvert[i], boxYvert[i], boxXvert[i+1], boxYvert[i+1]);
            testP = boxEdge.intersection(testEdge);
            if (testP.polyrec_index_is_valid()) intersectCount++; 
        }
        if (intersectCount == 2) return true;
        else return false;
    }
}

bool box::edgeIsInBox(edge & e) {
	correct();
	point tmpp1(e.get_fm_vertex());
	point tmpp2(e.get_to_vertex());
	return edgeIsInBox(tmpp1, tmpp2);
}

void box::boxIntersection(box & b) {
	correct();
	b.correct();
	// DO THE BOXES INTERSECT AT ALL?
	if((ll_corner.get_x() >= b.get_xMax())||(ur_corner.get_x() <= b.get_xMin())||
	(ll_corner.get_y() >= b.get_yMax())||(ur_corner.get_y() <= b.get_yMin())) {
		reset(0.0, 0.0, 0.0, 0.0); // make a box with zero area
	}
	else {
        //
		if(b.get_xMin() > ll_corner.get_x()) ll_corner.reset_x(b.get_xMin());
		if(b.get_xMax() < ur_corner.get_x()) ur_corner.reset_x(b.get_xMax());
		if(b.get_yMin() > ll_corner.get_y()) ll_corner.reset_y(b.get_yMin());
		if(b.get_yMax() < ur_corner.get_y()) ur_corner.reset_y(b.get_yMax());
	}
}

void box::boxUnion(box & b) {
	correct();
	b.correct();
	if(b.get_xMin() < ll_corner.get_x()) ll_corner.reset_x(b.get_xMin());
	if(b.get_xMax() > ur_corner.get_x()) ur_corner.reset_x(b.get_xMax());
	if(b.get_yMin() < ll_corner.get_y()) ll_corner.reset_y(b.get_yMin());
	if(b.get_yMax() > ur_corner.get_y()) ur_corner.reset_y(b.get_yMax());
}
