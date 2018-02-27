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
* Definitions of Member functions for point class
-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
* Definitions of Member functions for polygonPoint class
-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
* Definitions of Member functions for edge class
-----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
* Find the nearest point on an edge to a query point using orthagonal projection
* q - query point
------------------------------------------------------------------------------*/
polygonPoint edge::nearestPointOnEdge(point &q) {
	point n;                             // return variable
	double alpha_numer, alpha_denom, alpha_val;
	double dx1, dy1, dx2, dy2;
	
	point fmv(get_fm_vertex());
	point tov(get_to_vertex());
	
	// translate so from node of edge is the origin
	dx1 = q.get_x() - fmv.get_x();
	dy1 = q.get_y() - fmv.get_y();
	dx2 = tov.get_x() - fmv.get_x();
	dy2 = tov.get_y() - fmv.get_y();
	
	if (dx2 == 0) {
		if (dy2 == 0) { // edge is zero length-nearest point is either edge vertice
			n.reset(fmv);
		}
		else {
			// edge is vertical
			n.reset_x(fmv.get_x());
			if (dy2 > 0) {
				if (dy1 < 0) {
					n.reset_y(fmv.get_y());
				}
				else if (q.get_y() > tov.get_y()) {
					n.reset_y(tov.get_y());
				}
				else {
					n.reset_y(q.get_y());
				}
			}
			else {
				if (q.get_y() < tov.get_y()) {
					n.reset_y(tov.get_y());
				}
				else if (dy1 > 0) {
					n.reset_y(fmv.get_y());
				}
				else {
					n.reset_y(q.get_y());
				}
			}
		}
	}
	else { // edge is not vertical or zero length
		alpha_numer = dx1*dx2 + dy1*dy2; // might be negative
		if (alpha_numer < 0) {
			n.reset(fmv);
		}
		else {
			alpha_denom = dx2*dx2 + dy2*dy2; // must be positive
			alpha_val = alpha_numer/alpha_denom;
			if (alpha_val >= 1) {
				point tmpv(get_to_vertex());
				n.reset(tmpv);
			}
			else {
				double newx = alpha_val*dx2 + fmv.get_x();
				double newy = alpha_val*dy2 + fmv.get_y();
				n.reset(newx, newy);
			}
		}
	}
	return polygonPoint(n, fm_rec_index);
}

//
polygonPoint edge::nearestPointOnEdge(polygonPoint & q) {
	point tmpp(q.get_point());
	return nearestPointOnEdge(tmpp);
}

//
bool edge::edgeIsInBox(box & b) {
	int p1_xpos, p1_ypos, p2_xpos, p2_ypos;
	
	// -1 below both, 0 between, 1 above both
	if(fm_vertex.get_x() < b.get_xMin()) {
		p1_xpos = -1;
	}
	else if(fm_vertex.get_x() > b.get_xMax()) {
		p1_xpos = 1;
	}
	else {
		p1_xpos = 0;
	}
	if(to_vertex.get_x() < b.get_xMin()) {
		p2_xpos = -1;
	}
	else if(to_vertex.get_x() > b.get_xMax()) {
		p2_xpos = 1;
	}
	else {
		p2_xpos = 0;
	}
	if(fm_vertex.get_y() < b.get_yMin()) {
		p1_ypos = -1;
	}
	else if(fm_vertex.get_y() > b.get_yMax()) {
		p1_ypos = 1;
	}
	else {
		p1_ypos = 0;
	}
	if(to_vertex.get_y() < b.get_yMin()) {
		p2_ypos = -1;
	}
	else if(to_vertex.get_y() > b.get_yMax()) {
		p2_ypos = 1;
	}
	else {
		p2_ypos = 0;
	}
	
	if(((p1_xpos == 0)&&(p1_ypos == 0))||((p2_xpos == 0)&&(p2_ypos == 0))) {
		return true;
	}
	else if((p1_xpos != 0)&&(p1_xpos == p2_xpos)) {
		return false;
	}
	else if((p1_ypos != 0)&&(p1_ypos == p2_ypos)) {
		return false;
	}
	else if((p1_xpos != p2_xpos)&&(p1_ypos == 0)&&(p2_ypos == 0)) {
		return true;
	}
	else if((p1_ypos != p2_ypos)&&(p1_xpos == 0)&&(p2_ypos == 0)) {
		return true;
	}
	else {
		// do triangle method test for intersection
		double m; // slope of edge
		double a1, a2;  // areas of triangles
		m = (to_vertex.get_y() - fm_vertex.get_y())/(to_vertex.get_x() - fm_vertex.get_x());
		if(m > 0) {
			a1 = fm_vertex.get_x()*to_vertex.get_y() - fm_vertex.get_y()*to_vertex.get_x() + fm_vertex.get_y()*b.get_xMin() -
				fm_vertex.get_x()*b.get_yMax() + to_vertex.get_x()*b.get_yMax() - 
				b.get_xMin()*to_vertex.get_y();
			a2 = fm_vertex.get_x()*to_vertex.get_y() - fm_vertex.get_y()*to_vertex.get_x() + fm_vertex.get_y()*b.get_xMax() -
				fm_vertex.get_x()*b.get_yMin() + to_vertex.get_x()*b.get_yMin() - 
				b.get_xMax()*to_vertex.get_y();
			if (((a1 > 0)&&(a2 > 0))||((a1 < 0)&&(a2 < 0))) {
				return false;
			}
			else {
				return true;
			}
		}
		else {
			a1 = fm_vertex.get_x()*to_vertex.get_y() - fm_vertex.get_y()*to_vertex.get_x() + fm_vertex.get_y()*b.get_xMax() -
				fm_vertex.get_x()*b.get_yMax() + to_vertex.get_x()*b.get_yMax() - 
				b.get_xMax()*to_vertex.get_y();
			a2 = fm_vertex.get_x()*to_vertex.get_y() - fm_vertex.get_y()*to_vertex.get_x() + fm_vertex.get_y()*b.get_xMin() -
				fm_vertex.get_x()*b.get_yMin() + to_vertex.get_x()*b.get_yMin() - 
				b.get_xMin()*to_vertex.get_y();
			if (((a1 > 0)&&(a2 > 0))||((a1 < 0)&&(a2 < 0))) {
				return false;
			}
			else {
				return true;
			}
		}
	}
}


// Here, the interesection of 2 line segments is calculated, and a polygonPoint
// object is returned.  If the intersection point falls on the line segments
// (including the from), the polyrec_index is set to 1.  Otherwise, the
// the polyrec_index is set to VECT_INVALID_INDEX.  The validity of polyrec_index
// can be tested with polyrec_index_is_valid(), and the interesection point can
// be extracted with get_point().
polygonPoint edge::intersection(edge & e) {
	double Ax, Ay, Bx, By, Cx, Cy, Dx, Dy;
	double m1, m2, resX, resY;
    long resID;
	Ax = fm_vertex.get_x();
	Ay = fm_vertex.get_y();
	Bx = to_vertex.get_x();
	By = to_vertex.get_y();
	Cx = e.fm_vertex.get_x();
	Cy = e.fm_vertex.get_y();
	Dx = e.to_vertex.get_x();
	Dy = e.to_vertex.get_y();
	
	box thisBox(Ax, Ay, Bx, By);
	box eBox(Cx, Cy, Dx, Dy);

    if (((Ax == Bx)&&(Cx == Dx))||((Ay == By)&&(Cy == Dy))) {
        if ((Ax == Cx)&&(Ay == Cy)) { // both edges are the same point
            resX = Ax;
            resY = Ay;
            resID = 1;
        }
        else {
            resX = NAN;
            resY = NAN;
            resID = VECT_INVALID_INDEX;
        }
        
    }
    else if ((Ax == Bx)&&(Cy == Dy)) { // this ver, e horiz
        resX = Ax;
        resY = Cy;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else if ((Ay == By)&&(Cx == Dx)) { // this horz, e vert
        resX = Cx;
        resY = Ay;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else if (Ax == Bx) { // this is vert, e is not vert or horiz
        resX = Ax;
        // solve for resY
        m2 = (Dy - Cy)/(Dx - Cx);
        resY = m2*(resX - Cx) + Cy;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else if (Cx == Dx) { // e is vert, this is not vert or horiz
        resX = Cx;
        // solve for resY
        m1 = (By - Ay)/(Bx - Ax);
        resY = m1*(resX - Ax) + Ay;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else if (Ay == By) { // this is horiz, e is not vert or horiz
        resY = Ay;
        // solve for resX
        m2 = (Dy - Cy)/(Dx - Cx);
        resX = (resY - Cy)/m2 + Cx;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else if (Cy == Dy) { // e is horiz, this is not vert or horiz
        resY = Cy;
        // solve for resX
        m1 = (By - Ay)/(Bx - Ax);
        resX = (By - Ay)/m1 +Ax;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    else { // both this and e are not vert or horiz
        // solve for resX and resY
        m1 = (By - Ay)/(Bx - Ax);
        m2 = (Dy - Cy)/(Dx - Cx);
        resX = (Cy - Ay + m1*Ax - m2*Cx)/(m1 - m2);
        resY = m1*resX + Ay - m1*Ax;
        resID = 2; // must check again to see if point in both bounding boxes
    }
    if (resID == 2) {
        if ((thisBox.pointIsInBox(resX, resY))&&(eBox.pointIsInBox(resX, resY))) {
            resID = 1;
        }
        else { 
            resID = VECT_INVALID_INDEX;
        }
    }
    return polygonPoint(resX, resY, resID);
}

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
