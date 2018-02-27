/******************************************************************************
* GeometricPrimatives.h                                                       *
*                                                                             *
* Definitions for basic geometric data types.                                 *
*                                                                             *
* Author: Jeff A. Tracey                                                      *
* Created: 29 November 2002.                                                  *
* Copyright Jeff A. Tracey 2002. All Rights Reserved.                         *
******************************************************************************/

#ifndef GEO_PRIMATIVES_H
#define GEO_PRIMATIVES_H

// include files
#include "/home/jeff/projects/myCPPlib/basicDefinitions.h"
#include "/home/jeff/projects/myCPPlib/basicFunctions.h"

#ifdef _ANSI_STD_H // using ANSI/ISO standard header file names
#include <iostream>
#else              // NOT using ANSI/ISO standard header file names
#include <iostream.h>
#endif

#include <cmath>
#include <cfloat>
#include <cstdio>

// NOTES:
// const long VECT_INVALID_INDEX is defined in ibmm_defs.h
// const double PI is defined in ibmm_defs.h
// LATER:
// It may be useful to provide functions for different transformations
//		such as translate x, translate y, sheer, rotate, etc.

/*-----------------------------------------------------------------------------
* Some general use geometric functions
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
* Forward declarations
-----------------------------------------------------------------------------*/
class polygonPoint;
class vect;
class box;

/*-----------------------------------------------------------------------------
* Declaration Point
-----------------------------------------------------------------------------*/
class point {
	private:
		double x; // x-coordinate
		double y; // y-coordinate
		// 16 bytes of data
	public:
		point(); // default constructor
		point(double xx, double yy);      // constructor
		point(polygonPoint & old_p);      // constructor
		point(const point &old_p);        // copy constructor
		double get_x();                   // 
		double get_y();                   // 
		void reset_x(double xx);          // 
		void reset_y(double yy);          // 
		void reset(point & p);            // 
		void reset(polygonPoint & p);     // 
		void reset(double xx, double yy); // 
		double distance(point & q);       // get distance from the point to point q
		double distance(polygonPoint & q); // get distance from the point to point q
		double angle(point & q);          // get angle from the point to point q on (-pi, pi]
		double angle(polygonPoint & q);   // get angle from the point to point q on (-pi, pi]
		vect get_vector_to(point & q);    // 
		vect get_vector_to(polygonPoint & q); // 
		bool equals(point & p);           //
		void print();                     // 
};

/*-----------------------------------------------------------------------------
* Declaration for polygonPoint class (contains an index to a polygon record)
-----------------------------------------------------------------------------*/
class polygonPoint {
	private:
		double x;           // x-coordinate
		double y;           // y-coordinate
		long polyrec_index; // index to the polygon record across the boundary
		// 20 bytes of data
	public:
		polygonPoint();                               // default constructor
		polygonPoint(double xx, double yy);           // constructor
		polygonPoint(double xx, double yy, long pind); // constructor
		polygonPoint(point &old_p);                   // constructor
		polygonPoint(point &old_p, long pind);        // constructor
		polygonPoint(const polygonPoint &old_p);      // copy constructor
		double get_x();                               //
		double get_y();                               //
		point get_point();                            //
		long get_polyrec_index();                     //
		bool polyrec_index_is_valid();                //
		void reset_x(double xx);                      //
		void reset_y(double yy);                      // 
		void reset_polyrec_index(long pind);          //
		void reset(point & p);                        //
		void reset(polygonPoint & p);                 //
		void reset(point & p, long pind);             //
		void reset(double xx, double yy);             //
		void reset(double xx, double yy, long pind);  //
		double distance(point & q);                   //
		double distance(polygonPoint & q);            //
		double angle(point & q);                      //
		double angle(polygonPoint & q);               //
		vect get_vector_to(point & q);                //
		vect get_vector_to(polygonPoint & q);         //
		bool equals(polygonPoint & p);                //
		void print();                                 //
};

/*-----------------------------------------------------------------------------
* Declaration for edge class
-----------------------------------------------------------------------------*/
class edge {
	private:
		point fm_vertex;                    // coordinates of the from vertex
		point to_vertex;                    // coordinates of the to vertex
		long fm_rec_index;                  // index to polygon on left (across boundary polygon)
		long to_rec_index;                  // 
		// 40 bytes of data
	public:
		edge();                                     // default constructor
		edge(point & f, point & t);                 // constructor
		edge(polygonPoint & f, polygonPoint & t);   // constructor
		edge(point & f, vect & v);                  // constructor
		edge(polygonPoint & f, vect & v);           // constructor
		edge(double fx, double fy, double tx, double ty); // constructor
		edge(double fx, double fy, long f_pind, double tx, double ty, long t_pind); // constructor
		edge(const edge & e); // copy constructor
		void reset_fm_vertex(point & f);            //
		void reset_fm_vertex(polygonPoint & f);     //
		void reset_to_vertex(point & t);            //
		void reset_to_vertex(polygonPoint & t);     //
		void reset(double fx, double fy, double tx, double ty); //
		void reset(double fx, double fy, long f_pind, double tx, double ty, long t_pind); //
		void reset(point & fp, point & tp);         //
		void reset(point & fp, long fi, point & tp, long ti); //
		void reset(polygonPoint & fp, polygonPoint & tp); //
		void reset(edge & e);                       //
		void reset_fm_rec_index(long f_pind);       //
		void reset_to_rec_index(long t_pind);       //
		void flip();                                //
		point get_fm_vertex();                      //
		point get_to_vertex();                      //
		long get_fm_rec_index();                    //
		long get_to_rec_index();                    //
		polygonPoint nearestPointOnEdge(point & q); //
		polygonPoint nearestPointOnEdge(polygonPoint & q); //
		vect vectorToNearPointOnEdge(point & q);    //
		vect vectorToNearPointOnEdge(polygonPoint & q);   //
		double lengthOfEdge();                      //
		double angleOfEdge();                       //
		vect edgeToVector();                        //
		bool edgeIsInBox(box & b);                  //
		bool pointLocationIntersect(point & q);     //
		polygonPoint intersection(edge & e);        //
		bool equals(edge & e);                      //
		bool compare(edge & e, double snap_dist);   //
		void print();                               //
};

/*-----------------------------------------------------------------------------
* Class definition for a vector defined by an angle and length
* NOTE: Any input for a length < 0 is interpreted to mean "in the opposite
* direction, so the absolute value of the length is taken and pi radians is
* added to the angle.  Any input with length == 0 is not a valid vector, so
* ang and leng are set to 0, and is_null is set to true
-----------------------------------------------------------------------------*/
class vect {
	private:
		double ang;  // angle of vector
		double leng; // length of vector
		bool is_null; // true if vector is null
		// 17 bytes of data (not sure what a bool uses)
	public:
		vect();                                     // default constructor
		vect(double l, double a);                   // constructor
		vect(double l, double a, bool n);           // constructor
		vect(point & fmv, point & tov);             // constructor with 2 points as input
		vect(polygonPoint & fmv, polygonPoint & tov); // constructor with 2 polygonPoints as input
		vect(edge & e);                             // constructor with edge as input
		vect(const vect & v);                       // copy constructor
		double get_angle() { return ang; }          //
		double get_compass();                       //
		double get_length();                        //
		double get_deggrad();                       //
		void reset_angle(double a);                 //
		void reset_length(double l);                //
		void reset(double l, double a);             //
		void reset(vect & v);                       //
		double get_dx();                            //
		double get_dy();                            //
		bool is_null_vect() {return is_null; }      //
		vect operator+ (const vect & v1);           // overloaded "+"
		vect operator- (const vect & v1);           // overloaded "-"
		bool equals(vect & v);                      //
		void print();                               //
};

/*-----------------------------------------------------------------------------
* Definition for box class -- a minimum bounding box (mbb)
-----------------------------------------------------------------------------*/
class box {
	private:
		point ll_corner;
		point ur_corner;
		void correct();                                       //
		// 16 bytes of data
	public:
		box();                                                // default constructor
		box(point & ll, point & ur);                          // constructor
		box(polygonPoint & ll, polygonPoint & ur);            // constructor
		box(edge & e);                                        // constructor
		box(double xMin, double yMin, double xMax, double yMax); // constructor
		box(const box & b);                                   // copy constructor
		void reset(box & b);                                  //
		void reset(edge & e);                                 //
		void reset(point & p1, point & p2);                   //
		void reset(polygonPoint & p1, polygonPoint & p2);     //
		void reset(double xMin, double yMin, double xMax, double yMax);  //
		void reset_xMin(double xx);                           //
		void reset_xMax(double xx);                           //
		void reset_yMin(double yy);                           //
		void reset_yMax(double yy);                           //
		void expand(point & p);                               //
		void expand(polygonPoint & p);                        //
		bool pointIsInBox(point & p);                         //
		bool pointIsInBox(polygonPoint & p);                  //
        bool pointIsInBox(double x, double y);                //
        bool edgeIsInBox(edge & e);                           //
		bool edgeIsInBox(point & p1, point & p2);             //
		bool boxIntersectsBox(box & b);                       //
		bool contains(box & b);                               //
		bool isContainedBy(box & b);                          //
		void boxIntersection(box & b);                        //
		void boxUnion(box & b);                               //
		double get_xMin();                                    //
		double get_yMin();                                    //
		double get_xMax();                                    //
		double get_yMax();                                    //
		point get_ll();                                       //
		point get_ur();                                       //
		point get_midpoint();                                 //
		box get_NWquad();                                     //
		box get_NEquad();                                     //
		box get_SWquad();                                     //
		box get_SEquad();                                     //
		double area();                                        //
		bool equals(box & b);                                 //
		void print();                                         //
};

/*-----------------------------------------------------------------------------
* Definition for cirlce class
-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
* Member Functions for Point
-----------------------------------------------------------------------------*/
inline point::point() {
	x = 0.0;
	y = 0.0;
}

inline point::point(double xx, double yy) {
	x = xx;
	y = yy;
}

inline point::point(polygonPoint & old_p) {
	x = old_p.get_x();
	y = old_p.get_y();
}

inline point::point(const point &old_p) {
	x = old_p.x;
	y = old_p.y;
}

inline double point::get_x() {
	return x;
}

inline double point::get_y() {
	return y;
}

inline void point::reset_x(double xx) {
	x = xx;
}

inline void point::reset_y(double yy) {
	y = yy;
}

inline void point::reset(point & p) {
	x = p.x;
	y = p.y;
}

inline void point::reset(polygonPoint & p) {
	x = p.get_x();
	y = p.get_y();
}

inline void point::reset(double xx, double yy) {
	x = xx;
	y = yy;
}

inline double point::distance(point & q) {
	double dx, dy;
	dx = q.x - x;
	dy = q.y - y;
	return hypot(dx, dy);
}

inline double point::distance(polygonPoint & q) {
	double dx, dy;
	dx = q.get_x() - x;
	dy = q.get_y() - y;
	return hypot(dx, dy);
}

inline double point::angle(point & q) {
	double delta_x = q.x - x;
	double delta_y = q.y - y;
	
	if ((delta_x == 0) && (delta_y == 0)) return NAN;
    else return atan2(delta_y, delta_x);
}

inline double point::angle(polygonPoint & q) {
	double delta_x = q.get_x() - x;
	double delta_y = q.get_y() - y;
		
	if ((delta_x == 0) && (delta_y == 0)) return NAN;
	return atan2(delta_y, delta_x);
}

inline vect point::get_vector_to(point & q) {
	return vect(distance(q), angle(q));
}
 
inline vect point::get_vector_to(polygonPoint & q) {
	return vect(distance(q), angle(q));
}

inline bool point::equals(point & p) {
	if((x == p.x)&&(y==p.y)) return true;
	else return false;
}

inline void point::print() {
	std::cout << "POINT: (" << x << ", " << y << ")" << std::endl;
}

/*-----------------------------------------------------------------------------
* Member Functions for polygonPoint class
-----------------------------------------------------------------------------*/
inline polygonPoint::polygonPoint() {
	x = 0;
	y = 0;
	polyrec_index = VECT_INVALID_INDEX;
}

inline polygonPoint::polygonPoint(double xx, double yy) {
	x = xx;
	y = yy;
	polyrec_index = VECT_INVALID_INDEX;
}

inline polygonPoint::polygonPoint(double xx, double yy, long pind) {
	x = xx;
	y = yy;
	if(pind >= 0) polyrec_index = pind;
	else polyrec_index = VECT_INVALID_INDEX;
}

inline polygonPoint::polygonPoint(point &old_p) {
	x = old_p.get_x();
	y = old_p.get_y();
	polyrec_index = VECT_INVALID_INDEX;
}

inline polygonPoint::polygonPoint(point &old_p, long pind) {
	x = old_p.get_x();
	y = old_p.get_y();
	if(pind >= 0) polyrec_index = pind;
	else polyrec_index = VECT_INVALID_INDEX;
}

inline polygonPoint::polygonPoint(const polygonPoint &old_p) {
	x = old_p.x;
	y = old_p.y;
	polyrec_index = old_p.polyrec_index;
}

inline double polygonPoint::get_x() {
	return x;
}

inline double polygonPoint::get_y() {
	return y;
}

inline point polygonPoint::get_point() {
	return point(x, y);
}

inline long polygonPoint::get_polyrec_index() {
	return polyrec_index;
}

inline bool polygonPoint::polyrec_index_is_valid() {
	if(polyrec_index == VECT_INVALID_INDEX) return false;
	else return true;
}

// here, nothing is changed but the x-coord (the
// polyrec_index is not set to VECT_INVALID_INDEX)
inline void polygonPoint::reset_x(double xx) {
	x = xx;
}

// here, nothing is changed but the x-coord (the
// polyrec_index is not set to VECT_INVALID_INDEX)
inline void polygonPoint::reset_y(double yy) {
	y = yy;
}

inline void polygonPoint::reset_polyrec_index(long pind) {
	if(pind >= 0) polyrec_index = pind;
	else polyrec_index = VECT_INVALID_INDEX;
}

inline void polygonPoint::reset(point & p) {
	x = p.get_x();
	y = p.get_y();
	polyrec_index = VECT_INVALID_INDEX;
}

inline void polygonPoint::reset(polygonPoint & p) {
	x = p.get_x();
	y = p.get_y();
	polyrec_index = p.get_polyrec_index();
}

inline void polygonPoint::reset(point & p, long pind) {
	x = p.get_x();
	y = p.get_y();
	if(pind >= 0) polyrec_index = pind;
	else polyrec_index = VECT_INVALID_INDEX;
}


inline void polygonPoint::reset(double xx, double yy) {
	x = xx;
	y = yy;
	polyrec_index = VECT_INVALID_INDEX;
}

inline void polygonPoint::reset(double xx, double yy, long pind) {
	x = xx;
	y = yy;
	if(pind >= 0) polyrec_index = pind;
	else polyrec_index = VECT_INVALID_INDEX;
}

inline double polygonPoint::distance(point & q) {
	double out_dist;
	double dx, dy;
	dx = q.get_x() - x;
	dy = q.get_y() - y;
	out_dist = sqrt(dx*dx + dy*dy);
	return out_dist;
}

inline double polygonPoint::distance(polygonPoint & q) {
	double out_dist;
	double dx, dy;
	dx = q.x - x;
	dy = q.y - y;
	out_dist = sqrt(dx*dx + dy*dy);
	return out_dist;
}

inline double polygonPoint::angle(point & q) {
	double delta_x;
	double delta_y;
	double direction;
	
	delta_x = q.get_x() - x;
	delta_y = q.get_y() - y;
		
	if ((delta_x == 0) && (delta_y == 0)) {
		std::cerr << "Error in point::angle(): from point and to point are same." << std::endl;
		std::exit(1);
	}
	else {
		direction = atan2(delta_y, delta_x);
	}
	return direction;
}

inline double polygonPoint::angle(polygonPoint & q) {
	double delta_x;
	double delta_y;
	double direction;
	
	delta_x = q.x - x;
	delta_y = q.y - y;
		
	if ((delta_x == 0) && (delta_y == 0)) {
		std::cerr << "Error in point::angle(): from point and to point are same." << std::endl;
		std::exit(1);
	}
	else {
		direction = atan2(delta_y, delta_x);
	}
	return direction;
}

inline vect polygonPoint::get_vector_to(point & q) {
	return vect(distance(q), angle(q));
}

inline vect polygonPoint::get_vector_to(polygonPoint & q) {
	return vect(distance(q), angle(q));
}

inline bool polygonPoint::equals(polygonPoint & p) {
	if((x==p.x)&&(y==p.y)&&(polyrec_index==p.polyrec_index)) return true;
	else return false;
}

inline void polygonPoint::print() {
	std::cout << "POLYGONPOINT: (" << x << ", " << y << ") POLYRECINDEX: " << polyrec_index << std::endl;
}

/*-----------------------------------------------------------------------------
* Member Functions for edge class
-----------------------------------------------------------------------------*/
inline edge::edge() : fm_vertex(), to_vertex() {
	fm_rec_index = VECT_INVALID_INDEX;
	to_rec_index = VECT_INVALID_INDEX;
}

inline edge::edge(point & f, point & t) : fm_vertex(f), to_vertex(t) {
	fm_rec_index = VECT_INVALID_INDEX;
	to_rec_index = VECT_INVALID_INDEX;
}

inline edge::edge(polygonPoint & f, polygonPoint & t) : fm_vertex(f), to_vertex(t) {
	fm_rec_index = f.get_polyrec_index();
	to_rec_index = t.get_polyrec_index();
}

inline edge::edge(point & f, vect & v) : fm_vertex(f), to_vertex(f.get_x() + v.get_length()*cos(v.get_angle()), f.get_y() + v.get_length()*sin(v.get_angle())) {
}

inline edge::edge(polygonPoint & f, vect & v) : fm_vertex(f), to_vertex(f.get_x() + v.get_length()*cos(v.get_angle()), f.get_y() + v.get_length()*sin(v.get_angle())) {
}

inline edge::edge(double fx, double fy, double tx, double ty) : fm_vertex(fx, fy), to_vertex(tx, ty){
	fm_rec_index = VECT_INVALID_INDEX;
	to_rec_index = VECT_INVALID_INDEX;
}

inline edge::edge(double fx, double fy, long f_pind, double tx, double ty, long t_pind) : fm_vertex(fx, fy), to_vertex(tx, ty){
	if(f_pind >= 0) fm_rec_index = f_pind;
	else fm_rec_index = VECT_INVALID_INDEX;
	if(t_pind >= 0) to_rec_index = t_pind;
	else to_rec_index = VECT_INVALID_INDEX;
}

inline edge::edge(const edge & e) : fm_vertex(e.fm_vertex), to_vertex(e.to_vertex) {
	fm_rec_index = e.fm_rec_index;
	to_rec_index = e.to_rec_index;
}

inline void edge::reset_fm_vertex(point & f) {
	fm_vertex.reset(f);
	fm_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset_fm_vertex(polygonPoint & f) {
	fm_vertex.reset(f);
	fm_rec_index = f.get_polyrec_index();
}

inline void edge::reset_to_vertex(point & t) {
	to_vertex.reset(t);
	to_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset_to_vertex(polygonPoint & t) {
	to_vertex.reset(t);
	to_rec_index = t.get_polyrec_index();
}

inline void edge::reset(double fx, double fy, double tx, double ty) {
	fm_vertex.reset(fx, fy);
	fm_rec_index = VECT_INVALID_INDEX;
	to_vertex.reset(tx, ty);
	to_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset(double fx, double fy, long f_pind, double tx, double ty, long t_pind) {
	fm_vertex.reset(fx, fy);
	if(f_pind < 0) {
		fm_rec_index = VECT_INVALID_INDEX;
	}
	else {
		fm_rec_index = f_pind;
	}
	to_vertex.reset(tx, ty);
	if(t_pind < 0) {
		to_rec_index = VECT_INVALID_INDEX;
	}
	else {
		to_rec_index = t_pind;
	}
}

inline void edge::reset(point & fp, point & tp) {
	fm_vertex.reset(fp);
	fm_rec_index = VECT_INVALID_INDEX;
	to_vertex.reset(tp);
	to_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset(point & fp, long fi, point & tp, long ti) {
	fm_vertex.reset(fp);
	if(fi >= 0) fm_rec_index = fi;
	else fm_rec_index = VECT_INVALID_INDEX;
	to_vertex.reset(tp);
	if(ti >= 0) to_rec_index = ti;
	else to_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset(polygonPoint & fp, polygonPoint & tp) {
	fm_vertex.reset(fp);
	fm_rec_index = fp.get_polyrec_index();
	to_vertex.reset(tp);
	to_rec_index = tp.get_polyrec_index();
}

inline void edge::reset(edge & e) {
	fm_vertex.reset(e.fm_vertex);
	to_vertex.reset(e.to_vertex);
	fm_rec_index = e.fm_rec_index;
	to_rec_index = e.to_rec_index;
}

inline void edge::reset_fm_rec_index(long f_pind) {
	if(f_pind >= 0) fm_rec_index = f_pind;
	else fm_rec_index = VECT_INVALID_INDEX;
}

inline void edge::reset_to_rec_index(long t_pind) {
	if(t_pind >= 0) to_rec_index = t_pind;
	else to_rec_index = VECT_INVALID_INDEX;
}

inline void edge::flip() {
	point tmp_point = fm_vertex;
	long tmp_index = fm_rec_index;
	fm_vertex = to_vertex;
	fm_rec_index = to_rec_index;
	to_vertex = tmp_point;
	to_rec_index = tmp_index;
}

inline point edge::get_fm_vertex() {
	return fm_vertex;
}

inline point edge::get_to_vertex() {
	return to_vertex;
}

inline long edge::get_fm_rec_index() {
	return fm_rec_index;
}

inline long edge::get_to_rec_index() {
	return to_rec_index;
}

inline vect edge::vectorToNearPointOnEdge(point & q) {
	polygonPoint tmp_point(nearestPointOnEdge(q));
	return vect(q.distance(tmp_point), q.angle(tmp_point));
}

inline vect edge::vectorToNearPointOnEdge(polygonPoint & q) {
	polygonPoint tmp_point(nearestPointOnEdge(q));
	return vect(q.distance(tmp_point), q.angle(tmp_point));
}

inline double edge::lengthOfEdge() {
	return fm_vertex.distance(to_vertex);
}

inline double edge::angleOfEdge() {
	return fm_vertex.angle(to_vertex);
}

inline vect edge::edgeToVector() {
	return vect(fm_vertex.distance(to_vertex), fm_vertex.angle(to_vertex));
}

inline bool edge::pointLocationIntersect(point & q) {
	if(((q.get_y() <= fm_vertex.get_y())&&(q.get_y() < to_vertex.get_y()))||
	((q.get_y() >= fm_vertex.get_y())&&(q.get_y() > to_vertex.get_y()))) {
		return false;
	}
	else if((q.get_x() < fm_vertex.get_x())&&(q.get_x() < to_vertex.get_x())) {
		return false;
	}
	else if((q.get_x() > fm_vertex.get_x())&&(q.get_x() > to_vertex.get_x())) {
		return true;
	}
	else {
		double c, x_star;
		double dy = to_vertex.get_y() - fm_vertex.get_y();
		if(dy != 0.0) {
			c = (q.get_y() - fm_vertex.get_y())/dy;
			x_star = fm_vertex.get_x() + c*(to_vertex.get_x() - fm_vertex.get_x());
			if(q.get_x() <= x_star) {
				return false;
			}
			else {
				return true;
			}
		}
		else {
			if(q.get_y() >= to_vertex.get_y()) {
				return true;
			}
			else {
				return false;
			}
		}
	}
}

inline bool edge::equals(edge & e) {
	if((fm_vertex.equals(e.fm_vertex))&&(to_vertex.equals(e.to_vertex))&&
	(fm_rec_index==e.fm_rec_index)&&(to_rec_index==e.to_rec_index)) return true;
	else return false;
}

inline bool edge::compare(edge & e, double snap_dist) {
	if((fm_vertex.distance(e.fm_vertex) <= snap_dist)&&(to_vertex.distance(e.to_vertex) <= snap_dist)) {
		return true;
	}
	else {
		return false;
	}
}

inline void edge::print() {
	std::cout << "EDGE: FM: (" << fm_vertex.get_x() << ", " << fm_vertex.get_y() << ")";
	std::cout << ",{" << fm_rec_index << "}";
	std::cout << " TO: (" << to_vertex.get_x() << ", " << to_vertex.get_y() << ")";
	std::cout << ",{" << to_rec_index << "}" << std::endl;
}

/*-----------------------------------------------------------------------------
* Member Functions for vect class
-----------------------------------------------------------------------------*/
inline vect::vect() {
	leng = 0;
	ang = 0;
	is_null = true;
}

inline vect::vect(double l, double a) {
	if(l == 0.0) {
		leng = 0;
		ang = 0;
		is_null = true;
	}
	else if(l < 0.0) {
		leng = fabs(l);
		ang = angleModulus(a+PI);
		is_null = false;
	}
	else {
		leng = l;
		ang = angleModulus(a);
		is_null = false;
	}
}

inline vect::vect(double l, double a, bool n) {
	if(l == 0.0) {
		leng = 0;
		ang = 0;
		is_null = true;
	}
	else if(l < 0.0) {
		leng = fabs(l);
		ang = angleModulus(a+PI);
		is_null = n;
	}
	else {
		leng = l;
		ang = angleModulus(a);
		is_null = n;
	}
}

inline vect::vect(point & fmv, point & tov) {
	leng = fmv.distance(tov);
	if(leng == 0.0) {
		ang = 0.0;
		is_null = true;
	}
	else {
		ang = fmv.angle(tov);
		is_null = false;
	}
}

inline vect::vect(polygonPoint & fmv, polygonPoint & tov) {
	leng = fmv.distance(tov);
	if(leng == 0.0) {
		ang = 0.0;
		is_null = true;
	}
	else {
		ang = fmv.angle(tov);
		is_null = false;
	}
}

inline vect::vect(edge & e) {
	leng = e.lengthOfEdge();
	if(leng == 0.0) {
		ang = 0.0;
		is_null = true;
	}
	else {
		ang = e.angleOfEdge();
		is_null = false;
	}
}

inline vect::vect(const vect & v) {
	leng = v.leng;
	ang = v.ang;
	is_null = v.is_null;
}

// return angle in degrees and geographic convention
inline double vect::get_compass() {
	double tmp_deg = (PI/2 - ang)*(180.0/PI);
	if(tmp_deg < 0.0) {
		tmp_deg = 360.0 + tmp_deg;
	}
	return tmp_deg;
}

inline double vect::get_length() {
    if (is_null) return 0.0;
    else return leng;
}

// return length as a gradient in degrees
inline double vect::get_deggrad() {
	return atan2(leng, 1.0)*(180/PI); // check
}

inline void vect::reset_angle(double a) {
	ang = angleModulus(a);
}


inline void vect::reset_length(double l) {
	if(l == 0.0) {
		leng = 0.0;
		ang = 0.0;
		is_null = true;
	}
	else if(l < 0.0) {
		leng = fabs(l);
		ang = angleModulus(ang+PI);
        is_null = false;
	}
	else {
		leng = l;
        is_null = false;
	}
}


inline void vect::reset(double l, double a) {
	if(l == 0.0) {
		leng = 0.0;
		ang = 0.0;
		is_null = true;
	}
	else if(l < 0.0) {
		leng = fabs(l);
		ang = angleModulus(a+PI);
        is_null = false;
	}
	else {
		leng = l;
		ang = angleModulus(a);
        is_null = false;
	}
}


inline void vect::reset(vect & v) {
	leng = v.leng;
	ang = v.ang;
    is_null = v.is_null;
}

inline double vect::get_dx() {
    if (is_null) return 0.0;
    else return leng*cos(ang);
}

inline double vect::get_dy() {
    if (is_null) return 0.0;
    else return leng*sin(ang);
}

// non-member overloaded "+" operator
inline vect vect::operator+ (const vect & v1) {
	if(!is_null && !v1.is_null) {
		double dx = leng*cos(ang) + v1.leng*cos(v1.ang);
		double dy = leng*sin(ang) + v1.leng*sin(v1.ang);
		point p1(0.0, 0.0);
		point p2(dx, dy);
		return vect(p1, p2);
	}
	else {
		return vect(0.0, 0.0, true);
	}
}

// non-member overloaded "-" operator
inline vect vect::operator- (const vect & v1) {
	if(!is_null && !v1.is_null) {
		double dx = leng*cos(ang) - v1.leng*cos(v1.ang);
		double dy = leng*sin(ang) - v1.leng*sin(v1.ang);
		point p1(0.0, 0.0);
		point p2(dx, dy);
		return vect(p1, p2);
	}
	else {
		return vect(0.0, 0.0, true);
	}
}

inline bool vect::equals(vect & v) {
	if((ang==v.ang)&&(leng==v.leng)&&(is_null==v.is_null)) return true;
	else return false;
}

inline void vect::print() {
	std::cout << "VECT: " << "(length = " << leng;
	std::cout << ", angle = " << ang << ") is_null = ";
	if(is_null) {
		std::cout << "true";
	}
	else {
		std::cout << "false";
	}
	std::cout << std::endl;
}

/*-----------------------------------------------------------------------------
* Member Functions for box class
-----------------------------------------------------------------------------*/
inline box::box() : ll_corner(), ur_corner() {
}

inline box::box(point & ll, point & ur) : ll_corner(ll), ur_corner(ur) {
	correct();
}

inline box::box(polygonPoint & ll, polygonPoint & ur) : ll_corner(ll), ur_corner(ur) {
	correct();
}

inline box::box(edge & e) : ll_corner(e.get_fm_vertex()), ur_corner(e.get_to_vertex()) {
	correct();
}

inline box::box(double xMin, double yMin, double xMax, double yMax): ll_corner(xMin, yMin), ur_corner(xMax, yMax) {
	correct();
}

inline box::box(const box & b) : ll_corner(b.ll_corner), ur_corner(b.ur_corner) {
	correct();
}

inline void box::reset(box & b) {
	ll_corner.reset(b.ll_corner);
	ur_corner.reset(b.ur_corner);
	correct();
}

inline void box::reset(edge & e) {
	point tmp1(e.get_fm_vertex());
	point tmp2(e.get_to_vertex());
	ll_corner.reset(tmp1);
	ur_corner.reset(tmp2);
	correct();
}

inline void box::reset(point & p1, point & p2) {
	ll_corner.reset(p1);
	ur_corner.reset(p2);
	correct();
}

inline void box::reset(polygonPoint & p1, polygonPoint & p2) {
	ll_corner.reset(p1);
	ur_corner.reset(p2);
	correct();
}

inline void box::reset(double xMin, double yMin, double xMax, double yMax) {
	ll_corner.reset(xMin, yMin);
	ur_corner.reset(xMax, yMax);
	correct();
}

inline void box::reset_xMin(double xx) {
	ll_corner.reset_x(xx);
	correct();
}

inline void box::reset_xMax(double xx) {
	ur_corner.reset_x(xx);
	correct();
}

inline void box::reset_yMin(double yy) {
	ll_corner.reset_y(yy);
	correct();
}

inline void box::reset_yMax(double yy) {
	ur_corner.reset_y(yy);
	correct();
}

// fix the bounding box if it is incorrect.  Many of the functions
// for this class rely on xMin < xMax and yMin < yMax.  Calling this
// function after any changes will ensure this condition is true
inline void box::correct() {
	if (ll_corner.get_x() > ur_corner.get_x()) {
		double tmpx = ur_corner.get_x();
		ur_corner.reset_x(ll_corner.get_x());
		ll_corner.reset_x(tmpx);
	}
	if (ll_corner.get_y() > ur_corner.get_y()) {
		double tmpy = ur_corner.get_y();
		ur_corner.reset_y(ll_corner.get_y());
		ll_corner.reset_y(tmpy);
	}
}

// expand the bounding box to include a point
inline void box::expand(point & p) {
	if (p.get_x() < ll_corner.get_x()) ll_corner.reset_x(p.get_x());
	if (p.get_y() < ll_corner.get_y()) ll_corner.reset_y(p.get_y());
	if (p.get_x() > ur_corner.get_x()) ur_corner.reset_x(p.get_x());
	if (p.get_y() > ur_corner.get_y()) ur_corner.reset_y(p.get_y());
	correct();
}

inline void box::expand(polygonPoint & p) {
	if (p.get_x() < ll_corner.get_x()) ll_corner.reset_x(p.get_x());
	if (p.get_y() < ll_corner.get_y()) ll_corner.reset_y(p.get_y());
	if (p.get_x() > ur_corner.get_x()) ur_corner.reset_x(p.get_x());
	if (p.get_y() > ur_corner.get_y()) ur_corner.reset_y(p.get_y());
	correct();
}

inline bool box::pointIsInBox(point & p) {
	bool is_in;
	if((p.get_x() < ll_corner.get_x())||(p.get_x() > ur_corner.get_x())||
		(p.get_y() < ll_corner.get_y())||(p.get_y() > ur_corner.get_y())) {
		is_in = false;
	}
	else {
		is_in =  true;
	}
	return is_in;
}

inline bool box::pointIsInBox(polygonPoint & p) {
	bool is_in;
	if((p.get_x() < ll_corner.get_x())||(p.get_x() > ur_corner.get_x())||
		(p.get_y() < ll_corner.get_y())||(p.get_y() > ur_corner.get_y())) {
		is_in = false;
	}
	else {
		is_in =  true;
	}
	return is_in;
}

inline bool box::pointIsInBox(double x, double y) {
	if((x < ll_corner.get_x())||(x > ur_corner.get_x())||
		(y < ll_corner.get_y())||(y > ur_corner.get_y())) {
		return false;
	}
	else {
		return true;
	}
}


inline bool box::boxIntersectsBox(box & b) {
	correct();
	b.correct();
	if(b.get_xMax() < ll_corner.get_x()) {
		return false;
	}
	else if(b.get_xMin() > ur_corner.get_x()) {
		return false;
	}
	else if(b.get_yMax() < ll_corner.get_y()) {
		return false;
	}
	else if(b.get_yMin() > ur_corner.get_y()) {
		return false;
	}
	else {
		return true;
	}
}

// test if box b is in this box
inline bool box::contains(box & b) {
	if((get_xMin() <= b.get_xMin())&&
		(get_yMin() <= b.get_yMin())&&
		(get_xMax() >= b.get_xMax())&&
		(get_yMax() >= b.get_yMax())) {
		return true;
	}
	else {
		return false;
	}
}

// test if this box is in box b
inline bool box::isContainedBy(box & b) {
	if((get_xMin() >= b.get_xMin())&&
		(get_yMin() >= b.get_yMin())&&
		(get_xMax() <= b.get_xMax())&&
		(get_yMax() <= b.get_yMax())) {
		return true;
	}
	else {
		return false;
	}
}

inline double box::get_xMin() {
	return ll_corner.get_x();
}

inline double box::get_xMax() {
	return ur_corner.get_x();
}

inline double box::get_yMin() {
	return ll_corner.get_y();
}

inline double box::get_yMax() {
	return ur_corner.get_y();
}

inline point box::get_ll() {
	return ll_corner;
}

inline point box::get_ur() {
	return ur_corner;
}

inline point box::get_midpoint(){
	double xMid = 0.50*(ur_corner.get_x() + ll_corner.get_x());
	double yMid = 0.50*(ur_corner.get_y() + ll_corner.get_y());
	return point(xMid, yMid);
}

inline box box::get_NWquad() {
	double xMid = 0.50*(ur_corner.get_x() + ll_corner.get_x());
	double yMid = 0.50*(ur_corner.get_y() + ll_corner.get_y());
	return box(ll_corner.get_x(), yMid, xMid, ur_corner.get_y());
}

inline box box::get_NEquad() {
	double xMid = 0.50*(ur_corner.get_x() + ll_corner.get_x());
	double yMid = 0.50*(ur_corner.get_y() + ll_corner.get_y());
	return box(xMid, yMid, ur_corner.get_x(), ur_corner.get_y());
}

inline box box::get_SWquad() {
	double xMid = 0.50*(ur_corner.get_x() + ll_corner.get_x());
	double yMid = 0.50*(ur_corner.get_y() + ll_corner.get_y());
	return box(ll_corner.get_x(), ll_corner.get_y(), xMid, yMid);
}

inline box box::get_SEquad() {
	double xMid = 0.50*(ur_corner.get_x() + ll_corner.get_x());
	double yMid = 0.50*(ur_corner.get_y() + ll_corner.get_y());
	return box(xMid, ll_corner.get_y(), ur_corner.get_x(), yMid);
}

inline double box::area() {
	double a = (ur_corner.get_x() - ll_corner.get_x())*(ur_corner.get_y() - ll_corner.get_y());
	return a;
}

inline bool box::equals(box & b) {
	if((ll_corner.equals(b.ll_corner))&&(ur_corner.equals(b.ur_corner))) return true;
	else return false;
}

inline void box::print() {
	std::cout << "BOX: LL: (" << ll_corner.get_x() << ", " << ll_corner.get_y() << ")";
	std::cout << " UR: (" << ur_corner.get_x() << ", " << ur_corner.get_y() << ")" << std::endl;
}

#endif

