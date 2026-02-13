#pragma once

#include <CGAL/Real_timer.h>

// Correct and slow
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;

#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
 


//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Lazy_exact_nt.h>
//typedef CGAL::Lazy_exact_nt<double> NT;
//typedef CGAL::Simple_cartesian<NT> K;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

//#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>

struct HalfEdgeData
{
    int id;
    bool isPseudoSingular;
};

typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arr_extended_dcel<Traits_2, bool, HalfEdgeData, int> Dcel;
typedef CGAL::Arrangement_with_history_2<Traits_2, Dcel> Arrangement_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;


typedef CGAL::Arr_segment_2<K> Curve_2;

typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Face_handle Face_handle;
typedef Arrangement_2::Vertex_const_handle Vertex_const_handle;
typedef Arrangement_2::Vertex_iterator Vertex_iterator;

// Const iterators and handles
typedef Arrangement_2::Face_const_iterator Face_const_iterator;

//typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle Face_const_handle;




struct MySegment_2
{
    Segment_2 seg;
    std::optional<Halfedge_handle> originatingHalfEdge;
    int originatingRegularEdgeId;

    MySegment_2(const Segment_2& s, std::optional<Halfedge_handle> he = std::nullopt, int originatingEdge = -1)
        : seg(s), originatingHalfEdge(he), originatingRegularEdgeId(originatingEdge) {}
};


//typedef CGAL::Box_intersection_d::Box_d<double, 2> Box;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 2, const MySegment_2*> Box;
typedef CGAL::Bbox_2 Bbox;


// Used to find which face a point is in 
//#include <CGAL/Arr_trapezoid_ric_point_location.h>
//typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> Point_location;

//#include <CGAL/Arr_naive_point_location.h>
//typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

#include <CGAL/Arr_walk_along_line_point_location.h>
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;

#include <CGAL/Intersections_2/Segment_2_Segment_2.h>  // Correct header for CGAL 6




//
// Used for computing Barycentric coordinates for fibers
//
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Kernel/global_functions_2.h>  // For bounded_side_2
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> CartesianKernel;
typedef CartesianKernel::Point_2 CartesianPoint;
typedef CartesianKernel::Segment_2 CartesianSegment;
typedef CartesianKernel::Line_2 CartesianLine;
typedef CGAL::Polygon_2<CartesianKernel> CartesianPolygon_2;

#include <CGAL/centroid.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>


#include <CGAL/Union_find.h>






#include <CGAL/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>
#include <CGAL/algorithm.h>

// Custom interval that stores an ID
struct IntervalWithID : public CGAL::Interval_skip_list_interval<K::FT>
{
    int id; // your segment index

    IntervalWithID() : CGAL::Interval_skip_list_interval<K::FT>(), id(-1) {}

    IntervalWithID(const K::FT& low, const K::FT& high, int id_)
        : CGAL::Interval_skip_list_interval<K::FT>(low, high), id(id_) {}
};

typedef CGAL::Interval_skip_list<IntervalWithID> Interval_skip_listId;
