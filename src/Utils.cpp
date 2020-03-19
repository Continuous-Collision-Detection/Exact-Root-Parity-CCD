#include "Utils.hpp"

#include <array>
#include <fstream>

namespace ccd {

int orient3d(
    const Vector3r& a, const Vector3r& b, const Vector3r& c, const Vector3r& d)
{
    const Rational det = (a - d).dot(cross(b - d, c - d));
    return det.get_sign();
}
int orient2d(
    const Vector3r& a, const Vector3r& b, const Vector3r& c, const int axis)
{
    Vector3r v1, v2;
    v1 = a - b;
    v2 = c - b;
    int i1, i2;
    if (axis == 0) {
        i1 = 1;
        i2 = 2;
    }
    if (axis == 1) {
        i1 = 0;
        i2 = 2;
    }
    if (axis == 2) {
        i1 = 0;
        i2 = 1;
    }
    const Rational det = v1[i1] * v2[i2] - v1[i2] * v2[i1];
    return det.get_sign();
}

int origin_ray_triangle_inter(
    const Vector3d& dirf,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3)
{
    const Vector3r dir(dirf[0], dirf[1], dirf[2]);
    const Rational denom = dir[0] * t1[1] * t2[2] - dir[0] * t1[1] * t3[2]
        - dir[0] * t1[2] * t2[1] + dir[0] * t1[2] * t3[1]
        + dir[0] * t2[1] * t3[2] - dir[0] * t2[2] * t3[1]
        - dir[1] * t1[0] * t2[2] + dir[1] * t1[0] * t3[2]
        + dir[1] * t1[2] * t2[0] - dir[1] * t1[2] * t3[0]
        - dir[1] * t2[0] * t3[2] + dir[1] * t2[2] * t3[0]
        + dir[2] * t1[0] * t2[1] - dir[2] * t1[0] * t3[1]
        - dir[2] * t1[1] * t2[0] + dir[2] * t1[1] * t3[0]
        + dir[2] * t2[0] * t3[1] - dir[2] * t2[1] * t3[0];

    // const auto n = cross(t1 - t2, t3 - t2);
    // print(n);
    // std::cout<<denom<<std::endl;

    // infinite intersections
    if (denom.get_sign() == 0)
        return -1;

    // assert(denom.get_sign() > 0);

    // std::ofstream os("blaa.obj");
    // os << "v " << t1[0] << " " << t1[1] << " " << t1[2] << "\n";
    // os << "v " << t2[0] << " " << t2[1] << " " << t2[2] << "\n";
    // os << "v " << t3[0] << " " << t3[1] << " " << t3[2] << "\n";
    // os << "f 1 2 3\n";
    // os.close();

    const Rational u = (-1 * dir[0] * t1[1] * t3[2] + dir[0] * t1[2] * t3[1]
                        + dir[1] * t1[0] * t3[2] - dir[1] * t1[2] * t3[0]
                        - dir[2] * t1[0] * t3[1] + dir[2] * t1[1] * t3[0])
        / denom;
    const Rational v = (dir[0] * t1[1] * t2[2] - dir[0] * t1[2] * t2[1]
                        - dir[1] * t1[0] * t2[2] + dir[1] * t1[2] * t2[0]
                        + dir[2] * t1[0] * t2[1] - dir[2] * t1[1] * t2[0])
        / denom;
    const Rational t = (t1[0] * t2[1] * t3[2] - t1[0] * t2[2] * t3[1]
                        - t1[1] * t2[0] * t3[2] + t1[1] * t2[2] * t3[0]
                        + t1[2] * t2[0] * t3[1] - t1[2] * t2[1] * t3[0])
        / denom;

    // std::cout<<t<<std::endl;
    // std::cout << u << std::endl;
    // std::cout << v << std::endl;

    if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u + v <= 1 && t >= 0) {
        if (t.get_sign() == 0)
            return 2;
        // on a corner
        if (u.get_sign() == 0 || v.get_sign() == 0)
            return -1;

        return 1;
    }

    return 0;
}

bool segment_segment_inter(
    const Vector3r& s0,
    const Vector3r& e0,
    const Vector3r& s1,
    const Vector3r& e1,
    Vector3r& res,
    int axis)
{
    const int i1 = (axis + 1) % 3;
    const int i2 = (axis + 2) % 3;

    const Rational dd = e0[i1] * e1[i2] - e0[i1] * s1[i2] - e0[i2] * e1[i1]
        + e0[i2] * s1[i1] + e1[i1] * s0[i2] - e1[i2] * s0[i1] + s0[i1] * s1[i2]
        - s0[i2] * s1[i1];

    if (dd.get_sign() == 0) {
        return false;
    }

    const Rational t0 = (e1[i1] * s0[i2] - e1[i1] * s1[i2] - e1[i2] * s0[i1]
                         + e1[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1])
        / dd;
    const Rational t1 = (e0[i1] * s0[i2] - e0[i1] * s1[i2] - e0[i2] * s0[i1]
                         + e0[i2] * s1[i1] + s0[i1] * s1[i2] - s0[i2] * s1[i1])
        / dd;

    // we exclude intersection on corners
    if (t0 <= 0 || t0 >= 1 || t1 <= 0 || t1 >= 1) {
        return false;
    }

    res = (1 - t0) * s0 + t0 * e0;
#ifndef NDEBUG
    const Vector3r p1 = (1 - t1) * s1 + t1 * e1;

    assert(res[0] == p1[0] && res[1] == p1[1] && res[2] == p1[2]);
#endif
    return true;
}

bool is_segment_degenerated(const Vector3r& s0, const Vector3r& s1)
{
    if (s0[0] == s1[0] && s0[1] == s1[1] && s0[2] == s1[2])
        return true;
    return false;
}

// point and seg are colinear, now check if point is on the segment(can deal with segment degeneration)
bool colinear_point_on_segment(
    const Vector3r& pt, const Vector3r& s0, const Vector3r& s1)
{
    Vector3r dir1 = pt - s0;
    Vector3r dir2 = s1 - pt;
    if (dir1.dot(dir2) >= 0)
        return true;
    return false;
}
bool point_on_segment(const Vector3r& pt, const Vector3r& s0, const Vector3r& s1) {
	Vector3r dir1 = pt - s0;
	Vector3r dir2 = s1 - pt;
	Vector3r norm = cross(dir1, dir2);
	if (same_point(norm, ORIGIN))
		return colinear_point_on_segment(pt, s0, s1);
	else
		return false;
}
// check if two parallel segments (or any segment is degenerated as a point) have intersection
bool parallel_segments_inter(
    const Vector3r& s0,
    const Vector3r& e0,
    const Vector3r& s1,
    const Vector3r& e1,
    const bool dege0,
    const bool dege1)
{
    if (dege0)
        return colinear_point_on_segment(s0, s1, e1);
    if (dege1)
        return colinear_point_on_segment(s1, s0, e0);
    if (colinear_point_on_segment(s0, s1, e1))
        return true;
    if (colinear_point_on_segment(e0, s1, e1))
        return true;
    if (colinear_point_on_segment(s1, s0, e0))
        return true;
    if (colinear_point_on_segment(e1, s0, e0))
        return true;
    return false;
}
    // segment segment intersection, can deal with degenerated segments
bool segment_segment_intersection(
    const Vector3r& s0,
    const Vector3r& e0,
    const Vector3r& s1,
    const Vector3r& e1)
{
	// if four points not coplanar, not intersected
    if (orient3d(s0, e0, s1, e1) != 0)
        return false;

    bool dege0 = is_segment_degenerated(s0, e0);
    bool dege1 = is_segment_degenerated(s1, e1);
    if (!dege0 && !dege1) {// two segments not degenerated
        if (same_point(
                cross(e0 - s0, e1 - s1),
                ORIGIN)) // if two segments not degenerated are parallel
            return parallel_segments_inter(s0, e0, s1, e1, dege0, dege1);
        
		Vector3r norm = tri_norm(s0, s1, e0);
        if (same_point(norm, ORIGIN))
            return colinear_point_on_segment(s1, s0, e0);

        Vector3r norm1 = tri_norm(s0, s1, e1);
        if (same_point(norm1, ORIGIN))
            return colinear_point_on_segment(s0, s1, e1);

        Vector3r norm2 = tri_norm(e1, s1, e0);
        if (same_point(norm2, ORIGIN))
            return colinear_point_on_segment(e0, e1, s1);
        Vector3r norm3 = tri_norm(e1, s0, e0);
        if (same_point(norm3, ORIGIN))
            return colinear_point_on_segment(e1, s0, e0);
		// if norm is same direction with norm1, norm2, norm3, then intersected
		if (norm.dot(norm1).get_sign() > 0 && norm.dot(norm2).get_sign() > 0 && norm.dot(norm3).get_sign() > 0)
			return true;
		return false;
	}
	if (dege0&&dege1) {
		if (same_point(s0, s1))
			return true;
		else
			return false;
	}
	if (dege0) {
		return colinear_point_on_segment(s0, s1, e1);
	}
	if (dege1) {
		return colinear_point_on_segment(s1, s0, e0);
	}
	std::cout << " impossible to go here" << std::endl;
	return false;
}
// 0 not intersected; 1 intersected; 2 pt on s0
int point_on_ray(const Vector3r& s0,
	const Vector3r& dir0, const Vector3r& pt) {
	if (same_point(s0, pt))
		return 2;
	Vector3r dir1 = pt - s0;
	Vector3r norm = cross(dir0, dir1);
	if (!same_point(norm, ORIGIN))
		return 0;// if not coplanar, false
	if (dir0.dot(dir1) > 0)
		return 1;
	return 0;
}
// 0 not intersected, 1 intersected, 2 s0 on segment
// can deal with degenerated cases
int ray_segment_intersection(
	const Vector3r& s0,
	const Vector3r& dir0,
	const Vector3r& s1,
	const Vector3r& e1) {
	
	if (same_point(e1, s1))//degenerated case
		return point_on_ray(s0, dir0, s1);
	/////////////////////////////////////
	Vector3r norm = cross(s0-s1, e1-s0);
	if (same_point(norm, ORIGIN))
	{
		if (colinear_point_on_segment(s0, s1, e1))
			return 2;
		else return 0;
	}
	else {

		Vector3r norm1 = cross(s0 - s1, dir0);
		if (same_point(norm1, ORIGIN))
			return point_on_ray(s0, dir0, s1);

		Vector3r norm2 = cross(e1 - s0, dir0);
		if (same_point(norm2, ORIGIN)) {
			return point_on_ray(s0, dir0, e1);
		}

		if (norm.dot(norm1) > 0 && norm.dot(norm2) > 0) {
			return 1;
		}
		return 0;
	}


	/////////////////////////////////////
	//Vector3r dir1 = e1 - s1;
	//Vector3r norm = cross(dir1, dir0);
	//if (same_point(norm, ORIGIN))//parallel
	//{
	//	int inter1 = point_on_ray(s0, dir0, s1);
	//	int inter2 = point_on_ray(s0, dir0, e1);
	//	if (inter1 == 0 && inter2 == 0) return 0;
	//	if (inter1 == 2 && inter2 == 2) return 2;
	//	if (inter1 > 0 && inter2 == 0) return 2;
	//	if (inter2 > 0 && inter1 == 0) return 2;
	//	return 1;
	//}
	
}
int line_segment_intersection(
	const Vector3r& s0,
	const Vector3r& dir0,
	const Vector3r& s1,
	const Vector3r& e1) {

	if (same_point(e1, s1))//degenerated case
		return point_on_ray(s0, dir0, s1);
	/////////////////////////////////////
	Vector3r norm = cross(s0 - s1, e1 - s0);
	if (same_point(norm, ORIGIN))
	{
		if (colinear_point_on_segment(s0, s1, e1))
			return 2;
		else return 0;
	}
	else {

		Vector3r norm1 = cross(s0 - s1, dir0);
		if (same_point(norm1, ORIGIN))
			return point_on_ray(s0, dir0, s1);

		Vector3r norm2 = cross(e1 - s0, dir0);
		if (same_point(norm2, ORIGIN)) {
			return point_on_ray(s0, dir0, e1);
		}

		if (norm.dot(norm1) > 0 && norm.dot(norm2) > 0) {
			return 1;
		}
		return 0;
	}

}


void write(const Vector3d& v, std::ostream& out)
{
    out.write(reinterpret_cast<const char*>(&v[0]), sizeof(v[0]));
    out.write(reinterpret_cast<const char*>(&v[1]), sizeof(v[1]));
    out.write(reinterpret_cast<const char*>(&v[2]), sizeof(v[2]));
}

Vector3d read(std::istream& in)
{
    Vector3d res;
    double tmp;
    in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    res[0] = tmp;

    in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    res[1] = tmp;

    in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    res[2] = tmp;

    return res;
}
// 2 on edge, 1 interior, 0 not intersect, 3 intersect OPEN edge t2-t3
// norm follows right hand law
int point_inter_triangle(
	const Vector3r&pt, 
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	const bool& dege, const Vector3r& norm, const bool halfopen) {
	if (dege) {// check 2 edges are enough
		if (point_on_segment(pt, t1, t2))
			return 2;
		if (point_on_segment(pt, t1, t3))
			return 2;
		return 0;
	}
	else {
		if (orient3d(pt, t1, t2, t3) != 0)
			return false;
		/*if (point_on_segment(pt, t1, t2))
			return 2;
		if (point_on_segment(pt, t1, t3))
			return 2;
		if (point_on_segment(pt, t2, t3))
			return 2;*///no need to do above
		Vector3r np = t1 + norm;
		int o1 = orient3d(pt, np, t1, t2);
		int o2 = orient3d(pt, np, t2, t3);// this edge 
		int o3 = orient3d(pt, np, t3, t1);
		if (halfopen) {
			if (o2 == 0 && o1 == o3)
				return 3;// on open edge t2-t3
		}
		
		if (o1 == o2 && o1 == o3)
			return 1;
		if (o1 == 0 && o2 * o3 >= 0)
			return 2;
		if (o2 == 0 && o1 * o3 >= 0)
			return 2;
		if (o3 == 0 && o2* o1 >= 0)
			return 2;
		
		return 0;
	}
}
// we check if triangle intersect segment,
// this function is used in cube edge--prism tri and cube edge--bilinear tri
// if halfopen= true, can tell us if intersect the edge t2-t3
// return 0, 1, 2, 3
// this function is only used to check if seg intersect no degenerated bilinear
//and if seg intersect opposite facets of tet. 
int segment_triangle_intersection(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen) {
	
	Vector3r norm = cross(t2 - t1, t3 - t2);
	if (same_point(norm, ORIGIN))// triangle degeneration
	{
		return 0;//we already checked triangle (at least two edges)edge against cube
		//TODO need to make it a general one?
	}

	int o1 = orient3d(e0, t1, t2, t3);
	int o2 = orient3d(e1, t1, t2, t3);
	if (o1 > 0 && o2 > 0)
		return 0;
	if (o1 < 0 && o2 < 0)
		return 0;

	if (o1 == 0 && o2 != 0)// e0 on the plane
		return point_inter_triangle(e0, t1, t2, t3, false, norm, halfopen);
	if (o2 == 0 && o1 != 0)
		return point_inter_triangle(e1, t1, t2, t3, false, norm, halfopen);
	if (o1 == 0 && o2 == 0)// two points are all on the plane. we already checked seg-two edges before
	{
		if (same_point(e0,e1))
			return point_inter_triangle(e0, t1, t2, t3, false, norm, halfopen);
		else {
			int pinter0 = point_inter_triangle(e0, t1, t2, t3, false, norm, halfopen);//2 is impossible
			int pinter1 = point_inter_triangle(e1, t1, t2, t3, false, norm, halfopen);
			// two points all inside, inside; two outside, outside; others, intersect t2-t3
			if (pinter0 == 1 && pinter1 == 1) return 1; // interior
			if (pinter0 == 0 && pinter1 == 0) return 0;//  out
			if (halfopen)
				return 3;// intersect t2-t3
			else
				return 2;
		}
	}
	return is_line_cut_triangle(e0, e1, t1, t2, t3, halfopen, norm);
}
// 0 no intersection, 1 intersect, 2 point on triangle(including two edges), 3 point on t2-t3 edge, -1 shoot on border
int ray_triangle_intersection(
	const Vector3r& pt,
	const Vector3r& dir,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen) {

	Vector3r norm = cross(t2 - t1, t3 - t2);
	if (same_point(norm, ORIGIN))// triangle degeneration
	{
		int inter1 = ray_segment_intersection(pt, dir, t1, t2);
		if (inter1 == 1) return -1;
		if (inter1 == 2) return 2;

		int inter2 = ray_segment_intersection(pt, dir, t1, t3);// check two segs is enough
		if (inter2 == 1) return -1;
		if (inter2 == 2) return 2;

		return 0;
	}

	int o1 = orient3d(pt, t1, t2, t3);
	if (o1 == 0) {//point is on the plane
		int inter = point_inter_triangle(pt, t1, t2, t3, false, norm, halfopen);
		if (inter == 1 || inter == 2) return 2;
		if (inter == 3) return 3;
		//if (inter == 0) 
		else {// pt on the plane but not intersect triangle.
			if (norm.dot(dir) == 0) {// if ray is on the plane
				int inter1 = ray_segment_intersection(pt, dir, t1, t2);
				if (inter1 == 1) return -1;
				if (inter1 == 2) return 2;// acutally, cannot be 2 because already checked by pt_inter_tri

				int inter2 = ray_segment_intersection(pt, dir, t1, t3);
				if (inter2 == 1) return -1;
				if (inter2 == 2) return 2;
				//actually since point do not intersect triangle, check two segs are enough
				//int inter3 = ray_segment_intersection(pt, dir, t2, t3);
				//// ray overlaps t2-t3, shoot another ray, ray intersect it, shoot another one 
				//if (inter3 == 1) return -1;
				//if (inter3 == 2) return 2;

				return 0;
			}
			
			
		}
		return 0;
	}
	Rational dt = norm.dot(dir);// >0, dir is same direction with norm, pointing to +1 orientation
	if (dt.get_sign() > 0 && o1 >= 0) return 0;
	if (dt.get_sign() < 0 && o1 <= 0) return 0;
	Vector3r np = pt + dir;
	// if ray go across the plane, then get lpi and 3 orientations
	int inter = is_line_cut_triangle(pt, np, t1, t2, t3, halfopen, norm);
	if (inter == 0)return 0;
	if (inter == 1)return 1;
	if (inter == 2)return -1;//shoot on edge
	if (inter == 3) return 3;
	
	return 0;
}

// if a line (going across pt, pt+dir) intersects triangle
// triangle is not degenerated
// use ray_triangle_intersection twice, TODO modify this
// 0 not intersected, 1 intersected, 3 intersected t2-t3 edge 
//TODO make it to check open triangle, coplanar is not intersected
int line_triangle_intersection(
	const Vector3r& pt,
	const Vector3r& dir,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen) {

	int inter1 = ray_triangle_intersection(pt, dir, t1, t2, t3, halfopen);
	if (inter1 == 1) return 1;
	if (inter1 == 2) return 1;
	if (inter1 == 3) return 3;

	int inter2 = ray_triangle_intersection(pt, -dir, t1, t2, t3, halfopen);
	if (inter2 == 1) return 1;
	if (inter2 == 2) return 1;
	if (inter2 == 3) return 3;
	
	return 0;
}




} // namespace ccd
