#include <predicates/indirect_predicates.h>
#include <doubleCCD/double_subfunctions.h>
//#include <exact_subtraction.hpp>

//#include <ray_parity.h>
namespace ccd {

// cube
cube::cube(double eps)
{
    vr[0] = Vector3d(-eps, -eps, eps), vr[1] = Vector3d(eps, -eps, eps),
    vr[2] = Vector3d(eps, eps, eps), vr[3] = Vector3d(-eps, eps, eps),
    vr[4] = Vector3d(-eps, -eps, -eps), vr[5] = Vector3d(eps, -eps, -eps),
    vr[6] = Vector3d(eps, eps, -eps), vr[7] = Vector3d(-eps, eps, -eps);
    edgeid[0] = { { 0, 1 } };
    edgeid[1] = { { 1, 2 } };
    edgeid[2] = { { 2, 3 } };
    edgeid[3] = { { 3, 0 } };
    edgeid[4] = { { 4, 5 } };
    edgeid[5] = { { 5, 6 } };
    edgeid[6] = { { 6, 7 } };
    edgeid[7] = { { 7, 4 } };
    edgeid[8] = { { 0, 4 } };
    edgeid[9] = { { 1, 5 } };
    edgeid[10] = { { 2, 6 } };
    edgeid[11] = { { 3, 7 } };
    faceid[0] = { { 0, 1, 2, 3 } };
    faceid[1] = { { 4, 7, 6, 5 } };
    faceid[2] = { { 0, 4, 5, 1 } };
    faceid[3] = { { 1, 5, 6, 2 } };
    faceid[4] = { { 3, 2, 6, 7 } };
    faceid[5] = { { 0, 3, 7, 4 } }; // orientation out
    bmax = vr[2];
    bmin = vr[4];
    epsilon = eps;
}

// get aabb corners
void get__corners(const std::vector<Vector3d>& p, Vector3d& min, Vector3d& max)
{
    min = p[0];
    max = p[0];
    for (int i = 0; i < p.size(); i++) {
        if (min[0] > p[i][0])
            min[0] = p[i][0];
        if (min[1] > p[i][1])
            min[1] = p[i][1];
        if (min[2] > p[i][2])
            min[2] = p[i][2];

        if (max[0] < p[i][0])
            max[0] = p[i][0];
        if (max[1] < p[i][1])
            max[1] = p[i][1];
        if (max[2] < p[i][2])
            max[2] = p[i][2];
    }
}
/*
std::array<Vector3d, 6> get_prism_vertices_double(
    const Vector3d& x0,
    const Vector3d& x1,
    const Vector3d& x2,
    const Vector3d& x3,
    const Vector3d& x0b,
    const Vector3d& x1b,
    const Vector3d& x2b,
    const Vector3d& x3b,
    double& k,
    bool& correct,
    double& maxerror)
{
    std::vector<std::pair<double, double>> sub;
    sub.reserve(18);
    std::pair<double, double> temp;
    correct = true;
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x1[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x3[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x2[i];
        sub.push_back(temp);
    }
    //
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x1b[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x3b[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x2b[i];
        sub.push_back(temp);
    }
    std::vector<std::pair<double, double>> sub_record
        = sub; // record the origin data
    k = displaceSubtractions_double(sub);
    std::array<Vector3d, 6> result;
    int ct = 0;
    maxerror = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = sub[ct].first - sub[ct].second;

            if (Rational(sub[ct].first) - Rational(sub[ct].second)
                != result[i][j]) {
                correct = false;
            }
            double mns = sub_record[ct].first - sub_record[ct].second;
            Rational m = Rational(result[i][j]) - Rational(mns);
            double md = m.to_double();
            if (fabs(md) > maxerror) {
                maxerror = fabs(md);
            }

            ct++;
        }
    }
    return result;
}
*/
Vector3d get_prism_corner_double(
    const Vector3d& vertex_start,       // x0
    const Vector3d& face_vertex0_start, // x1
    const Vector3d& face_vertex1_start, // x2
    const Vector3d& face_vertex2_start, // x3
    const Vector3d& vertex_end,
    const Vector3d& face_vertex0_end,
    const Vector3d& face_vertex1_end,
    const Vector3d& face_vertex2_end,
    int i)
{
    Vector3d x0 = vertex_start, x1 = face_vertex0_start,
             x2 = face_vertex1_start, x3 = face_vertex2_start, x0b = vertex_end,
             x1b = face_vertex0_end, x2b = face_vertex1_end,
             x3b = face_vertex2_end;
    if (i == 0)
        return x0 - x1;
    if (i == 1)
        return x0 - x3;
    if (i == 2)
        return x0 - x2;
    if (i == 3)
        return x0b - x1b;
    if (i == 4)
        return x0b - x3b;
    if (i == 5)
        return x0b - x2b;

    else
        return Vector3d();
}

bool is_seg_intersect_cube(
    const double& eps, const Vector3d& e0, const Vector3d& e1)
{
    if (is_point_intersect_cube(eps, e0))
        return true;
    if (is_point_intersect_cube(eps, e1))
        return true;
    if (same_point(e0,e1))
        return false; // degenerate case: the segment is degenerated as a point
    // if intersected, must be coplanar with the edge, or intersect edge or face
    if (is_seg_intersect_cube_2d(eps, e0, e1, 0)
        && is_seg_intersect_cube_2d(eps, e0, e1, 1)
        && is_seg_intersect_cube_2d(eps, e0, e1, 2))
        return true;
    return false;
}
// check if a 2d segment intersects 2d cube
bool is_seg_intersect_cube_2d(
    const double eps, const Vector3d& e0, const Vector3d& e1, int axis)
{
	Vector2d p0, p1, p2, p3,
		e0p, e1p;// e0 and e1 projected to 2d
    projected_cube_edges(eps, p0, p1, p2, p3);// TODO move this out

    const int i1 = (axis + 1) % 3;
    const int i2 = (axis + 2) % 3;

    if (e0[i1] <= eps && e0[i1] >= -eps && e0[i2] <= eps && e0[i2] >= -eps)
        return true;
    if (e1[i1] <= eps && e1[i1] >= -eps && e1[i2] <= eps && e1[i2] >= -eps)
        return true;
	e0p = Vector2d(e0[i1], e0[i2]);
	e1p = Vector2d(e1[i1], e1[i2]);
    if (segment_segment_intersection_2d(e0p, e1p, p0, p1))
        return true; // check if segments has intersection, or if cube points
                     // p0, p1 on e0-e1
    if (segment_segment_intersection_2d(e0p, e1p, p1, p2))
        return true;
    if (segment_segment_intersection_2d(e0p, e1p, p2, p3))
        return true;
    if (segment_segment_intersection_2d(e0p, e1p, p3, p0))
        return true;

    return false;
}
void projected_cube_edges(
    const double eps,
    Vector2d& e0,
    Vector2d& e1,
    Vector2d& e2,
    Vector2d& e3)
{
	const int i1 = 0;
	const int i2 = 1;

    e0[i1] = -eps;
    e0[i2] = eps;

    e1[i1] = eps;
    e1[i2] = eps;

    e2[i1] = eps;
    e2[i2] = -eps;

    e3[i1] = -eps;
    e3[i2] = -eps;
}
bool is_point_intersect_cube(const double eps, const Vector3d& p)
{
    if (p[0] <= eps && p[0] >= -eps) {
        if (p[1] <= eps && p[1] >= -eps) {
            if (p[2] <= eps && p[2] >= -eps) {
                return true;
            }
        }
    }
    return false;
}

//int get_triangle_project_axis(
//    const Vector3r& t0, const Vector3r& t1, const Vector3r& t2)
//{
//    Vector3r normal = cross(t0 - t1, t0 - t2);
//    if (normal[0] == 0 && normal[1] == 0 && normal[2] == 0) {
//        return 3; // if triangle degenerated as a segment or point, then no
//                  // intersection, because before here we already check that
//    }
//
//    if (normal[0] * normal[0].get_sign() >= normal[1] * normal[1].get_sign())
//        if (normal[0] * normal[0].get_sign()
//            >= normal[2] * normal[2].get_sign())
//            return 0;
//    if (normal[1] * normal[1].get_sign() >= normal[0] * normal[0].get_sign())
//        if (normal[1] * normal[1].get_sign()
//            >= normal[2] * normal[2].get_sign())
//            return 1;
//    if (normal[2] * normal[2].get_sign() >= normal[1] * normal[1].get_sign())
//        if (normal[2] * normal[2].get_sign()
//            >= normal[0] * normal[0].get_sign())
//            return 2;
//}
bool is_cube_edges_intersect_triangle(
    const ccd::cube& cb,
    const Vector3d& t0,
    const Vector3d& t1,
    const Vector3d& t2)
{
    // the vertices of triangle are checked before going here, the edges are
    // also checked. so, only need to check if cube edge has intersection with
    // the open triangle.

    //int axis = get_triangle_project_axis(t0, t1, t2);
    //if (axis == 3)
    //    return false; // if triangle degenerated as a segment or point, then no
    //                  // intersection, because before here we already check that
    Vector3d s0, s1;
    for (int i = 0; i < 12; i++) {
        s0 = cb.vr[cb.edgeid[i][0]];
        s1 = cb.vr[cb.edgeid[i][1]];
		if (segment_triangle_intersection(s0, s1, t0, t1, t2, false) > 0)// return 0,1,2
            return true;
    }
    return false;
}

prism::prism(
    const Vector3d& vs,
    const Vector3d& fs0,
    const Vector3d& fs1,
    const Vector3d& fs2,
    const Vector3d& ve,
    const Vector3d& fe0,
    const Vector3d& fe1,
    const Vector3d& fe2)
{
    for (int i = 0; i < 3; i++) {
        vsr[i] = vs[i];
        ver[i] = ve[i];
        fs0r[i] = fs0[i];
        fs1r[i] = fs1[i];
        fs2r[i] = fs2[i];
        fe0r[i] = fe0[i];
        fe1r[i] = fe1[i];
        fe2r[i] = fe2[i];
    }
    p_vertices[0] = get_prism_corner(0, 0, 0);
    p_vertices[1] = get_prism_corner(0, 1, 0);
    p_vertices[2] = get_prism_corner(1, 0, 0);
    p_vertices[3] = get_prism_corner(0, 0, 1);
    p_vertices[4] = get_prism_corner(0, 1, 1);
    p_vertices[5] = get_prism_corner(1, 0, 1);
    // these are the 6 vertices of the prism,right hand law
    std::array<int, 2> eid;

    eid[0] = 0;
    eid[1] = 1;
    prism_edge_id[0] = eid;
    eid[0] = 1;
    eid[1] = 2;
    prism_edge_id[1] = eid;
    eid[0] = 2;
    eid[1] = 0;
    prism_edge_id[2] = eid;

    eid[0] = 3;
    eid[1] = 4;
    prism_edge_id[3] = eid;
    eid[0] = 4;
    eid[1] = 5;
    prism_edge_id[4] = eid;
    eid[0] = 5;
    eid[1] = 3;
    prism_edge_id[5] = eid;

    eid[0] = 0;
    eid[1] = 3;
    prism_edge_id[6] = eid;
    eid[0] = 1;
    eid[1] = 4;
    prism_edge_id[7] = eid;
    eid[0] = 2;
    eid[1] = 5;
    prism_edge_id[8] = eid;

}
bool prism::is_triangle_degenerated(const int up_or_bottom)
{
	int pid = up_or_bottom == 0 ? 0 : 3;
	const auto to_2d = [](const Vector3d &p, int t)
	{
		return Vector2d(p[(t + 1) % 3], p[(t + 2) % 3]);
	};
	double r = ((p_vertices[pid] - p_vertices[pid + 1]).cross(p_vertices[pid] - p_vertices[pid + 2])).norm();
	if (fabs(r) > 1e-8) return false;
	int ori;
	std::array<Vector2d, 3> p;
	for (int j = 0; j < 3; j++)
	{

		p[0] = to_2d(p_vertices[pid], j);
		p[1] = to_2d(p_vertices[pid+1], j);
		p[2] = to_2d(p_vertices[pid+2], j);

		ori = orient2d(p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1]);
		if (ori != 0)
		{
			return false;
		}
	}
	return true;
}
Vector3d prism::get_prism_corner(int u, int v, int t)
{

    return (1 - t) * vsr + t * ver
        - ((1 - t) * fs0r + t * fe0r) * (1 - u - v)
        - ((1 - t) * fs1r + t * fe1r) * u - ((1 - t) * fs2r + t * fe2r) * v;
}


// the facets of the tet are all oriented to outside. check if p is inside of
// OPEN tet
bool is_point_inside_tet(const bilinear& bl, const Vector3d& p)
{

    for (int i = 0; i < 4; i++) { // facets.size()==4
		Vector3d pt1= bl.v[bl.facets[i][0]], pt2= bl.v[bl.facets[i][1]], pt3= bl.v[bl.facets[i][2]];
		if (orient_3d(
			p, pt1, pt2, pt3) >= 0) {
            return false;
        }
    }
    return true; // all the orientations are -1, then point inside
}


// we already know the bilinear is degenerated, next check which kind
int bilinear_degeneration(const bilinear& bl)
{
	bool dege1 = is_triangle_degenerated(bl.v[0], bl.v[1], bl.v[2]);
	if (dege1) return BI_DEGE_PLANE;
	bool dege2 = is_triangle_degenerated(bl.v[0], bl.v[2], bl.v[3]);
	if (dege2) return BI_DEGE_PLANE;
	Vector3d np = Vector3d::Random();
	int ori = orient_3d(np, bl.v[0], bl.v[1], bl.v[2]);
	while (ori == 0) {// if coplanar, random
		np = Vector3d::Random();
		ori = orient_3d(np, bl.v[0], bl.v[1], bl.v[2]);
	}
	int ori1= orient_3d(np, bl.v[0], bl.v[2], bl.v[3]);
	if (ori*ori1 < 0) {
		return BI_DEGE_XOR_02;
	}
    // CAUTION: no need BI_DEGE_XOR_13 anymore, because it is the same with BI_DEGE_XOR_02
    return BI_DEGE_PLANE;
}

bool is_cube_intersect_degenerated_bilinear(
    const bilinear& bl, const cube& cube)
{
    int dege = bilinear_degeneration(bl);
    //int axis;
    bool res;
    if (dege == BI_DEGE_PLANE) {

        if (is_cube_edges_intersect_triangle(cube, bl.v[0], bl.v[1], bl.v[3]))
            return true;
        if (is_cube_edges_intersect_triangle(cube, bl.v[3], bl.v[1], bl.v[2]))
            return true;
        return false;
    } else {
		//int  axis = get_triangle_project_axis(bl.v[0], bl.v[1], bl.v[3]);
  //      if (axis == 3)
  //          axis = get_triangle_project_axis(bl.v[3], bl.v[1], bl.v[2]);
  //      if (axis == 3)
  //          return false; // both 2 triangles are all degenerated as a segment
        if (dege == BI_DEGE_XOR_02) { // triangle 0-1-2 and 0-2-3
            for (int i = 0; i < 12; i++) {
                res = int_seg_XOR(
					segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[1], bl.v[2], bl.v[0], true),// CAUTION: need to be careful for the order here
					segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[3], bl.v[2], bl.v[0], true));
                if (res == true)
                    return true;
            }
            return false;
        }

     //   if (dege == BI_DEGE_XOR_13) { // triangle 0-1-3 and 3-1-2
     //       for (int i = 0; i < 12; i++) {
     //           res = int_seg_XOR(
					//segment_triangle_intersection(
     //                   cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
     //                   bl.v[0], bl.v[1], bl.v[3], true),
					//segment_triangle_intersection(
     //                   cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
     //                   bl.v[2], bl.v[1], bl.v[3], true));
     //           if (res == true)
     //               return true;
     //       }
     //       return false;
     //   }
    }
    std::cout << "!! THIS CANNOT HAPPEN" << std::endl;
    return false;
}
// phisign gives the pair we want to check
bool line_shoot_same_pair_tet(
    const Vector3d& p0, const Vector3d& p1, const int phisign, bilinear& bl)
{
    int fid;
    if (bl.phi_f[0] == 2)
        get_tet_phi(bl);

    if (phisign > 0) {
        if (bl.phi_f[0] > 0)
            fid = 0;
        else
            fid = 2;
    } else {
        if (bl.phi_f[1] > 0)
            fid = 0;
        else
            fid = 2;
    }// get which pair of facets to check
	//if line is parallel to the triangle, return false
	int inter0 = is_line_cut_triangle(p0, p1, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
		bl.v[bl.facets[fid][2]], false);

	int inter1 = is_line_cut_triangle(p0, p1, bl.v[bl.facets[fid + 1][0]], bl.v[bl.facets[fid + 1][1]],
		bl.v[bl.facets[fid + 1][2]], false);

	if (inter0 == 1 && inter1 == 1)
        return true;
    return false;
}

Rational quadratic_function_value(
    const Rational& a, const Rational& b, const Rational& c, const Rational& t)
{

    if (t.get_sign() == 0) {
        return c;
    } else {
        return a * t * t + b * t + c;
    }
}
// f(t)=at^2+bt+c, when t is [t0, t1]
bool quadratic_function_rootfinder(
    const Rational& a,
    const Rational& b,
    const Rational& c,
    const Rational t0,
    const Rational t1)// t0 t1 no matter which is bigger
{
    Rational ft0, ft1;
    ft0 = quadratic_function_value(a, b, c, t0);
    ft1 = quadratic_function_value(a, b, c, t1);
    if (ft0.get_sign() == 0 || ft1.get_sign() == 0)
        return true;
    if (ft0.get_sign() != ft1.get_sign())
        return true;
    // the following are cases which signs of endpoints are same
    if (a.get_sign() == 0)
        return false;
    Rational t = -b / (2 * a);
    if (t < t1 && t > t0) {
        Rational ft
            = (4 * a * c - b * b); // actually it should be /4a. we only
                                   // care about the signs so doesnt matter
        int ftsign = a.get_sign() > 0 ? ft.get_sign() : (-1) * ft.get_sign();
        if (ft0.get_sign() != ftsign)
            return true;
    }
    return false;
}
// v0 is one point, dir is the direction,
// x0 - x3 are four points defined the shape
void get_quadratic_function(
    const Rational& v00,
    const Rational& v01,
    const Rational& v02,
    const Rational& dir0,
    const Rational& dir1,
    const Rational& dir2,
    const Vector3d x0d,
    const Vector3d x1d,
    const Vector3d x2d,
    const Vector3d x3d,
    Rational& a,
    Rational& b,
    Rational& c)
{
	Vector3r x0(x0d[0], x0d[1], x0d[2]), x1(x1d[0], x1d[1], x1d[2]), x2(x2d[0], x2d[1], x2d[2]),
		x3(x3d[0], x3d[1], x3d[2]);
    Rational x00 = x0[0], x01 = x0[1], x02 = x0[2];
    Rational x10 = x1[0], x11 = x1[1], x12 = x1[2];
    Rational x20 = x2[0], x21 = x2[1], x22 = x2[2];
    Rational x30 = x3[0], x31 = x3[1], x32 = x3[2];

    Rational x101 = x11 - x01, x100 = x10 - x00, x102 = x12 - x02,
             x201 = x21 - x01, x200 = x20 - x00, x300 = x30 - x00,
             x301 = x31 - x01, x310 = x30 - x10, x211 = x21 - x11,
             x311 = x31 - x11, x210 = x20 - x10, x202 = x22 - x02,
             x302 = x32 - x02, x212 = x22 - x12, x312 = x32 - x12;
    a = (dir0 * (x101 * x202 - x102 * x201)
         + dir1 * (-x100 * x202 + x102 * x200)
         + dir2 * (x100 * x201 - x101 * x200))
            * (dir0 * (-x201 * x302 + x202 * x301)
               + dir1 * (x200 * x302 - x202 * x300)
               + dir2 * (-x200 * x301 + x201 * x300))
        - (dir0 * (-x211 * x312 + x212 * x311)
           + dir1 * (x210 * x312 - x212 * x310)
           + dir2 * (-x210 * x311 + x211 * x310))
            * (dir0 * (x101 * x302 - x102 * x301)
               + dir1 * (-x100 * x302 + x102 * x300)
               + dir2 * (x100 * x301 - x101 * x300));
    b = ((v00 - x00) * (x101 * x202 - x102 * x201)
         + (v01 - x01) * (-x100 * x202 + x102 * x200)
         + (v02 - x02) * (x100 * x201 - x101 * x200))
            * (dir0 * (-x201 * x302 + x202 * x301)
               + dir1 * (x200 * x302 - x202 * x300)
               + dir2 * (-x200 * x301 + x201 * x300))
        + (dir0 * (x101 * x202 - x102 * x201)
           + dir1 * (-x100 * x202 + x102 * x200)
           + dir2 * (x100 * x201 - x101 * x200))
            * ((v00 - x00) * (-x201 * x302 + x202 * x301)
               + (v01 - x01) * (x200 * x302 - x202 * x300)
               + (v02 - x02) * (-x200 * x301 + x201 * x300))
        - ((v00 - x10) * (-x211 * x312 + x212 * x311)
           + (v01 - x11) * (x210 * x312 - x212 * x310)
           + (v02 - x12) * (-x210 * x311 + x211 * x310))
            * (dir0 * (x101 * x302 - x102 * x301)
               + dir1 * (-x100 * x302 + x102 * x300)
               + dir2 * (x100 * x301 - x101 * x300))
        - (dir0 * (-x211 * x312 + x212 * x311)
           + dir1 * (x210 * x312 - x212 * x310)
           + dir2 * (-x210 * x311 + x211 * x310))
            * ((v00 - x00) * (x101 * x302 - x102 * x301)
               + (v01 - x01) * (-x100 * x302 + x102 * x300)
               + (v02 - x02) * (x100 * x301 - x101 * x300));
    c = ((v00 - x00) * (x101 * x202 - x102 * x201)
         + (v01 - x01) * (-x100 * x202 + x102 * x200)
         + (v02 - x02) * (x100 * x201 - x101 * x200))
            * ((v00 - x00) * (-x201 * x302 + x202 * x301)
               + (v01 - x01) * (x200 * x302 - x202 * x300)
               + (v02 - x02) * (-x200 * x301 + x201 * x300))
        - ((v00 - x10) * (-x211 * x312 + x212 * x311)
           + (v01 - x11) * (x210 * x312 - x212 * x310)
           + (v02 - x12) * (-x210 * x311 + x211 * x310))
            * ((v00 - x00) * (x101 * x302 - x102 * x301)
               + (v01 - x01) * (-x100 * x302 + x102 * x300)
               + (v02 - x02) * (x100 * x301 - x101 * x300));
}
bool get_function_find_root(
    const bilinear& bl, const Vector3r& p0, const Vector3r& p1, const Rational &t0, const Rational &t1)
{
	if (t0 > 1 || t0 < 0) std::cout << "t is not right: exceed the limit" << std::endl;
	if (t1 > 1 || t1 < 0) std::cout << "t is not right: exceed the limit" << std::endl;
    Rational a, b, c;
	//Vector3r p0(p0d[0], p0d[1], p0d[2]), p1(p1d[0], p1d[1], p1d[2]);
    Vector3r dir = p1 - p0;
    get_quadratic_function(
        Rational(p0[0]), Rational(p0[1]), Rational(p0[2]), Rational(dir[0]), Rational(dir[1]), Rational(dir[2])
		, bl.v[0], bl.v[1], bl.v[2],
        bl.v[3], a, b, c);
    return quadratic_function_rootfinder(a, b, c, t0, t1);
}
bool rootfinder(
    const bilinear& bl,
    const Vector3d& p0d,
    const Vector3d& p1d,
    const bool p0in,
    const bool p1in,
    const int pairid)
{
	Vector3r p0(p0d[0], p0d[1], p0d[2]), p1(p1d[0], p1d[1], p1d[2]);
	//std::cout << "we use root finder here" << std::endl;
    if (p0in && p1in) {
		// t0=0, t1=1
        return get_function_find_root(bl, p0, p1, Rational(0), Rational(1));
    }
    int fid = 2 * pairid;// it should be 0 or 2
    Rational t;
    if (p0in) {
        bool res1=seg_triangle_inter_return_t(
            p0d, p1d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
            bl.v[bl.facets[fid][2]], t);
        if (!res1) {
            res1=seg_triangle_inter_return_t(
                p0d, p1d, bl.v[bl.facets[fid+1][0]], bl.v[bl.facets[fid+1][1]],
                bl.v[bl.facets[fid+1][2]], t);
        }
		if (!res1)
			return false; // means not really intersected
        // we got t
        return get_function_find_root(bl, p0, p1, 0, t);
    }
    if (p1in) { // change the order of input to get t just because we want
                // domain to be [0, t]
        bool res1 = seg_triangle_inter_return_t(//TODO here get n1, d1, n2, d2
            p1d, p0d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
            bl.v[bl.facets[fid][2]], t);
        if (!res1) {
			res1 = seg_triangle_inter_return_t(
                p1d, p0d, bl.v[bl.facets[fid + 1][0]],
                bl.v[bl.facets[fid + 1][1]], bl.v[bl.facets[fid + 1][2]], t);
        }
		if (!res1)
			return false; // means not really intersected
        // we got t
        return get_function_find_root(bl, p1, p0, 0, t);
    }

    Rational t1;
	bool res1 = seg_triangle_inter_return_t(
        p0d, p1d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
        bl.v[bl.facets[fid][2]], t);
	if (!res1)
		return false; // means not really intersected
	bool res2 = seg_triangle_inter_return_t(
        p0d, p1d, bl.v[bl.facets[fid + 1][0]], bl.v[bl.facets[fid + 1][1]],
        bl.v[bl.facets[fid + 1][2]], t1);
	if (!res2)
		return false; // means not really intersected
    // we got t, t1
    return get_function_find_root(bl, p0, p1, t, t1);
}
// segment intersect two opposite faces not included. compare phi, then use root
// finder
bool is_seg_intersect_not_degenerated_bilinear(
    bilinear& bl,
    const Vector3d& p0,
    const Vector3d& p1,
    const bool pin0,
    const bool pin1)
{
	//first compare phi, if phis are different, intersected;
	//then check if the line intersect two opposite facets of bilinear, if so, use
	//rootfinder, else, not intersected
    if (pin0 && pin1) { // two points are all inside

        Rational phi0 = phi(p0, bl.v);
        Rational phi1 = phi(p1, bl.v);
		if (phi0 == 0 || phi1 == 0 || phi0.get_sign() != phi1.get_sign())
            return true;
        if (line_shoot_same_pair_tet(p0, p1, phi1.get_sign(), bl)) {
            if (phi1.get_sign() == bl.phi_f[0])
                return rootfinder(bl, p0, p1, pin0, pin1, 0);
            else
                return rootfinder(bl, p0, p1, pin0, pin1, 1);
		}

        else
            return false; // if the phis are the same, and shoot same pair, need
                          // to use rootfinder
    }
    if (pin0) {
        Rational phi0 = phi(p0, bl.v);
        if (phi0 == 0)
            return true;
        int hitpair = -1;

        for (int i = 0; i < 4; i++) {
            if (segment_triangle_intersection(// 0,1,2,3. 1,2,3 are all intersected
                    p0, p1, bl.v[bl.facets[i][0]], bl.v[bl.facets[i][1]],
				bl.v[bl.facets[i][2]], false)
                > 0) {
                if (i < 2)
                    hitpair = 0;
                else
                    hitpair = 1;
                break;
            }
        }
        if (bl.phi_f[0] == 2)
            get_tet_phi(bl);
        if (hitpair == -1)
            return false;
        if (phi0.get_sign()
            != bl.phi_f[hitpair]) { // if different, intersected, if same,
                                    // extend; if shoot same, rootfinder, if
                                    // shoot diff, false
            return true;
        }

        else {
            if (line_shoot_same_pair_tet(p0, p1, phi0.get_sign(), bl)) {
                if (phi0.get_sign() == bl.phi_f[0])
                    return rootfinder(bl, p0, p1, pin0, pin1, 0);
                else
                    return rootfinder(bl, p0, p1, pin0, pin1, 1);
			}

            else
                return false; // if the phis are the same, and shoot same pair,
                              // need to use rootfinder
        }
    }
    if (pin1) {
        Rational phi1 = phi(p1, bl.v);
        if (phi1 == 0)
            return true;
        int hitpair = -1;
        for (int i = 0; i < 4; i++) {
            if (segment_triangle_intersection(// 0,1,2,3. 1,2,3 are all intersected
                    p0, p1, bl.v[bl.facets[i][0]], bl.v[bl.facets[i][1]],
				bl.v[bl.facets[i][2]], false)
                > 0) {
                if (i < 2)
                    hitpair = 0;
                else
                    hitpair = 1;
                break;
            }
        }
        if (bl.phi_f[0] == 2)
            get_tet_phi(bl);
        if (hitpair == -1)
            return false; // parallel , should be impossible
        if (phi1.get_sign() != bl.phi_f[hitpair]) {
            return true;
        }

        else {
            if (line_shoot_same_pair_tet(p0, p1, phi1.get_sign(), bl)) {
                if (phi1.get_sign() == bl.phi_f[0])
                    return rootfinder(bl, p0, p1, pin0, pin1, 0);
                else
                    return rootfinder(bl, p0, p1, pin0, pin1, 1);
            }
            else
                return false; // if the phis are the same, and shoot same pair,
                              // need to use rootfinder
        }
    }
    if (!pin0
        && !pin1) { // not intersect tet (false), or intersect same side(root
                    // finder) or intersect diff side(checked before)

        if (line_shoot_same_pair_tet(p0, p1, 1, bl)) {
            if (1 == bl.phi_f[0])
                return rootfinder(bl, p0, p1, pin0, pin1, 0);
            else
                return rootfinder(bl, p0, p1, pin0, pin1, 1);
		}

        else if (line_shoot_same_pair_tet(p0, p1, -1, bl)) {
            if (-1 == bl.phi_f[0])
                return rootfinder(bl, p0, p1, pin0, pin1, 0);
            else
                return rootfinder(bl, p0, p1, pin0, pin1, 1);
			}
        return false; // if the phis are the same, and shoot same pair,
                      // need to use rootfinder
    }
    std::cout
        << " it cannot happen here in is_seg_intersect_not_degenerated_bilinear"
        << std::endl;
    return false;
}

bool is_cube_edge_intersect_bilinear(
    bilinear& bl, const cube& cb, const std::array<bool, 8>&pin)
{
	if (bl.is_degenerated) return false;// we already checked degenerated cases
    for (int i = 0; i < 12; i++) {
        if (is_seg_intersect_not_degenerated_bilinear(
                bl, cb.vr[cb.edgeid[i][0]], cb.vr[cb.edgeid[i][1]],
                pin[cb.edgeid[i][0]], pin[cb.edgeid[i][1]]))
            return true;
	}
    return false;
}
    // vin is true, this vertex has intersection with open tet
// if tet is degenerated, just tell us if cube is intersected with the shape
bool is_cube_intersect_tet_opposite_faces(
    const bilinear& bl,
    const cube& cube,
    std::array<bool, 8>& vin,
    bool &cube_inter_tet)
{
	cube_inter_tet = false;
    if (!bl.is_degenerated) {
        for (int i = 0; i < 8; i++) {
			vin[i] = false;
            if (is_point_inside_tet(bl, cube.vr[i])) {
				cube_inter_tet = true;
                vin[i] = true;
            }
        }
    } else {

        return is_cube_intersect_degenerated_bilinear(bl, cube);
    }

    bool side1 = false;
    bool side2 = false;

    for (int i = 0; i < 12; i++) {
        if (vin[cube.edgeid[i][0]]
            && vin[cube.edgeid[i][1]]) { // if two vertices are all inside, it can not cut any edge
			cube_inter_tet = true;
            continue;
        }
        for (int j = 0; j < 4; j++) {
            if (segment_triangle_intersection(
                    cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                    bl.v[bl.facets[j][0]], bl.v[bl.facets[j][1]],
				bl.v[bl.facets[j][2]], false)
                > 0) {
				cube_inter_tet = true;
                if (j == 0 || j == 1)
                    side1 = true;
                if (j == 2 || j == 3)
                    side2 = true;
                if (side1 && side2)
                    return true;
            }
        }
    }
    if (side1 && side2)
        return true;
    return false;
}

} // namespace ccd