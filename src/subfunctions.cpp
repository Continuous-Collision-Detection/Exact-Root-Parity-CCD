#include "subfunctions.h"
namespace ccd {

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
//bool is_seg_intersect_cube(
//    const cube& cb, const Vector3r& e0, const Vector3r& e1)
//{
//	//if (e0[0] == e1[0] && e0[1] == e1[1] && e0[2] == e1[2])
//	//{
// //       for (int i = 0; i < 3; i++) {
// //               if (e0[i] < cb.bmin[i] || e0[i] > cb.bmax[i])
// //                   return false;
//	//	}
// //       return true;
//	//}
// //   Vector3r t0, t1, t2;
// //   int o1[6], o2[6], inter[6];
// //   bool inter_palne = false;
// //   for (int i = 0; i < 6; i++) {
// //       t0 = cb.vr[cb.faceid[i][0]];
// //       t1 = cb.vr[cb.faceid[i][1]];
// //       t2 = cb.vr[cb.faceid[i][2]];
// //       inter[i] = seg_cut_plane(e0, e1, t0, t1, t2);
// //       if (inter[i] == INTERSECTED) {
// //           if(segment_triangle_inter(e0, e1, t0, t1, t2)==1)
// //               return true;
// //           if (segment_triangle_inter(
// //                   e0, e1, cb.vr[cb.faceid[i][0]], cb.vr[cb.faceid[i][2]],
// //                   cb.vr[cb.faceid[i][3]])
// //               == 1)
// //               return true;
// //           inter_palne = true; 
//	//	}    
//	//	//
//	//}// the rest cases: totally inside, totally outside, totally in one (closed) rectangle
// //   if (inter_palne)
// //       return false;
//	//
//
//
//
//	return 0;
//
//}
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
    p_vertices[2] = get_prism_corner(0, 0, 1);
    p_vertices[3] = get_prism_corner(1, 0, 0);
    p_vertices[4] = get_prism_corner(1, 1, 0);
    p_vertices[5] = get_prism_corner(1, 0, 1);//these are the 6 vertices of the prism,right hand law
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
Vector3r prism::get_prism_corner(int u, int v, int t)
{
    const Rational ur(u);
    const Rational vr(v);
    const Rational tr(t);
    return (1 - tr) * vsr + tr * ver
        - ((1 - tr) * fs0r + t * fe0r) * (1 - ur - vr)
        - ((1 - tr) * fs1r + t * fe1r) * ur - ((1 - tr) * fs2r + t * fe2r) * vr;
}
cube::cube(double eps)
{
    vr[0] = Vector3r(-eps, -eps, eps), vr[1] = Vector3r(eps, -eps, eps),
    vr[2] = Vector3r(eps, eps, eps), vr[3] = Vector3r(-eps, eps, eps),
    vr[4] = Vector3r(-eps, -eps, -eps), vr[5] = Vector3r(eps, -eps, -eps),
    vr[6] = Vector3r(eps, eps, -eps), vr[7] = Vector3r(-eps, eps, -eps);
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
    faceid[5] = { { 0, 3, 7, 4 } };//orientation out
    bmax = vr[2];
    bmin = vr[4];
}

} // namespace ccd