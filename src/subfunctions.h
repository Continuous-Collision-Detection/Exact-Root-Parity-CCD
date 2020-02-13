#pragma once
#include "Utils.hpp"
#include <array>
#include <vector>
namespace ccd {
template <typename T, typename Y>
bool box_box_intersection(
    const T& min1, const T& max1, const Y& min2, const Y& max2)
{
    if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
        return 0;
    if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
        return 0;
    return 1;
}

bool is_seg_intersect_cube(
    const double& eps, const Vector3r& e0, const Vector3r& e1);
bool is_seg_intersect_cube_2d(
    const double eps, const Vector3r& e0, const Vector3r& e1, int axis);
void projected_cube_edges(
    const double eps,
    const int axis,
    Vector3r& e0,
    Vector3r& e1,
    Vector3r& e2,
    Vector3r& e3);
    // cube is centered at origin, and corner is (eps,eps,eps)
bool is_point_intersect_cube(
    const double eps, const Vector3r& p);
int seg_cut_plane(
    const Vector3r& seg0,
    const Vector3r& seg1,
    const Vector3r& t0,
    const Vector3r& t1,
    const Vector3r& t2);

//Causion: open triangle!
bool is_cube_intersect_triangle(
    const double eps,
    const Vector3r& t0,
    const Vector3r& t1,
    const Vector3r& t2);
bool is_seg_intersect_triangle(
    const Vector3r& s0,
    const Vector3r& s1,
    const Vector3r& t0,
    const Vector3r& t1,
    const Vector3r& t2);
//segment and triangle are coplanar, check intersection
bool is_coplanar_seg_intersect_triangle(
    const Vector3r& s0,
    const Vector3r& s1,
    const Vector3r& t0,
    const Vector3r& t1,
    const Vector3r& t2,
    const int axis);
void get__corners(const std::vector<Vector3d>& p, Vector3d& min, Vector3d& max);
template <typename T>
void get__corners(const std::array<T, 6>& p, T& min, T& max)
{

 min = p[0];
    max = p[0];
    for (int i = 0; i < 6; i++) {
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

class prism {
public:
    prism(
        const Vector3d& vertex_start,
        const Vector3d& face_vertex0_start,
        const Vector3d& face_vertex1_start,
        const Vector3d& face_vertex2_start,
        const Vector3d& vertex_end,
        const Vector3d& face_vertex0_end,
        const Vector3d& face_vertex1_end,
        const Vector3d& face_vertex2_end);

    template <typename T>
    bool is_prism_bbox_cut_bbox(const T& min, const T& max) const
    {
        Vector3r pmin, pmax;
        get__corners(p_vertices, pmin, pmax);
        return box_box_intersection(pmin, pmax, min, max);
    }
    std::array<std::array<int, 2>, 9> prism_edge_id;
    std::array<Vector3r, 6> p_vertices;

private:
    Vector3r vsr;
    Vector3r ver;
    Vector3r fs0r;
    Vector3r fs1r;
    Vector3r fs2r;
    Vector3r fe0r;
    Vector3r fe1r;
    Vector3r fe2r;

    

    Vector3r get_prism_corner(int u, int v, int t);
};

class cube {
public:
    cube(double eps);
    //std::array<Vector3d, 8> vd;
    std::array<Vector3r, 8> vr;
    std::array<std::array<int, 2>, 12> edgeid;
    std::array<std::array<int, 4>, 6> faceid;
    Vector3r bmin;
    Vector3r bmax;
};


} // namespace ccd