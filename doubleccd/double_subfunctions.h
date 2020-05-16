#pragma once

#include <array>
#include <vector>

#include <doubleCCD/double_Utils.hpp>

namespace doubleccd {

// x0 is the point, x1, x2, x3 is the triangle
struct vf_pair {
    Vector3d x0;
    Vector3d x1;
    Vector3d x2;
    Vector3d x3;
    Vector3d x0b;
    Vector3d x1b;
    Vector3d x2b;
    Vector3d x3b;

    vf_pair() {}

    vf_pair(
        const Vector3d& x0,
        const Vector3d& x1,
        const Vector3d& x2,
        const Vector3d& x3,
        const Vector3d& x0b,
        const Vector3d& x1b,
        const Vector3d& x2b,
        const Vector3d& x3b)
        : x0(x0)
        , x1(x1)
        , x2(x2)
        , x3(x3)
        , x0b(x0b)
        , x1b(x1b)
        , x2b(x2b)
        , x3b(x3b)
    {
    }
};
struct ee_pair {
    Vector3d a0;
    Vector3d a1;
    Vector3d b0;
    Vector3d b1;
    Vector3d a0b;
    Vector3d a1b;
    Vector3d b0b;
    Vector3d b1b;

    ee_pair() {}

    ee_pair(
        const Vector3d& a0,
        const Vector3d& a1,
        const Vector3d& b0,
        const Vector3d& b1,
        const Vector3d& a0b,
        const Vector3d& a1b,
        const Vector3d& b0b,
        const Vector3d& b1b)
        : a0(a0)
        , a1(a1)
        , b0(b0)
        , b1(b1)
        , a0b(a0b)
        , a1b(a1b)
        , b0b(b0b)
        , b1b(b1b)
    {
    }
};

double get_whole_mesh_shifted(
    const std::vector<vf_pair>& data1,
    const std::vector<ee_pair>& data2,
    std::vector<vf_pair>& shift_back1,
    std::vector<ee_pair>& shift_back2,
    Eigen::MatrixX3d& vertices);
// x0 is the point, x1, x2, x3 is the triangle
double get_whole_mesh_shifted(
    const std::vector<vf_pair>& data1,
    const std::vector<ee_pair>& data2,
    Eigen::MatrixX3d& vertices);

// Convenience function that just wrap get_whole_mesh_shifted(...)
double
shift_vertex_face(const vf_pair& input_vf_pair, vf_pair& shifted_vf_pair);
double shift_edge_edge(const ee_pair& input_ee_pair, ee_pair& shifted_ee_pair);

class cube {
public:
    cube(double eps);
    // std::array<Vector3d, 8> vd;
    std::array<Vector3d, 8> vr;
    std::array<std::array<int, 2>, 12> edgeid;
    std::array<std::array<int, 4>, 6> faceid;
    Vector3d bmin;
    Vector3d bmax;
    double epsilon;
};
void print_sub();

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
    const double& eps, const Vector3d& e0, const Vector3d& e1);
bool is_seg_intersect_cube_2d(
    const double eps, const Vector3d& e0, const Vector3d& e1, int axis);
void projected_cube_edges(
    const double eps, Vector2d& e0, Vector2d& e1, Vector2d& e2, Vector2d& e3);
// cube is centered at origin, and corner is (eps,eps,eps)
bool is_point_intersect_cube(const double eps, const Vector3d& p);

// Causion: open triangle!
bool is_cube_edges_intersect_triangle(
    const cube& cb, const Vector3d& t0, const Vector3d& t1, const Vector3d& t2);

// segment and triangle are coplanar, check intersection

void get_corners(const Eigen::MatrixX3d& p, Vector3d& min, Vector3d& max);
void get_tet_corners(const std::array<Vector3d, 4>& p, Vector3d& min, Vector3d& max);
void get_edge_coners(const Vector3d& e0, const Vector3d& e1, Vector3d &emin,Vector3d &emax);
template <typename T>
void get_bbd_corners(const std::array<T, 6>& p, T& min, T& max)
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
template <typename T>
void get_bbd_corners(const std::array<T, 8>& p, T& min, T& max)
{

    min = p[0];
    max = p[0];
    for (int i = 0; i < 8; i++) {
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
Vector3d get_prism_corner_double(
    const Vector3d& vertex_start,
    const Vector3d& face_vertex0_start,
    const Vector3d& face_vertex1_start,
    const Vector3d& face_vertex2_start,
    const Vector3d& vertex_end,
    const Vector3d& face_vertex0_end,
    const Vector3d& face_vertex1_end,
    const Vector3d& face_vertex2_end,
    int i);

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
        Vector3d pmin, pmax;
        get_bbd_corners(p_vertices, pmin, pmax);
        return box_box_intersection(pmin, pmax, min, max);
    }
    // 0 means up, 1 means bottom
    bool is_triangle_degenerated(const int up_or_bottom);
    std::array<std::array<int, 2>, 9> prism_edge_id;
    std::array<Vector3d, 6> p_vertices;

private:
    void get_prism_vertices(
        const Vector3d& x0,
        const Vector3d& x1,
        const Vector3d& x2,
        const Vector3d& x3,
        const Vector3d& x0b,
        const Vector3d& x1b,
        const Vector3d& x2b,
        const Vector3d& x3b,
        std::array<Vector3d, 6>& p_vertices);
};

class hex {
public:
    hex(const Vector3d& a0,
        const Vector3d& a1,
        const Vector3d& b0,
        const Vector3d& b1,
        const Vector3d& a0b,
        const Vector3d& a1b,
        const Vector3d& b0b,
        const Vector3d& b1b);
    //(1 - tar) * ((1 - tr) * a0rs_ + tr * a0re_) + tar * ((1 - tr) * a1rs_
    //+ tr
    //* a1re_) - ((1 - tbr) * ((1 - tr) * b0rs_ + tr * b0re_) + tbr * ((1 -
    // tr)
    //* b1rs_ + tr * b1re_));

    template <typename T>
    bool is_hex_bbox_cut_bbox(const T& min, const T& max) const
    {
        Vector3d pmin, pmax;
        get_bbd_corners(h_vertices, pmin, pmax);
        return box_box_intersection(pmin, pmax, min, max);
    }

    std::array<std::array<int, 2>, 12> hex_edge_id;
    std::array<Vector3d, 8> h_vertices;

private:
    void get_hex_vertices(
        const Vector3d& a0,
        const Vector3d& a1,
        const Vector3d& b0,
        const Vector3d& b1,
        const Vector3d& a0b,
        const Vector3d& a1b,
        const Vector3d& b0b,
        const Vector3d& b1b,
        std::array<Vector3d, 8>& h_vertices);
};
//pmin and pmax are the bounding box corners of bilinear bl
bool is_cube_intersect_tet_opposite_faces(
    const bilinear& bl,const Vector3d &pmin, const Vector3d &pmax,
    const cube& cube,
    std::array<bool, 8>& vin,
    bool& cube_inter_tet);
int bilinear_degeneration(const bilinear& bl);
// int get_triangle_project_axis(
//    const Vector3r& t0, const Vector3r& t1, const Vector3r& t2);
bool is_cube_edge_intersect_bilinear(
    bilinear& bl, const cube& cb, const std::array<bool, 8>& pin);
// phisign gives the pair we want to check
bool line_shoot_same_pair_tet(
    const Vector3d& p0, const Vector3d& p1, const int phisign, bilinear& bl);
bool rootfinder(
    const bilinear& bl,
    const Vector3d& p0d,
    const Vector3d& p1d,
    const bool p0in,
    const bool p1in,
    const int pairid);
bool cube_discrete_bilinear_intersection(
    const cube& cb, const bilinear& bl, int n);
bool is_point_inside_tet(const bilinear& bl, const Vector3d& p);

} // namespace doubleccd
