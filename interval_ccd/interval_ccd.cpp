#include <iostream>

#include <igl/Timer.h>

#include <interval_ccd/interval_ccd.hpp>
#include <interval_ccd/interval_root_finder.hpp>

namespace intervalccd {

double time0 = 0, time1 = 0, time2 = 0;
static double CCD_LENGTH_TOL = 1e-6;
double tol0 = 0, tol1 = 0, tol2 = 0, tol0n = 0, tol1n = 0, tol2n = 0;
double c00, c01, c10, c11, c20, c21;

double compute_face_vertex_tolerance_1d(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end)
{
    double dl = std::max(
        std::max(
            (vertex_end - vertex_start).norm(),
            (face_vertex0_end - face_vertex0_start).norm()),
        std::max(
            (face_vertex1_end - face_vertex1_start).norm(),
            (face_vertex2_end - face_vertex2_start).norm()));

    return CCD_LENGTH_TOL / dl;
}

Eigen::Vector3I cross(const Eigen::Vector3I& a, const Eigen::Vector3I& b)
{
    Eigen::Vector3I c;
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

bool are_points_on_same_side_of_edge(
    const Eigen::Vector3I& p1,
    const Eigen::Vector3I& p2,
    const Eigen::Vector3I& a,
    const Eigen::Vector3I& b)
{
    Eigen::Vector3I cp1 = cross(b - a, p1 - a); //(b - a).cross(p1 - a);
    Eigen::Vector3I cp2 = cross(b - a, p2 - a); //(b - a).cross(p2 - a);
    return cp1.dot(cp2).upper() >= 0;
}

bool is_point_inside_triangle_(
    const Eigen::Vector3I& point,
    const Eigen::Vector3I& triangle_vertex0,
    const Eigen::Vector3I& triangle_vertex1,
    const Eigen::Vector3I& triangle_vertex2)
{
    return are_points_on_same_side_of_edge(
               point, triangle_vertex0, triangle_vertex1, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex1, triangle_vertex0, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex2, triangle_vertex0, triangle_vertex1);
}

inline Eigen::Vector3I triangle_normal(
    const Eigen::Vector3I& face_vertex0,
    const Eigen::Vector3I& face_vertex1,
    const Eigen::Vector3I& face_vertex2)
{
    return cross(
        face_vertex1 - face_vertex0,
        face_vertex2
            - face_vertex0); //(face_vertex1 - face_vertex0).cross(face_vertex2
                             //- face_vertex0);
}

bool vertexFaceCCD_Redon(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi)
{
    const auto distance = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex
            = (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t
            + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0 = (face_vertex0_end.cast<Interval>()
                                        - face_vertex0_start.cast<Interval>())
                * t
            + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1 = (face_vertex1_end.cast<Interval>()
                                        - face_vertex1_start.cast<Interval>())
                * t
            + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2 = (face_vertex2_end.cast<Interval>()
                                        - face_vertex2_start.cast<Interval>())
                * t
            + face_vertex2_start.cast<Interval>();

        return (vertex - face_vertex0)
            .dot(triangle_normal(face_vertex0, face_vertex1, face_vertex2));
    };

    const auto is_point_inside_triangle = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex
            = (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t
            + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0 = (face_vertex0_end.cast<Interval>()
                                        - face_vertex0_start.cast<Interval>())
                * t
            + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1 = (face_vertex1_end.cast<Interval>()
                                        - face_vertex1_start.cast<Interval>())
                * t
            + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2 = (face_vertex2_end.cast<Interval>()
                                        - face_vertex2_start.cast<Interval>())
                * t
            + face_vertex2_start.cast<Interval>();

        return is_point_inside_triangle_(
            vertex, face_vertex0, face_vertex1, face_vertex2);
    };

    double tol = compute_face_vertex_tolerance_1d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_point_inside_triangle, Interval(0, 1), tol, toi_interval);

    // Return a conservative time-of-impact
    toi = toi_interval.lower();

    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

double compute_edge_edge_tolerance_1d(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max(
            (edge0_vertex0_end - edge0_vertex0_start).norm(),
            (edge0_vertex1_end - edge0_vertex1_start).norm()),
        std::max(
            (edge1_vertex0_end - edge1_vertex0_start).norm(),
            (edge1_vertex1_end - edge1_vertex1_start).norm()));

    return CCD_LENGTH_TOL / dl;
}

bool are_segments_intersecting(
    const Eigen::Vector3I& edge0_vertex0,
    const Eigen::Vector3I& edge0_vertex1,
    const Eigen::Vector3I& edge1_vertex0,
    const Eigen::Vector3I& edge1_vertex1)
{
    const Eigen::Vector3I& a = edge0_vertex0;
    const Eigen::Vector3I& b = edge0_vertex1;
    const Eigen::Vector3I& c = edge1_vertex0;
    const Eigen::Vector3I& d = edge1_vertex1;

    Eigen::Vector3I n = cross(b - a, c - a);

    return cross(b - d, c - d).dot(n).upper() >= 0
        && cross(c - d, a - d).dot(n).upper() >= 0
        && cross(a - d, b - d).dot(n).lower() <= 0;
}

bool edgeEdgeCCD_Redon(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi)
{
    const auto distance = [&](const Interval& t) {
        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();

        return cross(
                   edge0_vertex1 - edge0_vertex0, edge1_vertex1 - edge1_vertex0)
            .dot(edge1_vertex0 - edge0_vertex0);
    };

    const auto is_intersection_inside_edges = [&](const Interval& t) -> bool {
        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();

        return are_segments_intersecting(
            edge0_vertex0, edge0_vertex1, edge1_vertex0, edge1_vertex1);
    };

    double tol = compute_edge_edge_tolerance_1d(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_intersection_inside_edges, Interval(0, 1), tol,
        toi_interval);

    // Return a conservative time-of-impact
    toi = toi_interval.lower();

    return is_impacting;
}

Eigen::Vector3d compute_face_vertex_tolerance_3d(
    const Eigen::Vector3d& vs,
    const Eigen::Vector3d& f0s,
    const Eigen::Vector3d& f1s,
    const Eigen::Vector3d& f2s,
    const Eigen::Vector3d& ve,
    const Eigen::Vector3d& f0e,
    const Eigen::Vector3d& f1e,
    const Eigen::Vector3d& f2e)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max((ve - vs).norm(), (f0e - f0s).norm()),
        std::max((f1e - f1s).norm(), (f2e - f2s).norm()));

    double edge0_length = std::max((f1s - f0s).norm(), (f1e - f0e).norm());
    double edge1_length = std::max((f2s - f0s).norm(), (f2e - f0e).norm());
    // double edge_length = std::max(edge0_length, edge1_length);
    c00 = 0;
    c10 = 0;
    c20 = 0;
    c00 = dl;
    c10 = edge0_length;
    c20 = edge1_length;
    return Eigen::Vector3d(
        CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge0_length,
        CCD_LENGTH_TOL / edge1_length);
}

std::array<Eigen::Vector3d, 2>
bbd_4_pts_new(const std::array<Eigen::Vector3d, 4>& pts)
{
    Eigen::Vector3d min, max;
    min = pts[0];
    max = pts[0];
    for (int i = 1; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            if (min[j] > pts[i][j]) {
                min[j] = pts[i][j];
            }
            if (max[j] < pts[i][j]) {
                max[j] = pts[i][j];
            }
        }
    }
    std::array<Eigen::Vector3d, 2> rst;
    rst[0] = min;
    rst[1] = max;
    return rst;
}

double max_diff(
    const double b1min,
    const double b1max,
    const double b2min,
    const double b2max)
{
    double r = 0;
    if (r < b1max - b1min)
        r = b1max - b1min;
    if (r < b2max - b2min)
        r = b2max - b2min;
    if (r < fabs(b2min - b1max))
        r = fabs(b2min - b1max);
    if (r < fabs(b1min - b2max))
        r = fabs(b1min - b2max);
    return r;
}

// calculate maximum x, y and z diff
double get_max_axis_diff(
    const std::array<Eigen::Vector3d, 2>& b1,
    const std::array<Eigen::Vector3d, 2>& b2)
{

    double x = max_diff(b1[0][0], b1[1][0], b2[0][0], b2[1][0]);
    double y = max_diff(b1[0][1], b1[1][1], b2[0][1], b2[1][1]);
    double z = max_diff(b1[0][2], b1[1][2], b2[0][2], b2[1][2]);
    return std::max(std::max(x, y), z);
}

Eigen::Vector3d compute_face_vertex_tolerance_3d_new_old_version(
    const Eigen::Vector3d& vs,
    const Eigen::Vector3d& f0s,
    const Eigen::Vector3d& f1s,
    const Eigen::Vector3d& f2s,
    const Eigen::Vector3d& ve,
    const Eigen::Vector3d& f0e,
    const Eigen::Vector3d& f1e,
    const Eigen::Vector3d& f2e,
    const double tolerance)
{
    Eigen::Vector3d p000 = vs - f0s, p001 = vs - f2s,
                    p011 = vs - (f1s + f2s - f0s), p010 = vs - f1s;
    Eigen::Vector3d p100 = ve - f0e, p101 = ve - f2e,
                    p111 = ve - (f1e + f2e - f0e), p110 = ve - f1e;
    double dl = 0;
    double edge0_length = 0;
    double edge1_length = 0;
    std::array<Eigen::Vector3d, 2> t0box
        = bbd_4_pts_new({ { p000, p001, p010, p011 } });
    std::array<Eigen::Vector3d, 2> t1box
        = bbd_4_pts_new({ { p100, p101, p110, p111 } });
    std::array<Eigen::Vector3d, 2> u0box
        = bbd_4_pts_new({ { p000, p001, p100, p101 } });
    std::array<Eigen::Vector3d, 2> u1box
        = bbd_4_pts_new({ { p010, p011, p110, p111 } });
    std::array<Eigen::Vector3d, 2> v0box
        = bbd_4_pts_new({ { p000, p100, p010, p110 } });
    std::array<Eigen::Vector3d, 2> v1box
        = bbd_4_pts_new({ { p001, p101, p011, p111 } });
    dl = get_max_axis_diff(t0box, t1box);
    edge0_length = get_max_axis_diff(u0box, u1box);
    edge1_length = get_max_axis_diff(v0box, v1box);
    //     for(int i=0;i<3;i++){
    //         if(dl<fabs(ve[i]-vs[i]))
    //             dl=fabs(ve[i]-vs[i]);

    //         if(dl<fabs(f0e[i]-f0s[i]))
    //             dl=fabs(f0e[i]-f0s[i]);

    //         if(dl<fabs(f1e[i]-f1s[i]))
    //             dl=fabs(f1e[i]-f1s[i]);

    //         if(dl<fabs(f2e[i]-f2s[i]))
    //             dl=fabs(f2e[i]-f2s[i]);

    //         if(edge0_length<fabs(f1s[i] - f0s[i]))
    //             edge0_length=fabs(f1s[i] - f0s[i]);

    //         if(edge0_length<fabs(f1e[i] - f0e[i]))
    //             edge0_length=fabs(f1e[i] - f0e[i]);

    //         if(edge1_length<fabs(f2s[i]-f0s[i]))
    //             edge1_length=fabs(f2s[i]-f0s[i]);

    //         if(edge1_length<fabs(f2e[i]-f0e[i]))
    //             edge1_length=fabs(f2e[i]-f0e[i]);
    //     }
    //    //double edge_length = std::max(edge0_length, edge1_length);
    //     c00=0;c10=0;c20=0;
    //     c00=dl;c10=edge0_length;c20=edge1_length;
    return Eigen::Vector3d(
        tolerance / dl, tolerance / edge0_length, tolerance / edge1_length);
}

double max_linf_dist(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2)
{
    double r = 0;
    for (int i = 0; i < 3; i++) {
        if (r < fabs(p1[i] - p2[i])) {
            r = fabs(p1[i] - p2[i]);
        }
    }
    return r;
}

double max_linf_4(
    const Eigen::Vector3d& p1,
    const Eigen::Vector3d& p2,
    const Eigen::Vector3d& p3,
    const Eigen::Vector3d& p4,
    const Eigen::Vector3d& p1e,
    const Eigen::Vector3d& p2e,
    const Eigen::Vector3d& p3e,
    const Eigen::Vector3d& p4e)
{
    double r = 0, temp = 0;
    temp = max_linf_dist(p1e, p1);
    if (r < temp)
        r = temp;
    temp = max_linf_dist(p2e, p2);
    if (r < temp)
        r = temp;
    temp = max_linf_dist(p3e, p3);
    if (r < temp)
        r = temp;
    temp = max_linf_dist(p4e, p4);
    if (r < temp)
        r = temp;
    return r;
}

Eigen::Vector3d compute_face_vertex_tolerance_3d_new(
    const Eigen::Vector3d& vs,
    const Eigen::Vector3d& f0s,
    const Eigen::Vector3d& f1s,
    const Eigen::Vector3d& f2s,
    const Eigen::Vector3d& ve,
    const Eigen::Vector3d& f0e,
    const Eigen::Vector3d& f1e,
    const Eigen::Vector3d& f2e,
    const double tolerance)
{
    Eigen::Vector3d p000 = vs - f0s, p001 = vs - f2s,
                    p011 = vs - (f1s + f2s - f0s), p010 = vs - f1s;
    Eigen::Vector3d p100 = ve - f0e, p101 = ve - f2e,
                    p111 = ve - (f1e + f2e - f0e), p110 = ve - f1e;
    double dl = 0;
    double edge0_length = 0;
    double edge1_length = 0;
    dl = 3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
    edge0_length
        = 3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
    edge1_length
        = 3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);
    // double diag=max_linf_4(
    // p000,p100,p110,p010,
    // p111,p011,p001,p101);

    // std::array<Eigen::Vector3d,2>
    // t0box=bbd_4_pts_new({{p000,p001,p010,p011}});
    // std::array<Eigen::Vector3d,2>
    // t1box=bbd_4_pts_new({{p100,p101,p110,p111}});
    // std::array<Eigen::Vector3d,2>
    // u0box=bbd_4_pts_new({{p000,p001,p100,p101}});
    // std::array<Eigen::Vector3d,2>
    // u1box=bbd_4_pts_new({{p010,p011,p110,p111}});
    // std::array<Eigen::Vector3d,2>
    // v0box=bbd_4_pts_new({{p000,p100,p010,p110}});
    // std::array<Eigen::Vector3d,2>
    // v1box=bbd_4_pts_new({{p001,p101,p011,p111}});
    // dl=get_max_axis_diff(t0box,t1box);
    // edge0_length=get_max_axis_diff(u0box,u1box);
    // edge1_length=get_max_axis_diff(v0box,v1box);
    //     for(int i=0;i<3;i++){
    //         if(dl<fabs(ve[i]-vs[i]))
    //             dl=fabs(ve[i]-vs[i]);

    //         if(dl<fabs(f0e[i]-f0s[i]))
    //             dl=fabs(f0e[i]-f0s[i]);

    //         if(dl<fabs(f1e[i]-f1s[i]))
    //             dl=fabs(f1e[i]-f1s[i]);

    //         if(dl<fabs(f2e[i]-f2s[i]))
    //             dl=fabs(f2e[i]-f2s[i]);

    //         if(edge0_length<fabs(f1s[i] - f0s[i]))
    //             edge0_length=fabs(f1s[i] - f0s[i]);

    //         if(edge0_length<fabs(f1e[i] - f0e[i]))
    //             edge0_length=fabs(f1e[i] - f0e[i]);

    //         if(edge1_length<fabs(f2s[i]-f0s[i]))
    //             edge1_length=fabs(f2s[i]-f0s[i]);

    //         if(edge1_length<fabs(f2e[i]-f0e[i]))
    //             edge1_length=fabs(f2e[i]-f0e[i]);
    //     }
    //    //double edge_length = std::max(edge0_length, edge1_length);
    //     c00=0;c10=0;c20=0;
    //     c00=dl;c10=edge0_length;c20=edge1_length;
    return Eigen::Vector3d(
        tolerance / dl, tolerance / edge0_length, tolerance / edge1_length);
}

//
bool vertexFaceCCD(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi)
{
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 3);
        Interval t = params(0);
        Interval u = params(1);
        Interval v = params(2);

        Eigen::Vector3I vertex = vertex_start.cast<Interval>()
            + (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t;

        Eigen::Vector3I face0 = face_vertex0_start.cast<Interval>()
            + (face_vertex0_end.cast<Interval>()
               - face_vertex0_start.cast<Interval>())
                * t;

        Eigen::Vector3I face1 = face_vertex1_start.cast<Interval>()
            + (face_vertex1_end.cast<Interval>()
               - face_vertex1_start.cast<Interval>())
                * t;

        Eigen::Vector3I face2 = face_vertex2_start.cast<Interval>()
            + (face_vertex2_end.cast<Interval>()
               - face_vertex2_start.cast<Interval>())
                * t;

        Eigen::Vector3I face_vertex = face0.cast<Interval>()
            + u * (face1.cast<Interval>() - face0.cast<Interval>())
            + v * (face2.cast<Interval>() - face0.cast<Interval>());

        return (vertex - face_vertex).eval();
    };

    Eigen::Vector3d tol = compute_face_vertex_tolerance_3d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,
        true);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

Eigen::Vector3d compute_edge_edge_tolerance(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max(
            (edge0_vertex0_end - edge0_vertex0_start).norm(),
            (edge0_vertex1_end - edge0_vertex1_start).norm()),
        std::max(
            (edge1_vertex0_end - edge1_vertex0_start).norm(),
            (edge1_vertex1_end - edge1_vertex1_start).norm()));

    double edge0_length = std::max(
        (edge0_vertex1_start - edge0_vertex0_start).norm(),
        (edge0_vertex1_end - edge0_vertex0_end).norm());
    double edge1_length = std::max(
        (edge1_vertex1_start - edge1_vertex0_start).norm(),
        (edge1_vertex1_end - edge1_vertex0_end).norm());
    // double edge_length = std::max(edge0_length, edge1_length);

    return Eigen::Vector3d(
        CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge0_length,
        CCD_LENGTH_TOL / edge1_length);
}

Eigen::Vector3d compute_edge_edge_tolerance_new(
    const Eigen::Vector3d& edge0_vertex0_start, // a0s
    const Eigen::Vector3d& edge0_vertex1_start, // a1s
    const Eigen::Vector3d& edge1_vertex0_start, // b0s
    const Eigen::Vector3d& edge1_vertex1_start, // b1s
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    const double tolerance)
{

    Eigen::Vector3d p000 = edge0_vertex0_start - edge1_vertex0_start,
                    p001 = edge0_vertex0_start - edge1_vertex1_start,
                    p011 = edge0_vertex1_start - edge1_vertex1_start,
                    p010 = edge0_vertex1_start - edge1_vertex0_start;
    Eigen::Vector3d p100 = edge0_vertex0_end - edge1_vertex0_end,
                    p101 = edge0_vertex0_end - edge1_vertex1_end,
                    p111 = edge0_vertex1_end - edge1_vertex1_end,
                    p110 = edge0_vertex1_end - edge1_vertex0_end;
    double dl = 0;
    double edge0_length = 0;
    double edge1_length = 0;
    // {
    //     std::array<Eigen::Vector3d,2>
    //     t0box=bbd_4_pts_new({{p000,p001,p010,p011}});
    // std::array<Eigen::Vector3d,2>
    // t1box=bbd_4_pts_new({{p100,p101,p110,p111}});
    // std::array<Eigen::Vector3d,2>
    // u0box=bbd_4_pts_new({{p000,p001,p100,p101}});
    // std::array<Eigen::Vector3d,2>
    // u1box=bbd_4_pts_new({{p010,p011,p110,p111}});
    // std::array<Eigen::Vector3d,2>
    // v0box=bbd_4_pts_new({{p000,p100,p010,p110}});
    // std::array<Eigen::Vector3d,2>
    // v1box=bbd_4_pts_new({{p001,p101,p011,p111}});
    // dl=get_max_axis_diff(t0box,t1box);
    // edge0_length=get_max_axis_diff(u0box,u1box);
    // edge1_length=get_max_axis_diff(v0box,v1box);
    // }

    dl = 3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
    edge0_length
        = 3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
    edge1_length
        = 3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);
    // for(int i=0;i<3;i++){
    //     if(fabs(p000[i]-p100[i])>dl)
    //         dl=fabs(p000[i]-p100[i]);
    //     if(fabs(p001[i]-p101[i])>dl)
    //         dl=fabs(p001[i]-p101[i]);
    //     if(fabs(p010[i]-p110[i])>dl)
    //         dl=fabs(p010[i]-p110[i]);
    //     if(fabs(p011[i]-p111[i])>dl)
    //         dl=fabs(p011[i]-p111[i]);
    // }

    //    double dl=0;
    //    double edge0_length=0;
    //    double edge1_length=0;
    //     for(int i=0;i<3;i++){
    //         if(dl<fabs(edge0_vertex0_end[i]-edge0_vertex0_start[i]))
    //             dl=fabs(edge0_vertex0_end[i]-edge0_vertex0_start[i]);

    //         if(dl<fabs(edge0_vertex1_end[i]-edge0_vertex1_start[i]))
    //             dl=fabs(edge0_vertex1_end[i]-edge0_vertex1_start[i]);

    //         if(dl<fabs(edge1_vertex0_end[i]-edge1_vertex0_start[i]))
    //             dl=fabs(edge1_vertex0_end[i]-edge1_vertex0_start[i]);

    //         if(dl<fabs(edge1_vertex1_end[i]-edge1_vertex1_start[i]))
    //             dl=fabs(edge1_vertex1_end[i]-edge1_vertex1_start[i]);

    //         if(edge0_length<fabs(edge0_vertex1_start[i] -
    //         edge0_vertex0_start[i]))
    //             edge0_length=fabs(edge0_vertex1_start[i] -
    //             edge0_vertex0_start[i]);

    //         if(edge0_length<fabs(edge0_vertex1_end[i] -
    //         edge0_vertex0_end[i]))
    //             edge0_length=fabs(edge0_vertex1_end[i] -
    //             edge0_vertex0_end[i]);

    //         if(edge1_length<fabs(edge1_vertex1_start[i] -
    //         edge1_vertex0_start[i]))
    //             edge1_length=fabs(edge1_vertex1_start[i] -
    //             edge1_vertex0_start[i]);

    //         if(edge1_length<fabs(edge1_vertex1_end[i] -
    //         edge1_vertex0_end[i]))
    //             edge1_length=fabs(edge1_vertex1_end[i] -
    //             edge1_vertex0_end[i]);
    // }
    return Eigen::Vector3d(
        tolerance / dl, tolerance / edge0_length, tolerance / edge1_length);
}

bool edgeEdgeCCD(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi)
{
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 3);
        Interval t = params(0);
        Interval edge0_alpha = params(1);
        Interval edge1_alpha = params(2);

        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex
            = (edge0_vertex1 - edge0_vertex0) * edge0_alpha + edge0_vertex0;

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex
            = (edge1_vertex1 - edge1_vertex0) * edge1_alpha + edge1_vertex0;

        return (edge1_vertex - edge0_vertex).eval();
    };

    Eigen::Vector3d tol = compute_edge_edge_tolerance(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end);

    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,
        false);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

bool edgeEdgeCCD_opt(
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    double& toi)
{
    const auto distance
        = [&](const Numccd& tpara, const Numccd& upara, const Numccd& vpara) {
              int tu = tpara.first, td = tpara.second; // t=tu/(2^td)
              int uu = upara.first, ud = upara.second;
              int vu = vpara.first, vd = vpara.second;

              Eigen::Vector3I edge0_vertex0
                  = (a0e.cast<Interval>() - a0s.cast<Interval>()) * tu
                      / int(pow(2, td))
                  + a0s.cast<Interval>();
              Eigen::Vector3I edge0_vertex1
                  = (a1e.cast<Interval>() - a1s.cast<Interval>()) * tu
                      / int(pow(2, td))
                  + a1s.cast<Interval>();
              Eigen::Vector3I edge0_vertex
                  = (edge0_vertex1 - edge0_vertex0) * uu / int(pow(2, ud))
                  + edge0_vertex0;

              Eigen::Vector3I edge1_vertex0
                  = (b0e.cast<Interval>() - b0s.cast<Interval>()) * tu
                      / int(pow(2, td))
                  + b0s.cast<Interval>();
              Eigen::Vector3I edge1_vertex1
                  = (b1e.cast<Interval>() - b1s.cast<Interval>()) * tu
                      / int(pow(2, td))
                  + b1s.cast<Interval>();
              Eigen::Vector3I edge1_vertex
                  = (edge1_vertex1 - edge1_vertex0) * vu / int(pow(2, vd))
                  + edge1_vertex0;

              return (edge1_vertex - edge0_vertex).eval();
          };
    Eigen::Vector3d tol
        = compute_edge_edge_tolerance(a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

    Interval3 toi_interval;
    igl::Timer timer;
    timer.start();
    bool is_impacting
        = interval_root_finder_opt(distance, tol, toi_interval, false);
    time0 += timer.getElapsedTimeInMicroSec();
    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = double(toi_interval[0].first.first)
            / pow(2, toi_interval[0].first.second);
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
    return false;
}

bool edgeEdgeCCD_double(
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const std::array<double, 3>& err,
    const double ms,
    double& toi,
    const double tolerance,
    const double pre_check_t)
{

    Eigen::Vector3d tol = compute_edge_edge_tolerance_new(
        a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tolerance);
    
    //this should be the error of the whole mesh
    std::array<double, 3> err1;
    if (err[0] < 0) { // if error[0]<0, means we need to calculate error here
        std::vector<Eigen::Vector3d> vlist;
        vlist.emplace_back(a0s);
        vlist.emplace_back(a1s);
        vlist.emplace_back(b0s);
        vlist.emplace_back(b1s);

        vlist.emplace_back(a0e);
        vlist.emplace_back(a1e);
        vlist.emplace_back(b0e);
        vlist.emplace_back(b1e);

        err1 = get_numerical_error(vlist, false);
    } else {
        err1 = err;
    }

    //////////////////////////////////////////////////////////

    igl::Timer timer;
    timer.start();
    bool is_impacting;

    // 0 is normal ccd without pre-check, 1 is ccd with pre-check, 2 is ccd with pre-check and horizontal tree
    int CCD_TYPE=2;
    if(CCD_TYPE==0){
        is_impacting = interval_root_finder_double_normalCCD(
        tol, toi, false, err1, ms,a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
    }
    if(CCD_TYPE==1){
        is_impacting = interval_root_finder_double_pre_check(
        tol, toi, false, err1, ms,a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, pre_check_t);
    }
    if(CCD_TYPE==2){
        is_impacting = interval_root_finder_double_horizontal_tree(
        tol,tolerance, toi, false, err1, ms,a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, pre_check_t);
    }
    
    time0 += timer.getElapsedTimeInMicroSec();
    // Return a conservative time-of-impact
    //    if (is_impacting) {
    //        toi =
    //        double(toi_interval[0].first.first)/pow(2,toi_interval[0].first.second);
    //    }
    // This time of impact is very dangerous for convergence
    assert(!is_impacting || toi >= 0);
    return is_impacting;
    return false;
}

bool vertexFaceCCD_double(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    const std::array<double, 3>& err,
    const double ms,
    double& toi,
    const double tolerance,
    const double pre_check_t)
{
    Eigen::Vector3d tol = compute_face_vertex_tolerance_3d_new(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end, tolerance);
    
    
    
    //////////////////////////////////////////////////////////
    // this is the error of the whole mesh
    std::array<double, 3> err1;
    if (err[0] < 0) { // if error[0]<0, means we need to calculate error here
        std::vector<Eigen::Vector3d> vlist;
        vlist.emplace_back(vertex_start);
        vlist.emplace_back(face_vertex0_start);
        vlist.emplace_back(face_vertex1_start);
        vlist.emplace_back(face_vertex2_start);

        vlist.emplace_back(vertex_end);
        vlist.emplace_back(face_vertex0_end);
        vlist.emplace_back(face_vertex1_end);
        vlist.emplace_back(face_vertex2_end);

        err1 = get_numerical_error(vlist, false);

    } else {
        err1 = err;
    }
    //////////////////////////////////////////////////////////

    // Interval3 toi_interval;
    igl::Timer timer;
    timer.start();
    bool is_impacting;

    // 0 is normal ccd without pre-check, 1 is ccd with pre-check, 2 is ccd with pre-check and horizontal tree
    int CCD_TYPE=2;
    if(CCD_TYPE==0){
        is_impacting = interval_root_finder_double_normalCCD(
        tol, toi, true, err1, ms,vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    }
    if(CCD_TYPE==1){
        is_impacting = interval_root_finder_double_pre_check(
        tol, toi, true, err1, ms,vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end, pre_check_t);
    }
    if(CCD_TYPE==2){
        is_impacting = interval_root_finder_double_horizontal_tree(
        tol,tolerance, toi, true, err1, ms,vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end, pre_check_t);
    }
    // Return a conservative time-of-impact
    //    if (is_impacting) {
    //        toi =
    //        double(toi_interval[0].first.first)/pow(2,toi_interval[0].first.second);
    //    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
    return false;
}

bool edgeEdgeCCD_rational(
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const std::array<double, 3>& err,
    const double ms,
    double& toi)
{

    Eigen::Vector3d tol
        = compute_edge_edge_tolerance(a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

    //////////////////////////////////////////////////////////
    // TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(a0s);
    vlist.emplace_back(a1s);
    vlist.emplace_back(b0s);
    vlist.emplace_back(b1s);

    vlist.emplace_back(a0e);
    vlist.emplace_back(a1e);
    vlist.emplace_back(b0e);
    vlist.emplace_back(b1e);

    std::array<double, 3> err1;
    err1 = get_numerical_error(vlist, false);
    //////////////////////////////////////////////////////////

    std::array<std::pair<Rational, Rational>, 3> toi_interval;
    igl::Timer timer;
    timer.start();
    bool is_impacting = interval_root_finder_Rational(
        tol, toi_interval, false, err1, ms, a0s, a1s, b0s, b1s, a0e, a1e, b0e,
        b1e);
    time0 += timer.getElapsedTimeInMicroSec();
    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval[0].first.to_double();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
    return false;
}

bool vertexFaceCCD_rational(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    const std::array<double, 3>& err,
    const double ms,
    double& toi)
{

    Eigen::Vector3d tol = compute_face_vertex_tolerance_3d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    // std::cout<<"get tolerance successfully"<<std::endl;
    //////////////////////////////////////////////////////////
    // TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(vertex_start);
    vlist.emplace_back(face_vertex0_start);
    vlist.emplace_back(face_vertex1_start);
    vlist.emplace_back(face_vertex2_start);

    vlist.emplace_back(vertex_end);
    vlist.emplace_back(face_vertex0_end);
    vlist.emplace_back(face_vertex1_end);
    vlist.emplace_back(face_vertex2_end);

    std::array<double, 3> err1;
    err1 = get_numerical_error(vlist, false);
    // std::cout<<"get error successfully"<<std::endl;
    //////////////////////////////////////////////////////////

    std::array<std::pair<Rational, Rational>, 3> toi_interval;
    igl::Timer timer;
    timer.start();
    bool is_impacting = interval_root_finder_Rational(
        tol, toi_interval, true, err1, ms, vertex_start, face_vertex0_start,
        face_vertex1_start, face_vertex2_start, vertex_end, face_vertex0_end,
        face_vertex1_end, face_vertex2_end);
    time0 += timer.getElapsedTimeInMicroSec();
    // std::cout<<"get result successfully"<<std::endl;
    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval[0].first.to_double();
    }
    // std::cout<<"get time successfully"<<std::endl;
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
    return false;
}

void print_time_1()
{
    std::cout << "time of root finder, " << time0 << std::endl;
}

void print_tol()
{
    std::cout << "tol, " << tol0 << " " << tol1 << " " << tol2 << std::endl;
    std::cout << "toln, " << tol0n << " " << tol1n << " " << tol2n << std::endl;
    if (tol0 == 0 || tol1 == 0 || tol2 == 0) {
        std::cout << "***** tolerance is 0 for some cases" << std::endl;
    }
    std::cout << "lengths, " << c00 << ", " << c10 << ", " << c20 << std::endl;
}
} // namespace intervalccd
