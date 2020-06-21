#include "interval_ccd.hpp"

#include "interval_root_finder.hpp"

namespace ccd {

static double CCD_LENGTH_TOL = 1e-6;

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

Eigen::Vector2d compute_edge_vertex_tolerance(
    const Eigen::Vector2d& vertex_start,
    const Eigen::Vector2d& edge_vertex0_start,
    const Eigen::Vector2d& edge_vertex1_start,
    const Eigen::Vector2d& vertex_end,
    const Eigen::Vector2d& edge_vertex0_end,
    const Eigen::Vector2d& edge_vertex1_end)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        (vertex_end - vertex_start).norm(),
        std::max(
            (edge_vertex0_end - edge_vertex0_start).norm(),
            (edge_vertex1_end - edge_vertex1_start).norm()));

    double edge_length = std::max(
        (edge_vertex1_start - edge_vertex0_start).norm(),
        (edge_vertex1_end - edge_vertex0_end).norm());

    return Eigen::Vector2d(CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge_length);
}

bool edgeVertexCCD(
    const Eigen::Vector2d& vertex_start,
    const Eigen::Vector2d& edge_vertex0_start,
    const Eigen::Vector2d& edge_vertex1_start,
    const Eigen::Vector2d& vertex_end,
    const Eigen::Vector2d& edge_vertex0_end,
    const Eigen::Vector2d& edge_vertex1_end,
    double& toi)
{
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 2);
        Interval t = params(0);
        Interval alpha = params(1);

        // Get the vertex at time t
        Eigen::Vector2I vertex = (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t + vertex_start.cast<Interval>();

        // Get the vertex of the edge at time t
        Eigen::Vector2I edge_vertex0
            = (edge_vertex0_end.cast<Interval>() - edge_vertex0_start.cast<Interval>()) * t + edge_vertex0_start.cast<Interval>();
        Eigen::Vector2I edge_vertex1
            = (edge_vertex1_end - edge_vertex1_start).cast<Interval>() * t + edge_vertex1_start.cast<Interval>();
        Eigen::Vector2I edge_vertex
            = (edge_vertex1 - edge_vertex0).cast<Interval>() * alpha + edge_vertex0.cast<Interval>();

        return (vertex - edge_vertex).eval();
    };

    Eigen::Vector2d tol = compute_edge_vertex_tolerance(
        vertex_start, edge_vertex0_start, edge_vertex1_start, //
        vertex_end, edge_vertex0_end, edge_vertex1_end);

    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::VectorX3I::Constant(2, Interval(0, 1)),tol, toi_interval);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }

    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Face-Vertex

double compute_face_vertex_tolerance(
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

inline Eigen::Vector3I triangle_normal(
    const Eigen::Vector3I& face_vertex0,
    const Eigen::Vector3I& face_vertex1,
    const Eigen::Vector3I& face_vertex2)
{
    return (face_vertex1 - face_vertex0).cross(face_vertex2 - face_vertex0);
}

bool is_point_inside_triangle(
    const Eigen::Vector3I& point,
    const Eigen::Vector3I& triangle_vertex0,
    const Eigen::Vector3I& triangle_vertex1,
    const Eigen::Vector3I& triangle_vertex2)
{
    Eigen::Vector3I normal0
        = triangle_normal(triangle_vertex0, triangle_vertex1, point);
    Eigen::Vector3I normal1
        = triangle_normal(triangle_vertex0, point, triangle_vertex2);
    Eigen::Vector3I normal2
        = triangle_normal(triangle_vertex1, triangle_vertex2, point);
    for (int i = 0; i < normal0.size(); i++) {
        if (!overlap(intersect(normal0(i), normal1(i)), normal2(i))) {
            return false;
        }
    }
    return true;
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
    const auto distance = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex = (vertex_end - vertex_start).cast<Interval>() * t + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0
            = (face_vertex0_end - face_vertex0_start).cast<Interval>() * t + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1
            = (face_vertex1_end - face_vertex1_start).cast<Interval>() * t + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2
            = (face_vertex2_end - face_vertex2_start).cast<Interval>() * t + face_vertex2_start.cast<Interval>();

        const auto result=(vertex - face_vertex0)
            .dot(triangle_normal(face_vertex0, face_vertex1, face_vertex2));
        return result;
       // return ccd::Interval(0,1);
    };

    const auto is_point_in_triangle = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex = (vertex_end - vertex_start).cast<Interval>() * t + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0
            = (face_vertex0_end - face_vertex0_start).cast<Interval>() * t + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1
            = (face_vertex1_end - face_vertex1_start).cast<Interval>() * t + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2
            = (face_vertex2_end - face_vertex2_start).cast<Interval>() * t + face_vertex2_start.cast<Interval>();

        return is_point_inside_triangle(
            vertex, face_vertex0, face_vertex1, face_vertex2);
    };

    double tol = compute_face_vertex_tolerance(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_point_in_triangle, Interval(0, 1), tol, toi_interval);

    // Return a conservative time-of-impact
    toi = toi_interval.lower();

    return is_impacting;

    //return 0;
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
   double edge_length = std::max(edge0_length, edge1_length);

   return Eigen::Vector3d(CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge0_length, CCD_LENGTH_TOL / edge1_length);
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

       Eigen::Vector3I edge0_vertex0
           = (edge0_vertex0_end - edge0_vertex0_start).cast<Interval>() * t
           + edge0_vertex0_start.cast<Interval>();
       Eigen::Vector3I edge0_vertex1
           = (edge0_vertex1_end - edge0_vertex1_start).cast<Interval>() * t
           + edge0_vertex1_start.cast<Interval>();
       Eigen::Vector3I edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * edge0_alpha + edge0_vertex0;

       Eigen::Vector3I edge1_vertex0
           = (edge1_vertex0_end - edge1_vertex0_start).cast<Interval>() * t
           + edge1_vertex0_start.cast<Interval>();
       Eigen::Vector3I edge1_vertex1
           = (edge1_vertex1_end - edge1_vertex1_start).cast<Interval>() * t
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
       distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval);

   // Return a conservative time-of-impact
   if (is_impacting) {
       toi = toi_interval(0).lower();
   }
   // This time of impact is very dangerous for convergence
   // assert(!is_impacting || toi > 0);
   return is_impacting;
}

} // namespace ccd
