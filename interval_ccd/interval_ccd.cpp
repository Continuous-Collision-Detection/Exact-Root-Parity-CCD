#include "interval_ccd.hpp"

#include "interval_root_finder.hpp"

namespace intervalccd {

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
        distance, Eigen::VectorX3I::Constant(2, Interval(0, 1)),tol, toi_interval,false);

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
       std::max(
           (ve - vs).norm(),
           (f0e - f0s).norm()),
       std::max(
           (f1e - f1s).norm(),
           (f2e - f2s).norm()));

   double edge0_length = std::max(
       (f1s - f0s).norm(),
       (f1e - f0e).norm());
   double edge1_length = std::max(
       (f2s - f0s).norm(),
       (f2e - f0e).norm());
   //double edge_length = std::max(edge0_length, edge1_length);

   return Eigen::Vector3d(CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge0_length, CCD_LENGTH_TOL / edge1_length);
}
// Eigen::Vector3I cross(const Eigen::Vector3I&v1, const Eigen::Vector3I&v2){
//     Eigen::Vector3I res;
//     res[0] = v1[1] * v2[2] - v1[2] * v2[1];
//     res[1] = v1[2] * v2[0] - v1[0] * v2[2];
//     res[2] = v1[0] * v2[1] - v1[1] * v2[0];

//     return res;
// }
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
    // Eigen::Vector3I normal0
    //     = triangle_normal(triangle_vertex0, triangle_vertex1, point);
    // Eigen::Vector3I normal1
    //     = triangle_normal(triangle_vertex0, point, triangle_vertex2);
    // Eigen::Vector3I normal2
    //     = triangle_normal(triangle_vertex1, triangle_vertex2, point);
    // for (int i = 0; i < normal0.size(); i++) {
    //     if (!overlap(intersect(normal0(i), normal1(i)), normal2(i))) {
    //         return false;
    //     }
    // }
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
    const auto distance = [&](const Eigen::VectorX3I& params) {
       assert(params.size() == 3);
       Interval t = params(0);
       Interval u = params(1);
       Interval v = params(2);

        Eigen::Vector3I vertex=vertex_start.cast<Interval>()+(vertex_end.cast<Interval>()-vertex_start.cast<Interval>())*t;

        Eigen::Vector3I face0=face_vertex0_start.cast<Interval>()+
        (face_vertex0_end.cast<Interval>()-face_vertex0_start.cast<Interval>())*t;

        Eigen::Vector3I face1=face_vertex1_start.cast<Interval>()+
        (face_vertex1_end.cast<Interval>()-face_vertex1_start.cast<Interval>())*t;

        Eigen::Vector3I face2=face_vertex2_start.cast<Interval>()+
        (face_vertex2_end.cast<Interval>()-face_vertex2_start.cast<Interval>())*t;

        Eigen::Vector3I face_vertex=face0.cast<Interval>()+u*(face1.cast<Interval>()-face0.cast<Interval>())
        +v*(face2.cast<Interval>()-face0.cast<Interval>());

        return (vertex-face_vertex).eval();
       
   };

   Eigen::Vector3d tol = compute_face_vertex_tolerance_3d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    Eigen::VectorX3I toi_interval;
   bool is_impacting = interval_root_finder(
       distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,true);

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
   //double edge_length = std::max(edge0_length, edge1_length);

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
           = (edge0_vertex0_end.cast<Interval>() - edge0_vertex0_start.cast<Interval>()) * t
           + edge0_vertex0_start.cast<Interval>();
       Eigen::Vector3I edge0_vertex1
           = (edge0_vertex1_end.cast<Interval>() - edge0_vertex1_start.cast<Interval>()) * t
           + edge0_vertex1_start.cast<Interval>();
       Eigen::Vector3I edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * edge0_alpha + edge0_vertex0;

       Eigen::Vector3I edge1_vertex0
           = (edge1_vertex0_end.cast<Interval>() - edge1_vertex0_start.cast<Interval>()) * t
           + edge1_vertex0_start.cast<Interval>();
       Eigen::Vector3I edge1_vertex1
           = (edge1_vertex1_end.cast<Interval>() - edge1_vertex1_start.cast<Interval>()) * t
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
       distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,false);

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
    const Eigen::Vector3d& b1e, double &toi){
    const auto distance = [&](const Paraccd& params) {
       
       int tu = params[0].first, td=params[0].second;// t=tu/(2^td)
       int uu = params[1].first, ud=params[1].second;
       int vu = params[2].first, vd=params[2].second;

       Eigen::Vector3I edge0_vertex0
           = (a0e.cast<Interval>() - a0s.cast<Interval>()) * tu/pow(2,td)
           + a0s.cast<Interval>();
       Eigen::Vector3I edge0_vertex1
           = (a1e.cast<Interval>() - a1s.cast<Interval>()) * tu/pow(2,td)
           + a1s.cast<Interval>();
       Eigen::Vector3I edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/pow(2,ud) + edge0_vertex0;

       Eigen::Vector3I edge1_vertex0
           = (b0e.cast<Interval>() - b0s.cast<Interval>()) * tu/pow(2,td)
           + b0s.cast<Interval>();
       Eigen::Vector3I edge1_vertex1
           = (b1e.cast<Interval>() - b1s.cast<Interval>()) * tu/pow(2,td)
           + b1s.cast<Interval>();
       Eigen::Vector3I edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/pow(2,vd) + edge1_vertex0;

       return (edge1_vertex - edge0_vertex).eval();
   };
    Eigen::Vector3d tol = compute_edge_edge_tolerance(a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);

    Eigen::VectorX3I toi_interval;
//    bool is_impacting = interval_root_finder_opt(
//        distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,false);

//    // Return a conservative time-of-impact
//    if (is_impacting) {
//        toi = toi_interval(0).lower();
//    }
//    // This time of impact is very dangerous for convergence
//    // assert(!is_impacting || toi > 0);
//    return is_impacting;
   return false;
}
} // namespace ccd
