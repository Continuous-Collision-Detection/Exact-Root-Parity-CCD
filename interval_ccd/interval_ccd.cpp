#include "interval_ccd.hpp"

#include "interval_root_finder.hpp"
#include<igl/Timer.h>
#include<iostream>
namespace intervalccd {
double time0=0,time1=0,time2=0;
static double CCD_LENGTH_TOL = 1e-6;
double tol0=0,tol1=0,tol2=0,tol0n=0,tol1n=0,tol2n=0;
double c00,c01,c10,c11,c20,c21;

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
    c00=0;c10=0;c20=0;
    c00=dl;c10=edge0_length;c20=edge1_length;
   return Eigen::Vector3d(CCD_LENGTH_TOL / dl, CCD_LENGTH_TOL / edge0_length, CCD_LENGTH_TOL / edge1_length);
}
Eigen::Vector3d compute_face_vertex_tolerance_3d_new(
   const Eigen::Vector3d& vs,
   const Eigen::Vector3d& f0s,
   const Eigen::Vector3d& f1s,
   const Eigen::Vector3d& f2s,
   const Eigen::Vector3d& ve,
   const Eigen::Vector3d& f0e,
   const Eigen::Vector3d& f1e,
   const Eigen::Vector3d& f2e)
{
   
    double dl=0;
    double edge0_length=0;
    double edge1_length=0;
    for(int i=0;i<3;i++){
        if(dl<fabs(ve[i]-vs[i]))
            dl=fabs(ve[i]-vs[i]);

        if(dl<fabs(f0e[i]-f0s[i]))
            dl=fabs(f0e[i]-f0s[i]);

        if(dl<fabs(f1e[i]-f1s[i]))
            dl=fabs(f1e[i]-f1s[i]);

        if(dl<fabs(f2e[i]-f2s[i]))
            dl=fabs(f2e[i]-f2s[i]);

        if(edge0_length<fabs(f1s[i] - f0s[i]))
            edge0_length=fabs(f1s[i] - f0s[i]);

        if(edge0_length<fabs(f1e[i] - f0e[i]))
            edge0_length=fabs(f1e[i] - f0e[i]);

        if(edge1_length<fabs(f2s[i]-f0s[i]))
            edge1_length=fabs(f2s[i]-f0s[i]);

        if(edge1_length<fabs(f2e[i]-f0e[i]))
            edge1_length=fabs(f2e[i]-f0e[i]);
    }
   //double edge_length = std::max(edge0_length, edge1_length);
    c00=0;c10=0;c20=0;
    c00=dl;c10=edge0_length;c20=edge1_length;
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
Eigen::Vector3d compute_edge_edge_tolerance_new(
   const Eigen::Vector3d& edge0_vertex0_start,
   const Eigen::Vector3d& edge0_vertex1_start,
   const Eigen::Vector3d& edge1_vertex0_start,
   const Eigen::Vector3d& edge1_vertex1_start,
   const Eigen::Vector3d& edge0_vertex0_end,
   const Eigen::Vector3d& edge0_vertex1_end,
   const Eigen::Vector3d& edge1_vertex0_end,
   const Eigen::Vector3d& edge1_vertex1_end)
{
   
   double dl=0;
   double edge0_length=0;
   double edge1_length=0;
    for(int i=0;i<3;i++){
        if(dl<fabs(edge0_vertex0_end[i]-edge0_vertex0_start[i]))
            dl=fabs(edge0_vertex0_end[i]-edge0_vertex0_start[i]);

        if(dl<fabs(edge0_vertex1_end[i]-edge0_vertex1_start[i]))
            dl=fabs(edge0_vertex1_end[i]-edge0_vertex1_start[i]);

        if(dl<fabs(edge1_vertex0_end[i]-edge1_vertex0_start[i]))
            dl=fabs(edge1_vertex0_end[i]-edge1_vertex0_start[i]);

        if(dl<fabs(edge1_vertex1_end[i]-edge1_vertex1_start[i]))
            dl=fabs(edge1_vertex1_end[i]-edge1_vertex1_start[i]);

        if(edge0_length<fabs(edge0_vertex1_start[i] - edge0_vertex0_start[i]))
            edge0_length=fabs(edge0_vertex1_start[i] - edge0_vertex0_start[i]);

        if(edge0_length<fabs(edge0_vertex1_end[i] - edge0_vertex0_end[i]))
            edge0_length=fabs(edge0_vertex1_end[i] - edge0_vertex0_end[i]);


        if(edge1_length<fabs(edge1_vertex1_start[i] - edge1_vertex0_start[i]))
            edge1_length=fabs(edge1_vertex1_start[i] - edge1_vertex0_start[i]);

        if(edge1_length<fabs(edge1_vertex1_end[i] - edge1_vertex0_end[i]))
            edge1_length=fabs(edge1_vertex1_end[i] - edge1_vertex0_end[i]);
    }
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
    const auto distance = [&](const Numccd&tpara, const Numccd&upara, const Numccd&vpara) {
       
       int tu = tpara.first, td=tpara.second;// t=tu/(2^td)
       int uu = upara.first, ud=upara.second;
       int vu = vpara.first, vd=vpara.second;

       Eigen::Vector3I edge0_vertex0
           = (a0e.cast<Interval>() - a0s.cast<Interval>()) * tu/int(pow(2,td))
           + a0s.cast<Interval>();
       Eigen::Vector3I edge0_vertex1
           = (a1e.cast<Interval>() - a1s.cast<Interval>()) * tu/int(pow(2,td))
           + a1s.cast<Interval>();
       Eigen::Vector3I edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/int(pow(2,ud)) + edge0_vertex0;

       Eigen::Vector3I edge1_vertex0
           = (b0e.cast<Interval>() - b0s.cast<Interval>()) * tu/int(pow(2,td))
           + b0s.cast<Interval>();
       Eigen::Vector3I edge1_vertex1
           = (b1e.cast<Interval>() - b1s.cast<Interval>()) * tu/int(pow(2,td))
           + b1s.cast<Interval>();
       Eigen::Vector3I edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/int(pow(2,vd)) + edge1_vertex0;

       return (edge1_vertex - edge0_vertex).eval();
   };
    Eigen::Vector3d tol = compute_edge_edge_tolerance(a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);

    Interval3 toi_interval;
    igl::Timer timer;
    timer.start();
   bool is_impacting = interval_root_finder_opt(
       distance, tol,toi_interval,false);
    time0+=timer.getElapsedTimeInMicroSec();
   // Return a conservative time-of-impact
   if (is_impacting) {
       toi = double(toi_interval[0].first.first)/pow(2,toi_interval[0].first.second);
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
    const std::array<double,3>& err,
    const double ms,
    double &toi){
   //tol0=0;tol1=0;tol2=0;
    Eigen::Vector3d tol = compute_edge_edge_tolerance_new(a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
    // tol0=tol[0];
    // tol1=tol[1];
    // tol2=tol[2];
    //bool flag=false;
    // for(int i=0;i<3;i++){
    //     if (tol[i]==0)
    //     flag=true;
    // }
    // if(flag==true){

    // }
    //////////////////////////////////////////////////////////
    //TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(a0s);
    vlist.emplace_back(a1s);
    vlist.emplace_back(b0s);
    vlist.emplace_back(b1s);

    vlist.emplace_back(a0e);
    vlist.emplace_back(a1e);
    vlist.emplace_back(b0e);
    vlist.emplace_back(b1e);

    std::array<double,3> err1;
    err1=get_numerical_error(vlist,false);
    //////////////////////////////////////////////////////////

    Interval3 toi_interval;
    igl::Timer timer;
    timer.start();
   bool is_impacting = interval_root_finder_double(
        tol,toi_interval,false,err1,ms,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
    time0+=timer.getElapsedTimeInMicroSec();
   // Return a conservative time-of-impact
   if (is_impacting) {
       toi = double(toi_interval[0].first.first)/pow(2,toi_interval[0].first.second);
   }
   // This time of impact is very dangerous for convergence
   // assert(!is_impacting || toi > 0);
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
    const std::array<double,3>& err,
    const double ms,
    double& toi)
{
    
   Eigen::Vector3d tol = compute_face_vertex_tolerance_3d_new(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    tol0=tol[0];
    tol1=tol[1];
    tol2=tol[2];
    // Eigen::Vector3d toln=compute_face_vertex_tolerance_3d_new(
    //     vertex_start, face_vertex0_start, face_vertex1_start,
    //     face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
    //     face_vertex2_end);
    // tol0n=toln[0];
    // tol1n=toln[1];
    // tol2n=toln[2];
    //std::cout<<"get tolerance successfully"<<std::endl;
     //////////////////////////////////////////////////////////
    //TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(vertex_start);
    vlist.emplace_back(face_vertex0_start);
    vlist.emplace_back(face_vertex1_start);
    vlist.emplace_back(face_vertex2_start);

    vlist.emplace_back(vertex_end);
    vlist.emplace_back(face_vertex0_end);
    vlist.emplace_back(face_vertex1_end);
    vlist.emplace_back(face_vertex2_end);

    std::array<double,3> err1;
    err1=get_numerical_error(vlist,false);
    //////////////////////////////////////////////////////////

    Interval3 toi_interval;
    igl::Timer timer;
    timer.start();
   bool is_impacting = interval_root_finder_double(
        tol,toi_interval,true,err1,ms,vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    time0+=timer.getElapsedTimeInMicroSec();
   // Return a conservative time-of-impact
   if (is_impacting) {
       toi = double(toi_interval[0].first.first)/pow(2,toi_interval[0].first.second);
   }
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
    const std::array<double,3>& err,
    const double ms,
    double &toi){
   
    Eigen::Vector3d tol = compute_edge_edge_tolerance(a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);

    //////////////////////////////////////////////////////////
    //TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(a0s);
    vlist.emplace_back(a1s);
    vlist.emplace_back(b0s);
    vlist.emplace_back(b1s);

    vlist.emplace_back(a0e);
    vlist.emplace_back(a1e);
    vlist.emplace_back(b0e);
    vlist.emplace_back(b1e);

    std::array<double,3> err1;
    err1=get_numerical_error(vlist,false);
    //////////////////////////////////////////////////////////

    std::array<std::pair<Rational,Rational>, 3> toi_interval;
    igl::Timer timer;
    timer.start();
   bool is_impacting = interval_root_finder_Rational(
        tol,toi_interval,false,err1,ms,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
    time0+=timer.getElapsedTimeInMicroSec();
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
    const std::array<double,3>& err,
    const double ms,
    double& toi)
{
    
   Eigen::Vector3d tol = compute_face_vertex_tolerance_3d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
        //std::cout<<"get tolerance successfully"<<std::endl;
     //////////////////////////////////////////////////////////
    //TODO this should be the error of the whole mesh
    std::vector<Eigen::Vector3d> vlist;
    vlist.emplace_back(vertex_start);
    vlist.emplace_back(face_vertex0_start);
    vlist.emplace_back(face_vertex1_start);
    vlist.emplace_back(face_vertex2_start);

    vlist.emplace_back(vertex_end);
    vlist.emplace_back(face_vertex0_end);
    vlist.emplace_back(face_vertex1_end);
    vlist.emplace_back(face_vertex2_end);

    std::array<double,3> err1;
    err1=get_numerical_error(vlist,false);
    //std::cout<<"get error successfully"<<std::endl;
    //////////////////////////////////////////////////////////

    std::array<std::pair<Rational,Rational>, 3> toi_interval;
    igl::Timer timer;
    timer.start();
   bool is_impacting = interval_root_finder_Rational(
        tol,toi_interval,true,err1,ms,vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
    time0+=timer.getElapsedTimeInMicroSec();
    //std::cout<<"get result successfully"<<std::endl;
   // Return a conservative time-of-impact
   if (is_impacting) {
       toi = toi_interval[0].first.to_double();
   }
   //std::cout<<"get time successfully"<<std::endl;
   // This time of impact is very dangerous for convergence
   // assert(!is_impacting || toi > 0);
   return is_impacting;
   return false;
    
}
void print_time_1(){
    std::cout<<"time of root finder, "<<time0<<std::endl;
}
void print_tol(){
    std::cout<<"tol, "<<tol0<<" "<<tol1<<" "<<tol2<<std::endl;
    std::cout<<"toln, "<<tol0n<<" "<<tol1n<<" "<<tol2n<<std::endl;
    if(tol0==0||tol1==0||tol2==0){
        std::cout<<"***** tolerance is 0 for some cases"<<std::endl;
    }
    std::cout<<"lengths, "<<c00<<", "<<c10<<", "<<c20<<std::endl;
}
} // namespace ccd
