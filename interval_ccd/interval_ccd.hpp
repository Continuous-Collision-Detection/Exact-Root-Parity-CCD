// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <interval_ccd/interval.hpp>

#include <array>

/**
 * @namespace ccd:
 * @brief
 */
namespace intervalccd {

///////////////////////////////////////////////////////////////////////////////
// 2D

// This can be used in 3D too
bool edgeVertexCCD(
    const Eigen::Vector2d& vertex_start,
    const Eigen::Vector2d& edge_vertex0_start,
    const Eigen::Vector2d& edge_vertex1_start,
    const Eigen::Vector2d& vertex_end,
    const Eigen::Vector2d& edge_vertex0_end,
    const Eigen::Vector2d& edge_vertex1_end,
    double& toi);

///////////////////////////////////////////////////////////////////////////////
// 3D

bool vertexFaceCCD(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi);

bool edgeEdgeCCD(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi);

bool edgeEdgeCCD_opt(
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    double& toi);

// This function can give you the answer of continous collision detection with minimum 
// seperation, and the earlist collision time if collision happens.
// err is the filters calculated using the bounding box of the simulation scene.
// If you are checking a single query without a scene, please set it as [-1,-1,-1].
// ms is the minimum seperation. should set: ms < max(abs(x)), ms < max(abs(y)), ms < max(abs(z)) of the scene.
// toi is the earlist time of collision if collision happens. If there is no collision, toi will be infinate.
// tolerance is a user - input solving precision. we suggest to use 1e-6.
// pre_check_t is a number that allows you to check collision for time =[0, 1+2*tol_x+pre],
// to avoid collision time = 0 for next time step, which is useful in simulations using line-search
// max_itr is a user-defined value to terminate the algorithm earlier, and return a result under current 
// precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
// earlier and the precision will be user-defined precision -- tolerance.
// output_tolerance is the precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
// CCD_TYPE is a switch to choose root-finding methods. 0 is the normal CCD root finding which cannot be used 
// for line - search; 1 is the un-optimized pre-check method which can be used for line - search; 2 is the method 
// can be used for line - search, and has a user - input max_itr.
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
    const double pre_check_t,
    const int max_itr,
    double &output_tolerance,
    const int CCD_TYPE);

// This function can give you the answer of continous collision detection with minimum 
// seperation, and the earlist collision time if collision happens.
// err is the filters calculated using the bounding box of the simulation scene.
// If you are checking a single query without a scene, please set it as [-1,-1,-1].
// ms is the minimum seperation. should set: ms < max(abs(x)), ms < max(abs(y)), ms < max(abs(z)) of the scene.
// toi is the earlist time of collision if collision happens. If there is no collision, toi will be infinate.
// tolerance is a user - input solving precision. we suggest to use 1e-6.
// pre_check_t is a number that allows you to check collision for time =[0, 1+2*tol_x+pre],
// to avoid collision time = 0 for next time step, which is useful in simulations using line-search
// max_itr is a user-defined value to terminate the algorithm earlier, and return a result under current 
// precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
// earlier and the precision will be user-defined precision -- tolerance.
// output_tolerance is the precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
// CCD_TYPE is a switch to choose root-finding methods. 0 is the normal CCD root finding which cannot be used 
// for line - search; 1 is the un-optimized pre-check method which can be used for line - search; 2 is the method 
// can be used for line - search, and has a user - input max_itr.
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
    const double pre_check_t,
    const int max_itr,
    double &output_tolerance,
    const int CCD_TYPE);

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
    double& toi);

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
    double& toi);

bool vertexFaceCCD_Redon(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi);

bool edgeEdgeCCD_Redon(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi);

void print_time_1();

void print_tol();

} // namespace intervalccd
