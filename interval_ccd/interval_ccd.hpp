// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include "interval.hpp"

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
    const Eigen::Vector3d& b1e, double &toi);
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
    double &toi);
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
    double& toi);

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
    double &toi);

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
    double& toi);
void print_time_1();

} // namespace ccd
