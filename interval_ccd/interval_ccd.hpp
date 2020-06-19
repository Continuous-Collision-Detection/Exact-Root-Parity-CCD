// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <physics/rigid_body.hpp>

/**
 * @namespace ccd:
 * @brief
 */
namespace ccd {

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

} // namespace ccd
