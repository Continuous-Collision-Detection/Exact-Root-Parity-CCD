/// Our exact CCD method
#include "ccd.hpp"
#include "Rational.hpp"
#include "Utils.hpp"
#include "subfunctions.h"
namespace ccd {

// Detect collisions between a vertex and a triangular face.
bool vertexFaceCCD(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end)
{
    throw "not implemented";
    
}

bool vertexFaceCCD(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    const double eps)
{

    prism vfprism(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);

    Vector3d bmin(-eps, -eps, -eps), bmax(eps, eps, eps);
    bool bbox_intersection = vfprism.is_prism_bbox_cut_bbox(bmin, bmax);
    if (!bbox_intersection)
        return false;// if bounding box not intersected, then not intersected
    for (int i = 0; i < 9; i++) {
        i;
	}
	



	return false;
}

// Detect collisions between two edges as they move.
bool edgeEdgeCCD(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end)
{
    throw "not implemented";
}
void test() { std::cout << "compiles correct " << std::endl; }
} // namespace ccd
