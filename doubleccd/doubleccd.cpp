/// Our exact CCD method
#include <doubleCCD/doubleccd.hpp>

namespace doubleccd {

// Detect collisions between a vertex and a triangular face.
bool vertexFaceCCD(
    const Vector3d& vertex_start,
    const Vector3d& face_vertex0_start,
    const Vector3d& face_vertex1_start,
    const Vector3d& face_vertex2_start,
    const Vector3d& vertex_end,
    const Vector3d& face_vertex0_end,
    const Vector3d& face_vertex1_end,
    const Vector3d& face_vertex2_end)
{
    // Rational rtn=minimum_distance;

    // std::cout << "eps "<<rtn.to_double() << std::endl;
    prism vfprism(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);

    Vector3d bmin(0, 0, 0), bmax(0, 0, 0);
    bool intersection = vfprism.is_prism_bbox_cut_bbox(bmin, bmax);
    if (!intersection) {

        return false; // if bounding box not intersected, then not intersected
    }
    bilinear bl0(
        vfprism.p_vertices[0], vfprism.p_vertices[1], vfprism.p_vertices[4],
        vfprism.p_vertices[3]);
    bilinear bl1(
        vfprism.p_vertices[1], vfprism.p_vertices[2], vfprism.p_vertices[5],
        vfprism.p_vertices[4]);
    bilinear bl2(
        vfprism.p_vertices[0], vfprism.p_vertices[2], vfprism.p_vertices[5],
        vfprism.p_vertices[3]);
    std::array<bilinear, 3> bls = { { bl0, bl1, bl2 } };
    bool oin = shoot_origin_ray_prism(vfprism, bls);
    return oin;
}

// Detect collisions between two edges as they move.
bool edgeEdgeCCD(
    const Vector3d& edge0_vertex0_start,
    const Vector3d& edge0_vertex1_start,
    const Vector3d& edge1_vertex0_start,
    const Vector3d& edge1_vertex1_start,
    const Vector3d& edge0_vertex0_end,
    const Vector3d& edge0_vertex1_end,
    const Vector3d& edge1_vertex0_end,
    const Vector3d& edge1_vertex1_end)
{

    hex hx(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end);

    // step 1. bounding box checking

    Vector3d bmin(0, 0, 0), bmax(0, 0, 0);
    bool intersection = hx.is_hex_bbox_cut_bbox(bmin, bmax);

    if (!intersection)
        return false; // if bounding box not intersected, then not intersected
    bilinear bl0(
        hx.h_vertices[0], hx.h_vertices[1], hx.h_vertices[2], hx.h_vertices[3]);
    bilinear bl1(
        hx.h_vertices[4], hx.h_vertices[5], hx.h_vertices[6], hx.h_vertices[7]);
    bilinear bl2(
        hx.h_vertices[0], hx.h_vertices[1], hx.h_vertices[5], hx.h_vertices[4]);
    bilinear bl3(
        hx.h_vertices[1], hx.h_vertices[2], hx.h_vertices[6], hx.h_vertices[5]);
    bilinear bl4(
        hx.h_vertices[2], hx.h_vertices[3], hx.h_vertices[7], hx.h_vertices[6]);
    bilinear bl5(
        hx.h_vertices[0], hx.h_vertices[3], hx.h_vertices[7], hx.h_vertices[4]);
    std::array<bilinear, 6> bls = { { bl0, bl1, bl2, bl3, bl4, bl5 } };
    bool oin = shoot_origin_ray_hex(bls);
    return oin;
}

void test()
{
#ifdef CCD_ROUND_INPUTS
    std::cout << "rounding can be enabled" << std::endl;
#endif
}
} // namespace doubleccd
