/// Our exact CCD method
#include <ccd.hpp>
namespace ccd {

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
    throw "not implemented";
}

bool vertexFaceCCD(
    const Vector3d& vertex_start,
    const Vector3d& face_vertex0_start,
    const Vector3d& face_vertex1_start,
    const Vector3d& face_vertex2_start,
    const Vector3d& vertex_end,
    const Vector3d& face_vertex0_end,
    const Vector3d& face_vertex1_end,
    const Vector3d& face_vertex2_end,
    const double eps)
{

    prism vfprism(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);

    // step 1. bounding box checking
    Vector3r bmin(-eps, -eps, -eps), bmax(eps, eps, eps);
    bool intersection = vfprism.is_prism_bbox_cut_bbox(bmin, bmax);
    if (!intersection)
        return false; // if bounding box not intersected, then not intersected

    // step 2. prism edges & prism bottom triangles check
    // prism edges test, segment degenerate cases already handled
    for (int i = 0; i < 9; i++) {
        if (is_seg_intersect_cube(
                eps, vfprism.p_vertices[vfprism.prism_edge_id[i][0]],
                vfprism.p_vertices[vfprism.prism_edge_id[i][1]]))
            return true;
    }

    // prism top/bottom triangle test
    cube cb(eps);
    if (is_cube_edges_intersect_triangle(
            cb, vfprism.p_vertices[0], vfprism.p_vertices[1],
            vfprism.p_vertices[2]))
        return true;
    if (is_cube_edges_intersect_triangle(
            cb, vfprism.p_vertices[3], vfprism.p_vertices[4],
            vfprism.p_vertices[5]))
        return true;

    // step 3 tet facets- cube edges
    std::array<std::array<bool, 8>, 3> v_tet;//cube vertices - tets positions
	//TODO one solution for v_tet is to  make a boolean which can show if pt is on the border
    bilinear bl0(
        vfprism.p_vertices[0], vfprism.p_vertices[1], vfprism.p_vertices[4],
        vfprism.p_vertices[3]);
    bilinear bl1(
        vfprism.p_vertices[1], vfprism.p_vertices[2], vfprism.p_vertices[5],
        vfprism.p_vertices[4]);
    bilinear bl2(
        vfprism.p_vertices[0], vfprism.p_vertices[2], vfprism.p_vertices[5],
        vfprism.p_vertices[3]);
	std::array<bilinear, 3> bls = { {bl0,bl1,bl2} };
    bool cube_inter_tet[3];
	if (is_cube_intersect_tet_opposite_faces(bl0, cb, v_tet[0], cube_inter_tet[0]))
		return true;
	if (is_cube_intersect_tet_opposite_faces(bl1, cb, v_tet[1], cube_inter_tet[1]))
		return true;
	if (is_cube_intersect_tet_opposite_faces(bl2, cb, v_tet[2], cube_inter_tet[2]))
        return true;

	// if cube intersect any tet, need check if intersect bilinear;
		// if cube not intersect any tet, shoot a ray
	//TODO we can also have some information about the edge-face intersection above
	if (cube_inter_tet[0]) {
		if (is_cube_edge_intersect_bilinear(bl0, cb, v_tet[0]))
			return true;
	}
	if (cube_inter_tet[1]) {
		if (is_cube_edge_intersect_bilinear(bl1, cb, v_tet[1]))
			return true;
	}
	if (cube_inter_tet[2]) {
		if (is_cube_edge_intersect_bilinear(bl2, cb, v_tet[2]))
			return true;
	}
	//down here is the last part of the algorithm
	
	int min_v = 3;
    int curr_v;
    int target = 0;
    for (int i = 0; i < 8; i++) {
        curr_v = v_tet[0][i] + v_tet[1][i] + v_tet[2][i];
        if (curr_v < min_v) {
            min_v = curr_v;
            target = i;
        }
    }

    std::vector<bool> p_tet;
    p_tet.resize(3);
    p_tet[0] = v_tet[0][target];
    p_tet[1] = v_tet[1][target];
    p_tet[2] = v_tet[2][target];
	return retrial_ccd(vfprism, bls, cb.vr[target], p_tet);
	return 0;
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
    throw "not implemented";
}
void test() { std::cout << "compiles correct " << std::endl; }
} // namespace ccd
