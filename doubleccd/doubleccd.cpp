/// Our exact CCD method
#include <doubleCCD/doubleccd.hpp>
#include <igl/Timer.h>
igl::Timer timer1, timer2;
double time1 = 0, time2 = 0, time3 = 0, time4=0, time5=0, time6=0, time7=0, time8=0;
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
    const Vector3d& face_vertex2_end,
    const double minimum_distance)
{
    //Rational rtn=minimum_distance;

    //std::cout << "eps "<<rtn.to_double() << std::endl;
    prism vfprism(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end);
//std::cout<<"ori 0 "<<orient_3d(Vector3d(0,0,0),vfprism.p_vertices[0],vfprism.p_vertices[1],vfprism.p_vertices[2])<<std::endl;
    //std::cout<<"ori 1 "<<orient_3d(Vector3d(0,0,0),vfprism.p_vertices[3],vfprism.p_vertices[4],vfprism.p_vertices[5])<<std::endl;
    // step 1. bounding box checking
    Vector3d bmin(-minimum_distance, -minimum_distance, -minimum_distance),
        bmax(minimum_distance, minimum_distance, minimum_distance);
    bool intersection = vfprism.is_prism_bbox_cut_bbox(bmin, bmax);
    if (!intersection)
        return false; // if bounding box not intersected, then not intersected

    // step 2. prism edges & prism bottom triangles check
    // prism edges test, segment degenerate cases already handled
    // std::cout << "before seg- intersect cube in double" << std::endl;
    for (int i = 0; i < 9; i++) {
        if (is_seg_intersect_cube(
                minimum_distance,
                vfprism.p_vertices[vfprism.prism_edge_id[i][0]],
                vfprism.p_vertices[vfprism.prism_edge_id[i][1]])) {
             // std::cout << "which seg intersect cube "<<i << std::endl;
            return true;
        }
    }
    //std::cout << "after seg_inter_cube " << std::endl;
    // std::cout << "before cube edge intersect 2 triangles in double" <<
    // std::endl;
    // prism top/bottom triangle test
    cube cb(minimum_distance);
    if (!vfprism.is_triangle_degenerated(0)) {
        if (is_cube_edges_intersect_triangle(
                cb, vfprism.p_vertices[0], vfprism.p_vertices[1],
                vfprism.p_vertices[2])) // if this triangle is not degenerated
            return true;
    }
    if (!vfprism.is_triangle_degenerated(1)) {
        if (is_cube_edges_intersect_triangle(
                cb, vfprism.p_vertices[3], vfprism.p_vertices[4],
                vfprism.p_vertices[5]))
            return true;
    }
    //std::cout<<"after cube triangle"<<std::endl;
    // step 3 tet facets- cube edges
    std::array<std::array<bool, 8>, 3>
        v_tet; // cube vertices - tets positions
               // TODO one solution for v_tet is to  make a boolean which can
               // show if pt is on the border
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
    bool cube_inter_tet[3];
    // std::cout << "before cube - opposite faces in double" << std::endl;
    if (is_cube_intersect_tet_opposite_faces(
            bl0, cb, v_tet[0], cube_inter_tet[0]))
        return true;
    if (is_cube_intersect_tet_opposite_faces(
            bl1, cb, v_tet[1], cube_inter_tet[1]))
        return true;
    if (is_cube_intersect_tet_opposite_faces(
            bl2, cb, v_tet[2], cube_inter_tet[2]))
        return true;

    // if cube intersect any tet, need check if intersect bilinear;
    // if cube not intersect any tet, shoot a ray
    // TODO we can also have some information about the edge-face intersection
    // above std::cout << "before cube inter bilinear in double" << std::endl;
    if (cube_inter_tet[0]) {
        timer1.start();
        bool cit0 = is_cube_edge_intersect_bilinear(bl0, cb, v_tet[0]);
        time1 += timer1.getElapsedTimeInSec();
        if (cit0)
            return true;
    }
    if (cube_inter_tet[1]) {
        timer1.start();
        bool cit1 = is_cube_edge_intersect_bilinear(bl1, cb, v_tet[1]);
        time1 += timer1.getElapsedTimeInSec();
        if (cit1)
            return true;
    }
    if (cube_inter_tet[2]) {
        timer1.start();
        bool cit2 = is_cube_edge_intersect_bilinear(bl2, cb, v_tet[2]);
        time1 += timer1.getElapsedTimeInSec();
        if (cit2)
            return true;
    }
    // down here is the last part of the algorithm

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
    // std::cout << "shoot a ray in double" << std::endl;
    std::vector<bool> p_tet;
    p_tet.resize(3);
    p_tet[0] = v_tet[0][target];
    p_tet[1] = v_tet[1][target];
    p_tet[2] = v_tet[2][target];
    timer2.start();
    bool rtccd = retrial_ccd(vfprism, bls, cb.vr[target], p_tet);
    time2 += timer2.getElapsedTimeInSec();
    return rtccd;
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
    const Vector3d& edge1_vertex1_end,
    const double minimum_distance)
{
    //Rational rtn=minimum_distance;

    //std::cout << "eps "<<rtn.to_double() << std::endl;
    timer1.start();
    hex hx(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end);
    time3+=timer1.getElapsedTimeInSec();
    // step 1. bounding box checking
    timer1.start();
    Vector3d bmin(-minimum_distance, -minimum_distance, -minimum_distance),
        bmax(minimum_distance, minimum_distance, minimum_distance);
    bool intersection = hx.is_hex_bbox_cut_bbox(bmin, bmax);
    time3+=timer1.getElapsedTimeInSec();
    if (!intersection)
        return false; // if bounding box not intersected, then not intersected

    // step 2. prism edges & prism bottom triangles check
    // prism edges test, segment degenerate cases already handled
    // std::cout << "before seg- intersect cube in double" << std::endl;
    bool rt=false;
    timer1.start();
    for (int i = 0; i < 12; i++) {
        if (is_seg_intersect_cube(
                minimum_distance, hx.h_vertices[hx.hex_edge_id[i][0]],
                hx.h_vertices[hx.hex_edge_id[i][1]])) {
            // std::cout << "which seg intersect cube "<<i << std::endl;
            rt=true;
            break;
        }
    }
    
    time4+=timer1.getElapsedTimeInSec();
    if (rt) return true;
    //std::cout<<"go after seg_cube"<<std::endl;
    cube cb(minimum_distance);

    // step 3 tet facets- cube edges
    std::array<std::array<bool, 8>, 6> v_tet; // cube vertices - tets positions
    timer1.start();
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
    time5+=timer1.getElapsedTimeInSec();
    bool cube_inter_tet[6];
    // std::cout << "before cube - opposite faces in double" << std::endl;
    timer1.start();
    int discrete=5;
    for(int i=0;i<6;i++){
        //std::cout<<"start check opp"<<std::endl;
        if(is_cube_intersect_tet_opposite_faces(
            bls[i], cb, v_tet[i], cube_inter_tet[i])){

            // bool rr=cube_discrete_bilinear_intersection(cb,bls[i],5);
            // if(!rr){
            //     std::cout<<"result do not match in opposite check, ori vs discrete "<<1<<" "<<rr<<std::endl;
            // }
            //std::cout<<"ith iteration break "<<i<<std::endl;
            rt=true;
            break;
        }
            
    }
    time8+=timer1.getElapsedTimeInSec();
    if(rt) return true;
    //std::cout<<"go after oppo"<<std::endl;
    for(int i=0;i<6;i++){
        if (cube_inter_tet[i]) {
            timer1.start();
            bool cit0 = is_cube_edge_intersect_bilinear(bls[i], cb, v_tet[i]);
            // bool rr=cube_discrete_bilinear_intersection(cb,bls[i],5);
            // if(cit0!=rr){
            //     std::cout<<"result do not match, ori vs discrete "<<cit0<<" "<<rr<<std::endl;
            // }
            time6+=timer1.getElapsedTimeInSec();
            if (cit0)
            //std::cout<<"return intersect bilinear"<<std::endl;
                return true;
        }
    }

    timer1.start();
    int min_v = 6;
    int curr_v;
    int target = 0;
    for (int i = 0; i < 8; i++) {
        curr_v = v_tet[0][i] + v_tet[1][i] + v_tet[2][i] + v_tet[3][i]
            + v_tet[4][i] + v_tet[5][i];
        if (curr_v < min_v) {
            min_v = curr_v;
            target = i;
        }
    }
    // std::cout << "shoot a ray in double" << std::endl;
    std::vector<bool> p_tet;
    p_tet.resize(6);
    p_tet[0] = v_tet[0][target];
    p_tet[1] = v_tet[1][target];
    p_tet[2] = v_tet[2][target];
    p_tet[3] = v_tet[3][target];
    p_tet[4] = v_tet[4][target];
    p_tet[5] = v_tet[5][target];

    bool rtccd = retrial_ccd_hex(bls, cb.vr[target], p_tet);
    time7+=timer1.getElapsedTimeInSec();
    return rtccd;
    return false;
}

void test()
{
    //std::cout << "time for cube edge - bilinear " << time1 << std::endl;
    //std::cout << "time for retrail ccd  " << time2 << std::endl;
    std::cout << "edge time for initial " << time3 << std::endl;
    std::cout << "edge time for seg-cube " << time4 << std::endl;
    std::cout << "edge time for build bilinears " << time5 << std::endl;
    std::cout << "edge time for bilinear opposite faces " << time8 << std::endl;
    std::cout << "edge time for not degenerated bilinears " << time6 << std::endl;
    std::cout << "edge time for shooting ray " << time7 << std::endl;
    
}
} // namespace doubleccd
