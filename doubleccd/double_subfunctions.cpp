#include <Eigen/Dense>
#include <doubleCCD/double_subfunctions.h>
#include <doubleCCD/doubleccd.hpp>
#include <doubleCCD/exact_subtraction.hpp>




namespace doubleccd {

// cube
cube::cube(double eps)
{
    vr[0] = Vector3d(-eps, -eps, eps), vr[1] = Vector3d(eps, -eps, eps),
    vr[2] = Vector3d(eps, eps, eps), vr[3] = Vector3d(-eps, eps, eps),
    vr[4] = Vector3d(-eps, -eps, -eps), vr[5] = Vector3d(eps, -eps, -eps),
    vr[6] = Vector3d(eps, eps, -eps), vr[7] = Vector3d(-eps, eps, -eps);
    edgeid[0] = { { 0, 1 } };
    edgeid[1] = { { 1, 2 } };
    edgeid[2] = { { 2, 3 } };
    edgeid[3] = { { 3, 0 } };
    edgeid[4] = { { 4, 5 } };
    edgeid[5] = { { 5, 6 } };
    edgeid[6] = { { 6, 7 } };
    edgeid[7] = { { 7, 4 } };
    edgeid[8] = { { 0, 4 } };
    edgeid[9] = { { 1, 5 } };
    edgeid[10] = { { 2, 6 } };
    edgeid[11] = { { 3, 7 } };
    faceid[0] = { { 0, 1, 2, 3 } };
    faceid[1] = { { 4, 7, 6, 5 } };
    faceid[2] = { { 0, 4, 5, 1 } };
    faceid[3] = { { 1, 5, 6, 2 } };
    faceid[4] = { { 3, 2, 6, 7 } };
    faceid[5] = { { 0, 3, 7, 4 } }; // orientation out
    bmax = vr[2];
    bmin = vr[4];
    epsilon = eps;
}

// get aabb corners
void get_corners(const Eigen::MatrixX3d& p, Vector3d& min, Vector3d& max)
{
    min = p.colwise().minCoeff();
    max = p.colwise().maxCoeff();
}
void get_tet_corners(const std::array<Vector3d, 4>& p, Vector3d& min, Vector3d& max)
{

    min = p[0];
    max = p[0];
    for (int i = 0; i < 4; i++) {
        if (min[0] > p[i][0])
            min[0] = p[i][0];
        if (min[1] > p[i][1])
            min[1] = p[i][1];
        if (min[2] > p[i][2])
            min[2] = p[i][2];

        if (max[0] < p[i][0])
            max[0] = p[i][0];
        if (max[1] < p[i][1])
            max[1] = p[i][1];
        if (max[2] < p[i][2])
            max[2] = p[i][2];
    }
}
void get_edge_coners(const Vector3d& e0, const Vector3d& e1, Vector3d &emin,Vector3d &emax){
    for(int i=0;i<3;i++){
    if(e0[i]>e1[i]){
        emin[i]=e1[i];
        emax[i]=e0[i];
    }
    else{
        emin[i]=e0[i];
        emax[i]=e1[i];
    }
    }
    
}

Vector3d get_prism_corner_double(
    const Vector3d& vertex_start,       // x0
    const Vector3d& face_vertex0_start, // x1
    const Vector3d& face_vertex1_start, // x2
    const Vector3d& face_vertex2_start, // x3
    const Vector3d& vertex_end,
    const Vector3d& face_vertex0_end,
    const Vector3d& face_vertex1_end,
    const Vector3d& face_vertex2_end,
    int i)
{
    Vector3d x0 = vertex_start, x1 = face_vertex0_start,
             x2 = face_vertex1_start, x3 = face_vertex2_start, x0b = vertex_end,
             x1b = face_vertex0_end, x2b = face_vertex1_end,
             x3b = face_vertex2_end;
    if (i == 0)
        return x0 - x1;
    if (i == 1)
        return x0 - x3;
    if (i == 2)
        return x0 - x2;
    if (i == 3)
        return x0b - x1b;
    if (i == 4)
        return x0b - x3b;
    if (i == 5)
        return x0b - x2b;

    else
        return Vector3d();
}

bool is_seg_intersect_cube(
    const double& eps, const Vector3d& e0, const Vector3d& e1)
{
    if (is_point_intersect_cube(eps, e0))
        return true;
    if (is_point_intersect_cube(eps, e1))
        return true;
    if (same_point(e0, e1))
        return false; // degenerate case: the segment is degenerated as a point
    // if intersected, must be coplanar with the edge, or intersect edge or face
    if (is_seg_intersect_cube_2d(eps, e0, e1, 0)
        && is_seg_intersect_cube_2d(eps, e0, e1, 1)
        && is_seg_intersect_cube_2d(eps, e0, e1, 2)) {

        return true;
    }

    return false;
}
// check if a 2d segment intersects 2d cube
bool is_seg_intersect_cube_2d(
    const double eps, const Vector3d& e0, const Vector3d& e1, int axis)
{
    Vector2d p0, p1, p2, p3, e0p, e1p;         // e0 and e1 projected to 2d
    projected_cube_edges(eps, p0, p1, p2, p3); // TODO move this out

    const int i1 = (axis + 1) % 3;
    const int i2 = (axis + 2) % 3;

    if (e0[i1] <= eps && e0[i1] >= -eps && e0[i2] <= eps && e0[i2] >= -eps)
        return true;
    if (e1[i1] <= eps && e1[i1] >= -eps && e1[i2] <= eps && e1[i2] >= -eps)
        return true;
    e0p = Vector2d(e0[i1], e0[i2]);
    e1p = Vector2d(e1[i1], e1[i2]);
    if (segment_segment_intersection_2d(e0p, e1p, p0, p1))
        return true; // check if segments has intersection, or if cube points
                     // p0, p1 on e0-e1
    if (segment_segment_intersection_2d(e0p, e1p, p1, p2))
        return true;
    if (segment_segment_intersection_2d(e0p, e1p, p2, p3))
        return true;
    if (segment_segment_intersection_2d(e0p, e1p, p3, p0))
        return true;

    return false;
}
void projected_cube_edges(
    const double eps, Vector2d& e0, Vector2d& e1, Vector2d& e2, Vector2d& e3)
{
    const int i1 = 0;
    const int i2 = 1;

    e0[i1] = -eps;
    e0[i2] = eps;

    e1[i1] = eps;
    e1[i2] = eps;

    e2[i1] = eps;
    e2[i2] = -eps;

    e3[i1] = -eps;
    e3[i2] = -eps;
}
bool is_point_intersect_cube(const double eps, const Vector3d& p)
{
    if (p[0] <= eps && p[0] >= -eps) {
        if (p[1] <= eps && p[1] >= -eps) {
            if (p[2] <= eps && p[2] >= -eps) {
                return true;
            }
        }
    }
    return false;
}


bool is_cube_edges_intersect_triangle(
    const cube& cb, const Vector3d& t0, const Vector3d& t1, const Vector3d& t2)
{
    // the vertices of triangle are checked before going here, the edges are
    // also checked. so, only need to check if cube edge has intersection with
    // the open triangle.

    /// if triangle degenerated as a segment or point, then
    //    no
    //                  // intersection, because before here we already check
    //                  that
    Vector3d s0, s1;
    for (int i = 0; i < 12; i++) {
        s0 = cb.vr[cb.edgeid[i][0]];
        s1 = cb.vr[cb.edgeid[i][1]];
        if (segment_triangle_intersection(s0, s1, t0, t1, t2, false)
            > 0) // return 0,1,2

            return true;
    }
    return false;
}

// convert a array of subtraction pair to vertices
void convert_to_shifted_v(
    const std::array<std::pair<double, double>, 18>& dt, vf_pair& vs)
{
    vs.x0[0] = dt[0].first;
    vs.x0[1] = dt[1].first;
    vs.x0[2] = dt[2].first;

    vs.x0b[0] = dt[15].first;
    vs.x0b[1] = dt[16].first;
    vs.x0b[2] = dt[17].first;

    vs.x1[0] = dt[0].second;
    vs.x1[1] = dt[1].second;
    vs.x1[2] = dt[2].second;

    vs.x3[0] = dt[3].second;
    vs.x3[1] = dt[4].second;
    vs.x3[2] = dt[5].second;

    vs.x2[0] = dt[6].second;
    vs.x2[1] = dt[7].second;
    vs.x2[2] = dt[8].second;

    vs.x1b[0] = dt[9].second;
    vs.x1b[1] = dt[10].second;
    vs.x1b[2] = dt[11].second;

    vs.x3b[0] = dt[12].second;
    vs.x3b[1] = dt[13].second;
    vs.x3b[2] = dt[14].second;

    vs.x2b[0] = dt[15].second;
    vs.x2b[1] = dt[16].second;
    vs.x2b[2] = dt[17].second;
}
void convert_to_shifted_v(
    const std::array<std::pair<double, double>, 24>& dt, ee_pair& vs)
{
    vs.a0[0] = dt[0].first;
    vs.a0[1] = dt[1].first;
    vs.a0[2] = dt[2].first;

    vs.a1[0] = dt[3].first;
    vs.a1[1] = dt[4].first;
    vs.a1[2] = dt[5].first;

    vs.b0[0] = dt[0].second;
    vs.b0[1] = dt[1].second;
    vs.b0[2] = dt[2].second;

    vs.b1[0] = dt[6].second;
    vs.b1[1] = dt[7].second;
    vs.b1[2] = dt[8].second;
    //////
    vs.a0b[0] = dt[12].first;
    vs.a0b[1] = dt[13].first;
    vs.a0b[2] = dt[14].first;

    vs.a1b[0] = dt[15].first;
    vs.a1b[1] = dt[16].first;
    vs.a1b[2] = dt[17].first;

    vs.b0b[0] = dt[12].second;
    vs.b0b[1] = dt[13].second;
    vs.b0b[2] = dt[14].second;

    vs.b1b[0] = dt[18].second;
    vs.b1b[1] = dt[19].second;
    vs.b1b[2] = dt[20].second;
}
void push_vers_into_subtract_pair(
    const std::vector<vf_pair>& data1,
    const std::vector<ee_pair>& data2,
    std::vector<std::pair<double, double>>& sub)
{
    sub.clear();
    sub.reserve(
        18 * data1.size() + 24 * data2.size()); // each vf_pair has 6*3, ee_pair
                                                // has 8*3 subtractions
    std::pair<double, double> temp;

    for (int j = 0; j < data1.size(); j++) {
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0[i];
            temp.second = data1[j].x1[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0[i];
            temp.second = data1[j].x3[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0[i];
            temp.second = data1[j].x2[i];
            sub.push_back(temp);
        }
        //
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0b[i];
            temp.second = data1[j].x1b[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0b[i];
            temp.second = data1[j].x3b[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data1[j].x0b[i];
            temp.second = data1[j].x2b[i];
            sub.push_back(temp);
        }
    }
    // ee
    for (int j = 0; j < data2.size(); j++) {
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a0[i];
            temp.second = data2[j].b0[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a1[i];
            temp.second = data2[j].b0[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a1[i];
            temp.second = data2[j].b1[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a0[i];
            temp.second = data2[j].b1[i];
            sub.push_back(temp);
        }
        //
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a0b[i];
            temp.second = data2[j].b0b[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a1b[i];
            temp.second = data2[j].b0b[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a1b[i];
            temp.second = data2[j].b1b[i];
            sub.push_back(temp);
        }
        for (int i = 0; i < 3; i++) {
            temp.first = data2[j].a0b[i];
            temp.second = data2[j].b1b[i];
            sub.push_back(temp);
        }
    }
}
// void push_single_subtraction(const Vector3d& p, const Vector3d &q,
// std::vector<std::pair<double, double>>& sub_x,
//     std::vector<std::pair<double, double>>& sub_y,
//     std::vector<std::pair<double, double>>& sub_z ){
//         std::pair<double, double> temp;
//         temp.first=p[0];
//         temp.second=q[0];
//         sub_x.emplace_back(temp);

//         temp.first=p[1];
//         temp.second=q[1];
//         sub_y.emplace_back(temp);

//         temp.first=p[2];
//         temp.second=q[2];
//         sub_z.emplace_back(temp);
//     }
// void push_VF_vers_into_subtract_pair(
//     const std::vector<vf_pair>& data1,
    
//     std::vector<std::pair<double, double>>& sub_x,
//     std::vector<std::pair<double, double>>& sub_y,
//     std::vector<std::pair<double, double>>& sub_z)
// {
//     sub_x.clear();sub_y.clear();sub_z.clear();
//     sub_x.reserve(data1.size()*6);
//     sub_y.reserve(data1.size()*6);
//     sub_z.reserve(data1.size()*6);
    
//     std::pair<double, double> temp;

//     for (int j = 0; j < data1.size(); j++) {
//         push_single_subtraction(data1[j].x0,data1[j].x1,sub_x,sub_y,sub_z);
//         push_single_subtraction(data1[j].x0,data1[j].x3,sub_x,sub_y,sub_z);
//         push_single_subtraction(data1[j].x0,data1[j].x2,sub_x,sub_y,sub_z);

//         push_single_subtraction(data1[j].x0b,data1[j].x1b,sub_x,sub_y,sub_z);
//         push_single_subtraction(data1[j].x0b,data1[j].x3b,sub_x,sub_y,sub_z);
//         push_single_subtraction(data1[j].x0b,data1[j].x2b,sub_x,sub_y,sub_z);
        
        
//     }
// }
double vf_shift_error(const vf_pair& d1, const vf_pair& d2)
{
    double err = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(d1.x0[i] - d2.x0[i]) > err)
            err = fabs(d1.x0[i] - d2.x0[i]);
        if (fabs(d1.x1[i] - d2.x1[i]) > err)
            err = fabs(d1.x1[i] - d2.x1[i]);
        if (fabs(d1.x2[i] - d2.x2[i]) > err)
            err = fabs(d1.x2[i] - d2.x2[i]);
        if (fabs(d1.x3[i] - d2.x3[i]) > err)
            err = fabs(d1.x3[i] - d2.x3[i]);

        if (fabs(d1.x0b[i] - d2.x0b[i]) > err)
            err = fabs(d1.x0b[i] - d2.x0b[i]);
        if (fabs(d1.x1b[i] - d2.x1b[i]) > err)
            err = fabs(d1.x1b[i] - d2.x1b[i]);
        if (fabs(d1.x2b[i] - d2.x2b[i]) > err)
            err = fabs(d1.x2b[i] - d2.x2b[i]);
        if (fabs(d1.x3b[i] - d2.x3b[i]) > err)
            err = fabs(d1.x3b[i] - d2.x3b[i]);
    }
    return err;
}
double ee_shift_error(const ee_pair& d1, const ee_pair& d2)
{
    double err = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(d1.a0[i] - d2.a0[i]) > err)
            err = fabs(d1.a0[i] - d2.a0[i]);
        if (fabs(d1.a1[i] - d2.a1[i]) > err)
            err = fabs(d1.a1[i] - d2.a1[i]);
        if (fabs(d1.b0[i] - d2.b0[i]) > err)
            err = fabs(d1.b0[i] - d2.b0[i]);
        if (fabs(d1.b1[i] - d2.b1[i]) > err)
            err = fabs(d1.b1[i] - d2.b1[i]);

        if (fabs(d1.a0b[i] - d2.a0b[i]) > err)
            err = fabs(d1.a0b[i] - d2.a0b[i]);
        if (fabs(d1.a1b[i] - d2.a1b[i]) > err)
            err = fabs(d1.a1b[i] - d2.a1b[i]);
        if (fabs(d1.b0b[i] - d2.b0b[i]) > err)
            err = fabs(d1.b0b[i] - d2.b0b[i]);
        if (fabs(d1.b1b[i] - d2.b1b[i]) > err)
            err = fabs(d1.b1b[i] - d2.b1b[i]);
    }
    return err;
}

void push_mesh_vers_into_sub_pair(
    const Eigen::MatrixX3d& V, std::vector<std::pair<double, double>>& sub)
{
    sub.resize(V.size());
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < V.cols(); j++) {
            sub[i * 3 + j].first = V(i, j);
            sub[i * 3 + j].second = 0;
        }
    }
}
void push_mesh_vers_into_sub_pair(
    const Eigen::MatrixX3d& V, 
    std::vector<std::pair<double, double>>& sub_x,
    std::vector<std::pair<double, double>>& sub_y,
    std::vector<std::pair<double, double>>& sub_z)
{   sub_x.clear();
    sub_y.clear();
    sub_z.clear();
    sub_x.resize(V.rows());
    sub_y.resize(V.rows());
    sub_z.resize(V.rows());
    for (int i = 0; i < V.rows(); i++) {
        sub_x[i].first=V(i,0);sub_x[i].second=0;
        sub_y[i].first=V(i,1);sub_y[i].second=0;
        sub_z[i].first=V(i,2);sub_z[i].second=0;  
    }
}


void convert_sub_pairs_to_mesh_vers(
    const std::vector<std::pair<double, double>>& sub, Eigen::MatrixX3d& V)
{
    //assert(sub.size() % 3 == 0 && V.size() == sub.size());
    V.resize(sub.size() / 3, 3);
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            V(i, j) = sub[i * 3 + j].first;
        }
    }
}
void convert_sub_pairs_to_mesh_vers(
    const std::vector<std::pair<double, double>>& sub_x,
    const std::vector<std::pair<double, double>>& sub_y,
    const std::vector<std::pair<double, double>>& sub_z, Eigen::MatrixX3d& V){
    V.resize(sub_x.size(),3);
    for(int i=0;i<sub_x.size();i++){
        V(i,0)=sub_x[i].first;
        V(i,1)=sub_y[i].first;
        V(i,2)=sub_z[i].first;
    }
    
}

void get_whole_mesh_shifted(
    Eigen::MatrixX3d& vertices,
    const Vector3d &pmin, const Vector3d& pmax){
    std::vector<std::pair<double, double>> box_x(1), box_y(1), box_z(1),whole_x,whole_y,whole_z;
    box_x[0].first=pmax[0];box_x[0].second=pmin[0];
    box_y[0].first=pmax[1];box_y[0].second=pmin[1];
    box_z[0].first=pmax[2];box_z[0].second=pmin[2];
    push_mesh_vers_into_sub_pair(vertices,whole_x,whole_y,whole_z);
    perturbSubtractions(whole_x,box_x);
    perturbSubtractions(whole_y,box_y);
    perturbSubtractions(whole_z,box_z);
    convert_sub_pairs_to_mesh_vers(whole_x,whole_y,whole_z,vertices);
}
void compare_whole_mesh_err(const Eigen::MatrixX3d& vertices,const Eigen::MatrixX3d& vertices1){
    double ex=0,ey=0,ez=0;
    for(int i=0;i<vertices.rows();i++){
        if(fabs(vertices(i,0)-vertices1(i,0))>ex)
            ex=fabs(vertices(i,0)-vertices1(i,0));

        if(fabs(vertices(i,1)-vertices1(i,1))>ey)
            ey=fabs(vertices(i,1)-vertices1(i,1));

        if(fabs(vertices(i,2)-vertices1(i,2))>ez)
            ez=fabs(vertices(i,2)-vertices1(i,2));
    }
    std::cout<<"vertices diff x, "<<ex<<" y, "<<ey<<" z, "<<ez<<std::endl;
}

// x0 is the point, x1, x2, x3 is the triangle
double get_whole_mesh_shifted(
    const std::vector<vf_pair>& data1,
    const std::vector<ee_pair>& data2,
    std::vector<vf_pair>& shift_back1,
    std::vector<ee_pair>& shift_back2,
    Eigen::MatrixX3d& vertices)
{
    std::vector<std::pair<double, double>> whole, suback;

    push_vers_into_subtract_pair(data1, data2, suback);
    int subsize = suback.size();
    push_mesh_vers_into_sub_pair(vertices, whole);
    // suback = sub;// this is for shift back
    // k = displaceSubtractions_double(sub);
    suback.insert(suback.end(), whole.begin(), whole.end());
    perturbSubtractions(suback); // get shifted back data

    shift_back1.resize(data1.size());

    shift_back2.resize(data2.size());
    int c = 0;
    // Vector3d kvec(k, k, k);
    int d1size = data1.size();
    int d2size = data2.size();
    int datasize = d1size + d2size;
    std::array<std::pair<double, double>, 18> dtback1;
    std::array<std::pair<double, double>, 24> dtback2;
    for (int r = 0; r < datasize; r++) {
        if (r < d1size) {
            for (int i = 0; i < 18; i++) {
                // dt1[i] = sub[r * 18 + i];
                dtback1[i] = suback[r * 18 + i];
            }

            convert_to_shifted_v(dtback1, shift_back1[r]);
        }      // r is d1size-1, sub has been read to d1size*18-1
        else { // r is from d1size to datasize-1, sub is from d1*18
            for (int i = 0; i < 24; i++) {
                // dt2[i] = sub[d1size * 18 + (r - d1size) * 24 + i];
                dtback2[i] = suback[d1size * 18 + (r - d1size) * 24 + i];
            }

            convert_to_shifted_v(dtback2, shift_back2[r - d1size]);
        }
    }
    std::vector<std::pair<double, double>> vernew;
    vernew.resize(vertices.size());
    int c1 = 0;
    for (int i = subsize; i < suback.size(); i++) {
        vernew[c1] = suback[i];
        c1++;
    }
    convert_sub_pairs_to_mesh_vers(vernew, vertices);
    double err = 0, temerr;
    for (int i = 0; i < d1size; i++) {
        temerr = vf_shift_error(data1[i], shift_back1[i]);
        if (temerr > err)
            err = temerr;
    }
    for (int i = 0; i < d2size; i++) {
        temerr = ee_shift_error(data2[i], shift_back2[i]);
        if (temerr > err)
            err = temerr;
    }
    return err;
}

// x0 is the point, x1, x2, x3 is the triangle
// double get_whole_mesh_shifted(
//     const std::vector<vf_pair>& data1,
//     const std::vector<ee_pair>& data2,
//     Eigen::MatrixX3d& vertices)
// {
//     // std::vector<std::pair<double, double>> whole, suback;

//     // push_vers_into_subtract_pair(data1, data2, suback);
//     // int subsize = suback.size();
//     // push_mesh_vers_into_sub_pair(vertices, whole);
//     // // suback = sub;// this is for shift back
//     // // k = displaceSubtractions_double(sub);
//     // suback.insert(suback.end(), whole.begin(), whole.end());
//     // perturbSubtractions(suback); // get shifted back data

//     // std::vector<std::pair<double, double>> vernew;
//     // vernew.resize(vertices.size());
//     // int c = 0;
//     // for (int i = subsize; i < suback.size(); i++) {
//     //     vernew[c] = suback[i];
//     //     c++;
//     // }
//     // assert(whole.size() == vernew.size());
//     double err = 0;
//     // for (int i = 0; i < whole.size(); i++) {
//     //     if (fabs(whole[i].first - vernew[i].first) > err) {
//     //         err = fabs(whole[i].first - vernew[i].first);
//     //     }
//     // }

//     // convert_sub_pairs_to_mesh_vers(vernew, vertices);

//     return err;
// }

double shift_vertex_face(const vf_pair& input_vf_pair, vf_pair& shifted_vf_pair)
{
    // std::vector<vf_pair> input_vf_pairs = { { input_vf_pair } };
    // std::vector<vf_pair> shifted_vf_pairs;
    // std::vector<ee_pair> input_ee_pairs, shifted_ee_pairs;
    // Eigen::MatrixX3d V;
    // double err = get_whole_mesh_shifted(
    //     input_vf_pairs, input_ee_pairs, shifted_vf_pairs, shifted_ee_pairs, V);
    // shifted_vf_pair = shifted_vf_pairs[0];
    // return err;
    Vector3d x0=input_vf_pair.x0, x1=input_vf_pair.x1,x2=input_vf_pair.x2,x3=input_vf_pair.x3,
    x0b=input_vf_pair.x0b, x1b=input_vf_pair.x1b,x2b=input_vf_pair.x2b,x3b=input_vf_pair.x3b;

    std::vector<std::pair<double, double>> subs;
    subs.resize(6*3);
    for(int i=0;i<3;i++){
        subs[3*0+i].first=x0[i];subs[3*0+i].second=x1[i];
        subs[3*1+i].first=x0[i];subs[3*1+i].second=x3[i];
        subs[3*2+i].first=x0[i];subs[3*2+i].second=x2[i];

        subs[3*3+i].first=x0b[i];subs[3*3+i].second=x1b[i];
        subs[3*4+i].first=x0b[i];subs[3*4+i].second=x3b[i];
        subs[3*5+i].first=x0b[i];subs[3*5+i].second=x2b[i];
    }
    perturbSubtractions(subs);
    for(int i=0;i<3;i++){
        x0[i]=subs[3*0+i].first;x1[i]=subs[3*0+i].second;
        x3[i]=subs[3*1+i].second;
        x2[i]=subs[3*2+i].second;

        x0b[i]=subs[3*3+i].first;x1b[i]=subs[3*3+i].second;
        x3b[i]=subs[3*4+i].second;
        x2b[i]=subs[3*5+i].second;
    }
    shifted_vf_pair.x0=x0; 
    shifted_vf_pair.x1=x1;
    shifted_vf_pair.x2=x2;
    shifted_vf_pair.x3=x3;

    shifted_vf_pair.x0b=x0b;
    shifted_vf_pair.x1b=x1b;
    shifted_vf_pair.x2b=x2b;
    shifted_vf_pair.x3b=x3b;
    // double err = 0, temerr;
    for(int i=0;i<subs.size();i++){
        Rational a=subs[i].first;
        Rational b=subs[i].second;
        Rational rst=a-b;
        double rd=subs[i].first-subs[i].second;
        if(rst==rd){

        }
        else{
            std::cout<<"diff, "<<(rst-rd)<<", "<<rst<<", "<<rd<<std::endl;

        }
    }

    return vf_shift_error(input_vf_pair, shifted_vf_pair);
 

    // for (int i = 0; i < d2size; i++) {
    //     temerr = ee_shift_error(data2[i], shift_back2[i]);
    //     if (temerr > err)
    //         err = temerr;
    // }
    // return err;
}

double shift_edge_edge(const ee_pair& input_ee_pair, ee_pair& shifted_ee_pair)
{
    // std::vector<vf_pair> input_vf_pairs, shifted_vf_pairs;
    // std::vector<ee_pair> input_ee_pairs = { { input_ee_pair } };
    // std::vector<ee_pair> shifted_ee_pairs;
    // Eigen::MatrixX3d V;
    // double err = get_whole_mesh_shifted(
    //     input_vf_pairs, input_ee_pairs, shifted_vf_pairs, shifted_ee_pairs, V);
    // shifted_ee_pair = shifted_ee_pairs[0];
    // return err;
    Vector3d a0=input_ee_pair.a0, a1=input_ee_pair.a1,b0=input_ee_pair.b0,b1=input_ee_pair.b1,
    a0b=input_ee_pair.a0b, a1b=input_ee_pair.a1b,b0b=input_ee_pair.b0b,b1b=input_ee_pair.b1b;

    std::vector<std::pair<double, double>> subs;
    subs.resize(8*3);
    for(int i=0;i<3;i++){
        subs[3*0+i].first=a0[i];subs[3*0+i].second=b0[i];
        subs[3*1+i].first=a1[i];subs[3*1+i].second=b0[i];
        subs[3*2+i].first=a1[i];subs[3*2+i].second=b1[i];
        subs[3*3+i].first=a0[i];subs[3*3+i].second=b1[i];

        subs[3*4+i].first=a0b[i];subs[3*4+i].second=b0b[i];
        subs[3*5+i].first=a1b[i];subs[3*5+i].second=b0b[i];
        subs[3*6+i].first=a1b[i];subs[3*6+i].second=b1b[i];
        subs[3*7+i].first=a0b[i];subs[3*7+i].second=b1b[i];
    }
    perturbSubtractions(subs);
    for(int i=0;i<3;i++){
        a0[i]=subs[3*0+i].first;b0[i]=subs[3*0+i].second;
        a1[i]=subs[3*1+i].first;
        b1[i]=subs[3*2+i].second;
        

        a0b[i]=subs[3*4+i].first;b0b[i]=subs[3*4+i].second;
        a1b[i]=subs[3*5+i].first;
        b1b[i]=subs[3*6+i].second;
        
    }
    shifted_ee_pair.a0=a0; 
    shifted_ee_pair.a1=a1;
    shifted_ee_pair.b0=b0;
    shifted_ee_pair.b1=b1;

    shifted_ee_pair.a0b=a0b;
    shifted_ee_pair.a1b=a1b;
    shifted_ee_pair.b0b=b0b;
    shifted_ee_pair.b1b=b1b;
    // double err = 0, temerr;
    for(int i=0;i<subs.size();i++){
        Rational a=subs[i].first;
        Rational b=subs[i].second;
        Rational rst=a-b;
        double rd=subs[i].first-subs[i].second;
        if(rst==rd){

        }
        else{
            std::cout<<"diff, "<<(rst-rd)<<", "<<rst<<", "<<rd<<std::endl;

        }
    }
    return ee_shift_error(input_ee_pair, shifted_ee_pair);
}

// x0 is the point, x1, x2, x3 is the triangle
void get_prism_shifted_vertices_double(
    const Vector3d& x0,
    const Vector3d& x1,
    const Vector3d& x2,
    const Vector3d& x3,
    const Vector3d& x0b,
    const Vector3d& x1b,
    const Vector3d& x2b,
    const Vector3d& x3b,
    double& k,
    std::array<Vector3d, 6>& p_vertices)
{
    std::vector<std::pair<double, double>> sub;

    sub.clear();
    sub.reserve(18);
    std::pair<double, double> temp;
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x1[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x3[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0[i];
        temp.second = x2[i];
        sub.push_back(temp);
    }
    //
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x1b[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x3b[i];
        sub.push_back(temp);
    }
    for (int i = 0; i < 3; i++) {
        temp.first = x0b[i];
        temp.second = x2b[i];
        sub.push_back(temp);
    }

    k = displaceSubtractions_double(sub);
    int c = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            p_vertices[i][j] = sub[c].first - sub[c].second;
            c++;
        }
    }
}
// check if p1-p2 has truncation
bool have_no_truncation(const Vector3d&p1,const Vector3d&p2){
    for (int i=0;i<3;i++){
        double x=p1[i]-p2[i];
        Rational xr=Rational(p1[i])-Rational(p2[i]);
        if(xr>x||xr<x) {
            std::cout<<"double, "<<x<<std::endl;
            std::cout<<"rational, "<<xr<<std::endl;
            if(xr>x){
                std::cout<<"larger"<<std::endl;
            }
            if(xr<x){
                std::cout<<"smaller"<<std::endl;
            }
            std::cout<<"diff,"<<xr-x<<std::endl;
            return false;}
    }
    return true;


}

// x0 is the point, x1, x2, x3 is the triangle
void prism::get_prism_vertices(
    const Vector3d& x0,
    const Vector3d& x1,
    const Vector3d& x2,
    const Vector3d& x3,
    const Vector3d& x0b,
    const Vector3d& x1b,
    const Vector3d& x2b,
    const Vector3d& x3b,
    std::array<Vector3d, 6>& p_vertices)
{
    p_vertices[0] = x0 - x1;
    p_vertices[1] = x0 - x3;
    p_vertices[2] = x0 - x2;
    p_vertices[3] = x0b - x1b;
    p_vertices[4] = x0b - x3b;
    p_vertices[5] = x0b - x2b;
    assert(have_no_truncation(x0,x1));
    assert(have_no_truncation(x0,x3));
    assert(have_no_truncation(x0,x2));
    assert(have_no_truncation(x0b , x1b));
    assert(have_no_truncation(x0b , x3b));
    assert(have_no_truncation(x0b , x2b));

}
prism::prism(
    const Vector3d& vs,
    const Vector3d& fs0,
    const Vector3d& fs1,
    const Vector3d& fs2,
    const Vector3d& ve,
    const Vector3d& fe0,
    const Vector3d& fe1,
    const Vector3d& fe2)
{

   // double k;
    // these are the 6 vertices of the prism,right hand law

    get_prism_vertices(
        vs, fs0, fs1, fs2, ve, fe0, fe1, fe2,
        p_vertices); //  before use this we need to shift all the vertices
    std::array<int, 2> eid;

    eid[0] = 0;
    eid[1] = 1;
    prism_edge_id[0] = eid;
    eid[0] = 1;
    eid[1] = 2;
    prism_edge_id[1] = eid;
    eid[0] = 2;
    eid[1] = 0;
    prism_edge_id[2] = eid;

    eid[0] = 3;
    eid[1] = 4;
    prism_edge_id[3] = eid;
    eid[0] = 4;
    eid[1] = 5;
    prism_edge_id[4] = eid;
    eid[0] = 5;
    eid[1] = 3;
    prism_edge_id[5] = eid;

    eid[0] = 0;
    eid[1] = 3;
    prism_edge_id[6] = eid;
    eid[0] = 1;
    eid[1] = 4;
    prism_edge_id[7] = eid;
    eid[0] = 2;
    eid[1] = 5;
    prism_edge_id[8] = eid;
}
bool prism::is_triangle_degenerated(const int up_or_bottom)
{
    int pid = up_or_bottom == 0 ? 0 : 3;
    const auto to_2d = [](const Vector3d& p, int t) {
        return Vector2d(p[(t + 1) % 3], p[(t + 2) % 3]);
    };
    double r = ((p_vertices[pid] - p_vertices[pid + 1])
                    .cross(p_vertices[pid] - p_vertices[pid + 2]))
                   .norm();
    if (fabs(r) > 1e-8)
        return false;
    int ori;
    std::array<Vector2d, 3> p;
    for (int j = 0; j < 3; j++) {

        p[0] = to_2d(p_vertices[pid], j);
        p[1] = to_2d(p_vertices[pid + 1], j);
        p[2] = to_2d(p_vertices[pid + 2], j);

        ori = orient_2d(p[0], p[1], p[2]);
        if (ori != 0) {
            return false;
        }
    }
    return true;
}

void hex::get_hex_vertices(
    const Vector3d& a0,
    const Vector3d& a1,
    const Vector3d& b0,
    const Vector3d& b1,
    const Vector3d& a0b,
    const Vector3d& a1b,
    const Vector3d& b0b,
    const Vector3d& b1b,
    std::array<Vector3d, 8>& h_vertices)
{
    h_vertices[0] = a0 - b0;
    h_vertices[1] = a1 - b0;
    h_vertices[2] = a1 - b1;
    h_vertices[3] = a0 - b1;
    h_vertices[4] = a0b - b0b;
    h_vertices[5] = a1b - b0b;
    h_vertices[6] = a1b - b1b;
    h_vertices[7] = a0b - b1b;
    assert(have_no_truncation(a0,b0));
    assert(have_no_truncation(a1,b0));
    assert(have_no_truncation(a1,b1));
    assert(have_no_truncation(a0,b1));

    assert(have_no_truncation(a0b,b0b));
    assert(have_no_truncation(a1b,b0b));
    assert(have_no_truncation(a1b,b1b));
    assert(have_no_truncation(a0b,b1b));
}
// a0, a1 is one edge, b0, b1 is another  edge
hex::hex(
    const Vector3d& a0,
    const Vector3d& a1,
    const Vector3d& b0,
    const Vector3d& b1,
    const Vector3d& a0b,
    const Vector3d& a1b,
    const Vector3d& b0b,
    const Vector3d& b1b)
{


    get_hex_vertices(
        a0, a1, b0, b1, a0b, a1b, b0b, b1b,
        h_vertices); //before use this we need to shift all the vertices
    std::array<int, 2> eid;

    eid[0] = 0;
    eid[1] = 1;
    hex_edge_id[0] = eid;
    eid[0] = 1;
    eid[1] = 2;
    hex_edge_id[1] = eid;
    eid[0] = 2;
    eid[1] = 3;
    hex_edge_id[2] = eid;
    eid[0] = 3;
    eid[1] = 0;
    hex_edge_id[3] = eid;

    eid[0] = 4;
    eid[1] = 5;
    hex_edge_id[4] = eid;
    eid[0] = 5;
    eid[1] = 6;
    hex_edge_id[5] = eid;
    eid[0] = 6;
    eid[1] = 7;
    hex_edge_id[6] = eid;
    eid[0] = 7;
    eid[1] = 4;
    hex_edge_id[7] = eid;

    eid[0] = 0;
    eid[1] = 4;
    hex_edge_id[8] = eid;
    eid[0] = 1;
    eid[1] = 5;
    hex_edge_id[9] = eid;
    eid[0] = 2;
    eid[1] = 6;
    hex_edge_id[10] = eid;
    eid[0] = 3;
    eid[1] = 7;
    hex_edge_id[11] = eid;
}

// the facets of the tet are all oriented to outside. check if p is inside of
// OPEN tet
bool is_point_inside_tet(const bilinear& bl, const Vector3d& p)
{

    for (int i = 0; i < 4; i++) { // facets.size()==4
        Vector3d pt1 = bl.v[bl.facets[i][0]], pt2 = bl.v[bl.facets[i][1]],
                 pt3 = bl.v[bl.facets[i][2]];
        if (orient_3d(p, pt1, pt2, pt3) >= 0) {
            return false;
        }
    }
    return true; // all the orientations are -1, then point inside
}

// we already know the bilinear is degenerated, next check which kind
int bilinear_degeneration(const bilinear& bl)
{
    bool dege1 = is_triangle_degenerated(bl.v[0], bl.v[1], bl.v[2]);
    bool dege2 = is_triangle_degenerated(bl.v[0], bl.v[2], bl.v[3]);

    if (dege1 && dege2) {
        return BI_DEGE_PLANE;
    }
    Vector3d p0, p1, p2;

    if (dege1) {
        p0 = bl.v[0];
        p1 = bl.v[2];
        p2 = bl.v[3];
    } else {
        p0 = bl.v[0];
        p1 = bl.v[1];
        p2 = bl.v[2];
    }

    Vector3d np = Vector3d::Random();
    int ori = orient_3d(np, p0, p1, p2);
    while (ori == 0) { // if coplanar, random
        np = Vector3d::Random();
        ori = orient_3d(np, p0, p1, p2);
    }
    int ori0 = orient_3d(np, bl.v[0], bl.v[1], bl.v[2]);
    int ori1 = orient_3d(np, bl.v[0], bl.v[2], bl.v[3]);
    if (ori0 * ori1 <= 0) {
        return BI_DEGE_XOR_02;
    }
    ori0 = orient_3d(np, bl.v[0], bl.v[1], bl.v[3]);
    ori1 = orient_3d(np, bl.v[3], bl.v[1], bl.v[2]);
    if (ori0 * ori1 <= 0) {
        return BI_DEGE_XOR_13;
    }
    return BI_DEGE_PLANE;
}

bool is_cube_intersect_degenerated_bilinear(
    const bilinear& bl, const cube& cube)
{
    int dege = bilinear_degeneration(bl);
    // int axis;
    bool res;
    if (dege == BI_DEGE_PLANE) {

        if (is_cube_edges_intersect_triangle(cube, bl.v[0], bl.v[1], bl.v[3]))
            return true;
        if (is_cube_edges_intersect_triangle(cube, bl.v[3], bl.v[1], bl.v[2]))
            return true;
        return false;
    } else {

        if (dege == BI_DEGE_XOR_02) { // triangle 0-1-2 and 0-2-3
            for (int i = 0; i < 12; i++) {
                res = int_seg_XOR(
                    segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[1], bl.v[2], bl.v[0],
                        true), // CAUTION: need to be careful for the order here
                    segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[3], bl.v[2], bl.v[0], true));
                if (res == true)
                    return true;
            }
            return false;
        }
        if (dege == BI_DEGE_XOR_13) { // triangle 0-1-2 and 0-2-3
            for (int i = 0; i < 12; i++) {
                res = int_seg_XOR(
                    segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[0], bl.v[1], bl.v[3],
                        true), // CAUTION: need to be careful for the order here
                    segment_triangle_intersection(
                        cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                        bl.v[2], bl.v[1], bl.v[3], true));
                if (res == true)
                    return true;
            }
            return false;
        }

    }
    std::cout << "!! THIS CANNOT HAPPEN" << std::endl;
    return false;
}
// phisign gives the pair we want to check
bool line_shoot_same_pair_tet(
    const Vector3d& p0, const Vector3d& p1, const int phisign, bilinear& bl)
{
    int fid;
    if (bl.phi_f[0] == 2)
        get_tet_phi(bl);

    if (phisign > 0) {
        if (bl.phi_f[0] > 0)
            fid = 0;
        else
            fid = 2;
    } else {
        if (bl.phi_f[1] > 0)
            fid = 0;
        else
            fid = 2;
    } // get which pair of facets to check
    // if line is parallel to the triangle, return false
    int inter0 = is_line_cut_triangle(
        p0, p1, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
        bl.v[bl.facets[fid][2]], false);

    int inter1 = is_line_cut_triangle(
        p0, p1, bl.v[bl.facets[fid + 1][0]], bl.v[bl.facets[fid + 1][1]],
        bl.v[bl.facets[fid + 1][2]], false);

    if (inter0 == 1 && inter1 == 1)
        return true;
    return false;
}

Rational quadratic_function_value(
    const Rational& a, const Rational& b, const Rational& c, const Rational& t)
{

    if (t.get_sign() == 0) {
        return c;
    } else {
        return a * t * t + b * t + c;
    }
}
// f(t)=at^2+bt+c, when t is [t0, t1]
bool quadratic_function_rootfinder(
    const Rational& a,
    const Rational& b,
    const Rational& c,
    const Rational t0,
    const Rational t1) // t0 t1 no matter which is bigger
{
    Rational ft0, ft1;
    ft0 = quadratic_function_value(a, b, c, t0);
    ft1 = quadratic_function_value(a, b, c, t1);
    if (ft0.get_sign() == 0 || ft1.get_sign() == 0)
        return true;
    if (ft0.get_sign() != ft1.get_sign())
        return true;
    // the following are cases which signs of endpoints are same
    if (a.get_sign() == 0)
        return false;
    Rational t = -b / (2 * a);
    if (t < t1 && t > t0) {
        Rational ft
            = (4 * a * c - b * b); // actually it should be /4a. we only
                                   // care about the signs so doesnt matter
        int ftsign = a.get_sign() > 0 ? ft.get_sign() : (-1) * ft.get_sign();
        if (ft0.get_sign() != ftsign)
            return true;
    }
    return false;
}
// v0 is one point, dir is the direction,
// x0 - x3 are four points defined the shape
void get_quadratic_function(
    const Rational& v00,
    const Rational& v01,
    const Rational& v02,
    const Rational& dir0,
    const Rational& dir1,
    const Rational& dir2,
    const Vector3d x0d,
    const Vector3d x1d,
    const Vector3d x2d,
    const Vector3d x3d,
    Rational& a,
    Rational& b,
    Rational& c)
{
    Vector3r x0(x0d[0], x0d[1], x0d[2]), x1(x1d[0], x1d[1], x1d[2]),
        x2(x2d[0], x2d[1], x2d[2]), x3(x3d[0], x3d[1], x3d[2]);
    Rational x00 = x0[0], x01 = x0[1], x02 = x0[2];
    Rational x10 = x1[0], x11 = x1[1], x12 = x1[2];
    Rational x20 = x2[0], x21 = x2[1], x22 = x2[2];
    Rational x30 = x3[0], x31 = x3[1], x32 = x3[2];

    Rational x101 = x11 - x01, x100 = x10 - x00, x102 = x12 - x02,
             x201 = x21 - x01, x200 = x20 - x00, x300 = x30 - x00,
             x301 = x31 - x01, x310 = x30 - x10, x211 = x21 - x11,
             x311 = x31 - x11, x210 = x20 - x10, x202 = x22 - x02,
             x302 = x32 - x02, x212 = x22 - x12, x312 = x32 - x12;
    a = (dir0 * (x101 * x202 - x102 * x201)
         + dir1 * (-x100 * x202 + x102 * x200)
         + dir2 * (x100 * x201 - x101 * x200))
            * (dir0 * (-x201 * x302 + x202 * x301)
               + dir1 * (x200 * x302 - x202 * x300)
               + dir2 * (-x200 * x301 + x201 * x300))
        - (dir0 * (-x211 * x312 + x212 * x311)
           + dir1 * (x210 * x312 - x212 * x310)
           + dir2 * (-x210 * x311 + x211 * x310))
            * (dir0 * (x101 * x302 - x102 * x301)
               + dir1 * (-x100 * x302 + x102 * x300)
               + dir2 * (x100 * x301 - x101 * x300));
    b = ((v00 - x00) * (x101 * x202 - x102 * x201)
         + (v01 - x01) * (-x100 * x202 + x102 * x200)
         + (v02 - x02) * (x100 * x201 - x101 * x200))
            * (dir0 * (-x201 * x302 + x202 * x301)
               + dir1 * (x200 * x302 - x202 * x300)
               + dir2 * (-x200 * x301 + x201 * x300))
        + (dir0 * (x101 * x202 - x102 * x201)
           + dir1 * (-x100 * x202 + x102 * x200)
           + dir2 * (x100 * x201 - x101 * x200))
            * ((v00 - x00) * (-x201 * x302 + x202 * x301)
               + (v01 - x01) * (x200 * x302 - x202 * x300)
               + (v02 - x02) * (-x200 * x301 + x201 * x300))
        - ((v00 - x10) * (-x211 * x312 + x212 * x311)
           + (v01 - x11) * (x210 * x312 - x212 * x310)
           + (v02 - x12) * (-x210 * x311 + x211 * x310))
            * (dir0 * (x101 * x302 - x102 * x301)
               + dir1 * (-x100 * x302 + x102 * x300)
               + dir2 * (x100 * x301 - x101 * x300))
        - (dir0 * (-x211 * x312 + x212 * x311)
           + dir1 * (x210 * x312 - x212 * x310)
           + dir2 * (-x210 * x311 + x211 * x310))
            * ((v00 - x00) * (x101 * x302 - x102 * x301)
               + (v01 - x01) * (-x100 * x302 + x102 * x300)
               + (v02 - x02) * (x100 * x301 - x101 * x300));
    c = ((v00 - x00) * (x101 * x202 - x102 * x201)
         + (v01 - x01) * (-x100 * x202 + x102 * x200)
         + (v02 - x02) * (x100 * x201 - x101 * x200))
            * ((v00 - x00) * (-x201 * x302 + x202 * x301)
               + (v01 - x01) * (x200 * x302 - x202 * x300)
               + (v02 - x02) * (-x200 * x301 + x201 * x300))
        - ((v00 - x10) * (-x211 * x312 + x212 * x311)
           + (v01 - x11) * (x210 * x312 - x212 * x310)
           + (v02 - x12) * (-x210 * x311 + x211 * x310))
            * ((v00 - x00) * (x101 * x302 - x102 * x301)
               + (v01 - x01) * (-x100 * x302 + x102 * x300)
               + (v02 - x02) * (x100 * x301 - x101 * x300));
}
bool get_function_find_root(
    const bilinear& bl,
    const Vector3r& p0,
    const Vector3r& p1,
    const Rational& t0,
    const Rational& t1)
{
    if (t0 > 1 || t0 < 0)
        std::cout << "t is not right: exceed the limit" << std::endl;
    if (t1 > 1 || t1 < 0)
        std::cout << "t is not right: exceed the limit" << std::endl;
    Rational a, b, c;

    Vector3r dir = p1 - p0;
    get_quadratic_function(
        Rational(p0[0]), Rational(p0[1]), Rational(p0[2]), Rational(dir[0]),
        Rational(dir[1]), Rational(dir[2]), bl.v[0], bl.v[1], bl.v[2], bl.v[3],
        a, b, c);
    bool res = quadratic_function_rootfinder(a, b, c, t0, t1);

    return res;
}
void print_sub()
{
    //std::cout<<"time of rootfinder "<<rftime<<std::endl;
}
double root_finder_time(){
    return 0;
}
bool rootfinder(
    const bilinear& bl,
    const Vector3d& p0d,
    const Vector3d& p1d,
    const bool p0in,
    const bool p1in,
    const int pairid)
{

    Vector3r p0(p0d[0], p0d[1], p0d[2]), p1(p1d[0], p1d[1], p1d[2]);

    if (p0in && p1in) {
        // t0=0, t1=1
        return get_function_find_root(bl, p0, p1, Rational(0), Rational(1));
    }
    int fid = 2 * pairid; // it should be 0 or 2
    Rational t;
    if (p0in) {
        bool res1 = seg_triangle_inter_return_t(
            p0d, p1d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
            bl.v[bl.facets[fid][2]], t);
        if (!res1) {
            res1 = seg_triangle_inter_return_t(
                p0d, p1d, bl.v[bl.facets[fid + 1][0]],
                bl.v[bl.facets[fid + 1][1]], bl.v[bl.facets[fid + 1][2]], t);
        }
        if (!res1)
            return false; // means not really intersected
        // we got t
        return get_function_find_root(bl, p0, p1, 0, t);
    }
    if (p1in) { // change the order of input to get t just because we want
                // domain to be [0, t]
        bool res1 = seg_triangle_inter_return_t( // here get n1, d1, n2, d2
            p1d, p0d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
            bl.v[bl.facets[fid][2]], t);
        if (!res1) {
            res1 = seg_triangle_inter_return_t(
                p1d, p0d, bl.v[bl.facets[fid + 1][0]],
                bl.v[bl.facets[fid + 1][1]], bl.v[bl.facets[fid + 1][2]], t);
        }
        if (!res1)
            return false; // means not really intersected
        // we got t
        return get_function_find_root(bl, p1, p0, 0, t);
    }

    Rational t1;
    bool res1 = seg_triangle_inter_return_t(
        p0d, p1d, bl.v[bl.facets[fid][0]], bl.v[bl.facets[fid][1]],
        bl.v[bl.facets[fid][2]], t);
    if (!res1)
        return false; // means not really intersected
    bool res2 = seg_triangle_inter_return_t(
        p0d, p1d, bl.v[bl.facets[fid + 1][0]], bl.v[bl.facets[fid + 1][1]],
        bl.v[bl.facets[fid + 1][2]], t1);
    if (!res2)
        return false; // means not really intersected
    // we got t, t1
    return get_function_find_root(bl, p0, p1, t, t1);
}
// segment intersect two opposite faces not included. compare phi, then use root
// finder
bool is_seg_intersect_not_degenerated_bilinear(
    bilinear& bl,
    const Vector3d& p0,
    const Vector3d& p1,
    const bool pin0,
    const bool pin1)
{
    // first compare phi, if phis are different, intersected;
    // then check if the line intersect two opposite facets of bilinear, if so,
    // use rootfinder, else, not intersected
    
    if (pin0 && pin1) { // two points are all inside

        Rational phi0 = phi(p0, bl.v);
        Rational phi1 = phi(p1, bl.v);
        if (phi0 == 0 || phi1 == 0 || phi0.get_sign() != phi1.get_sign())
            return true;
        if (line_shoot_same_pair_tet(p0, p1, phi1.get_sign(), bl)) {
            if (phi1.get_sign() == bl.phi_f[0]){
                
                bool rf=rootfinder(bl, p0, p1, pin0, pin1, 0);
                
                return rf;
            }
                
            else{
                
                bool rf=rootfinder(bl, p0, p1, pin0, pin1, 1);
                
                return rf;
            }
                
        }

        else
            return false; // if the phis are the same, and shoot same pair, need
                          // to use rootfinder
    }
    if (pin0) {
        Rational phi0 = phi(p0, bl.v);
        if (phi0 == 0)
            return true;
        int hitpair = -1;

        for (int i = 0; i < 4; i++) {
            if (segment_triangle_intersection( // 0,1,2,3. 1,2,3 are all
                                               // intersected
                    p0, p1, bl.v[bl.facets[i][0]], bl.v[bl.facets[i][1]],
                    bl.v[bl.facets[i][2]], false)
                > 0) {
                if (i < 2)
                    hitpair = 0;
                else
                    hitpair = 1;
                break;
            }
        }
        if (bl.phi_f[0] == 2)
            get_tet_phi(bl);
        if (hitpair == -1)
            return false;
        if (phi0.get_sign()
            != bl.phi_f[hitpair]) { // if different, intersected, if same,
                                    // extend; if shoot same, rootfinder, if
                                    // shoot diff, false
            return true;
        }

        else {
            if (line_shoot_same_pair_tet(p0, p1, phi0.get_sign(), bl)) {
                if (phi0.get_sign() == bl.phi_f[0]){
                    
                    bool rf= rootfinder(bl, p0, p1, pin0, pin1, 0);
                    
                    return rf;
                }
                    
                else{
                    
                    bool rf=rootfinder(bl, p0, p1, pin0, pin1, 1);
                    
                    return rf;
                }
                     
            }

            else
                return false; // if the phis are the same, and shoot same pair,
                              // need to use rootfinder
        }
    }
    if (pin1) {
        Rational phi1 = phi(p1, bl.v);
        if (phi1 == 0)
            return true;
        int hitpair = -1;
        for (int i = 0; i < 4; i++) {
            if (segment_triangle_intersection( // 0,1,2,3. 1,2,3 are all
                                               // intersected
                    p0, p1, bl.v[bl.facets[i][0]], bl.v[bl.facets[i][1]],
                    bl.v[bl.facets[i][2]], false)
                > 0) {
                if (i < 2)
                    hitpair = 0;
                else
                    hitpair = 1;
                break;
            }
        }
        if (bl.phi_f[0] == 2)
            get_tet_phi(bl);
        if (hitpair == -1)
            return false; // parallel , should be impossible
        if (phi1.get_sign() != bl.phi_f[hitpair]) {
            return true;
        }

        else {
            if (line_shoot_same_pair_tet(p0, p1, phi1.get_sign(), bl)) {
                if (phi1.get_sign() == bl.phi_f[0]){
                    
                    bool rf=rootfinder(bl, p0, p1, pin0, pin1, 0);
                    
                    return rf;
                }

                else{
                    
                    bool rf=  rootfinder(bl, p0, p1, pin0, pin1, 1);
                    
                    return rf;
                }
                   
            } else
                return false; // if the phis are the same, and shoot same pair,
                              // need to use rootfinder
        }
    }
    if (!pin0
        && !pin1) { // not intersect tet (false), or intersect same side(root
                    // finder) or intersect diff side(checked before)

        if (line_shoot_same_pair_tet(p0, p1, 1, bl)) {
            if (1 == bl.phi_f[0]){
                
                bool rf= rootfinder(bl, p0, p1, pin0, pin1, 0);
                
                return rf;
            }
                
            else{
                
                bool rf= rootfinder(bl, p0, p1, pin0, pin1, 1);
                
                return rf;
            }
                
        }

        else if (line_shoot_same_pair_tet(p0, p1, -1, bl)) {
            if (-1 == bl.phi_f[0]){
                
                bool rf= rootfinder(bl, p0, p1, pin0, pin1, 0);
                
                return rf;
            }
                
            else{
                
                bool rf= rootfinder(bl, p0, p1, pin0, pin1, 1);
                
                return rf;
            }
                
        }
        return false; // if the phis are the same, and shoot same pair,
                      // need to use rootfinder
    }
    std::cout
        << " it cannot happen here in is_seg_intersect_not_degenerated_bilinear"
        << std::endl;
    return false;
}

bool is_cube_edge_intersect_bilinear(
    bilinear& bl, const cube& cb, const std::array<bool, 8>& pin)
{
    if (bl.is_degenerated)
        return false; // we already checked degenerated cases
    for (int i = 0; i < 12; i++) {
        if (is_seg_intersect_not_degenerated_bilinear(
                bl, cb.vr[cb.edgeid[i][0]], cb.vr[cb.edgeid[i][1]],
                pin[cb.edgeid[i][0]], pin[cb.edgeid[i][1]]))
            return true;
    }
    return false;
}
// vin is true, this vertex has intersection with open tet
// if tet is degenerated, just tell us if cube is intersected with the shape
bool is_cube_intersect_tet_opposite_faces(
    const bilinear& bl,const Vector3d &pmin, const Vector3d &pmax,
    const cube& cube,
    std::array<bool, 8>& vin,
    bool& cube_inter_tet)
{

    cube_inter_tet = false;
    if (!bl.is_degenerated) {
        for (int i = 0; i < 8; i++) {
            vin[i] = false;

            if (is_point_inside_tet(bl, cube.vr[i])) {
                cube_inter_tet = true;
                vin[i] = true;
            }

        }
    } else {

        bool rst = is_cube_intersect_degenerated_bilinear(bl, cube);

        return rst;
    }

    bool side1 = false;
    bool side2 = false;
    Vector3d vmin,vmax;
    for (int i = 0; i < 12; i++) {
       
        if (vin[cube.edgeid[i][0]]
            && vin[cube.edgeid[i][1]]) { // if two vertices are all inside, it
                                         // can not cut any edge
            cube_inter_tet = true;
            continue;
        }
        get_edge_coners(cube.vr[cube.edgeid[i][0]],cube.vr[cube.edgeid[i][1]],vmin,vmax);
        if(!box_box_intersection(pmin,pmax,vmin,vmax)){
            continue;
        }
        for (int j = 0; j < 4; j++) {

            int inter = segment_triangle_intersection(
                cube.vr[cube.edgeid[i][0]], cube.vr[cube.edgeid[i][1]],
                bl.v[bl.facets[j][0]], bl.v[bl.facets[j][1]],
                bl.v[bl.facets[j][2]], false);

            if (inter > 0) {
                cube_inter_tet = true;
                if (j == 0 || j == 1)
                    side1 = true;
                if (j == 2 || j == 3)
                    side2 = true;
                if (side1 && side2)
                    return true;
            }
        }
    }
    if (side1 && side2)
        return true;
    return false;
}
bool cube_discrete_bilinear_intersection(
    const cube& cb, const bilinear& bl, int n)
{
    Vector3d s0, s1;
    for (int i = 0; i < 12; i++) {
        s0 = cb.vr[cb.edgeid[i][0]];
        s1 = cb.vr[cb.edgeid[i][1]];
        if (seg_discrete_bilinear_intersection(bl, n, s0, s1)) // return 0,1,2

            return true;
    }
    return false;
}

} // namespace doubleccd
