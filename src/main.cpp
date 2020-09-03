#include <iostream>
//
//#include <ccd.hpp>
#include <vector>
//#include<Utils.hpp>
//#include<subfunctions.h>
#include <CCD/ccd.hpp>
#include <CCD/exact_subtraction.hpp>
#include <array>
#include <doubleCCD/doubleccd.hpp>
#include <doubleCCD/hack.h>
#include <fstream>
#include "read_collision_data.hpp"
#include<sstream>
#include<interval_ccd/interval_ccd.hpp>
#include<interval_ccd/interval_root_finder.hpp>

//#include <predicates/indirect_predicates.h>
//#include <exact_subtraction.hpp>
//#include<subfunctions.h>
using namespace doubleccd;
using namespace std;

std::string root_path = CCD_DATA_PATH;
#ifdef WIN32
std::string path_sep = "\\";
#else
std::string path_sep = "/";
#endif

struct sccd {
    Vector3d pts;
    Vector3d pte;
    Vector3d v1s;
    Vector3d v2s;
    Vector3d v3s;
    Vector3d v1e;
    Vector3d v2e;
    Vector3d v3e;
};
Vector3d construct1(double a, double b, double c) { return Vector3d(a, b, c); }
//
bool read_CSV(const string inputFileName, vector<sccd>& data)
{

    ifstream infile;
    infile.open(inputFileName);
    if (!infile.is_open()) {
        cout << "Path Wrong!!!!" << endl;
        return false;
    }

    int l = 0;
    while (infile) // there is input overload classfile
    {
        l++;
        string s;
        if (!getline(infile, s))
            break;
        if (s[0] != '#') {
            istringstream ss(s);
            array<double, 24> record;
            int c = 0;
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record[c] = stod(line);
                    c++;
                } catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line "
                         << l << endl;
                    e.what();
                }
            }
            sccd tmp;
            tmp.pts = construct1(record[0], record[1], record[2]);
            tmp.v1s = construct1(record[3], record[4], record[5]);
            tmp.v2s = construct1(record[6], record[7], record[8]);
            tmp.v3s = construct1(record[9], record[10], record[11]);
            tmp.pte = construct1(record[12], record[13], record[14]);
            tmp.v1e = construct1(record[15], record[16], record[17]);
            tmp.v2e = construct1(record[18], record[19], record[20]);
            tmp.v3e = construct1(record[21], record[22], record[23]);
            data.emplace_back(tmp);
        }
    }
    if (!infile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
    }
    cout << "data size " << data.size() << endl;
    return true;
}
bool read_result(const string inputFileName, vector<bool>& data)
{

    ifstream infile;
    infile.open(inputFileName);
    if (!infile.is_open()) {
        cout << "Path Wrong!!!!" << endl;
        return false;
    }

    int l = 0;
    while (infile) // there is input overload classfile
    {
        l++;
        string s;
        if (!getline(infile, s))
            break;
        if (s[0] != '#') {
            istringstream ss(s);
            bool record;
            int c = 0;
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record = stod(line);
                    c++;
                } catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line "
                         << l << endl;
                    e.what();
                }
            }

            data.push_back(record);
        }
    }
    if (!infile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
    }
    cout << "data size " << data.size() << endl;
    return true;
}

void test1()
{
    std::vector<sccd> data;
    read_CSV(root_path + path_sep + "cow-head-collisions.csv", data);
    std::vector<bool> results;
    int fn = data.size(); // 50000;
    results.resize(fn);
    double k;
    int count = 0;
    double maxvalue = 0, merror = 0, maxerror;

    bool correct;
    for (int i = 0; i < fn; i++) {
        if (i % 200 == 0)
            std::cout << "i " << i << std::endl;

        /*get_prism_vertices_double(
            data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s, data[i].pte,
            data[i].v1e, data[i].v2e, data[i].v3e, k, correct, maxerror);
        if (k != 0)
            count++;
        if (fabs(k) > maxvalue)
            maxvalue = fabs(k);
        if (fabs(k) > 0.9)
            std::cout << "big k " << k << " i " << i << std::endl;
        if (correct == false)
            cout << "wrong minus, " << i << std::endl;
        if (maxerror > merror) {
            merror = maxerror;
        }*/
    }
    // cout << "percentage of using k " << double(count) / double(fn) <<
    // std::endl; cout << "max of k " << maxvalue << std::endl; cout << "max
    // error " << merror << std::endl;

    // std::cout << "Exact CCD" << std::endl;
}
void test_compare()
{
    std::vector<sccd> data;
    read_CSV(root_path + path_sep + "cow-head-collisions.csv", data);
    vector<bool> rst;
    read_result(root_path + path_sep + "result_all.csv", rst);
    int rst_true = 0;
    for (int i = 0; i < rst.size(); i++) {
        if (rst[i])
            rst_true++;
    }
    std::cout << "original collision nbr, " << rst_true << std::endl;
    std::vector<bool> results;
    int fn = data.size(); // 50000;
    results.resize(fn);
    int inside = 0;
    for (int i = 244; i < 245; i++) {
        if (i % 200 == 0)
            std::cout << "i " << i << std::endl;
        // std::cout << "i " << std::endl;
        /*results[i] = vertexFaceCCD(
                data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
                data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-8);*/
        results[i] = ccd::vertexFaceCCD(
            data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s, data[i].pte,
            data[i].v1e, data[i].v2e, data[i].v3e, 1e-8);
        if (rst[i] == 1
            && results[i] == 0) { // when old method says yes but we missed it
            std::cout << "wrong! i= " << i << std::endl;
        }
        // std::cout << "result " << results[i] << std::endl;
        if (results[i] == true)
            inside++;
    }
    std::cout << int_seg_XOR(1, 0) << std::endl;
    cube cb(0.1);
    std::cout << "inside number " << inside << std::endl;
}

void test_shifted_compare()
{
    std::vector<sccd> data;
    read_CSV(root_path + path_sep + "cow-head-collisions.csv", data);
    vector<bool> rst;
    read_result(root_path + path_sep + "result_all.csv", rst);
    int rst_true = 0;
    for (int i = 0; i < rst.size(); i++) {
        if (rst[i])
            rst_true++;
    }
    std::cout << "original collision nbr, " << rst_true << std::endl;
    std::vector<bool> results, results1;
    int fn = data.size(); // 50000;
    results.resize(fn);
    results1.resize(fn);
    int inside = 0;
    igl::Timer timer;
    double time = 0;
    for (int i = 0; i < fn; i++) {
        if (i % 2000 == 0)
            std::cout << "i " << i << std::endl;
        // std::cout << "i " << std::endl;
        timer.start();
        results[i] = edgeEdgeCCD( // double
            data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s, data[i].pte,
            data[i].v1e, data[i].v2e, data[i].v3e, 1e-30);
            if(results[1]==0){
               // std::cout<<"this case not collision "<<i<<std::endl;
               // std::cout<<"x0,\n "<<data[i].pts<<std::endl;std::cout<<"x1,\n "<<data[i].v1s<<std::endl;std::cout<<"x2,\n "<<data[i].v2s<<std::endl;std::cout<<"x3,\n "<<data[i].v3s<<std::endl;
              //  std::cout<<"x0b,\n "<<data[i].pte<<std::endl;std::cout<<"x1b,\n "<<data[i].v1e<<std::endl;std::cout<<"x2b,\n "<<data[i].v2e<<std::endl;std::cout<<"x3b,\n "<<data[i].v3e<<std::endl;

            }
        time += timer.getElapsedTimeInSec();
        // results1[i] = ccd::vertexFaceCCD(//rational
        //	data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
        //	data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-3);

        // std::cout << "Rational vs double " << results1[i] << " " <<
        // results[i] << std::endl; if (results1[i] != results[i]) {//when old
        // method says yes but we missed it 	std::cout << "double don't match
        // rational! i= " << i << std::endl; 	std::cout << "Rational vs double
        // "
        //<< results1[i]<<" "<<results[i] << std::endl;
        //}
        // std::cout << "result " << results[i] << std::endl;
        if (results[i] == true)
            inside++;
    }

    std::cout << "total time " << time << std::endl;
    std::cout << "inside number " << inside << std::endl;
    test();

    print_sub();
    ray_time();
    double phit=print_phi_time();
    std::cout<<"percentage "<<phit/time<<"\n"<<std::endl;
}

void test_rootfinder()
{
    Vector3d p0(0, 0, 1);
    Vector3d p1(0, 1, 1);
    Vector3d p2(0, 0, 0);
    Vector3d p3(1, 1, 0);
    bilinear bl(p0, p1, p2, p3);
    if (bl.phi_f[0] == 2) {
        std::cout << "calculate bilinear face signs" << std::endl;
        get_tet_phi(bl);
    }

    Vector3d s0(0.01, 0.02, 2);
    Vector3d s1(0.01, 0.02, -1);
    bool pin0 = false, pin1 = false, result;
    if (line_shoot_same_pair_tet(s0, s1, 1, bl)) {
        std::cout << "shoot first pair" << std::endl;
        if (1 == bl.phi_f[0])
            result = rootfinder(bl, s0, s1, pin0, pin1, 0);
        else
            result = rootfinder(bl, s0, s1, pin0, pin1, 1);
    }

    else if (line_shoot_same_pair_tet(s0, s1, -1, bl)) {
        std::cout << "shoot second pair" << std::endl;
        if (-1 == bl.phi_f[0])
            result = rootfinder(bl, s0, s1, pin0, pin1, 0);
        else
            result = rootfinder(bl, s0, s1, pin0, pin1, 1);
    }
    std::cout << "result of rootfinder is shoot or not " << result << std::endl;

    std::vector<std::array<Vector3r, 3>> tris;
    tri_bilinear(bl, 10, tris);
    ccd::Vector3r e0(0.01, 0.02, 2), e1(0.01, 0.02, -1);
    /*for (int i = 0; i < tris.size(); i++) {
            if (ccd::segment_triangle_intersection(Vector3r(0.01, 0.02, 2),
    Vector3r(0.01, 0.02, -1), tris[i][0], tris[i][1], tris[i][2], false)>0) {
                    std::cout << "here seg intersected bilinear" << std::endl;
                    break;
            }
    }*/
    save_obj("D:\\vs\\collision\\test_rootfinder.obj", tris);
}

void testing_difference(Vector3d p, Vector3d q){
    double x,y,z;    
    x=(q[0]-p[0]);
    y=(q[1]-p[1]);
    z=(q[2]-p[2]);
    Rational xr=Rational(q[0])-Rational(p[0]);
    Rational yr=Rational(q[1])-Rational(p[1]);
    Rational zr=Rational(q[2])-Rational(p[2]);
    

    if(xr!=x) std::cout<<"this x is not 0"<<std::endl;
    if(yr!=y) std::cout<<"this y is not 0"<<std::endl;
    if(zr!=z) std::cout<<"this z is not 0"<<std::endl;
}
void test_shift_maxerror()
{
    std::vector<sccd> data;
    Eigen::MatrixX3d vertices, vertices1, vertices2,vertices3;
    read_CSV(root_path + path_sep + "cow-head-collisions.csv", data);
    std::vector<vf_pair> vfs, shifted_vfs;
    std::vector<ee_pair> ees, shifted_ees;
    vfs.resize(data.size());
    vertices.resize(8 * data.size(), 3);
    for (int i = 0; i < data.size(); i++) {
        vfs[i].x0 = data[i].pts;
        vfs[i].x1 = data[i].v1s;
        vfs[i].x2 = data[i].v2s;
        vfs[i].x3 = data[i].v3s;

        vfs[i].x0b = data[i].pte;
        vfs[i].x1b = data[i].v1e;
        vfs[i].x2b = data[i].v2e;
        vfs[i].x3b = data[i].v3e;
        vertices.row(8 * i + 0) = data[i].pts;
        vertices.row(8 * i + 1) = data[i].pte;
        vertices.row(8 * i + 2) = data[i].v1s;
        vertices.row(8 * i + 3) = data[i].v1e;
        vertices.row(8 * i + 4) = data[i].v2s;
        vertices.row(8 * i + 5) = data[i].v2e;
        vertices.row(8 * i + 6) = data[i].v3s;
        vertices.row(8 * i + 7) = data[i].v3e;
    }
    vertices1 = vertices;
    vertices2 = vertices;
    vertices3=vertices;
    
    get_whole_mesh_shifted(vfs, ees, shifted_vfs, shifted_ees, vertices1);
    get_whole_mesh_shifted(vfs, ees, vertices2);

    // shift accroding to the bounding box
    Vector3d pmin, pmax;
    get_corners(vertices,pmin,pmax);
    get_whole_mesh_shifted(vertices3,pmin,pmax);
    compare_whole_mesh_err(vertices3,vertices);
    int count=0;
    for(int i=0;i<vertices3.rows();i++){
        for(int j=0;j<vertices3.rows();j++){
            count++;
            if (count>=1000000) return;
            Vector3d p,q;
            for(int k=0;k<3;k++){
                p[k]=vertices3(i,k);
                q[k]=vertices3(j,k);
            }
            testing_difference(p,q);
        }
    }

}
void read_vf_data(const string file, std::vector<vf_pair>& vfdata) {
	std::vector<Eigen::Matrix<double, 8, 3>> vertex_face_data;
	read_vertex_face_data(file, vertex_face_data);
	vfdata.resize(vertex_face_data.size());
	for (int i = 0; i < vertex_face_data.size(); i++) {
		for (int j = 0; j < 3; j++) {
			vfdata[i].x0[j] = vertex_face_data[i](0, j);
			vfdata[i].x1[j] = vertex_face_data[i](1, j);
			vfdata[i].x2[j] = vertex_face_data[i](2, j);
			vfdata[i].x3[j] = vertex_face_data[i](3, j);

			vfdata[i].x0b[j] = vertex_face_data[i](4, j);
			vfdata[i].x1b[j] = vertex_face_data[i](5, j);
			vfdata[i].x2b[j] = vertex_face_data[i](6, j);
			vfdata[i].x3b[j] = vertex_face_data[i](7, j);
		}

	}
}

void read_ee_data(const string file, std::vector<ee_pair>& eedata) {
	std::vector<Eigen::Matrix<double, 8, 3>> edge_edge_data;
	read_edge_edge_data(file, edge_edge_data);
	eedata.resize(edge_edge_data.size());
	for (int i = 0; i < edge_edge_data.size(); i++) {
		for (int j = 0; j < 3; j++) {
			eedata[i].a0[j] = edge_edge_data[i](0, j);
			eedata[i].a1[j] = edge_edge_data[i](1, j);
			eedata[i].b0[j] = edge_edge_data[i](2, j);
			eedata[i].b1[j] = edge_edge_data[i](3, j);

			eedata[i].a0b[j] = edge_edge_data[i](4, j);
			eedata[i].a1b[j] = edge_edge_data[i](5, j);
			eedata[i].b0b[j] = edge_edge_data[i](6, j);
			eedata[i].b1b[j] = edge_edge_data[i](7, j);
		}

	}
}
void test_more() {

}
void test_edge_edge(){
    string filename="/home/bw1760/scratch/cow-head/edge-edge/edge-edge-collisions-part010.hdf5";
    std::vector<ee_pair> eedata;
	read_ee_data(filename,eedata);

    std::vector<bool> results, results1;
    int fn = eedata.size(); // 50000;
    std::cout<<"fn "<<fn<<std::endl;
    results.resize(fn);
    results1.resize(fn);
    int inside = 0;
    igl::Timer timer;
    double time = 0;
    for (int i = 0; i < fn; i++) {
        if (i % 200 == 0)
            std::cout << "i " << i << std::endl;
        // std::cout << "i " << std::endl;
        timer.start();
        results[i] = edgeEdgeCCD( // double
            eedata[i].a0, eedata[i].a1, eedata[i].b0, eedata[i].b1,
            eedata[i].a0b, eedata[i].a1b, eedata[i].b0b, eedata[i].b1b, 1e-3);
        time += timer.getElapsedTimeInSec();
        // results1[i] = ccd::vertexFaceCCD(//rational
        //	data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
        //	data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-3);

        // std::cout << "Rational vs double " << results1[i] << " " <<
        // results[i] << std::endl; if (results1[i] != results[i]) {//when old
        // method says yes but we missed it 	std::cout << "double don't match
        // rational! i= " << i << std::endl; 	std::cout << "Rational vs double
        // "
        //<< results1[i]<<" "<<results[i] << std::endl;
        //}
        // std::cout << "result " << results[i] << std::endl;
        if (results[i] == true)
            inside++;
    }
    std::cout<<"edge edge total time "<<time<<std::endl;
    std::cout<<"collision number "<<inside<<" out of total number "<<fn<<std::endl;
    test();
    print_sub();
}

void check_false(){
    H5Easy::File file(root_path + path_sep +"vertex-face-collisions.hdf5");
    //H5Easy::File file("/home/zachary/Development/ccd-queries/erleben-cube-internal-edges/edge-edge/edge-edge-collisions.hdf5");
    Eigen::Matrix<double, 8, 3> vertex_face_data;
    string test_case= "/edge_edge_0000416/shifted/points";
    vertex_face_data=H5Easy::load<Eigen::Matrix<double, 8, 3>>(file,test_case);
    std::cout<<test_case<<std::endl;
    vf_pair dt;
    for (int j = 0; j < 3; j++) {
			dt.x0[j] = vertex_face_data(0, j);
			dt.x1[j] = vertex_face_data(1, j);
			dt.x2[j] = vertex_face_data(2, j);
			dt.x3[j] = vertex_face_data(3, j);

			dt.x0b[j] = vertex_face_data(4, j);
			dt.x1b[j] = vertex_face_data(5, j);
			dt.x2b[j] = vertex_face_data(6, j);
			dt.x3b[j] = vertex_face_data(7, j);
	}
    // dt.x0=Vector3d(0.1,0.1,1);
    // dt.x1=Vector3d(0,0,0);
    // dt.x2=Vector3d(0, -0.5, 0);
    // dt.x3=Vector3d(-0.5, -0.5, 0);
    // dt.x0b=Vector3d(0.5,0.5,0.5);
    // dt.x1b=Vector3d(0.4, 0.4, 0.5);
    // dt.x2b=Vector3d(0.4, -0.1, 0.5);
    // dt.x3b=Vector3d(-0.1,-0.1,0.5);

    // dt.x0=Vector3d(0.1,0.1,1);
    // dt.x1=Vector3d(-0.05,-0.05,-0.05);
    // dt.x2=Vector3d(0, -0.9, 0);
    // dt.x3=Vector3d(-0.8, -0.8, 0);
    // dt.x0b=Vector3d(0.5,0.5,0.5);
    // dt.x1b=Vector3d(0.3, 0.3, 0.5);
    // dt.x2b=Vector3d(0.4, -0.4, 0.5);
    // dt.x3b=Vector3d(-0.2,0,0.5);

    double ms=1e-8;
    std::cout<<"ms "<<ms<<std::endl;

    std::cout<<"\n*method double"<<std::endl;
    int r2=doubleccd::edgeEdgeCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,ms);
    //std::cout<<" exact dir "<<doubleccd::hack::getInstance().dir[0]<<" "<< doubleccd::hack::getInstance().dir[1]<<" "<< doubleccd::hack::getInstance().dir[2]<<std::endl;
    std::cout<<"*\n\nmethod rational"<<std::endl;
    int r1=ccd::edgeEdgeCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,ms);


    std::cout<<"the double vf ccd result is "<<r2<<std::endl;
    std::cout<<"the rational vf ccd result is "<<r1<<std::endl;
    //std::cout<<"ori1 "<< orient_3d(dt.x0,dt.x1,dt.x2,dt.x3)<<std::endl;
    //std::cout<<"ori2 "<< orient_3d(dt.x0b,dt.x1b,dt.x2b,dt.x3b)<<std::endl;
       // if(dt.x0[0]==dt.x0b[0]&&dt.x0[2]==dt.x0b[2]&&dt.x0[1]==dt.x0b[1]) std::cout<<"x0, x1 same point"<<std::endl;
    // if(dt.x1[0]==dt.x0[0]&&dt.x1[2]==dt.x0[2])
    //     if(dt.x1[0]==dt.x1b[0]&&dt.x1[2]==dt.x1b[2])
    //         if(dt.x1[1]>dt.x0[1]&&dt.x1b[1]<dt.x0[1])
    //             std::cout<<"point x1 hit x0 some time"<<std::endl;
    //std::cout<<"x0,\n "<<dt.x0<<std::endl;std::cout<<"x1,\n "<<dt.x1<<std::endl;std::cout<<"x2,\n "<<dt.x2<<std::endl;std::cout<<"x3,\n "<<dt.x3<<std::endl;
   //std::cout<<"x0b,\n "<<dt.x0b<<std::endl;std::cout<<"x1b,\n "<<dt.x1b<<std::endl;std::cout<<"x2b,\n "<<dt.x2b<<std::endl;std::cout<<"x3b,\n "<<dt.x3b<<std::endl;
    test();

    print_sub();
    ray_time();
}
void case_check(){
    vf_pair dt;
    
// const Eigen::Vector3d a0s(-30022200, 2362580, 165247);
//     const Eigen::Vector3d a1s(-32347850, 8312380, -1151003);
//     const Eigen::Vector3d a0e(-28995600, 345838, 638580);
//     const Eigen::Vector3d a1e(-31716930, 6104858, -713340);
//     const Eigen::Vector3d b0s(-30319900, 3148750, 0);
//     const Eigen::Vector3d b1s(-28548800, 900349, 0);
//     const Eigen::Vector3d b0e(-30319900, 3148750, 0);
//     const Eigen::Vector3d b1e(-28548800, 900349, 0);

    const Eigen::Vector3d a0s(5.55495e-11,    0.160014,    0.914204);
    const Eigen::Vector3d a1s(0.098607, 0.160014, 0.898586);
    const Eigen::Vector3d a0e(-1.02941e-07     ,0.470588          ,0.6);
    const Eigen::Vector3d a1e(0.19509, -0.529412,   1.58079 );
    const Eigen::Vector3d b0s(1.42931e-05 ,   0.159178,    0.914963);
    const Eigen::Vector3d b1s(0.098703, 0.159491, 0.899153);
    const Eigen::Vector3d b0e(-1.02941e-07,     0.470588,          0.6);
    const Eigen::Vector3d b1e(0.19509, -0.529412,   1.58079);
    // dt.x0 = Vector3d(0,0,1);
    // dt.x1 = Vector3d(0,1,1);
    // dt.x2 = Vector3d(0.1, 0.2, 2);
    // dt.x3 = Vector3d(0.1, 0.2, -1);

    // dt.x0b = Vector3d(1,1,0);
    // dt.x1b = Vector3d(0,0,0);
    // dt.x2b = Vector3d(0.1, 0.2, 2);
    // dt.x3b = Vector3d(0.1, 0.2, -1);
    double toi;
    bool res=intervalccd::edgeEdgeCCD_double(
        a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,{{4.9737991503207013e-11, 4.9737991503207013e-11, 4.9737991503207013e-11}},0,toi,1e-6
    );
    std::cout<<"the ee ccd result is "<<res<<", t, "<<toi<<std::endl;
}
void get_multiple(string data,string method, string part, double ms){
//double ms=1e-30;
    std::cout<<"ms "<<ms<<std::endl;
    std::cout<<"data "<<data<<std::endl;
     std::cout<<"method "<<method<<std::endl;
     std::cout<<"part "<<part<<std::endl;
    //string method="vertex-face";
    //string method="edge-edge";
    H5Easy::File file("/home/zachary/Development/ccd-queries/"+data+"/"+method+"/"+method+"-"+part+".hdf5");
const auto query_names = file.getGroup("/").listObjectNames();
std::cout<<"data size "<<query_names.size()<<std::endl;
igl::Timer timer;
    double time = 0;
for (size_t i = 0; i < query_names.size(); i++) {
    // if (i % 1000 == 0)
    //         std::cout << "i " << i << std::endl;
std::stringstream ss;
ss << query_names[i] << "/shifted/points";
Eigen::Matrix<double, 8, 3> vertex_face_data = H5Easy::load<Eigen::Matrix<double, 8, 3>>(
                    file, ss.str());
    vf_pair dt;
    for (int j = 0; j < 3; j++) {
			dt.x0[j] = vertex_face_data(0, j);
			dt.x1[j] = vertex_face_data(1, j);
			dt.x2[j] = vertex_face_data(2, j);
			dt.x3[j] = vertex_face_data(3, j);

			dt.x0b[j] = vertex_face_data(4, j);
			dt.x1b[j] = vertex_face_data(5, j);
			dt.x2b[j] = vertex_face_data(6, j);
			dt.x3b[j] = vertex_face_data(7, j);
	}
    if(method=="vertex-face"){
    timer.start();
    doubleccd::vertexFaceCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,ms);
    time += timer.getElapsedTimeInSec();

    }
    else{
        timer.start();
    doubleccd::edgeEdgeCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,ms);
    time += timer.getElapsedTimeInSec();

    }
    
}
std::cout << "total time " << time << std::endl;

    test();

    print_sub();
    ray_time();
    double phit=print_phi_time();
    std::cout<<"percentage "<<phit/time<<"\n"<<std::endl;
}
void call_multi(){
    std::vector<string> data,method,part;
    //data.push_back("chain");
    data.push_back("cow-heads");
    method.push_back("vertex-face");
    method.push_back("edge-edge");
    //part.push_back("collisions");
    for(int i=30;i<35;i++){
        part.push_back("collisions-part000"+to_string(i));
    }
    // part.push_back("collisions-part00000");
    // part.push_back("collisions-part00001");
    // part.push_back("collisions-part00002");
    // part.push_back("collisions-part00003");
    // part.push_back("collisions-part00004");
    std::vector<double> ms;
    ms.push_back(1e-8);
    ms.push_back(1e-30);
    for(int i=0;i<data.size();i++){
        for(int j=0;j<method.size();j++){
            for (int k=0;k<part.size();k++){
                for (int r=0;r<ms.size();r++){
                    get_multiple(data[i],method[j],part[k],ms[r]);
                    
                }
            }
        }
    }
   

}
// void test_ori(){
//     Vector3d p0,p1,p2,p3;
//     p0=Vector3d(0,0,0);
//     p1=Vector3d(1,1,0);
//     p2=Vector3d(2,0,1);
//     p3=Vector3d(3,4,0);
//     Vector2d q0(0,0),q1(1,0),q2(0,1);
//     std::cout<<"ori ind "<<ind_orient_3d(p0,p1,p2,p3)<<std::endl;
//     std::cout<<"ori geo "<<geo_orient_3d(p0,p1,p2,p3)<<std::endl;
//     std::cout<<"ori igl "<<igl_orient_3d(p0,p1,p2,p3)<<std::endl;
//     std::cout<<"ori2d ind "<<ind_orient_2d(q0,q1,q2)<<std::endl;
//     std::cout<<"ori2d geo "<<geo_orient_2d(q0,q1,q2)<<std::endl;
//     std::cout<<"ori2d igl "<<igl_orient_2d(q0,q1,q2)<<std::endl;


// }
void test_small(){
    using namespace intervalccd;
    long a;
    std::cout<<"reduction power "<<intervalccd::reduction(6,a)<<std::endl;
    std::cout<<"result, "<<a<<std::endl;

    intervalccd::Numccd low(1,3), up(1,2), mid;
    intervalccd::Singleinterval oi(low,up);
    std::pair<Singleinterval,Singleinterval> bi=bisect(oi);
    std::cout<<"new midpoint,"<<bi.first.second.first<<" "<<bi.first.second.second<<std::endl;
} 
void get_rand_paras_vf(double& a, double&b, double&c, double&d,
double&e, double&f, double&g, double&h ){
    const auto paras=Eigen::MatrixXd::Random(3,1);
    double t=fabs(paras(0));
    double u=fabs(paras(1));
    double v=fabs(paras(2))*(1-u);// make sure 1-u-v>=0
    //std::cout<<"tuv, "<<t<<", "<<u<<", "<<v<<std::endl;
     a=1-t, b=t, c=-(1-u-v)*(1-t),
    d=-t*(1-u-v), e=-u*(1-t), f=-u*(t),
    g=-v*(1-t),h=-v*t;
}
void get_rand_paras_ee(double& a, double&b, double&c, double&d,
double&e, double&f, double&g, double&h ){
    const auto paras=Eigen::MatrixXd::Random(3,1);
    double t=fabs(paras(0));
    double u=fabs(paras(1));
    double v=fabs(paras(2));
    //std::cout<<"tuv, "<<t<<", "<<u<<", "<<v<<std::endl;
     a=(1-t)*(1-u), b=t*(1-u), c=u*(1-t),
    d=u*t, e=-(1-v)*(1-t), f=-(1-v)*(t),
    g=-v*(1-t),h=-v*t;
}
void get_wrong_paras_ee(double& a, double&b, double&c, double&d,
double&e, double&f, double&g, double&h){
    const auto paras=Eigen::MatrixXd::Random(3,1);
    double t=fabs(paras(0));
    double u=fabs(paras(1));
    double v;
    int rand_switch=rand()%2;
    if(rand_switch==1){
        v =1+fabs(paras(2));
    }
    else{
        v=-1-fabs(paras(2));
    }
     //std::cout<<"tuv, "<<t<<", "<<u<<", "<<v<<std::endl;
     a=(1-t)*(1-u), b=t*(1-u), c=u*(1-t),
    d=u*t, e=-(1-v)*(1-t), f=-(1-v)*(t),
    g=-v*(1-t),h=-v*t;
    
}
void get_wrong_paras_vf(double& a, double&b, double&c, double&d,
double&e, double&f, double&g, double&h){

    int rand_switch=rand()%2;
    const auto paras=Eigen::MatrixXd::Random(3,1);
    double t=fabs(paras(0));
    double u=fabs(paras(1));
    double v;
    if(rand_switch==1){
        v =(1-u)+(1+fabs(paras(2)));// make sure u+v>1
    }
    else{
        v=-fabs(paras(2));
    }
   
    //std::cout<<"tuv, "<<t<<", "<<u<<", "<<v<<std::endl;
     a=1-t, b=t, c=-(1-u-v)*(1-t),
    d=-t*(1-u-v), e=-u*(1-t), f=-u*(t),
    g=-v*(1-t),h=-v*t;
}
Vector3d get_rand_p(){
    const auto paras=Eigen::MatrixXd::Random(3,1);
    return Vector3d(paras(0),paras(1),paras(2));
}
void print_V_to_string(const Eigen::MatrixXd &V, std::array<std::array<std::string,6>,8>&s){// V is 8x3
    
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){
            intervalccd::Rational r(V(i,j));
            std::cout<<r.get_numerator_str()<<", "<<r.get_denominator_str();
            s[i][j*2]=r.get_numerator_str();
            s[i][j*2+1]=r.get_denominator_str();
            
            if(j<2){
                std::cout<<",";
            }

        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}
void get_queries_from_funcs(bool vf,Eigen::MatrixXd& Vout,std::array<std::array<std::string,6>,8>&s, 
    bool two_roots){
    Vector3d x0s=get_rand_p(),x0e=get_rand_p(),
    x1s=get_rand_p(),x1e=get_rand_p(),x2s=get_rand_p();
    Vector3d x2e,x3s,x3e; //these are the points to be solved
    Vector3d a,b,c,d,e,f,g,h;
    if(vf){
     get_rand_paras_vf(a[0],b[0],c[0],d[0],e[0],f[0],g[0],h[0]);
     get_rand_paras_vf(a[1],b[1],c[1],d[1],e[1],f[1],g[1],h[1]);
     if(two_roots){
         get_wrong_paras_vf(a[2],b[2],c[2],d[2],e[2],f[2],g[2],h[2]);
     }
     else{
         get_rand_paras_vf(a[2],b[2],c[2],d[2],e[2],f[2],g[2],h[2]);
     }
     
    }
    else{
        get_rand_paras_ee(a[0],b[0],c[0],d[0],e[0],f[0],g[0],h[0]);
     get_rand_paras_ee(a[1],b[1],c[1],d[1],e[1],f[1],g[1],h[1]);
     
     if(two_roots){
         get_wrong_paras_ee(a[2],b[2],c[2],d[2],e[2],f[2],g[2],h[2]);
     }
     else{
         get_rand_paras_ee(a[2],b[2],c[2],d[2],e[2],f[2],g[2],h[2]);
     }
    }

    Vector3d B0,B1,B2;
    
    B0=a[0]*x0s+b[0]*x0e+c[0]*x1s+d[0]*x1e+e[0]*x2s;
    B1=a[1]*x0s+b[1]*x0e+c[1]*x1s+d[1]*x1e+e[1]*x2s;
    B2=a[2]*x0s+b[2]*x0e+c[2]*x1s+d[2]*x1e+e[2]*x2s;
    // A0X=-B0, A1X=-B1, A2X=-B2
    Vector3d nBx=-Vector3d(B0[0],B1[0],B2[0]),
    nBy=-Vector3d(B0[1],B1[1],B2[1]),nBz=-Vector3d(B0[2],B1[2],B2[2]);
    Eigen::Matrix3d A;
    A<<f[0],g[0],h[0],
    f[1],g[1],h[1],
    f[2],g[2],h[2];
    //AX_x=nBx. X_x=(x2ex, x3sx, x3ex)
    Vector3d x=A.lu().solve(nBx);
    Vector3d y=A.lu().solve(nBy);
    Vector3d z=A.lu().solve(nBz);
    x2e=Vector3d(x[0],y[0],z[0]);
    x3s=Vector3d(x[1],y[1],z[1]);
    x3e=Vector3d(x[2],y[2],z[2]);
    // std::cout<<"V\n"<<
    // x0s.transpose()<<"\n"<<
    // x1s.transpose()<<"\n"<<
    // x2s.transpose()<<"\n"<<
    // x3s.transpose()<<"\n"<<
    // x0e.transpose()<<"\n"<<
    // x1e.transpose()<<"\n"<<
    // x2e.transpose()<<"\n"<<
    // x3e.transpose()<<"\n"<<
    // std::endl<<std::endl;

    Eigen::MatrixXd V(8,3);
    for(int i=0;i<3;i++){
        V(0,i)=x0s[i];
        V(1,i)=x1s[i];
        V(2,i)=x2s[i];
        V(3,i)=x3s[i];
        V(4,i)=x0e[i];
        V(5,i)=x1e[i];
        V(6,i)=x2e[i];
        V(7,i)=x3e[i];
    }
    Vout=V;
    double max_coor=1;
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){
            if(max_coor<fabs(V(i,j)))
            max_coor=fabs(V(i,j));
        }
    }
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){
            V(i,j)/=max_coor;
        }
    }
    
    // std::cout<<"V\n"<<V<<std::endl;
    
    print_V_to_string(V,s);

}
void get_coplanar_queries(bool vf,Eigen::MatrixXd& Vout,std::array<std::array<std::string,6>,8>&s){
    // vf
    std::array<Vector3d,3> x;//trinagle
    x[0]=get_rand_p(); x[1]=get_rand_p(); x[2]=get_rand_p(); 
    Vector3d x3s=get_rand_p();// vertex
    int axis=rand()%3;
    x[0][axis]=0;x[1][axis]=0;x[2][axis]=0;x3s[axis]=0;
    int edge=rand()%3;
    Vector3d mid_p;
    if(vf){
        mid_p=(x[edge]+x[(edge+1)%3])/2;
    }
    else{
        mid_p=(x[1]+x[2])/2;
    }
    
    Vector3d rd=get_rand_p();
    double t=(1+fabs(rd[0]));
    Vector3d x3e=(mid_p-x3s)*t+x3s;
    Eigen::MatrixXd v(8,3);
    for(int i=0;i<3;i++){
        
        v(0,i)=x3s[i];
        v(1,i)=x[0][i];
        v(2,i)=x[1][i];
        v(3,i)=x[2][i];
        v(4,i)=x3e[i];
        v(5,i)=x[0][i];
        v(6,i)=x[1][i];
        v(7,i)=x[2][i];
    }

    Vout=v;
    std::cout<<v<<std::endl<<std::endl;
    print_V_to_string(v,s);
}

void get_coplanar_queries_no_collision(bool vf,Eigen::MatrixXd& Vout,std::array<std::array<std::string,6>,8>&s){
    // vf
    if(vf){
        double far=1;
    std::array<Vector3d,3> x;//trinagle
    x[0]=get_rand_p(); x[1]=get_rand_p(); x[2]=get_rand_p(); 
    Vector3d rd=get_rand_p();// vertex
    int axis=rand()%3;
    x[0][axis]=0;x[1][axis]=0;x[2][axis]=0;

    Vector3d bmin,bmax;
    for(int i=0;i<3;i++){
        bmin[i]=std::min(std::min(x[0][i],x[1][i]),x[2][i]);
        bmax[i]=std::max(std::max(x[0][i],x[1][i]),x[2][i]);
    }
    Vector3d x3s=bmax+far*Vector3d(fabs(rd[0]),fabs(rd[1]),fabs(rd[2]));
    x3s[axis]=0;
    rd=get_rand_p();
    Vector3d x3e=x3s+far*Vector3d(fabs(rd[0]),fabs(rd[1]),fabs(rd[2]));
    x3e[axis]=0;
    Eigen::MatrixXd v(8,3);
    for(int i=0;i<3;i++){
        
        v(0,i)=x3s[i];
        v(1,i)=x[0][i];
        v(2,i)=x[1][i];
        v(3,i)=x[2][i];
        v(4,i)=x3e[i];
        v(5,i)=x[0][i];
        v(6,i)=x[1][i];
        v(7,i)=x[2][i];
    }

    Vout=v;
    std::cout<<v<<std::endl<<std::endl;
    print_V_to_string(v,s);
    }
    else{  //  ee
        double far=1;
    Vector3d x0, x1,x2,x3s,x3e;//triangle
    x0=get_rand_p(); x1=get_rand_p(); 
    Vector3d rd=get_rand_p();// vertex
    int axis=rand()%3;
    x0[axis]=0;x1[axis]=0;
    Vector3d bmin,bmax;
    for(int i=0;i<3;i++){
        bmin[i]=std::min(x0[i],x1[i]);
        bmax[i]=std::max(x0[i],x1[i]);
    } 
    
    x2=bmax+far*Vector3d(fabs(rd[0]),fabs(rd[1]),fabs(rd[2]));
    x2[axis]=0;
    
    rd=get_rand_p();
    x3s=x2+far*Vector3d(fabs(rd[0]),fabs(rd[1]),fabs(rd[2]));
    
    rd=get_rand_p();
    x3e=x3s+far*Vector3d(fabs(rd[0]),fabs(rd[1]),fabs(rd[2]));
    
    Eigen::MatrixXd v(8,3);
    for(int i=0;i<3;i++){
        
        v(0,i)=x3s[i];
        v(1,i)=x2[i];
        v(2,i)=x1[i];
        v(3,i)=x0[i];
        v(4,i)=x3e[i];
        v(5,i)=x2[i];
        v(6,i)=x1[i];
        v(7,i)=x0[i];
    }

    Vout=v;
    std::cout<<v<<std::endl<<std::endl;
    print_V_to_string(v,s);
    }
}

void print_matrix_string(
    const std::vector<std::array<std::array<std::string,6>,8>>&s,std::string& filename){
    
    std::ofstream fout;

    
    fout.open(filename);
    for(int k=0;k<s.size();k++){
        for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){
            fout<<s[k][i][2*j]<<","<<s[k][i][2*j+1];
            if(j<2){
                fout<<",";
            }
            else{
                fout<<"\n";
            }
        }
        
        }
    }
    fout.close();
}
void generate_butterfly(bool vf,Eigen::MatrixXd& Vout,std::array<std::array<std::string,6>,8>&s){
    Eigen::MatrixXd v;
    if(vf){
        Vector3d x0s(0.1,0.1,0.1),
        x1s(0,0,1),x2s(1,0,1),x3s(0,1,1);
        
        Vector3d rd=get_rand_p();
        double pratio=fabs(rd[0]);
        x0s=x0s*pratio;

        double tratio=fabs(rd[1])+pratio;
        Vector3d x0e=x0s,x1e(0,0,0),x2e(0,1,0),x3e(1,0,0);
        x1s=x1s*tratio;
        x2s=x2s*tratio;
        x3s=x3s*tratio;

        x1e=x1e*tratio;
        x2e=x2e*tratio;
        x3e=x3e*tratio;


        Eigen::MatrixXd v(8,3);
        for(int i=0;i<3;i++){

            v(0,i)=x0s[i];
            v(1,i)=x1s[i];
            v(2,i)=x2s[i];
            v(3,i)=x3s[i];
            v(4,i)=x0e[i];
            v(5,i)=x1e[i];
            v(6,i)=x2e[i];
            v(7,i)=x3e[i];
        }
        Vout=v;
        std::cout<<v<<std::endl<<std::endl;
        print_V_to_string(v,s);

    }
    else{
        // Vector3d x0s(0.1,0.1,0.1),
        // x1s(0,0,1),x2s(1,0,1),x3s(0,1,1);
        
        // Vector3d rd=get_rand_p();
        // double pratio=fabs(rd[0]);
        // x0s=x0s*pratio;

        // double tratio=fabs(rd[1])+pratio;
        // Vector3d x0e=x0s,x1e(0,0,0),x2e(0,1,0),x3e(1,0,0);
        // x1s=x1s*tratio;
        // x2s=x2s*tratio;
        // x3s=x3s*tratio;

        // x1e=x1e*tratio;
        // x2e=x2e*tratio;
        // x3e=x3e*tratio;


        // Eigen::MatrixXd v(8,3);
        // for(int i=0;i<3;i++){

        //     v(0,i)=x0s[i];
        //     v(1,i)=x1s[i];
        //     v(2,i)=x2s[i];
        //     v(3,i)=x3s[i];
        //     v(4,i)=x0e[i];
        //     v(5,i)=x1e[i];
        //     v(6,i)=x2e[i];
        //     v(7,i)=x3e[i];
        // }
        // Vout=v;
        // std::cout<<v<<std::endl<<std::endl;
        // print_V_to_string(v,s);
    }
}

void generate_different_queries(){
    Eigen::MatrixXd Vout;
    std::array<std::array<std::string,6>,8> s;
    std::vector<std::array<std::array<std::string,6>,8>> svec;
   bool vf_flag=false;

   std::string filename;
   if(vf_flag){
        filename="/home/bolun1/interval/special_data_vf.csv";
    }
    else{
        filename="/home/bolun1/interval/special_data_ee.csv";
    }
    int length=5;
   // three roots
   for(int i=0;i<length;i++){
        get_queries_from_funcs(vf_flag,Vout,s,false);
        svec.emplace_back(s);
    }
    // two roots
    for(int i=0;i<length;i++){
        get_queries_from_funcs(vf_flag,Vout,s,true);
        svec.emplace_back(s);
    }
    //coplanar with roots
    for(int i=0;i<length;i++){
        get_coplanar_queries(vf_flag,Vout,s);
        svec.emplace_back(s);
    }
    //coplanar without roots
    for(int i=0;i<length;i++){
        get_coplanar_queries_no_collision(vf_flag,Vout,s);
        svec.emplace_back(s);
    }
    //butterfly
    
    if(vf_flag){
        std::cout<<"start butterfly"<<std::endl;
        for(int i=0;i<length;i++){
        generate_butterfly(vf_flag,Vout,s);
        svec.emplace_back(s);
    }
    }
    
    print_matrix_string(svec,filename);

}

int main(int argc, char* argv[])
{
    // generate_different_queries();
    // test_small();
//     Eigen::MatrixXd Vout;
//     std::array<std::array<std::string,6>,8> s;
//     // get_queries_from_funcs(false,Vout,s);
//     // get_queries_from_funcs(false,Vout,s);
//     // get_queries_from_funcs(false,Vout,s);
//     // get_queries_from_funcs(false,Vout,s);
//     //print_matrix_string(s);

//     std::vector<std::array<std::array<std::string,6>,8>> svec;
//    bool vf_flag=false;
    // for(int i=0;i<125;i++){
    //     get_queries_from_funcs(vf_flag,Vout,s);
    //     svec.emplace_back(s);
    // }
    // print_matrix_string(svec,vf_flag);

// for(int i=0;i<3;i++){
//         get_coplanar_queries(vf_flag,Vout,s);
//         svec.emplace_back(s);
//     }

    // for(int i=0;i<3;i++){
    //     get_coplanar_queries_no_collision(vf_flag,Vout,s);
    //     svec.emplace_back(s);
    // }
    //print_matrix_string(svec,vf_flag);

    
    // get_queries_from_funcs();
    // TODO: Put something more relevant here
    // ccd::test();
   // test_shifted_compare();
    // test_rootfinder();
    //test_shift_maxerror();
	//test_edge_edge();
    // const string filename="/home/bw1760/scratch/cow-head/edge-edge/edge-edge-collisions-part004.hdf5";
    // std::vector<Eigen::Matrix<double, 8, 3>> edge_edge_data;
	// read_edge_edge_data(filename, edge_edge_data);
    case_check();
    //check_false();
//     for(int i=0;i<20;i++){
//  doubleccd::Vector3d p0=Vector3d::Random(),p1=Vector3d::Random(),p2=Vector3d::Random(),p3=Vector3d::Random();
//     ccd::Vector3r p0r(p0[0],p0[1],p0[2]),p1r(p1[0],p1[1],p1[2]),p2r(p2[0],p2[1],p2[2]),p3r(p3[0],p3[1],p3[2]);
//     std::cout<<"orient "<<doubleccd::orient_3d(p0,p1,p2,p3)<<" "<<ccd::orient3d(p0r,p1r,p2r,p3r)<<std::endl;
//     }
   //compare_lpi_results();
   //call_multi();
//    test_ori();
//std::cout<<"pwr,"<<2^3<<std::endl;

    std::cout<<"done"<<std::endl;

    return 1;
}
