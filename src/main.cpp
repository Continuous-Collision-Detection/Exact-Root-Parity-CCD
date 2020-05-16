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
        if (i % 200 == 0)
            std::cout << "i " << i << std::endl;
        // std::cout << "i " << std::endl;
        timer.start();
        results[i] = vertexFaceCCD( // double
            data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s, data[i].pte,
            data[i].v1e, data[i].v2e, data[i].v3e, 1e-8);
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
void test_shift_maxerror()
{
    std::vector<sccd> data;
    Eigen::MatrixX3d vertices, vertices1, vertices2;
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
    double k;
    get_whole_mesh_shifted(vfs, ees, shifted_vfs, shifted_ees, vertices1);
    get_whole_mesh_shifted(vfs, ees, vertices2);
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
    //H5Easy::File file(root_path + path_sep +"vertex-face-collisions.hdf5");
    H5Easy::File file("/home/zachary/Development/ccd-queries/erleben-cube-internal-edges/edge-edge/edge-edge-collisions.hdf5");
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
    dt.x0 = Vector3d(5.06670675754481e-17,0.571825903784439,1);
    dt.x1 = Vector3d(0,0.5,1);
    dt.x2 = Vector3d(0,0.5,0);
    dt.x3 = Vector3d(-0.25,-0.5,-0.25);

    dt.x0b = Vector3d(1.26333711630941e-16,0.277626403784438,1);
    dt.x1b = Vector3d(0,0.5,1);
    dt.x2b = Vector3d(0,0.5,0);
    dt.x3b = Vector3d(-0.25,-0.5,-0.25);
    std::cout<<"the vf ccd result is "<<doubleccd::vertexFaceCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,1e-300)<<std::endl;
}
void get_multiple(){
double ms=1e-8;
    std::cout<<"ms "<<ms<<std::endl;

    H5Easy::File file("/home/zachary/Development/ccd-queries/erleben-cube-cliff-edges/edge-edge/edge-edge-collisions.hdf5");
const auto query_names = file.getGroup("/").listObjectNames();
std::cout<<"data size "<<query_names.size()<<std::endl;
igl::Timer timer;
    double time = 0;
for (size_t i = 0; i < query_names.size(); i++) {
    if (i % 1000 == 0)
            std::cout << "i " << i << std::endl;
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
    timer.start();
    doubleccd::edgeEdgeCCD(dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b,ms);
    time += timer.getElapsedTimeInSec();

}
std::cout << "total time " << time << std::endl;

    test();

    print_sub();
    ray_time();
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

int main(int argc, char* argv[])
{
    // TODO: Put something more relevant here
    // ccd::test();
   test_shifted_compare();
    // test_rootfinder();
    //test_shift_maxerror();
	//test_edge_edge();
    // const string filename="/home/bw1760/scratch/cow-head/edge-edge/edge-edge-collisions-part004.hdf5";
    // std::vector<Eigen::Matrix<double, 8, 3>> edge_edge_data;
	// read_edge_edge_data(filename, edge_edge_data);
   // case_check();
    //check_false();
//     for(int i=0;i<20;i++){
//  doubleccd::Vector3d p0=Vector3d::Random(),p1=Vector3d::Random(),p2=Vector3d::Random(),p3=Vector3d::Random();
//     ccd::Vector3r p0r(p0[0],p0[1],p0[2]),p1r(p1[0],p1[1],p1[2]),p2r(p2[0],p2[1],p2[2]),p3r(p3[0],p3[1],p3[2]);
//     std::cout<<"orient "<<doubleccd::orient_3d(p0,p1,p2,p3)<<" "<<ccd::orient3d(p0r,p1r,p2r,p3r)<<std::endl;
//     }
   //compare_lpi_results();
   //get_multiple();
//    test_ori();
    std::cout<<"done"<<std::endl;

    return 1;
}
