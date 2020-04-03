#include <iostream>
//
//#include <ccd.hpp>
#include <vector>
//#include<Utils.hpp>
//#include<subfunctions.h>
#include <array>
#include <doubleCCD/doubleccd.hpp>
#include <CCD/exact_subtraction.hpp>
#include <CCD/ccd.hpp>
#include <fstream>
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
Vector3d construct1(double a, double b, double c)
{
    return Vector3d(a, b, c);
}
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


void test1() {
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
void test_compare() {
	std::vector<sccd> data;
	read_CSV(root_path + path_sep+ "cow-head-collisions.csv", data);
	vector<bool> rst;
	read_result(root_path + path_sep+ "result_all.csv", rst);
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
		if (i % 200 == 0)std::cout << "i " << i << std::endl;
		//std::cout << "i " << std::endl;
		/*results[i] = vertexFaceCCD(
			data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
			data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-8);*/
		results[i] = ccd:: vertexFaceCCD(
			data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
			data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-8); 
		if (rst[i] == 1 && results[i] == 0) {//when old method says yes but we missed it
			std::cout << "wrong! i= " << i << std::endl;
		}
		//std::cout << "result " << results[i] << std::endl;
		if (results[i] == true) inside++;
	}
	std::cout << int_seg_XOR(1, 0) << std::endl;
	cube cb(0.1);
	std::cout << "inside number " << inside << std::endl;

}

void test_shifted_compare() {
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
	std::vector<bool> results,results1;
	int fn = data.size(); // 50000;
	results.resize(fn);
	results1.resize(fn);
	int inside = 0;
	for (int i = 0; i < fn; i++) {
		if (i % 200 == 0)std::cout << "i " << i << std::endl;
		//std::cout << "i " << std::endl;
		results[i] = vertexFaceCCD(//double 
			data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
			data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-3);
		results1[i] = ccd::vertexFaceCCD(//rational
			data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s,
			data[i].pte, data[i].v1e, data[i].v2e, data[i].v3e, 1e-3);
		//if (rst[i] == 1 && results[i] == 0) {//when old method says yes but we missed it
		//	std::cout << "wrong! i= " << i << std::endl;
		//}
		//std::cout << "Rational vs double " << results1[i] << " " << results[i] << std::endl;
		if (results1[i] != results[i]) {//when old method says yes but we missed it
			std::cout << "double don't match rational! i= " << i << std::endl;
			std::cout << "Rational vs double " << results1[i]<<" "<<results[i] << std::endl;
		}
		//std::cout << "result " << results[i] << std::endl;
		if (results[i] == true) inside++;
	}
	cube cb(0.1);
	std::cout << "inside number " << inside << std::endl;

}
int main(int argc, char* argv[])
{
    // TODO: Put something more relevant here
    //ccd::test();
	test_shifted_compare();
    return 1;
}
