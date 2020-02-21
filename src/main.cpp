#include <iostream>

#include <ccd.hpp>
#include <vector>
#include<Utils.hpp>
#include <fstream>
#include <array>
#include <exact_subtraction.hpp>
#include<subfunctions.h>
using namespace ccd;
using namespace std;
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

            data.emplace_back(record);
        }
    }
    if (!infile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
    }
    cout << "data size " << data.size() << endl;
    return true;
}


void test1() {

}

int main(int argc, char* argv[])
{
    // TODO: Put something more relevant here
    //ccd::test();
    std::vector<sccd> data;
    read_CSV("D:\\vs\\collision\\CCD\\data\\cow-head-collisions.csv", data);
    std::vector<bool> results;
    int fn = data.size(); // 50000;
    results.resize(fn);
    double k;
    int count = 0;
    double maxvalue = 0, merror = 0,maxerror;
    
    bool correct;
    for (int i = 0; i < fn; i++) {
        if (i % 200 == 0)
            std::cout << "i " << i << std::endl;

		get_prism_vertices_double(
            data[i].pts, data[i].v1s, data[i].v2s, data[i].v3s, data[i].pte,
            data[i].v1e, data[i].v2e, data[i].v3e,k,correct,maxerror);
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
		}
        
    }
    cout << "percentage of using k " << double(count) / double(fn) << std::endl;
    cout << "max of k " << maxvalue << std::endl;
    cout << "max error " << merror << std::endl;
	
    std::cout << "Exact CCD" << std::endl;
    return 1;
}
