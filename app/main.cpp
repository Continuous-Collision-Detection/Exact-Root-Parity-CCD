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
#include <iomanip>

using namespace doubleccd;
void case_check(){
    doubleccd::vf_pair dt,dtshift;
    dt.x0 = Vector3d(0,0,1);
    dt.x1 = Vector3d(0,1,1);
    dt.x2 = Vector3d(0.1, 0.2, 2);
    dt.x3 = Vector3d(0.1, 0.2, -1);

    dt.x0b = Vector3d(1,1,0);
    dt.x1b = Vector3d(0,0,0);
    dt.x2b = Vector3d(0.1, 0.2, 2);
    dt.x3b = Vector3d(0.1, 0.2, -1);

	double time;
	double err=shift_vertex_face(dt,dtshift,time);
    dt=dtshift;
    bool res=doubleccd::vertexFaceCCD(
        dt.x0,dt.x1,dt.x2,dt.x3,dt.x0b,dt.x1b,dt.x2b,dt.x3b
    );
    std::cout<<"the double ccd result is "<<res<<std::endl;
    std::cout<<"the rational ccd result is "<<ccd::vertexFaceCCD(
        ccd::Vector3d(dt.x0),ccd::Vector3d(dt.x1),ccd::Vector3d(dt.x2),ccd::Vector3d(dt.x3),
        ccd::Vector3d(dt.x0b),ccd::Vector3d(dt.x1b),ccd::Vector3d(dt.x2b),ccd::Vector3d(dt.x3b))<<std::endl;
    std::cout<<"shift err, "<<err<<std::endl;
}

void test_shift(){
    std::vector<doubleccd::vf_pair> data;
    std::vector<std::pair<double, double>> subs;
    subs.resize(1);
    // subs[0].first=0.49827393449335344; subs[0].second=0.22518885388748636;
    for(int i=0;i<100000000;i++){
        subs.emplace_back(double(rand())/double((RAND_MAX)),double(rand())/double((RAND_MAX)));
    }

    displaceSubtractions_double(subs);
    for(int i=0;i<subs.size();i++){
        if(i%1000000==0){
            std::cout<<subs[i].first<<", "<<subs[i].second<<std::endl;
        }
        Rational a=subs[i].first;
        Rational b=subs[i].second;
        Rational rst=a-b;
        double rd=subs[i].first-subs[i].second;
        if(rst==rd){
            // std::cout<<"no truncation err"<<std::endl;
        }
        else{
            std::cout<<"diff, "<<(rst-rd)<<", "<<rst<<", "<<rd<<std::endl;
            exit(0);

        }
    }

    // 0.49827393449335344,0.22518885388748636
}
std::vector<std::pair<double,double>> read_rational_CSV(const std::string inputFileName) {



    std::vector<std::pair<double,double>> vs;

    vs.clear();

	std::ifstream infile;

	infile.open(inputFileName);

    // std::array<double,3> v;

	if (!infile.is_open())

	{

		std::cout << "Path Wrong!!!!" << std::endl;

        return vs;

	}
	int l = 0;

	while (infile) // there is input overload classfile

	{
		l++;
		std::string s;
		if (!getline(infile, s)) break;
		if (s[0] != '#') {
			std::istringstream ss(s);
			std::array<std::string,4> record;
			int c = 0;
			while (ss) {
				std::string line;
				if (!getline(ss, line, ','))
					break;
				try {
					record[c] = line;
					c++;
				}
				catch (const std::invalid_argument e) {
					std::cout << "NaN found in file " << inputFileName << " line " << l
						<< std::endl;
					e.what();
				}
			}
            Rational rt;
            double x=rt.get_double(record[0],record[1]),y=rt.get_double(record[2],record[3]);
            std::pair<double,double> pr;pr.first=x;pr.second=y;
            vs.push_back(pr);
		}
	}
    return vs;
    // Eigen::MatrixXd all_v(vs.size(),3);
    // for(int i=0;i<vs.size();i++){
    //     all_v(i,0)=vs[i][0];
    //   all_v(i,1)=vs[i][1];
    //     all_v(i,2)=vs[i][2];
    // }
	// if (!infile.eof()) {
	// 	std::cerr << "Could not read file " << inputFileName << "\n";
	// }
	// return all_v;
}
bool check_subs_err(std::vector<std::pair<double, double>> subs, std::string discribe){
    bool have_err=false;
    for(int i=0;i<subs.size();i++){
        Rational a=subs[i].first;
        Rational b=subs[i].second;
        Rational rst=a-b;
        double rd=subs[i].first-subs[i].second;
        if(rst==rd){

        }
        else{
            std::cout<<discribe<<" diff, "<<std::setprecision(17)<<i <<", "<<(rst-rd)<<", "<<rst<<", "<<rd<<std::endl;
            have_err=true;
        }
    }
    return have_err;
}
void read_and_test_shift(){
    std::vector<std::pair<double,double>> readed = read_rational_CSV(
    "/home/bolun1/interval/Round/CCD-Wrapper/build/single_ee.csv");
    displaceSubtractions_double(readed);
    check_subs_err(readed,"we successfully reproduce this problem");
    
}

int main(int argc, char* argv[])
{
    
   case_check();
    // test_shift();
    // read_and_test_shift();
	doubleccd::test();

    std::cout<<"done"<<std::endl;

    return 1;
}
