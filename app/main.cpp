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
    
// const Eigen::Vector3d a0s(-30022200, 2362580, 165247);
//     const Eigen::Vector3d a1s(-32347850, 8312380, -1151003);
//     const Eigen::Vector3d a0e(-28995600, 345838, 638580);
//     const Eigen::Vector3d a1e(-31716930, 6104858, -713340);
//     const Eigen::Vector3d b0s(-30319900, 3148750, 0);
//     const Eigen::Vector3d b1s(-28548800, 900349, 0);
//     const Eigen::Vector3d b0e(-30319900, 3148750, 0);
//     const Eigen::Vector3d b1e(-28548800, 900349, 0);

    // const Eigen::Vector3d a0s(5.55495e-11,    0.160014,    0.914204);
    // const Eigen::Vector3d a1s(0.098607, 0.160014, 0.898586);
    // const Eigen::Vector3d a0e(-1.02941e-07     ,0.470588          ,0.6);
    // const Eigen::Vector3d a1e(0.19509, -0.529412,   1.58079 );
    // const Eigen::Vector3d b0s(1.42931e-05 ,   0.159178,    0.914963);
    // const Eigen::Vector3d b1s(0.098703, 0.159491, 0.899153);
    // const Eigen::Vector3d b0e(-1.02941e-07,     0.470588,          0.6);
    // const Eigen::Vector3d b1e(0.19509, -0.529412,   1.58079);
    dt.x0 = Vector3d(0,0,1);
    dt.x1 = Vector3d(0,1,1);
    dt.x2 = Vector3d(0.1, 0.2, 2);
    dt.x3 = Vector3d(0.1, 0.2, -1);

    dt.x0b = Vector3d(1,1,0);
    dt.x1b = Vector3d(0,0,0);
    dt.x2b = Vector3d(0.1, 0.2, 2);
    dt.x3b = Vector3d(0.1, 0.2, -1);
    // std::vector<vf_pair> data1;std::vector<ee_pair> data2;
    // Eigen::MatrixX3d vertices;
    // data1.push_back(dt);
    
    // get_whole_mesh_shifted(
    // data1,data2,vertices
    // );
    // std::cout<<"data1 size, "<<data1.size()<<std::endl;
    // std::cout<<"data2 size, "<<data2.size()<<std::endl;
    // dt=data1[0];
    // std::cout<<dt.x0[0]<<","<<dt.x0[1]<<","<<dt.x0[2]<<std::endl;
    // std::cout<<dt.x1[0]<<","<<dt.x1[1]<<","<<dt.x1[2]<<std::endl;
    // std::cout<<dt.x2[0]<<","<<dt.x2[1]<<","<<dt.x2[2]<<std::endl;
    // std::cout<<dt.x3[0]<<","<<dt.x3[1]<<","<<dt.x3[2]<<std::endl;
    // std::cout<<dt.x0b[0]<<","<<dt.x0b[1]<<","<<dt.x0b[2]<<std::endl;
    // std::cout<<dt.x1b[0]<<","<<dt.x1b[1]<<","<<dt.x1b[2]<<std::endl;
    // std::cout<<dt.x2b[0]<<","<<dt.x2b[1]<<","<<dt.x2b[2]<<std::endl;
    // std::cout<<dt.x3b[0]<<","<<dt.x3b[1]<<","<<dt.x3b[2]<<std::endl;
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

    // for(int i=0;i<subs.size();i++){
    //     if(i%1000000==0){
    //         std::cout<<subs[i].first<<", "<<subs[i].second<<std::endl;
    //     }
    //     Rational a=subs[i].first;
    //     Rational b=subs[i].second;
    //     Rational rst=a-b;
    //     double rd=subs[i].first-subs[i].second;
    //     if(rst==rd){
    //         // std::cout<<"no truncation err"<<std::endl;
    //     }
    //     else{
    //         std::cout<<"before perturb, diff, "<<(rst-rd)<<", "<<rst<<", "<<rd<<std::endl;
    //         // exit(0);

    //     }
    // }
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
    
   // case_check();
    // test_shift();
    // read_and_test_shift();
	doubleccd::test();

    std::cout<<"done"<<std::endl;

    return 1;
}
