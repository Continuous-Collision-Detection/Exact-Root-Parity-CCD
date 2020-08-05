#pragma once
#include<interval_ccd/interval.hpp>
namespace intervalccd{
    
    std::array<double,8> function_vf(
        const double& vs,
        const double& t0s,
        const double& t1s,const double& t2s,
    const double& ve,const double& t0e,
    const double& t1e,const double& t2e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw);

    std::array<double,4> function_vf(
    const double& vs,const double& t0s,
    const double& t1s,const double& t2s,
    const double& ve,const double& t0e,
    const double& t1e,const double& t2e,
    const std::array<double,4>& t_up,const std::array<double,4>& t_dw,
    const std::array<double,4>& u_up,const std::array<double,4>& u_dw,
    const std::array<double,4>& v_up,const std::array<double,4>& v_dw);

    std::array<double,8> function_ee(
    const double& a0s,const double& a1s,const double& b0s,const double& b1s,
    const double& a0e,const double& a1e,const double& b0e,const double& b1e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw);
    std::array<double,4> function_ee(
    const double& a0s,const double& a1s,const double& b0s,const double& b1s,
    const double& a0e,const double& a1e,const double& b0e,const double& b1e,
    const std::array<double,4>& t_up,const std::array<double,4>& t_dw,
    const std::array<double,4>& u_up,const std::array<double,4>& u_dw,
    const std::array<double,4>& v_up,const std::array<double,4>& v_dw);
    std::array<int,4> select_useful_pt_id(const int split_i);

    class Interval3_2_double{
        public:
        Interval3_2_double(const Interval3& itv);
        double t0_up;
        double t0_dw;
        double t1_up;
        double t1_dw;
        double u0_up;
        double u0_dw;
        double u1_up;
        double u1_dw;
        double v0_up;
        double v0_dw;
        double v1_up;
        double v1_dw;
    };
    void convert_tuv_to_array(const Interval3_2_double& itv, 
    std::array<double,8>& t_up,std::array<double,8>&t_dw,
    std::array<double,8>& u_up,std::array<double,8>&u_dw,
    std::array<double,8>& v_up,std::array<double,8>&v_dw);
    void convert_tuv_to_array_half(const Interval3_2_double& itv, 
    std::array<double,4>& t_up,std::array<double,4>&t_dw,
    std::array<double,4>& u_up,std::array<double,4>&u_dw,
    std::array<double,4>& v_up,std::array<double,4>&v_dw,const int split_i);
    
    
    // calculate a*(2^b)
int power(const int a, const int b);
    //class get_
}