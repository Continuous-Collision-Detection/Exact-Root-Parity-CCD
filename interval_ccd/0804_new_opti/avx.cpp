//#include <immintrin.h>
#include <stdio.h>
#include <array>
#include <interval_ccd/avx.h>

namespace intervalccd{
//#include <interval_ccd/avx.h>
//#include <iostream>
// __m512d function_f_ee_vector(
//     __m512d a0s,__m512d a1s,__m512d b0s,__m512d b1s,
//     __m512d a0e,__m512d a1e,__m512d b0e,__m512d b1e,
//     __m512d t_up,__m512d t_dw,__m512d u_up,__m512d u_dw,__m512d v_up,__m512d v_dw){

//         __m512d edge0_vertex0=_mm512_sub_pd(a0e,a0s);
//         edge0_vertex0=_mm512_mul_pd(edge0_vertex0,t_up);
//         edge0_vertex0=_mm512_div_pd(edge0_vertex0,t_dw);
//         edge0_vertex0=_mm512_add_pd(edge0_vertex0,a0s);

//         __m512d edge0_vertex1=_mm512_sub_pd(a1e,a1s);
//         edge0_vertex1=_mm512_mul_pd(edge0_vertex1,t_up);
//         edge0_vertex1=_mm512_div_pd(edge0_vertex1,t_dw);
//         edge0_vertex1=_mm512_add_pd(edge0_vertex1,a1s);

//         __m512d edge1_vertex0=_mm512_sub_pd(b0e,b0s);
//         edge1_vertex0=_mm512_mul_pd(edge1_vertex0,t_up);
//         edge1_vertex0=_mm512_div_pd(edge1_vertex0,t_dw);
//         edge1_vertex0=_mm512_add_pd(edge1_vertex0,b0s);

//         __m512d edge1_vertex1=_mm512_sub_pd(b1e,b1s);
//         edge1_vertex1=_mm512_mul_pd(edge1_vertex1,t_up);
//         edge1_vertex1=_mm512_div_pd(edge1_vertex1,t_dw);
//         edge1_vertex1=_mm512_add_pd(edge1_vertex1,b1s);

//         __m512d edge0_vertex=_mm512_sub_pd(edge0_vertex1,edge0_vertex0);
//         edge0_vertex=_mm512_mul_pd(edge0_vertex,u_up);
//         edge0_vertex=_mm512_div_pd(edge0_vertex,u_dw);
//         edge0_vertex=_mm512_add_pd(edge0_vertex,edge0_vertex0);

//         __m512d edge1_vertex=_mm512_sub_pd(edge1_vertex1,edge1_vertex0);
//         edge1_vertex=_mm512_mul_pd(edge1_vertex,v_up);
//         edge1_vertex=_mm512_div_pd(edge1_vertex,v_dw);
//         edge1_vertex=_mm512_add_pd(edge1_vertex,edge1_vertex0);

//         return _mm512_sub_pd(edge1_vertex,edge0_vertex);
//     }


// calculate a*(2^b)
int power(const int a, const int b){
  return a<<b;
}
std::array<double,8> function_ee(
    const double& a0s,const double& a1s,const double& b0s,const double& b1s,
    const double& a0e,const double& a1e,const double& b0e,const double& b1e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw){
        std::array<double,8> rst;
        for(int i=0;i<8;i++){
            double edge0_vertex0=(a0e-a0s)*t_up[i]/t_dw[i]+a0s;
            double edge0_vertex1=(a1e-a1s)*t_up[i]/t_dw[i]+a1s;
            double edge1_vertex0=(b0e-b0s)*t_up[i]/t_dw[i]+b0s;
            double edge1_vertex1=(b1e-b1s)*t_up[i]/t_dw[i]+b1s;

            double edge0_vertex=(edge0_vertex1-edge0_vertex0)*u_up[i]/u_dw[i]+edge0_vertex0;
            double edge1_vertex=(edge1_vertex1-edge1_vertex0)*v_up[i]/v_dw[i]+edge1_vertex0;
            rst[i]=edge0_vertex-edge1_vertex;
            
        }
        return rst;
    }
std::array<double,4> function_ee(
    const double& a0s,const double& a1s,const double& b0s,const double& b1s,
    const double& a0e,const double& a1e,const double& b0e,const double& b1e,
    const std::array<double,4>& t_up,const std::array<double,4>& t_dw,
    const std::array<double,4>& u_up,const std::array<double,4>& u_dw,
    const std::array<double,4>& v_up,const std::array<double,4>& v_dw){
        std::array<double,4> rst;
        for(int i=0;i<4;i++){
            double edge0_vertex0=(a0e-a0s)*t_up[i]/t_dw[i]+a0s;
            double edge0_vertex1=(a1e-a1s)*t_up[i]/t_dw[i]+a1s;
            double edge1_vertex0=(b0e-b0s)*t_up[i]/t_dw[i]+b0s;
            double edge1_vertex1=(b1e-b1s)*t_up[i]/t_dw[i]+b1s;

            double edge0_vertex=(edge0_vertex1-edge0_vertex0)*u_up[i]/u_dw[i]+edge0_vertex0;
            double edge1_vertex=(edge1_vertex1-edge1_vertex0)*v_up[i]/v_dw[i]+edge1_vertex0;
            rst[i]=edge0_vertex-edge1_vertex;
            
        }
        return rst;
    }
    std::array<double,8> function_vf(
    const double& vs,const double& t0s,
    const double& t1s,const double& t2s,
    const double& ve,const double& t0e,
    const double& t1e,const double& t2e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw){
        std::array<double,8> rst;
        for(int i=0;i<8;i++){
            double v=(ve-vs)*t_up[i]/t_dw[i]+vs;
            double t0=(t0e-t0s)*t_up[i]/t_dw[i]+t0s;
            double t1=(t1e-t1s)*t_up[i]/t_dw[i]+t1s;
            double t2=(t2e-t2s)*t_up[i]/t_dw[i]+t2s;
            double pt=(t1-t0)*u_up[i]/u_dw[i]+(t2-t0)*v_up[i]/v_dw[i]+t0;
            rst[i]=v-pt;
        }
        return rst;
    }
    std::array<double,4> function_vf(
    const double& vs,const double& t0s,
    const double& t1s,const double& t2s,
    const double& ve,const double& t0e,
    const double& t1e,const double& t2e,
    const std::array<double,4>& t_up,const std::array<double,4>& t_dw,
    const std::array<double,4>& u_up,const std::array<double,4>& u_dw,
    const std::array<double,4>& v_up,const std::array<double,4>& v_dw){
        std::array<double,4> rst;
        for(int i=0;i<4;i++){
            double v=(ve-vs)*t_up[i]/t_dw[i]+vs;
            double t0=(t0e-t0s)*t_up[i]/t_dw[i]+t0s;
            double t1=(t1e-t1s)*t_up[i]/t_dw[i]+t1s;
            double t2=(t2e-t2s)*t_up[i]/t_dw[i]+t2s;
            double pt=(t1-t0)*u_up[i]/u_dw[i]+(t2-t0)*v_up[i]/v_dw[i]+t0;
            rst[i]=v-pt;
        }
        return rst;
    }

    // std::array<double,8> double_to_array(const double& a){
    //     std::array<double,8> rst;
    //     for(int i=0;i<8;i++){

    //     }
    // }

    Interval3_2_double::Interval3_2_double(const Interval3& itv){
        t0_up=itv[0].first.first,t0_dw=power(1,itv[0].first.second),
        t1_up=itv[0].second.first,t1_dw=power(1,itv[0].second.second),

        u0_up=itv[1].first.first,u0_dw=power(1,itv[1].first.second),
        u1_up=itv[1].second.first,u1_dw=power(1,itv[1].second.second),

        v0_up=itv[2].first.first,v0_dw=power(1,itv[2].first.second),
        v1_up=itv[2].second.first,v1_dw=power(1,itv[2].second.second);
    }
    void convert_tuv_to_array(const Interval3_2_double& itv, 
    std::array<double,8>& t_up,std::array<double,8>&t_dw,
    std::array<double,8>& u_up,std::array<double,8>&u_dw,
    std::array<double,8>& v_up,std::array<double,8>&v_dw){
        // t order: 0,0,0,0,1,1,1,1
        // u order: 0,0,1,1,0,0,1,1
        // v order: 0,1,0,1,0,1,0,1
        
        t_up={{itv.t0_up,itv.t0_up,itv.t0_up,itv.t0_up,itv.t1_up,itv.t1_up,itv.t1_up,itv.t1_up}};
        t_dw={{itv.t0_dw,itv.t0_dw,itv.t0_dw,itv.t0_dw,itv.t1_dw,itv.t1_dw,itv.t1_dw,itv.t1_dw}};
        u_up={{itv.u0_up,itv.u0_up,itv.u1_up,itv.u1_up,itv.u0_up,itv.u0_up,itv.u1_up,itv.u1_up}};
        u_dw={{itv.u0_dw,itv.u0_dw,itv.u1_dw,itv.u1_dw,itv.u0_dw,itv.u0_dw,itv.u1_dw,itv.u1_dw}};
        v_up={{itv.v0_up,itv.v1_up,itv.v0_up,itv.v1_up,itv.v0_up,itv.v1_up,itv.v0_up,itv.v1_up}};
        v_dw={{itv.v0_dw,itv.v1_dw,itv.v0_dw,itv.v1_dw,itv.v0_dw,itv.v1_dw,itv.v0_dw,itv.v1_dw}};
    }
    void convert_tuv_to_array_half(const Interval3_2_double& itv, 
    std::array<double,4>& t_up,std::array<double,4>&t_dw,
    std::array<double,4>& u_up,std::array<double,4>&u_dw,
    std::array<double,4>& v_up,std::array<double,4>&v_dw,const int split_i){
        // t order: 0,0,0,0,1,1,1,1
        // u order: 0,0,1,1,0,0,1,1
        // v order: 0,1,0,1,0,1,0,1
        
        if(split_i==0){
        t_up={{itv.t0_up,itv.t0_up,itv.t0_up,itv.t0_up}};
        t_dw={{itv.t0_dw,itv.t0_dw,itv.t0_dw,itv.t0_dw}};
        u_up={{itv.u0_up,itv.u0_up,itv.u1_up,itv.u1_up}};
        u_dw={{itv.u0_dw,itv.u0_dw,itv.u1_dw,itv.u1_dw}};
        v_up={{itv.v0_up,itv.v1_up,itv.v0_up,itv.v1_up}};
        v_dw={{itv.v0_dw,itv.v1_dw,itv.v0_dw,itv.v1_dw}};
        }
        if(split_i==1){
        t_up={{itv.t0_up,itv.t0_up,itv.t1_up,itv.t1_up}};
        t_dw={{itv.t0_dw,itv.t0_dw,itv.t1_dw,itv.t1_dw}};
        u_up={{itv.u0_up,itv.u0_up,itv.u0_up,itv.u0_up}};
        u_dw={{itv.u0_dw,itv.u0_dw,itv.u0_dw,itv.u0_dw}};
        v_up={{itv.v0_up,itv.v1_up,itv.v0_up,itv.v1_up}};
        v_dw={{itv.v0_dw,itv.v1_dw,itv.v0_dw,itv.v1_dw}};
        }
        if(split_i==2){
        t_up={{itv.t0_up,itv.t1_up,itv.t0_up,itv.t1_up}};
        t_dw={{itv.t0_dw,itv.t1_dw,itv.t0_dw,itv.t1_dw}};
        u_up={{itv.u0_up,itv.u0_up,itv.u1_up,itv.u1_up}};
        u_dw={{itv.u0_dw,itv.u0_dw,itv.u1_dw,itv.u1_dw}};
        v_up={{itv.v0_up,itv.v0_up,itv.v0_up,itv.v0_up}};
        v_dw={{itv.v0_dw,itv.v0_dw,itv.v0_dw,itv.v0_dw}};
        }
    }
    std::array<int,4> select_useful_pt_id(const int split_i){
        
        if(split_i==0){
            return {{0,1,2,3}};
        }
        if(split_i==1){
            return {{0,1,4,5}};
        }
        else{
            return {{0,2,4,6}};
        }
    }
    
//     __m512d function_f_vf_vector(
//     __m512d vs,__m512d t0s,__m512d t1s,__m512d t2s,
//     __m512d ve,__m512d t0e,__m512d t1e,__m512d t2e,
//     __m512d t_up,__m512d t_dw,__m512d u_up,__m512d u_dw,__m512d v_up,__m512d v_dw){
//          __m512d v=_mm512_sub_pd(ve,vs);
//         v=_mm512_mul_pd(v,t_up);
//         v=_mm512_div_pd(v,t_dw);
//         v=_mm512_add_pd(v,vs);

//          __m512d t0=_mm512_sub_pd(t0e,t0s);
//         t0=_mm512_mul_pd(t0,t_up);
//         t0=_mm512_div_pd(t0,t_dw);
//         t0=_mm512_add_pd(t0,t0s);

//          __m512d t1=_mm512_sub_pd(t1e,t1s);
//         t1=_mm512_mul_pd(t1,t_up);
//         t1=_mm512_div_pd(t1,t_dw);
//         t1=_mm512_add_pd(t1,t1s);

//          __m512d t2=_mm512_sub_pd(t2e,t2s);
//         t2=_mm512_mul_pd(t2,t_up);
//         t2=_mm512_div_pd(t2,t_dw);
//         t2=_mm512_add_pd(t2,t2s);

//          __m512d t01=_mm512_sub_pd(t1,t0);
//         t01=_mm512_mul_pd(t01,u_up);
//         t01=_mm512_div_pd(t01,u_dw);

//         __m512d t02=_mm512_sub_pd(t2,t0);
//         t02=_mm512_mul_pd(t02,v_up);
//         t02=_mm512_div_pd(t02,v_dw);
       
//        __m512d pt=_mm512_add_pd(t01,t02);
//         pt=_mm512_add_pd(pt,t0);

//        return _mm512_sub_pd(v,pt);

// } 

// void convert_to_vector_pts(
//     const double asd,const double bsd,const double csd,const double dsd,
//     const double aed,const double bed,const double ced,const double ded,


//     __m512d& as, __m512d& bs, __m512d& cs, __m512d& ds,
//     __m512d& ae, __m512d& be, __m512d& ce, __m512d& de

//     ){
//         as=_mm512_setr_pd(asd,asd,asd,asd,asd,asd,asd,asd);
//         bs=_mm512_setr_pd(bsd,bsd,bsd,bsd,bsd,bsd,bsd,bsd);
//         cs=_mm512_setr_pd(csd,csd,csd,csd,csd,csd,csd,csd);
//         ds=_mm512_setr_pd(dsd,dsd,dsd,dsd,dsd,dsd,dsd,dsd);

//         ae=_mm512_setr_pd(aed,aed,aed,aed,aed,aed,aed,aed);
//         be=_mm512_setr_pd(bed,bed,bed,bed,bed,bed,bed,bed);
//         ce=_mm512_setr_pd(ced,ced,ced,ced,ced,ced,ced,ced);
//         de=_mm512_setr_pd(ded,ded,ded,ded,ded,ded,ded,ded);
//     }
// void convert_to_vector_pts_uvt(
//     const std::array<double,8> &t_up, const std::array<double,8> &t_dw,
//     const std::array<double,8> &u_up, const std::array<double,8> &u_dw,
//     const std::array<double,8> &v_up, const std::array<double,8> &v_dw, 

//     __m512d& tu, __m512d& td, __m512d& uu, __m512d& ud,__m512d& vu, __m512d& vd
//     ){
//         tu=_mm512_setr_pd(t_up[0],t_up[1],t_up[2],t_up[3],t_up[4],t_up[5],t_up[6],t_up[7]);
//         td=_mm512_setr_pd(t_dw[0],t_dw[1],t_dw[2],t_dw[3],t_dw[4],t_dw[5],t_dw[6],t_dw[7]);

//         uu=_mm512_setr_pd(u_up[0],u_up[1],u_up[2],u_up[3],u_up[4],u_up[5],u_up[6],u_up[7]);
//         ud=_mm512_setr_pd(u_dw[0],u_dw[1],u_dw[2],u_dw[3],u_dw[4],u_dw[5],u_dw[6],u_dw[7]);

//         vu=_mm512_setr_pd(v_up[0],v_up[1],v_up[2],v_up[3],v_up[4],v_up[5],v_up[6],v_up[7]);
//         vd=_mm512_setr_pd(v_dw[0],v_dw[1],v_dw[2],v_dw[3],v_dw[4],v_dw[5],v_dw[6],v_dw[7]);
//     }
// //(a-b)/c
// __m512d function_test(__m512d a,__m512d b, __m512d c){
//     __m512d v=_mm512_sub_pd(a,b);
//     v= _mm512_div_pd(v,c);
//     return v;
// }

// void test(){
//     __m512d a=_mm512_setr_pd(3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0);
//     __m512d b=_mm512_setr_pd(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0);
//     __m512d c=_mm512_setr_pd(4.0,4.0,4.0,4.0,8.0,8.0,8.0,8.0);
//     __m512d rst=function_test(a,b,c);
//     double *ptr = (double*)&rst;
//     for(int i=0;i<8;i++){
//         //printf("p%i %f\n", i, ptr[i]);
// double test_(){
//     std::array<double,8> a0s;
// std::array<double,8> a1s;
// std::array<double,8> b0s;
// std::array<double,8> b1s;
// std::array<double,8> a0e;
// std::array<double,8> a1e;
// std::array<double,8> b0e;
// std::array<double,8> b1e;
// std::array<double,8> t_up;
// std::array<double,8> t_dw;
// std::array<double,8> u_up;
// std::array<double,8> u_dw;
// std::array<double,8> v_up;
// std::array<double,8> v_dw;

// auto x=function_ee(a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,t_up,t_dw,u_up,u_dw,v_up,v_dw);
// return x[0]+x[7];
// }
        
//     }
//     printf("%f %f %f %f %f %f %f %f\n",
//     ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
// }
}
// int main(){
//     double sum=0;
//     // for(long i=0;i<1e9;i++){
//     //     sum+=test_();
//     // }
    

//     return sum;
// }