// A root finder using interval arithmetic.
#include <interval_ccd/interval_root_finder.hpp>

#include <stack>
#include<igl/Timer.h>
#include<iostream>
#include<interval_ccd/Rational.hpp>
#include <interval_ccd/avx.h>
#include <queue>
#include<fstream>
// #define COMPARE_WITH_RATIONAL
// #define USE_TIMER
#define DEBUGING
namespace intervalccd {
double time20=0,time21=0,time22=0, time23=0,time24=0,time25=0,time_rational=0;
int refine=0;
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x)
{
    return interval_root_finder(
        f, [](const Interval&) { return true; }, x0, tol, x);
}

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x)
{
    Eigen::VectorX3I x0_vec = Eigen::VectorX3I::Constant(1, x0), x_vec;
    Eigen::VectorX3d tol_vec = Eigen::VectorX3d::Constant(1, tol);
    bool found_root = interval_root_finder(
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return Eigen::VectorX3I::Constant(1, f(x(0)));
        },
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        x0_vec, tol_vec, x_vec);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
    return true;
}
bool interval_root_finder_Redon(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x)
{
    Eigen::VectorX3I x0_vec = Eigen::VectorX3I::Constant(1, x0), x_vec;
    Eigen::VectorX3d tol_vec = Eigen::VectorX3d::Constant(1, tol);
    bool found_root = interval_root_finder(
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return Eigen::VectorX3I::Constant(1, f(x(0)));
        },
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        x0_vec, tol_vec, x_vec);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
    
}
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x, const bool check_vf)
{
    return interval_root_finder(
        f, [](const Eigen::VectorX3I&) { return true; }, x0, tol, x,check_vf);
}

inline Eigen::VectorX3d width(const Eigen::VectorX3I& x)
{
    Eigen::VectorX3d w(x.size());
    for (int i = 0; i < x.size(); i++) {
        w(i) = width(x(i));
    }
    return w;
}
// convert Numccd to double number
double Numccd2double(const Numccd& n){
    double r=double(n.first)/power(1,n.second);
    return r;
}
Eigen::VectorX3d width(const Interval3&x){
    Eigen::VectorX3d w;
    w.resize(3);
    for(int i=0;i<3;i++){
        w[i]=Numccd2double(x[i].second)-Numccd2double(x[i].first);
        assert(w[i]>=0);
    }
    return w;
    
}

std::array<Rational,3> width(const std::array<std::pair<Rational,Rational>, 3>&x){
    std::array<Rational,3> w;
    
    for(int i=0;i<3;i++){
        Rational sub=x[i].first-x[i].second;
        w[i]=sub>=0?sub:-sub;
        assert(w[i]>=0);
    }
    return w;
    
}
template <int dim, int max_dim = dim>
inline bool zero_in(Eigen::Vector<Interval, dim, max_dim> X)
{
    // Check if the origin is in the n-dimensional interval
    for (int i = 0; i < X.size(); i++) {
        if (!boost::numeric::zero_in(X(i))) {
            return false;
        }
    }
    return true;
}
// check if (i1,i2) overlaps {(u,v)|u+v<1,u>=0,v>=0}
// by checking if i1.lower()+i2.lower()<=1
bool interval_satisfy_constrain(const Interval &i1, const Interval &i2){
    Interval l1(i1.lower(),i1.lower());
    Interval l2(i2.lower(),i2.lower());
    Interval sum=l1+l2;
    if(sum.lower()>1)
    return false;
    else return true;
    }


bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x,const bool check_vf)
{
    // Stack of intervals and the last split dimension
    std::stack<std::pair<Eigen::VectorX3I, int>> xs;
    xs.emplace(x0, -1);
    while (!xs.empty()) {
        x = xs.top().first;
        int last_split = xs.top().second;
        xs.pop();

        Eigen::VectorX3I y = f(x);

        if (!zero_in(y)) {
            continue;
        }

        Eigen::VectorX3d widths = width(x);
        if ((widths.array() <= tol.array()).all()) {
            if (constraint_predicate(x)) {
                return true;
            }
            continue;
        }

        // Bisect the next dimension that is greater than its tolerance
        int split_i;
        for (int i = 1; i <= x.size(); i++) {
            split_i = (last_split + i) % x.size();
            if (widths(split_i) > tol(split_i)) {
                break;
            }
        }
        std::pair<Interval, Interval> halves = bisect(x(split_i));
        Eigen::VectorX3I x1 = x;
        // Push the second half on first so it is examined after the first half
        if(check_vf){
            if(split_i==1){
                if(interval_satisfy_constrain(halves.second,x(2))){
                    x(split_i) = halves.second;
                    xs.emplace(x, split_i);
                }
                if(interval_satisfy_constrain(halves.first,x(2))){
                    x(split_i) = halves.first;
                    xs.emplace(x, split_i);
                }
            }
            if(split_i==2){
                if(interval_satisfy_constrain(halves.second,x(1))){
                    x(split_i) = halves.second;
                    xs.emplace(x, split_i);
                }
                if(interval_satisfy_constrain(halves.first,x(1))){
                    x(split_i) = halves.first;
                    xs.emplace(x, split_i);
                }
            }
            if(split_i==0){
                x(split_i) = halves.second;
                xs.emplace(x, split_i);
                x(split_i) = halves.first;
                xs.emplace(x, split_i);
            }
        }
        else{
            x(split_i) = halves.second;
            xs.emplace(x, split_i);
            x(split_i) = halves.first;
            xs.emplace(x, split_i);
        }
        
        
    }
    return false;
}
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x){
        return interval_root_finder(f,constraint_predicate,x0,tol,x,false);
    }

// check if 0 is in interval
bool interval_bounding_box_check(const Eigen::Vector3I&in, std::array<bool,6>& flag){
    
    for(int i=0;i<3;i++){
        if(!flag[2*i]){
            if(in(i).lower()<=0){
                flag[2*i]=true;
            }
        }
        if(!flag[2*i+1]){
            if(in(i).upper()>=0){
                flag[2*i+1]=true;
            }
        }
    }
    if(flag[0]&&flag[1]&&flag[2]&&flag[3]&&flag[4]&&flag[5]) 
        return true;
    else return false;
}

std::array<Eigen::Vector3d,2> bbd_4_pts(const std::array<Eigen::Vector3d,4>& pts){
    Eigen::Vector3d min,max;
    min=pts[0];max=pts[0];
    for(int i=1;i<4;i++){
        for(int j=0;j<3;j++){
            if(min[j]>pts[i][j]){
                min[j]=pts[i][j];
            }
            if(max[j]<pts[i][j]){
                max[j]=pts[i][j];
            }
        }
    }
    std::array<Eigen::Vector3d,2> rst;
    rst[0]=min;rst[1]=max;
    return rst;

}
std::array<Eigen::Vector3d,2> bbd_6_pts(const std::array<Eigen::Vector3d,6>& pts){
    Eigen::Vector3d min,max;
    min=pts[0];max=pts[0];
    for(int i=1;i<6;i++){
        for(int j=0;j<3;j++){
            if(min[j]>pts[i][j]){
                min[j]=pts[i][j];
            }
            if(max[j]<pts[i][j]){
                max[j]=pts[i][j];
            }
        }
    }
    std::array<Eigen::Vector3d,2> rst;
    rst[0]=min;rst[1]=max;
    return rst;

}

// eps is the interval [-eps,eps] we need to check
// if [-eps,eps] overlap, return true
template<typename T>
bool evaluate_bbox_one_dimension(
    const std::array<Numccd,2>& t,
    const std::array<Numccd,2>& u,
    const std::array<Numccd,2>& v,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const int dimension,T tp,const bool check_vf,const double eps){
#ifdef USE_TIMER
    igl::Timer timer;
#endif
    
/*   
{// smart way but no possibility of parallazation
    double eva;
    bool flag0=false, flag1=false;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                if(!check_vf){
#ifdef USE_TIMER
                    timer.start();
#endif
                    eva=function_f_ee(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif
                }
                else{
#ifdef USE_TIMER
                    timer.start();
#endif
                    eva=function_f_vf(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif
                }
                if(eva<=eps&&eva>=-eps){
                    return true;
                }
                if(eva<-eps){
                    flag0=true;
                }
                if(eva>eps){
                    flag1=true;
                }
                if(flag0&&flag1){
                    return true;
                }
            }
        }
    }
    if(flag0&&flag1)
    return true;
    return false;
}        
*/
if(check_vf){// test
    std::array<double,8> vs;
    int count=0;
#ifdef USE_TIMER
    timer.start();
#endif    
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){

                vs[count]=function_f_vf(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                
                count++;

            }
        }
    }
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif

    double minv=vs[0], maxv=vs[0];
    
    for(int i=1;i<8;i++){
        if(minv>vs[i]){
            minv=vs[i];
        }
        if(maxv<vs[i]){
            maxv=vs[i];
        }
    }
    if(minv>eps||maxv<-eps)
        return false;
    return true;

    }// test end
    
     else{// test
    std::array<double,8> vs;
    int count=0;
#ifdef USE_TIMER
    timer.start();
#endif    
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){

                vs[count]=function_f_ee(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                count++;

            }
        }
    }
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif

    double minv=vs[0], maxv=vs[0];
    
    for(int i=1;i<8;i++){
        if(minv>vs[i]){
            minv=vs[i];
        }
        if(maxv<vs[i]){
            maxv=vs[i];
        }
    }
    if(minv>eps||maxv<-eps)
        return false;
    return true;

    }// test end

}
// eps is the interval [-eps,eps] we need to check
// if [-eps,eps] overlap, return true
// bbox_in_eps tell us if the box is totally in eps box
bool evaluate_bbox_one_dimension_vector(
    std::array<double,8>& t_up,std::array<double,8>&t_dw,
    std::array<double,8>& u_up,std::array<double,8>&u_dw,
    std::array<double,8>& v_up,std::array<double,8>&v_dw,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const int dimension,const bool check_vf,const double eps, bool& bbox_in_eps){
#ifdef USE_TIMER
    igl::Timer timer;
#endif
    std::array<double,8> vs;
    int count=0;
    bbox_in_eps=false;
#ifdef USE_TIMER
    timer.start();
#endif  
    if(check_vf){
        vs=function_vf(a0s[dimension],a1s[dimension],b0s[dimension],b1s[dimension],a0e[dimension],a1e[dimension],b0e[dimension],b1e[dimension],t_up,t_dw,
        u_up,u_dw,v_up,v_dw);
    }
    else{
        vs=function_ee(a0s[dimension],a1s[dimension],b0s[dimension],b1s[dimension],a0e[dimension],a1e[dimension],b0e[dimension],b1e[dimension],t_up,t_dw,
        u_up,u_dw,v_up,v_dw);
    }
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
                    
#endif
    double minv=vs[0], maxv=vs[0];
    
    for(int i=1;i<8;i++){
        if(minv>vs[i]){
            minv=vs[i];
        }
        if(maxv<vs[i]){
            maxv=vs[i];
        }
    }
    if(minv>eps||maxv<-eps)
        return false;
    if(minv>=-eps&&maxv<=eps){
        bbox_in_eps=true;
    }
    return true;

}

// the bounding boxes generated are t0, t1, u, u1, v0, v1 boxes
template<typename T>
void evaluate_tuv_bboxes(
    const std::array<Numccd,2>& t,
    const std::array<Numccd,2>& u,
    const std::array<Numccd,2>& v,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    T tp,const bool check_vf,std::array<std::array<Eigen::Vector3d,2>,6>& bboxes){
    

    int count=0;
    std::array<Eigen::Vector3d,8> pts;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                if(!check_vf){
                    pts[count][0]=
                    function_f_ee(t[i],u[j],v[k],tp,0,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);

                    pts[count][1]=
                    function_f_ee(t[i],u[j],v[k],tp,1,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                    
                    pts[count][2]=
                    function_f_ee(t[i],u[j],v[k],tp,2,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                else{
                    pts[count][0]=
                    function_f_vf(t[i],u[j],v[k],tp,0,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);

                    pts[count][1]=
                    function_f_vf(t[i],u[j],v[k],tp,1,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                    
                    pts[count][2]=
                    function_f_vf(t[i],u[j],v[k],tp,2,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                count++;
                
            }
        }
    }
    // now the parameters of pts are:
    // 000, 001, 010, 011, 100, 101, 110, 111
    std::array<Eigen::Vector3d,4> bps;
    bps[0]=pts[0];
    bps[1]=pts[1];
    bps[2]=pts[2];
    bps[3]=pts[3];
    bboxes[0]=bbd_4_pts(bps);//t0

    bps[0]=pts[4];
    bps[1]=pts[5];
    bps[2]=pts[6];
    bps[3]=pts[7];
    bboxes[1]=bbd_4_pts(bps);//t1

    bps[0]=pts[0];
    bps[1]=pts[1];
    bps[2]=pts[4];
    bps[3]=pts[5];
    bboxes[2]=bbd_4_pts(bps);//u0

    bps[0]=pts[2];
    bps[1]=pts[3];
    bps[2]=pts[6];
    bps[3]=pts[7];
    bboxes[3]=bbd_4_pts(bps);//u1

    bps[0]=pts[0];
    bps[1]=pts[2];
    bps[2]=pts[4];
    bps[3]=pts[6];
    bboxes[4]=bbd_4_pts(bps);//v0

    bps[0]=pts[1];
    bps[1]=pts[3];
    bps[2]=pts[5];
    bps[3]=pts[7];
    bboxes[5]=bbd_4_pts(bps);//v1
        
}


Vector3r function_f_ee_Rational (
const Numccd&tpara, const Numccd&upara, const Numccd&vpara, 
const Eigen::Vector3d& a0sd,
    const Eigen::Vector3d& a1sd,
    const Eigen::Vector3d& b0sd,
    const Eigen::Vector3d& b1sd,
    const Eigen::Vector3d& a0ed,
    const Eigen::Vector3d& a1ed,
    const Eigen::Vector3d& b0ed,
    const Eigen::Vector3d& b1ed ) {
       
       long tu = tpara.first; int td=tpara.second;// t=tu/(2^td)
       long uu = upara.first; int ud=upara.second;
       long vu = vpara.first; int vd=vpara.second;
        Vector3r 
        a0s(a0sd[0],a0sd[1],a0sd[2]), 
        a1s(a1sd[0],a1sd[1],a1sd[2]),
        b0s(b0sd[0],b0sd[1],b0sd[2]),
        b1s(b1sd[0],b1sd[1],b1sd[2]),
        a0e(a0ed[0],a0ed[1],a0ed[2]),
        a1e(a1ed[0],a1ed[1],a1ed[2]),
        b0e(b0ed[0],b0ed[1],b0ed[2]),
        b1e(b1ed[0],b1ed[1],b1ed[2]);
       Vector3r edge0_vertex0
           = (a0e - a0s) * tu/power(1,td)
           + a0s;
       Vector3r edge0_vertex1
           = (a1e - a1s) * tu/power(1,td)
           + a1s;
       Vector3r edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/power(1,ud) + edge0_vertex0;

       Vector3r edge1_vertex0
           = (b0e - b0s) * tu/power(1,td)
           + b0s;
       Vector3r edge1_vertex1
           = (b1e - b1s) * tu/power(1,td)
           + b1s;
       Vector3r edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/power(1,vd) + edge1_vertex0;

       
       return edge1_vertex-edge0_vertex;
       
}
Vector3r function_f_ee_Rational (
const Rational&tpara, const Rational&upara, const Rational&vpara, 
const Eigen::Vector3d& a0sd,
    const Eigen::Vector3d& a1sd,
    const Eigen::Vector3d& b0sd,
    const Eigen::Vector3d& b1sd,
    const Eigen::Vector3d& a0ed,
    const Eigen::Vector3d& a1ed,
    const Eigen::Vector3d& b0ed,
    const Eigen::Vector3d& b1ed ) {
       
        Vector3r 
        a0s(a0sd[0],a0sd[1],a0sd[2]), 
        a1s(a1sd[0],a1sd[1],a1sd[2]),
        b0s(b0sd[0],b0sd[1],b0sd[2]),
        b1s(b1sd[0],b1sd[1],b1sd[2]),
        a0e(a0ed[0],a0ed[1],a0ed[2]),
        a1e(a1ed[0],a1ed[1],a1ed[2]),
        b0e(b0ed[0],b0ed[1],b0ed[2]),
        b1e(b1ed[0],b1ed[1],b1ed[2]);

        // Vector3r las = (1-upara)*a0s+upara*a1s;
        // //std::cout<<"las, "<<las[0]<<", "<<las[1]<<", "<<las[2]<<std::endl;
        // Vector3r lae = (1-upara)*a0e+upara*a1e;
        // Vector3r lbs = (1-vpara)*b0s+vpara*b1s;
        // Vector3r lbe = (1-vpara)*b0e+vpara*b1e;
        // Vector3r lla=(1-tpara)*las + tpara*lae;
        // Vector3r llb=(1-tpara)*lbs + tpara*lbe;

        // //std::cout<<"lae, "<<lae[0]<<", "<<lae[1]<<", "<<lae[2]<<std::endl;
        // //std::cout<<"lbs, "<<lbs[0]<<", "<<lbs[1]<<", "<<lbs[2]<<std::endl;
        // //std::cout<<"lbe, "<<lbe[0]<<", "<<lbe[1]<<", "<<lbe[2]<<std::endl;
        // //std::cout<<"lla, "<<lla[0]<<", "<<lla[1]<<", "<<lla[2]<<std::endl;
        // //std::cout<<"llb, "<<llb[0]<<", "<<llb[1]<<", "<<llb[2]<<std::endl<<std::endl;
        // return lla-llb;





       Vector3r edge0_vertex0
           = (a0e - a0s) * tpara
           + a0s;
       Vector3r edge0_vertex1
           = (a1e - a1s) * tpara
           + a1s;
       Vector3r edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * upara + edge0_vertex0;

       Vector3r edge1_vertex0
           = (b0e - b0s) * tpara
           + b0s;
       Vector3r edge1_vertex1
           = (b1e - b1s) * tpara
           + b1s;
       Vector3r edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vpara + edge1_vertex0;

       
       return edge1_vertex-edge0_vertex;
       
}
Vector3r function_f_vf_Rational (
const Numccd&tpara, const Numccd&upara, const Numccd&vpara, 
const Eigen::Vector3d& a0sd,
    const Eigen::Vector3d& a1sd,
    const Eigen::Vector3d& b0sd,
    const Eigen::Vector3d& b1sd,
    const Eigen::Vector3d& a0ed,
    const Eigen::Vector3d& a1ed,
    const Eigen::Vector3d& b0ed,
    const Eigen::Vector3d& b1ed ) {
       
       long tu = tpara.first; int td=tpara.second;// t=tu/(2^td)
       long uu = upara.first; int ud=upara.second;
       long vu = vpara.first; int vd=vpara.second;
        Vector3r 
        vs(a0sd[0],a0sd[1],a0sd[2]), 
        t0s(a1sd[0],a1sd[1],a1sd[2]),
        t1s(b0sd[0],b0sd[1],b0sd[2]),
        t2s(b1sd[0],b1sd[1],b1sd[2]),

        ve(a0ed[0],a0ed[1],a0ed[2]),
        t0e(a1ed[0],a1ed[1],a1ed[2]),
        t1e(b0ed[0],b0ed[1],b0ed[2]),
        t2e(b1ed[0],b1ed[1],b1ed[2]);

        Vector3r v=(ve-vs)*tu/power(1,td)+vs;

        Vector3r t0=(t0e-t0s)*tu/power(1,td)+t0s;
        Vector3r t1=(t1e-t1s)*tu/power(1,td)+t1s;
        Vector3r t2=(t2e-t2s)*tu/power(1,td)+t2s;
        Vector3r p=(t1-t0)*uu/power(1,ud)+(t2-t0)*vu/power(1,vd)+t0;
        return v-p;
       
}
Vector3r function_f_vf_Rational (
const Rational&tpara, const Rational&upara, const Rational&vpara, 
const Eigen::Vector3d& a0sd,
    const Eigen::Vector3d& a1sd,
    const Eigen::Vector3d& b0sd,
    const Eigen::Vector3d& b1sd,
    const Eigen::Vector3d& a0ed,
    const Eigen::Vector3d& a1ed,
    const Eigen::Vector3d& b0ed,
    const Eigen::Vector3d& b1ed ) {

        Vector3r 
        vs(a0sd[0],a0sd[1],a0sd[2]), 
        t0s(a1sd[0],a1sd[1],a1sd[2]),
        t1s(b0sd[0],b0sd[1],b0sd[2]),
        t2s(b1sd[0],b1sd[1],b1sd[2]),

        ve(a0ed[0],a0ed[1],a0ed[2]),
        t0e(a1ed[0],a1ed[1],a1ed[2]),
        t1e(b0ed[0],b0ed[1],b0ed[2]),
        t2e(b1ed[0],b1ed[1],b1ed[2]);

        Vector3r v=(ve-vs)*tpara+vs;

        Vector3r t0=(t0e-t0s)*tpara+t0s;
        Vector3r t1=(t1e-t1s)*tpara+t1s;
        Vector3r t2=(t2e-t2s)*tpara+t2s;
        Vector3r p=(t1-t0)*upara+(t2-t0)*vpara+t0;
        return v-p;
       
}
bool Origin_in_function_bounding_box(
    const Interval3& paras,
    const std::function<Eigen::VectorX3I(const Numccd&, const Numccd&, const Numccd&)>& f){
    igl::Timer timer;
    std::array<Numccd,2> t,u,v;
    
    t[0]=paras[0].first;
    t[1]=paras[0].second;
    u[0]=paras[1].first;
    u[1]=paras[1].second;
    v[0]=paras[2].first;
    v[1]=paras[2].second;
    Eigen::Vector3I result;
    std::array<bool,6> flag;// when all flags are true, return true;
    for(int i=0;i<6;i++){
        flag[i]=false;
    }
    
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                timer.start();
                result=f(t[i],u[j],v[k]);
                timer.stop();
                time23+=timer.getElapsedTimeInMicroSec();
                timer.start();
                bool check=interval_bounding_box_check(result,flag);
                timer.stop();
                //time24+=timer.getElapsedTimeInMicroSec();
                if(check)
                    return true;
            }
        }
    }
    return false;
}

bool Origin_in_function_bounding_box_double(
    const Interval3& paras,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const bool check_vf,
    const std::array<double,3>& box){
#ifdef USE_TIMER
    igl::Timer timer;
#endif
    std::array<Numccd,2> t,u,v;
    
    t[0]=paras[0].first;
    t[1]=paras[0].second;
    u[0]=paras[1].first;
    u[1]=paras[1].second;
    v[0]=paras[2].first;
    v[1]=paras[2].second;
    //bool zero_0=false, zer0_1=false, zero_2=false;
    double input_type;
    bool ck;
#ifdef USE_TIMER
    timer.start();
#endif
    ck=evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,0,input_type,check_vf,box[0]);
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
    if(!ck)
        return false;
#ifdef USE_TIMER
    timer.start();
#endif
    ck=evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,1,input_type,check_vf,box[1]);
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
    if(!ck)
        return false;
#ifdef USE_TIMER
    timer.start();
#endif
    ck=evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,2,input_type,check_vf,box[2]);
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
    if(!ck)
        return false;
    return true;
    
}
bool Origin_in_function_bounding_box_double_vector(
    const Interval3& paras,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const bool check_vf,
    const std::array<double,3>& box, bool &box_in_eps){
#ifdef USE_TIMER
    igl::Timer timer;
#endif  
    box_in_eps=false;
    std::array<double,8> t_up;std::array<double,8>t_dw;
    std::array<double,8> u_up;std::array<double,8>u_dw;
    std::array<double,8> v_up;std::array<double,8>v_dw;
#ifdef USE_TIMER
    timer.start();
#endif
    convert_tuv_to_array(paras,t_up,t_dw,u_up,u_dw,v_up,v_dw);
#ifdef USE_TIMER
    timer.stop();
    time24+=timer.getElapsedTimeInMicroSec();
#endif
    bool ck;
    bool box_in[3];
    for(int i=0;i<3;i++){
#ifdef USE_TIMER
    timer.start();
#endif
        ck=evaluate_bbox_one_dimension_vector(t_up,t_dw,u_up,u_dw,v_up,v_dw,
    a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,i,check_vf,box[i],box_in[i]);
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
     if(!ck)
        return false;
    }
    if(box_in[0]&&box_in[1]&&box_in[2]){
        box_in_eps=true;
    }
    return true;
    
    
}
bool bounding_box_intersection(const Eigen::Vector3d &pmin,const Eigen::Vector3d &pmax,
const Eigen::Vector3d &qmin,const Eigen::Vector3d &qmax){
    if(pmax[0]<qmin[0]||pmax[1]<qmin[1]||pmax[2]<qmin[2]){
        return false;
    }
    if(qmax[0]<pmin[0]||qmax[1]<pmin[1]||qmax[2]<pmin[2]){
        return false;
    }
    return true;
}
bool estimate_tuv_through_bbox(
    const Interval3& paras,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const bool check_vf,
    const std::array<double,3>& box){
    //igl::Timer timer;
    std::array<Numccd,2> t,u,v;
    
    t[0]=paras[0].first;
    t[1]=paras[0].second;
    u[0]=paras[1].first;
    u[1]=paras[1].second;
    v[0]=paras[2].first;
    v[1]=paras[2].second;
    //bool zero_0=false, zer0_1=false, zero_2=false;
    double input_type;
    std::array<std::array<Eigen::Vector3d,2>,6> bboxes;
    // get the bounding boxes 
    evaluate_tuv_bboxes(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,input_type,check_vf,bboxes);
     //TODO
    if(!evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,0,input_type,check_vf,box[0]))
        return false;
    if(!evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,1,input_type,check_vf,box[1]))
        return false;
    if(!evaluate_bbox_one_dimension(t,u,v,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,2,input_type,check_vf,box[2]))
        return false;
    return true;
    
}
bool Origin_in_function_bounding_box_Rational(
    const Interval3& paras,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,const bool check_vf){
    //igl::Timer timer;
    std::array<Numccd,2> t,u,v;
    
    t[0]=paras[0].first;
    t[1]=paras[0].second;
    u[0]=paras[1].first;
    u[1]=paras[1].second;
    v[0]=paras[2].first;
    v[1]=paras[2].second;
    //std::cout<<"t, ["<<t[0].first/pow(2,t[0].second)<<","<<t[1].first/pow(2,t[1].second)<<"]"<<std::endl;
    //std::cout<<"u, ["<<u[0].first/pow(2,u[0].second)<<","<<u[1].first/pow(2,u[1].second)<<"]"<<std::endl;
    //std::cout<<"v, ["<<v[0].first/pow(2,v[0].second)<<","<<v[1].first/pow(2,v[1].second)<<"]"<<std::endl<<std::endl;
    Vector3r minv, maxv;
    std::array<Vector3r,8> pts;
    int c=0;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                //std::cout<<"c="<<c<<std::endl;
                if(!check_vf){
                    //std::cout<<"ee"<<std::endl;
                    pts[c]=function_f_ee_Rational(t[i],u[j],v[k],a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                else{
                    pts[c]=function_f_vf_Rational(t[i],u[j],v[k],a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                
                c++;
            }
        }
    }
    minv=pts[0]; maxv=pts[0];
    for(int i=0;i<8;i++){
        if(minv[0]>pts[i][0]){
            minv[0]=pts[i][0];
        }
        if(minv[1]>pts[i][1]){
            minv[1]=pts[i][1];
        }
        if(minv[2]>pts[i][2]){
            minv[2]=pts[i][2];
        }
        if(maxv[0]<pts[i][0]){
            maxv[0]=pts[i][0];
        }
        if(maxv[1]<pts[i][1]){
            maxv[1]=pts[i][1];
        }
        if(maxv[2]<pts[i][2]){
            maxv[2]=pts[i][2];
        }
    }
    if(minv[0]<=0&&minv[1]<=0&&minv[2]<=0)
    {
        if(maxv[0]>=0&&maxv[1]>=0&&maxv[2]>=0){
            return true;
        }
    }
    return false;
    
    
}

bool Origin_in_function_bounding_box_Rational(
    const std::array<std::pair<Rational,Rational>, 3>& paras,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,const bool check_vf){
    //igl::Timer timer;
    std::array<Rational,2> t,u,v;
    
    t[0]=paras[0].first;
    t[1]=paras[0].second;
    u[0]=paras[1].first;
    u[1]=paras[1].second;
    v[0]=paras[2].first;
    v[1]=paras[2].second;
    Vector3r minv, maxv;
    std::array<Vector3r,8> pts;
    int c=0;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                
                if(!check_vf){
                    pts[c]=function_f_ee_Rational(t[i],u[j],v[k],a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                else{
                    pts[c]=function_f_vf_Rational(t[i],u[j],v[k],a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
                }
                
                c++;
            }
        }
    }
    minv=pts[0]; maxv=pts[0];
    for(int i=0;i<8;i++){
        if(minv[0]>pts[i][0]){
            minv[0]=pts[i][0];
        }
        if(minv[1]>pts[i][1]){
            minv[1]=pts[i][1];
        }
        if(minv[2]>pts[i][2]){
            minv[2]=pts[i][2];
        }
        if(maxv[0]<pts[i][0]){
            maxv[0]=pts[i][0];
        }
        if(maxv[1]<pts[i][1]){
            maxv[1]=pts[i][1];
        }
        if(maxv[2]<pts[i][2]){
            maxv[2]=pts[i][2];
        }
    }
    if(minv[0]<=0&&minv[1]<=0&&minv[2]<=0)
    {
        if(maxv[0]>=0&&maxv[1]>=0&&maxv[2]>=0){
            return true;
        }
    }
    return false;
    
    
}
// return power t. n=result*2^t
long reduction(const long n, long& result){
    int t=0;
    int newn=n;
    while (newn%2==0){
        newn=newn/2;
        t++;
    }
    result=newn;
    return t;
}
std::pair<Singleinterval, Singleinterval> bisect(const Singleinterval& inter){
    Numccd low=inter.first;
    Numccd up=inter.second;

    // interval is [k1/pow(2,n1), k2/pow(2,n2)], k1,k2,n1,n2 are all not negative
    long k1=low.first;
    int n1=low.second;
    long k2=up.first;
    int n2=up.second;

    assert(k1 >= 0 && k2 >= 0 && n1 >= 0 && n2 >= 0);

    std::pair<Singleinterval, Singleinterval> result;
    long  k;
    int n; 
    int p;
    long r;
    if(n2==n1){
        p=reduction(k1+k2,r);
        k=r;
        n=n2-p+1;
    }
    if(n2>n1){
        k=k1*power(1,n2-n1)+k2; assert(k%2==1);
        n=n2+1;
    }
    if(n2<n1){
        k=k1+k2*power(1,n1-n2); assert(k%2==1);
        n=n1+1;
    }
    Numccd newnum(k,n);
    Singleinterval i1(low,newnum),i2(newnum,up);
    // std::cout<<"low,"<<Numccd2double(low)<<",up,"<<Numccd2double(up)<<", num, "<<Numccd2double(newnum)<<std::endl;
    // std::cout<<"new, k1, "<<newnum.first<<", n1, "<<newnum.second<<std::endl;
    assert(Numccd2double(newnum)>Numccd2double(low)&&Numccd2double(newnum)<Numccd2double(up));
    result.first=i1;result.second=i2;
    return result;
}
std::pair<std::pair<Rational,Rational>, std::pair<Rational,Rational>> bisect(const std::pair<Rational,Rational>& inter){
    std::pair<std::pair<Rational,Rational>, std::pair<Rational,Rational>> result;
    std::pair<Rational,Rational> single;
    Rational mid=(inter.first+inter.second)/2;
    single.first=inter.first;single.second=mid;
    result.first=single;
    single.first=mid;single.second=inter.second;
    result.second=single;
    return result;
}
bool sum_no_larger_1(const Numccd& num1, const Numccd& num2){
    long k1=num1.first;
    int n1=num1.second;
    long k2=num2.first;
    int n2=num2.second;
    long k; int n;
    if(n1==n2){
        k=k1+k2;
        n=n1;
    }
    if(n1<n2){
        k=power(1,n2-n1)*k1+k2;
        n=n2;
    }
    if(n1>n2){
        k=power(1,n1-n2)*k2+k1;
        n=n1;
    }
    assert(k>=0&&n>=0);
    if(k>power(1,n)) return false;
    else return true;

}
//check if num1<num2
bool less_than(const Numccd& num1, const Numccd& num2){
    long k1=num1.first;
    int n1=num1.second;
    long k2=num2.first;
    int n2=num2.second;
    
    if(n1<n2){
        k1=power(1,n2-n1)*k1;
    }
    if(n1>n2){
        k2=power(1,n1-n2)*k2;
    }
    if(k1<k2) return true;
    return false;
}
bool sum_no_larger_1_Rational(const Numccd& num1, const Numccd& num2){
    long k1=num1.first;
    int n1=num1.second;
    long k2=num2.first;
    int n2=num2.second;
    Rational nbr1,nbr2;
    nbr1=Rational(k1)/Rational(power(1,n1));
    nbr2=Rational(k2)/Rational(power(1,n2));
    Rational rst=nbr1+nbr2-Rational(1);
    if(rst>0) return false;
    else return true;
    
}
bool interval_root_finder_opt(
    const std::function<Eigen::VectorX3I(const Numccd&, const Numccd&, const Numccd&)>& f,
    //const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    //const Eigen::VectorX3I& x0,// initial interval, must be [0,1]x[0,1]x[0,1]
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    Interval3& final,
    const bool check_vf){
    
    Numccd low_number; low_number.first=0; low_number.second=0;// low_number=0;
    Numccd up_number; up_number.first=1; up_number.second=0;// up_number=1;
    // initial interval [0,1]
    Singleinterval init_interval;init_interval.first=low_number;init_interval.second=up_number;
    //build interval set [0,1]x[0,1]x[0,1]
    Interval3 iset;
    iset[0]=init_interval;iset[1]=init_interval;iset[2]=init_interval;
    // Stack of intervals and the last split dimension
    std::stack<std::pair<Interval3,int>> istack;
    istack.emplace(iset,-1);
    
    // current intervals
    Interval3 current;
    if(check_vf){
        std::cout<<"NOT IMPLEMENTED, DO NOT USE THIS FUNCTION"<<std::endl;
        return false;
    }
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();
#ifdef USE_TIMER
        igl::Timer timer;
        timer.start();
#endif
        bool zero_in = Origin_in_function_bounding_box(current,f);
        // TODO need to add vf function check here 
#ifdef USE_TIMER
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
        if(!zero_in) continue;
#ifdef USE_TIMER
        timer.start();
#endif
        Eigen::VectorX3d widths = width(current);
#ifdef USE_TIMER
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        if ((widths.array() <= tol.array()).all()) {
            final=current;
                return true;
        }

        // Bisect the next dimension that is greater than its tolerance
        int split_i;
        for (int i = 1; i <= 3; i++) {
            split_i = (last_split + i) % 3;
            if (widths(split_i) > tol(split_i)) {
                break;
            }
        }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);

        if(check_vf){
            if(split_i==1){
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){

                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
        }
        else{
            current[split_i] = halves.second;
            istack.emplace(current, split_i);
            current[split_i] = halves.first;
            istack.emplace(current, split_i);
        }
    }
    return false;
    
}
// calculate the sign of f. dim is the dimension we are evaluating.
template<typename T>
T function_f_ee (
const Numccd&tpara, const Numccd&upara, const Numccd&vpara,const T& type, const int dim,
const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e ) {
       
       long tu = tpara.first; int td=tpara.second;// t=tu/(2^td)
       long uu = upara.first; int ud=upara.second;
       long vu = vpara.first; int vd=vpara.second;

       T edge0_vertex0
           = (T(a0e[dim]) - T(a0s[dim])) * tu/power(1,td)
           + T(a0s[dim]);
       T edge0_vertex1
           = (T(a1e[dim]) - T(a1s[dim])) * tu/power(1,td)
           + T(a1s[dim]);
       T edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/power(1,ud) + edge0_vertex0;

       T edge1_vertex0
           = (T(b0e[dim]) - T(b0s[dim])) * tu/power(1,td)
           + T(b0s[dim]);
       T edge1_vertex1
           = (T(b1e[dim]) - T(b1s[dim])) * tu/power(1,td)
           + T(b1s[dim]);
       T edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/power(1,vd) + edge1_vertex0;

       
       return edge1_vertex-edge0_vertex;
       
}

template<typename T>
T function_f_vf (
const Numccd&tpara, const Numccd&upara, const Numccd&vpara,const T& type, const int dim,
    const Eigen::Vector3d& vs,
    const Eigen::Vector3d& t0s,
    const Eigen::Vector3d& t1s,
    const Eigen::Vector3d& t2s,

    const Eigen::Vector3d& ve,
    const Eigen::Vector3d& t0e,
    const Eigen::Vector3d& t1e,
    const Eigen::Vector3d& t2e ) {
       
       long tu = tpara.first; int td=tpara.second;// t=tu/(2^td)
       long uu = upara.first; int ud=upara.second;
       long vu = vpara.first; int vd=vpara.second;

        T v= (T(ve[dim])-T(vs[dim]))* tu/power(1,td)+T(vs[dim]);
        T t0=(T(t0e[dim])-T(t0s[dim]))* tu/power(1,td)+T(t0s[dim]);
        T t1=(T(t1e[dim])-T(t1s[dim]))* tu/power(1,td)+T(t1s[dim]);
        T t2=(T(t2e[dim])-T(t2s[dim]))* tu/power(1,td)+T(t2s[dim]);
        T p=(t1-t0)*uu/power(1,ud)+(t2-t0)*vu/power(1,vd)+t0;
        return v-p;
       
}


bool interval_root_finder_double_(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    // Interval3& final,
    double& toi,
    const bool check_vf,
    const std::array<double,3> err,
    const double ms,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e){
    auto cmp = [](std::pair<Interval3,int> i1, std::pair<Interval3,int> i2) { 
        return !(less_than(i1.first[0].first, i2.first[0].first));
        };
    Numccd low_number; low_number.first=0; low_number.second=0;// low_number=0;
    Numccd up_number; up_number.first=1; up_number.second=0;// up_number=1;
    // initial interval [0,1]
    Singleinterval init_interval;init_interval.first=low_number;init_interval.second=up_number;
    //build interval set [0,1]x[0,1]x[0,1]
    Interval3 iset;
    iset[0]=init_interval;iset[1]=init_interval;iset[2]=init_interval;
    // Stack of intervals and the last split dimension
    // std::stack<std::pair<Interval3,int>> istack;
    std::priority_queue<std::pair<Interval3,int>, std::vector<std::pair<Interval3,int>>, decltype(cmp)> istack(cmp);
    istack.emplace(iset,-1);

    // current intervals
    Interval3 current;
    std::array<double,3> err_and_ms;
    err_and_ms[0]=err[0]+ms;
    err_and_ms[1]=err[1]+ms;
    err_and_ms[2]=err[2]+ms;
    refine=0;
    toi=std::numeric_limits<double>::infinity();
    Numccd TOI; TOI.first=1;TOI.second=0;
    //std::array<double,3> 
    bool collision=false;
    int rnbr=0;
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();
        // if(rnbr>0&&less_than( current[0].first,TOI)){
        //     std::cout<<"not the first"<<std::endl;
        //     // continue;
        // }
        if(!less_than(current[0].first,TOI)){
            // std::cout<<"not the first"<<std::endl;
            continue;
        }
        //TOI should always be no larger than current
        
            
        // if(Numccd2double(current[0].first)>=Numccd2double(TOI)){
        //     std::cout<<"here wrong, comparing"<<std::endl;
        // } 
#ifdef USE_TIMER
        igl::Timer timer;

        timer.start();
#endif
        refine++;
        bool zero_in;
       bool box_in;
        zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,box_in);
        
#ifdef USE_TIMER
    
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL// this is defined in the begining of this file
        
        zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);

#endif
        if(!zero_in) continue;
#ifdef USE_TIMER
        timer.start();
#endif
        Eigen::VectorX3d widths = width(current);
#ifdef USE_TIMER
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        if ((widths.array() <= tol.array()).all()) {
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI);
            return true;
        }
        if(box_in){
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI);
            return true;
        }

        std::array<bool , 3> check;
        Eigen::VectorX3d widthratio;
        widthratio.resize(3);
        check[0]=false;check[1]=false; check[2]=false;
        for(int i=0;i<3;i++){
            widthratio(i)=widths(i)/tol(i);
            if(widths(i) > tol(i))
                check[i]=true;// means this need to be checked
        }
        
        int split_i=-1;
        for(int i=0;i<3;i++){
            if(check[i]){
                if(check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)&&widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }    
                }
                if(check[(i+1)%3]&&!check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&!check[(i+2)%3]){
                   
                        split_i=i;
                        break;
                    
                }
            }
        }
        if(split_i<0){
            std::cout<<"ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT"<<std::endl;
        }
        // Bisect the next dimension that is greater than its tolerance
        // int split_i;
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
        if(!less_than(halves.first.first, halves.first.second)){
                std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
                return true;
            }
        if(!less_than(halves.second.first, halves.second.second)){
            std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
            return true;
        }
        if(check_vf){
            //std::cout<<"*** check_vf"<<std::endl;
            if(split_i==1){
               // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
               // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));
                
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){
                //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));
                
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
            
        }
        else{
            current[split_i] = halves.second;
            istack.emplace(current, split_i);
            current[split_i] = halves.first;
            istack.emplace(current, split_i);
        }

    }
    if(collision) toi=Numccd2double(TOI);
    // if(toi==0)std::cout<<"infinate roots, "<<std::endl;
    // if(rnbr>5) std::cout<<"nbr of roots, "<<rnbr<<", time, "<<toi<<std::endl;
    return collision;
    return false;
    
}


void sum_up(const Numccd& nbr1, const Numccd& nbr2, Numccd& result){
    long k1=nbr1.first;
    long k2=nbr2.first;
    int n1=nbr1.second;
    int n2=nbr2.second;
    long  k;
    int n; 
    int p;
    long r;
    if(n2==n1){
        p=reduction(k1+k2,r);
        k=r;
        n=n2-p;
    }
    if(n2>n1){
        k=k1*power(1,n2-n1)+k2; assert(k%2==1);
        n=n2;
    }
    if(n2<n1){
        k=k1+k2*power(1,n1-n2); assert(k%2==1);
        n=n1;
    }
    result.first=k;
    result.second=n;
    return;
}
// find a t value that t<tol
void t_tol_width( Numccd&x, const double tol){
    x.first=1;x.second=0;
    while (Numccd2double(x)>=tol){
        x.second+=1;
    }
  
}
bool interval_root_finder_double(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    // Interval3& final,
    double& toi,
    const bool check_vf,
    const std::array<double,3> err,
    const double ms,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e){
    auto cmp = [](std::pair<Interval3,int> i1, std::pair<Interval3,int> i2) { 
        return !(less_than(i1.first[0].first, i2.first[0].first));
        };
    Numccd low_number; low_number.first=0; low_number.second=0;// low_number=0;
    Numccd up_number; up_number.first=1; up_number.second=0;// up_number=1;
    // initial interval [0,1]
    Singleinterval init_interval;init_interval.first=low_number;init_interval.second=up_number;
    //build interval set [0,1]x[0,1]x[0,1]
    Interval3 iset;
    iset[0]=init_interval;iset[1]=init_interval;iset[2]=init_interval;
    // Stack of intervals and the last split dimension
    // std::stack<std::pair<Interval3,int>> istack;
    std::priority_queue<std::pair<Interval3,int>, std::vector<std::pair<Interval3,int>>, decltype(cmp)> istack(cmp);
    istack.emplace(iset,-1);

    // current intervals
    Interval3 current;
    std::array<double,3> err_and_ms;
    err_and_ms[0]=err[0]+ms;
    err_and_ms[1]=err[1]+ms;
    err_and_ms[2]=err[2]+ms;
    refine=0;
    double impact_ratio=0.8;
    toi=std::numeric_limits<double>::infinity();
    Numccd TOI; TOI.first=1;TOI.second=0;
    //std::array<double,3> 
    bool collision=false;
    int rnbr=0;
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();
        // if(rnbr>0&&less_than( current[0].first,TOI)){
        //     std::cout<<"not the first"<<std::endl;
        //     // continue;
        // }
        if(!less_than(current[0].first,TOI)){
            // std::cout<<"not the first"<<std::endl;
            continue;
        }
        //TOI should always be no larger than current
        
            
        // if(Numccd2double(current[0].first)>=Numccd2double(TOI)){
        //     std::cout<<"here wrong, comparing"<<std::endl;
        // } 
#ifdef USE_TIMER
        igl::Timer timer;

        timer.start();
#endif
        refine++;
        bool zero_in;
       bool box_in;
        zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,box_in);
        
#ifdef USE_TIMER
    
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL// this is defined in the begining of this file
        
        zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);

#endif
        if(!zero_in) continue;
#ifdef USE_TIMER
        timer.start();
#endif
        Eigen::VectorX3d widths = width(current);
#ifdef USE_TIMER
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        if ((widths.array() <= tol.array()).all()) {
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI)*impact_ratio;
            return true;
        }
        if(box_in){
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI)*impact_ratio;
            return true;
        }

        std::array<bool , 3> check;
        Eigen::VectorX3d widthratio;
        widthratio.resize(3);
        check[0]=false;check[1]=false; check[2]=false;
        for(int i=0;i<3;i++){
            widthratio(i)=widths(i)/tol(i);
            if(widths(i) > tol(i))
                check[i]=true;// means this need to be checked
        }
        
        int split_i=-1;
        for(int i=0;i<3;i++){
            if(check[i]){
                if(check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)&&widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }    
                }
                if(check[(i+1)%3]&&!check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&!check[(i+2)%3]){
                   
                        split_i=i;
                        break;
                    
                }
            }
        }
        if(split_i<0){
            std::cout<<"ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT"<<std::endl;
        }
        // Bisect the next dimension that is greater than its tolerance
        // int split_i;
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
        if(!less_than(halves.first.first, halves.first.second)){
                std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
                return true;
            }
        if(!less_than(halves.second.first, halves.second.second)){
            std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
            return true;
        }
        if(check_vf){
            //std::cout<<"*** check_vf"<<std::endl;
            if(split_i==1){
               // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
               // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));
                
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){
                //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));
                
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
            
        }
        else{
            current[split_i] = halves.second;
            istack.emplace(current, split_i);
            current[split_i] = halves.first;
            istack.emplace(current, split_i);
        }

    }
    // if(collision) toi=Numccd2double(TOI);
// std::cout<<"checking tail"<<std::endl;
    Numccd delta_t;
    t_tol_width(delta_t,2*tol(0)); // make delta larger to be conservative
    
    iset[0].first=up_number;//t0=1;
    Numccd tail;
    sum_up(up_number,delta_t,tail);// t1=1+delta
    iset[0].second=tail;
    TOI=tail;
    // now we get the interval. next 
    /////////////////////////////////////////////////////////////////////////////////////
    if(!istack.empty()) std::cout<<"ERROR HERE, STACK SHOULD BE EMPTY"<<std::endl;
    istack.emplace(iset,-1);

    // current intervals
    
    // toi=std::numeric_limits<double>::infinity();
    // Numccd TOI; TOI.first=1;TOI.second=0;
    //std::array<double,3> 
    // bool collision=false;
    // int rnbr=0;
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();
 

#ifdef USE_TIMER
        igl::Timer timer;

        timer.start();
#endif
        refine++;
        bool zero_in;
       bool box_in;
        zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,box_in);
        
#ifdef USE_TIMER
    
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL// this is defined in the begining of this file
        
        zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);

#endif
        if(!zero_in) continue;
#ifdef USE_TIMER
        timer.start();
#endif
        Eigen::VectorX3d widths = width(current);
#ifdef USE_TIMER
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        if ((widths.array() <= tol.array()).all()) {
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI)*impact_ratio;
            return true;
        }
        if(box_in){
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI)*impact_ratio;
            return true;
        }

        std::array<bool , 3> check;
        Eigen::VectorX3d widthratio;
        widthratio.resize(3);
        check[0]=false;check[1]=false; check[2]=false;
        for(int i=0;i<3;i++){
            widthratio(i)=widths(i)/tol(i);
            if(widths(i) > tol(i))
                check[i]=true;// means this need to be checked
        }
        
        int split_i=-1;
        for(int i=0;i<3;i++){
            if(check[i]){
                if(check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)&&widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }    
                }
                if(check[(i+1)%3]&&!check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&!check[(i+2)%3]){
                   
                        split_i=i;
                        break;
                    
                }
            }
        }
        if(split_i<0){
            std::cout<<"ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT"<<std::endl;
        }
        // Bisect the next dimension that is greater than its tolerance
        // int split_i;
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
        if(!less_than(halves.first.first, halves.first.second)){
                std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
                return true;
            }
        if(!less_than(halves.second.first, halves.second.second)){
            std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
            return true;
        }
        if(check_vf){
            //std::cout<<"*** check_vf"<<std::endl;
            if(split_i==1){
               // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
               // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));
                
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){
                //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));
                
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
            
        }
        else{
            current[split_i] = halves.second;
            istack.emplace(current, split_i);
            current[split_i] = halves.first;
            istack.emplace(current, split_i);
        }

    }
    return collision;
    return false;
    
}

bool distance_less_than(const Numccd& n1,const Numccd& n2, const double nbr){
    if(n1==n2) return false;//if distance ==0, we can continue to check
    double r=fabs(Numccd2double(n1)-Numccd2double(n2));
    return r<nbr;
}
bool interval_root_finder_double_two_stacks(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    // Interval3& final,
    double& toi,
    const bool check_vf,
    const std::array<double,3> err,
    const double ms,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e){
        auto cmp = [](std::pair<Interval3,int> i1, std::pair<Interval3,int> i2) { 
        return !(less_than(i1.first[0].first, i2.first[0].first));
        };
        auto cmp3 = [](std::pair<Interval3,int> i1, std::pair<Interval3,int> i2) { 
         if(i1.first[0].first!=i2.first[0].first) {
             return !(less_than(i1.first[0].first, i2.first[0].first));
         }  
         else{
             if(i1.first[1].first!=i2.first[1].first){
                 return !(less_than(i1.first[1].first, i2.first[1].first));
             }
             else{
                 return !(less_than(i1.first[2].first, i2.first[2].first));
             }
             
         }
        
        };
    Numccd low_number; low_number.first=0; low_number.second=0;// low_number=0;
    Numccd up_number; up_number.first=1; up_number.second=0;// up_number=1;
    // initial interval [0,1]
    Singleinterval init_interval;init_interval.first=low_number;init_interval.second=up_number;
    //build interval set [0,1]x[0,1]x[0,1]
    Interval3 iset;
    iset[0]=init_interval;iset[1]=init_interval;iset[2]=init_interval;
    // Stack of intervals and the last split dimension
    // std::stack<std::pair<Interval3,int>> istack, istack2; 
    std::priority_queue<std::pair<Interval3,int>, 
    std::vector<std::pair<Interval3,int>>, decltype(cmp3)> istack(cmp3), istack2(cmp3);
    istack.emplace(iset,-1);
    // std::cout<<"tol,"<<tol(0)<<","<<tol(1)<<","<<tol(2)<<std::endl;
    // current intervals
    Interval3 current;
    std::array<double,3> err_and_ms;
    err_and_ms[0]=err[0]+ms;
    err_and_ms[1]=err[1]+ms;
    err_and_ms[2]=err[2]+ms;
    refine=0;
    toi=std::numeric_limits<double>::infinity();
    Numccd TOI; TOI.first=1;TOI.second=0;
    //std::array<double,3> 
    bool collision=false;
    int rnbr=0;
    bool check_stack_2=false;
    double root_dis=tol(0)*10000;
    
    int last_split;
    #ifdef DEBUGING
    std::cout<<"root_dist, "<<root_dis<<std::endl;
    std::ofstream fout;
    
    fout.open("uvt.csv");
    #endif
    while(!(istack.empty()&&istack2.empty())){
        if(istack.empty()){
            check_stack_2=true;
        }
        if(istack2.empty()){
            check_stack_2=false;
        }
        if(check_stack_2){
        current=istack2.top().first;
        last_split=istack2.top().second;
        istack2.pop();
        }
        else{
            current=istack.top().first;
        last_split=istack.top().second;
        istack.pop();
        }
        

        if(!less_than(current[0].first,TOI)){
            continue;
        }
        // if(Numccd2double(current[0].first)>=Numccd2double(TOI)){
        //     std::cout<<"here wrong, comparing"<<std::endl;
        // } 
#ifdef USE_TIMER
        igl::Timer timer;

        timer.start();
#endif
        
        bool zero_in;
       bool bbox_in;
        zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,bbox_in);
        
#ifdef USE_TIMER
    
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL// this is defined in the begining of this file
        
        zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);

#endif
        if(!zero_in) continue;
#ifdef USE_TIMER
        timer.start();
#endif
        Eigen::VectorX3d widths = width(current);
        #ifdef DEBUGING
        fout<<Numccd2double(current[0].first)<<","<< Numccd2double(current[1].first)<<","<<
        Numccd2double(current[2].first)<<","<<widths(0)<<","<<widths(1)<<","<<widths(2)<<"\n";
        std::cout<<"t,"<<Numccd2double(current[0].first)<< "\r" << std::flush;
        #endif
#ifdef USE_TIMER
        
        // std::cout<<"width,"<<widths(0)<<","<<widths(1)<<","<<widths(1)<< "\r" << std::flush;
        // std::cout<<"tuv,"<<widths(0)<<","<<widths(1)<<","<<widths(1)<< "\r" << std::flush;
        
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        if ((widths.array() <= tol.array()).all()) {
            if(!less_than(current[0].first,TOI)){
                std::cout<<"find bigger?"<<std::endl;
            }
            TOI=current[0].first;
            collision=true;
            rnbr++;
            
            // continue;
            toi=Numccd2double(TOI);
            return true;
        }
        if(bbox_in){
            if(!less_than(current[0].first,TOI)){
                std::cout<<"find bigger?"<<std::endl;
            }
            TOI=current[0].first;
            collision=true;
            rnbr++;
            // continue;
            toi=Numccd2double(TOI);
            return true;
        }
        refine++;
        #ifdef DEBUGING
        if(refine>100000){
            fout.close();
        }
        #endif
        // std::cout << "roots"<<rnbr << "\r" << std::flush;
        // std::cout << "refine"<<refine << "\r" << std::flush;
        std::array<bool , 3> check;
        Eigen::VectorX3d widthratio;
        widthratio.resize(3);
        check[0]=false;check[1]=false; check[2]=false;
        for(int i=0;i<3;i++){
            widthratio(i)=widths(i)/tol(i);
            if(widths(i) > tol(i))
                check[i]=true;// means this need to be checked
        }
        
        int split_i=-1;
        for(int i=0;i<3;i++){
            if(check[i]){
                if(check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)&&widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }    
                }
                if(check[(i+1)%3]&&!check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+1)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&check[(i+2)%3]){
                    if(widthratio(i)>=widthratio((i+2)%3)){
                        split_i=i;
                        break;
                    }
                }
                if(!check[(i+1)%3]&&!check[(i+2)%3]){
                   
                        split_i=i;
                        break;
                    
                }
            }
        }
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        if(split_i<0){
            std::cout<<"ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT"<<std::endl;
        }
        // Bisect the next dimension that is greater than its tolerance
        // int split_i;
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
        if(!less_than(halves.first.first, halves.first.second)){
                std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
                return true;
            }
        if(!less_than(halves.second.first, halves.second.second)){
            std::cout<<"OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"<<std::endl;
            return true;
        }
        if(check_vf){
            //std::cout<<"*** check_vf"<<std::endl;
            Interval3 Ipush1=current;
            Interval3 Ipush2=current;
            bool push1=false;
            bool push2=false;
            if(split_i==1){
               // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
               // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));
                
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    Ipush1=current;
                    push1=true;
                    //istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    Ipush2=current;
                    push2=true;
                    //istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){
                //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));
                
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    Ipush1=current;
                    push1=true;
                    //istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    Ipush2=current;
                    push2=true;
                    //istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                Ipush1=current;
                push1=true;
                // istack.emplace(current, split_i);
                current[split_i] = halves.first;
                Ipush2=current;
                push2=true;
                // istack.emplace(current, split_i);
            }
            if(check_stack_2){
                
                if(push1){
                    if (distance_less_than(Ipush1[0].first,TOI,root_dis)){
                        istack.emplace(Ipush1,split_i);
                    }
                    else{
                        istack2.emplace(Ipush1,split_i);
                    }
                }
                if(push2){
                    if (distance_less_than(Ipush2[0].first,TOI,root_dis)){
                        istack.emplace(Ipush2,split_i);
                    }
                    else{
                        istack2.emplace(Ipush2,split_i);
                    }
                }
            }
            if(!check_stack_2){
                
                if(push1){
                    if (distance_less_than(Ipush1[0].first,TOI,root_dis)){
                        istack2.emplace(Ipush1,split_i);
                    }
                    else{
                        istack.emplace(Ipush1,split_i);
                    }
                }
                if(push2){
                    if (distance_less_than(Ipush2[0].first,TOI,root_dis)){
                        istack2.emplace(Ipush2,split_i);
                    }
                    else{
                        istack.emplace(Ipush2,split_i);
                    }
                }
            }
            
        }
        else{
            
            if(check_stack_2){
                current[split_i] = halves.second;
                if(less_than(current[0].first,TOI)){
                    if (distance_less_than(current[0].first,TOI,root_dis)){
                        istack.emplace(current, split_i);
                    }
                    else{
                         istack2.emplace(current, split_i);
                    }
                }
                
                current[split_i] = halves.first;
                if(less_than(current[0].first,TOI)){
                    if(distance_less_than(current[0].first,TOI,root_dis)){
                        istack.emplace(current, split_i);
                    }
                    else{
                         istack2.emplace(current, split_i);
                    }
                }
                
            }
            if(!check_stack_2){
                current[split_i] = halves.second;
                if(less_than(current[0].first,TOI)){
                    if (distance_less_than(current[0].first,TOI,root_dis)){
                        istack2.emplace(current, split_i);
                    }
                    else{
                         istack.emplace(current, split_i);
                    }
                }
                
                current[split_i] = halves.first;
                if(less_than(current[0].first,TOI)){
                    if(distance_less_than(current[0].first,TOI,root_dis)){
                        istack2.emplace(current, split_i);
                    }
                    else{
                         istack.emplace(current, split_i);
                    }
                }
                
            }
            
            
        }

    }
    if(collision) toi=Numccd2double(TOI);
    // if(toi==0)std::cout<<"infinate roots, "<<std::endl;
    // if(rnbr>5) std::cout<<"nbr of roots, "<<rnbr<<", time, "<<toi<<std::endl;
    return collision;
    return false;
    
}
int print_refine(){
    return refine;
}


bool interval_root_finder_Rational(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    std::array<std::pair<Rational,Rational>, 3>& final,
    const bool check_vf,
    const std::array<double,3> err,
    const double ms,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e){
    
    std::pair<Rational,Rational> interval01;interval01.first=0;interval01.second=1;
    std::array<std::pair<Rational,Rational>, 3> paracube;
    paracube[0]=interval01;paracube[1]=interval01;paracube[2]=interval01;
    
    
    // Stack of intervals and the last split dimension
    std::stack<std::pair<std::array<std::pair<Rational,Rational>, 3>,int>> istack;
    istack.emplace(paracube,-1);

    // current intervals
    std::array<std::pair<Rational,Rational>, 3> current;
    std::array<double,3> err_and_ms;
    err_and_ms[0]=err[0]+ms;
    err_and_ms[1]=err[1]+ms;
    err_and_ms[2]=err[2]+ms;
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();
        igl::Timer timer;

        timer.start();
        bool zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);
        timer.stop();
        time_rational+=timer.getElapsedTimeInMicroSec();

        if(!zero_in) continue;
        timer.start();
        std::array<Rational,3> widths = width(current);
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
        if (widths[0].to_double()<=tol(0)&&widths[1].to_double()<=tol(1)&&widths[2].to_double()<=tol(2)) {
            final=current;
                return true;
        }

        // Bisect the next dimension that is greater than its tolerance
        int split_i;
        for (int i = 1; i <= 3; i++) {
            split_i = (last_split + i) % 3;
            if (widths[split_i].to_double() > tol(split_i)) {
                break;
            }
        }
        std::pair<std::pair<Rational,Rational>, std::pair<Rational,Rational>> halves = bisect(current[split_i]);

        if(check_vf){
            if(split_i==1){
                if(halves.first.first+current[2].first<=1){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                if(halves.second.first+current[2].first<=1){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){

                if(halves.first.first+current[1].first<=1){
                    current[split_i]=halves.first;
                    istack.emplace(current, split_i);
                }
                if(halves.second.first+current[1].first<=1){
                    current[split_i]=halves.second;
                    istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
        }
        else{
            current[split_i] = halves.second;
            istack.emplace(current, split_i);
            current[split_i] = halves.first;
            istack.emplace(current, split_i);
        }

    }
    return false;
    
}


void print_time_2(){
    std::cout<<"origin predicates, "<<time20<<std::endl;
    std::cout<<"width, "<<time21<<std::endl;
    std::cout<<"bisect, "<<time22<<std::endl;
    std::cout<<"origin part1(evaluate 1 dimension), "<<time23<<std::endl;
    std::cout<<"origin part2(convert tuv), "<<time24<<std::endl;
    std::cout<<"time of call the vertex solving function, "<<time25<<std::endl;
    std::cout<<"how many times of interval check for this query, "<<refine<<std::endl;
}
double print_time_rational(){
    return time_rational;
}

std::array<double, 3> get_numerical_error(const std::vector<Eigen::Vector3d> &vertices,const bool& check_vf){
    double xmax=fabs(vertices[0][0]);
    double ymax=fabs(vertices[0][1]);
    double zmax=fabs(vertices[0][2]);
    for(int i=0;i<vertices.size();i++){
        if(xmax<fabs(vertices[i][0])){
            xmax=fabs(vertices[i][0]);
        }
        if(ymax<fabs(vertices[i][1])){
            ymax=fabs(vertices[i][1]);
        }
        if(zmax<fabs(vertices[i][2])){
            zmax=fabs(vertices[i][2]);
        }
    }
    double delta_x=xmax>1?xmax:1;
    double delta_y=ymax>1?ymax:1;
    double delta_z=zmax>1?zmax:1;
    std::array<double, 3> result;
    if(!check_vf){
        result[0]=delta_x*delta_x*delta_x*6.217248937900877e-15;
        result[1]=delta_y*delta_y*delta_y*6.217248937900877e-15;
        result[2]=delta_z*delta_z*delta_z*6.217248937900877e-15;
    }
    else{
        result[0]=delta_x*delta_x*delta_x*6.661338147750939e-15;
        result[1]=delta_y*delta_y*delta_y*6.661338147750939e-15;
        result[2]=delta_z*delta_z*delta_z*6.661338147750939e-15;
    }
    return result;
}
} // namespace ccd
