// A root finder using interval arithmetic.
#include <interval_ccd/interval_root_finder.hpp>

#include <stack>
#include<igl/Timer.h>
#include<iostream>
#include<interval_ccd/Rational.hpp>
#include <interval_ccd/avx.h>
// #define COMPARE_WITH_RATIONAL
// #define USE_TIMER

namespace intervalccd {
double time20=0,time21=0,time22=0, time23=0,time24=0,time25=0,time_rational=0,time_convert=0;;
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
    const int dimension,const bool check_vf,const double eps,std::array<double,8> &vs){
#ifdef USE_TIMER
    igl::Timer timer;
#endif
    
    int count=0;
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
    return true;

}
bool evaluate_bbox_one_dimension_vector_half(
    std::array<double,4>& t_up,std::array<double,4>&t_dw,
    std::array<double,4>& u_up,std::array<double,4>&u_dw,
    std::array<double,4>& v_up,std::array<double,4>&v_dw,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e,
    const int dimension,const bool check_vf,const double eps,
    const std::array<double,8> &prism,const int split_i){
#ifdef USE_TIMER
    igl::Timer timer;
#endif
    
    int count=0;
    std::array<double,4> vs;
#ifdef USE_TIMER
    timer.start();
#endif  
    if(check_vf){
        vs=function_vf(a0s[dimension],a1s[dimension],b0s[dimension],b1s[dimension],a0e[dimension],a1e[dimension],b0e[dimension],b1e[dimension],
        t_up,t_dw,u_up,u_dw,v_up,v_dw);
    }
    else{
        vs=function_ee(a0s[dimension],a1s[dimension],b0s[dimension],b1s[dimension],a0e[dimension],a1e[dimension],b0e[dimension],b1e[dimension],
        t_up,t_dw,u_up,u_dw,v_up,v_dw);
    }
#ifdef USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif
    double minv=vs[0], maxv=vs[0];
    
    for(int i=1;i<4;i++){
        if(minv>vs[i]){
            minv=vs[i];
        }
        if(maxv<vs[i]){
            maxv=vs[i];
        }
    }
    
    
    std::array<int,4>ids=select_useful_pt_id(split_i);
    for(int i=0;i<4;i++){
        if(minv>prism[ids[i]]){
            minv=prism[ids[i]];
        }
        if(maxv<prism[ids[i]]){
            maxv=prism[ids[i]];
        }
    }
    if(minv>eps||maxv<-eps)
        return false;
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
       
       int tu = tpara.first, td=tpara.second;// t=tu/(2^td)
       int uu = upara.first, ud=upara.second;
       int vu = vpara.first, vd=vpara.second;
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
           = (a0e - a0s) * tu/int(power(1,td))
           + a0s;
       Vector3r edge0_vertex1
           = (a1e - a1s) * tu/int(power(1,td))
           + a1s;
       Vector3r edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/int(power(1,ud)) + edge0_vertex0;

       Vector3r edge1_vertex0
           = (b0e - b0s) * tu/int(power(1,td))
           + b0s;
       Vector3r edge1_vertex1
           = (b1e - b1s) * tu/int(power(1,td))
           + b1s;
       Vector3r edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/int(power(1,vd)) + edge1_vertex0;

       
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
       
       int tu = tpara.first, td=tpara.second;// t=tu/(2^td)
       int uu = upara.first, ud=upara.second;
       int vu = vpara.first, vd=vpara.second;
        Vector3r 
        vs(a0sd[0],a0sd[1],a0sd[2]), 
        t0s(a1sd[0],a1sd[1],a1sd[2]),
        t1s(b0sd[0],b0sd[1],b0sd[2]),
        t2s(b1sd[0],b1sd[1],b1sd[2]),

        ve(a0ed[0],a0ed[1],a0ed[2]),
        t0e(a1ed[0],a1ed[1],a1ed[2]),
        t1e(b0ed[0],b0ed[1],b0ed[2]),
        t2e(b1ed[0],b1ed[1],b1ed[2]);

        Vector3r v=(ve-vs)*tu/int(power(1,td))+vs;

        Vector3r t0=(t0e-t0s)*tu/int(power(1,td))+t0s;
        Vector3r t1=(t1e-t1s)*tu/int(power(1,td))+t1s;
        Vector3r t2=(t2e-t2s)*tu/int(power(1,td))+t2s;
        Vector3r p=(t1-t0)*uu/int(power(1,ud))+(t2-t0)*vu/int(power(1,vd))+t0;
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
                time24+=timer.getElapsedTimeInMicroSec();
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
    const std::array<double,3>& box,std::array<std::array<double,8>,3> &prism,int&checked){
#ifdef USE_TIMER
    igl::Timer timer;
#endif  
    std::array<double,8> t_up;std::array<double,8>t_dw;
    std::array<double,8> u_up;std::array<double,8>u_dw;
    std::array<double,8> v_up;std::array<double,8>v_dw;
Interval3_2_double itv(paras);
#ifdef USE_TIMER
    timer.start();
#endif

    convert_tuv_to_array(itv,t_up,t_dw,u_up,u_dw,v_up,v_dw);
#ifdef USE_TIMER
    timer.stop();
    time_convert+=timer.getElapsedTimeInMicroSec();
#endif
    bool ck;
    checked=0;
    for(int i=0;i<3;i++){
#ifdef USE_TIMER
    timer.start();
#endif
        ck=evaluate_bbox_one_dimension_vector(t_up,t_dw,u_up,u_dw,v_up,v_dw,
    a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,i,check_vf,box[i],prism[i]);
    checked++;
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
     if(!ck)
        return false;
    }
    return true;
    
    
}
bool Origin_in_function_bounding_box_double_vector_half(//TODO write sth to deal with checked
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
    const std::array<double,3>& box,const std::array<std::array<double,8>,3> &prism,const int split_i){

Interval3_2_double itv(paras);
#ifdef USE_TIMER
    igl::Timer timer;
#endif  
    std::array<double,4> t_up;std::array<double,4>t_dw;
    std::array<double,4> u_up;std::array<double,4>u_dw;
    std::array<double,4> v_up;std::array<double,4>v_dw;
    
    convert_tuv_to_array_half(itv,t_up,t_dw,u_up,u_dw,v_up,v_dw,split_i);
    bool ck;
    for(int i=0;i<3;i++){
#ifdef USE_TIMER
    timer.start();
#endif
        ck=evaluate_bbox_one_dimension_vector_half(t_up,t_dw,u_up,u_dw,v_up,v_dw,
    a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,i,check_vf,box[i],prism[i],split_i);
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
     if(!ck)
        return false;
    }
    return true;
    
    
}
//split_i is the para splited, checked is the dimension already calculated
bool Origin_in_function_bounding_box_double_vector_half(//TODO write sth to deal with checked
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
    const std::array<double,3>& box,
    const std::array<std::array<double,8>,3> &prism,const int split_i,const int checked){
#ifdef USE_TIMER
    igl::Timer timer;
#endif  
    std::array<double,4> t_up;std::array<double,4>t_dw;
    std::array<double,4> u_up;std::array<double,4>u_dw;
    std::array<double,4> v_up;std::array<double,4>v_dw;

    std::array<double,8> t_up_whole;std::array<double,8>t_dw_whole;
    std::array<double,8> u_up_whole;std::array<double,8>u_dw_whole;
    std::array<double,8> v_up_whole;std::array<double,8>v_dw_whole;

Interval3_2_double itv(paras);
#ifdef USE_TIMER
    timer.start();
#endif
    
    convert_tuv_to_array_half(itv,t_up,t_dw,u_up,u_dw,v_up,v_dw,split_i);
    std::array<double,8> temv;
    if(checked<3){
        convert_tuv_to_array(itv,t_up_whole,t_dw_whole,u_up_whole,u_dw_whole,v_up_whole,v_dw_whole);
    }
#ifdef USE_TIMER
    timer.stop();
    time_convert+=timer.getElapsedTimeInMicroSec();
#endif    
    bool ck;
    for(int i=0;i<3;i++){
#ifdef USE_TIMER
    timer.start();
#endif

    if(i<checked){
        ck=evaluate_bbox_one_dimension_vector_half(t_up,t_dw,u_up,u_dw,v_up,v_dw,
    a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,i,check_vf,box[i],prism[i],split_i);
    }
        
    else{
        ck=evaluate_bbox_one_dimension_vector(t_up_whole,t_dw_whole,u_up_whole,u_dw_whole,v_up_whole,v_dw_whole,
    a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,i,check_vf,box[i],temv);
    }
#ifdef USE_TIMER
    timer.stop();
    time23+=timer.getElapsedTimeInMicroSec();
#endif
     if(!ck)
        return false;
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
    //std::cout<<"t, ["<<t[0].first/power(1,t[0].second)<<","<<t[1].first/power(1,t[1].second)<<"]"<<std::endl;
    //std::cout<<"u, ["<<u[0].first/power(1,u[0].second)<<","<<u[1].first/power(1,u[1].second)<<"]"<<std::endl;
    //std::cout<<"v, ["<<v[0].first/power(1,v[0].second)<<","<<v[1].first/power(1,v[1].second)<<"]"<<std::endl<<std::endl;
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
int reduction(const int n, int& result){
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

    // interval is [k1/power(1,n1), k2/power(1,n2)], k1,k2,n1,n2 are all not negative
    int k1=low.first;
    int n1=low.second;
    int k2=up.first;
    int n2=up.second;

    assert(k1 >= 0 && k2 >= 0 && n1 >= 0 && n2 >= 0);

    std::pair<Singleinterval, Singleinterval> result;
    int k,n, p,r;
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
    int k1=num1.first;
    int n1=num1.second;
    int k2=num2.first;
    int n2=num2.second;
    int k,n;
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
bool sum_no_larger_1_Rational(const Numccd& num1, const Numccd& num2){
    int k1=num1.first;
    int n1=num1.second;
    int k2=num2.first;
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
       
       int tu = tpara.first, td=tpara.second;// t=tu/(2^td)
       int uu = upara.first, ud=upara.second;
       int vu = vpara.first, vd=vpara.second;

       T edge0_vertex0
           = (T(a0e[dim]) - T(a0s[dim])) * tu/int(power(1,td))
           + T(a0s[dim]);
       T edge0_vertex1
           = (T(a1e[dim]) - T(a1s[dim])) * tu/int(power(1,td))
           + T(a1s[dim]);
       T edge0_vertex
           = (edge0_vertex1 - edge0_vertex0) * uu/int(power(1,ud)) + edge0_vertex0;

       T edge1_vertex0
           = (T(b0e[dim]) - T(b0s[dim])) * tu/int(power(1,td))
           + T(b0s[dim]);
       T edge1_vertex1
           = (T(b1e[dim]) - T(b1s[dim])) * tu/int(power(1,td))
           + T(b1s[dim]);
       T edge1_vertex
           = (edge1_vertex1 - edge1_vertex0) * vu/int(power(1,vd)) + edge1_vertex0;

       
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
       
       int tu = tpara.first, td=tpara.second;// t=tu/(2^td)
       int uu = upara.first, ud=upara.second;
       int vu = vpara.first, vd=vpara.second;

        T v= (T(ve[dim])-T(vs[dim]))* tu/int(power(1,td))+T(vs[dim]);
        T t0=(T(t0e[dim])-T(t0s[dim]))* tu/int(power(1,td))+T(t0s[dim]);
        T t1=(T(t1e[dim])-T(t1s[dim]))* tu/int(power(1,td))+T(t1s[dim]);
        T t2=(T(t2e[dim])-T(t2s[dim]))* tu/int(power(1,td))+T(t2s[dim]);
        T p=(t1-t0)*uu/int(power(1,ud))+(t2-t0)*vu/int(power(1,vd))+t0;
        return v-p;
       
}


bool interval_root_finder_double(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    Interval3& final,
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
        {//la = 0.2954570784440031002, lb = 1.1090868531538015859, t = 4635.2393745312786939
            //Vector3r value=function_f_ee_Rational(Rational(4635.2393745312786939),Rational( 0.2954570784440031002),
            //Rational(1.1090868531538015859),a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
            // std::cout<<"a0s,"<<a0s.transpose()<<std::endl;
            // std::cout<<"a1s,"<<a1s.transpose()<<std::endl;
            // std::cout<<"b0s,"<<b0s.transpose()<<std::endl;
            // std::cout<<"b1s,"<<b1s.transpose()<<std::endl;
            // std::cout<<"a0e,"<<a0e.transpose()<<std::endl;
            // std::cout<<"a1e,"<<a1e.transpose()<<std::endl;
            // std::cout<<"b0e,"<<b0e.transpose()<<std::endl;
            // std::cout<<"b1e,"<<b1e.transpose()<<std::endl;
            // std::cout<<"rational value, "<<value[0]<<", "<<value[1]<<", "<<value[2]<<std::endl;
        }
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
    Interval3 current, interval1,interval2;
    std::array<double,3> err_and_ms;
    err_and_ms[0]=err[0]+ms;
    err_and_ms[1]=err[1]+ms;
    err_and_ms[2]=err[2]+ms;
    refine=0;
    bool first_check=true;
    std::array<std::array<double,8>,3> prism;
    int checked;
    //std::array<double,3> 
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
        Eigen::VectorX3d widths = width(current);
        if(first_check){
            zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,checked);
            if(!zero_in){
                return false;
            }
            if ((widths.array() <= tol.array()).all()) {
                final=current;
                return true;
            }
            
        }
        first_check=false;
        
#ifdef USE_TIMER
        
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL// this is defined in the begining of this file
        
        zero_in=Origin_in_function_bounding_box_Rational(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf);

#endif
        
#ifdef USE_TIMER
        timer.start();
#endif
        
#ifdef USE_TIMER
        timer.stop();
        time21+=timer.getElapsedTimeInMicroSec();
#endif
        

        std::array<bool , 3> check;
        Eigen::VectorX3d widthratio;
        widthratio.resize(3);
        check[0]=false;check[1]=false; check[2]=false;
        for(int i=0;i<3;i++){
            widthratio(i)=widths(i)/tol(i);
            if(widths(i) > tol(i))
                check[i]=true;// means this need to be checked
        }
        // if(!check[0]&&!check[1]&&!check[2]){
        //     std::cout<<"all check false"<<std::endl;
        // }
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

        // if(split_i<0){
        //     std::cout<<"ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT"<<std::endl;
        //     std::cout<<"check,"<<check[0]<<","<<check[1]<<","<<check[2]<<std::endl;
        //     std::cout<<"ratio,"<<widthratio[0]<<","<<widthratio[1]<<","<<widthratio[2]<<std::endl;
        //     std::cout<<"widths,"<<widths[0]<<","<<widths[1]<<","<<widths[2]<<std::endl;
            
        // }
        assert(split_i>=0);
        //return false;
        // Bisect the next dimension that is greater than its tolerance
        // int split_i;
        // for (int i = 1; i <= 3; i++) {
        //     split_i = (last_split + i) % 3;
        //     if (widths(split_i) > tol(split_i)) {
        //         break;
        //     }
        // }
        std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);

        if(check_vf){
            bool f2=false, f1=false;
           
            
            if(split_i==1){
               // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
               // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));
                
                if(sum_no_larger_1(halves.second.first, current[2].first)){
                    current[split_i]=halves.second;
                    f2=true;
                    interval2=current;
                    
                    //istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[2].first)){
                    current[split_i]=halves.first;
                    f1=true;
                    interval1=current;
                    
                    //istack.emplace(current, split_i);
                }
                
            }

            if(split_i==2){
                //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));
                
                if(sum_no_larger_1(halves.second.first, current[1].first)){
                    current[split_i]=halves.second;
                    f2=true;
                    interval2=current;
                    
                    //istack.emplace(current, split_i);
                }
                if(sum_no_larger_1(halves.first.first, current[1].first)){
                    current[split_i]=halves.first;
                    f1=true;
                    interval1=current;
                    
                    //istack.emplace(current, split_i);
                }

            }
            if(split_i==0){
                f1=true;f2=true;
                current[split_i] = halves.second;
                interval2=current;
                //istack.emplace(current, split_i);
                current[split_i] = halves.first;
                interval1=current;
                
                //istack.emplace(current, split_i);
            }
            checked=0;
            widths(split_i) /= 2;
            if(f2){
#ifdef USE_TIMER      
        timer.start();
        
#endif
                zero_in= Origin_in_function_bounding_box_double_vector(interval2,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,checked);
#ifdef USE_TIMER
        
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif                
                if(zero_in){
                    

                    if ((widths.array() <= tol.array()).all()) {
                        final=interval2;
                            return true;
                    }
                    
                    istack.emplace(interval2,split_i);
                }
                
            }
            if(f1){//some vertices already calculated
#ifdef USE_TIMER
        
        timer.start();
#endif
                zero_in= Origin_in_function_bounding_box_double_vector_half(interval1,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,split_i,checked);
#ifdef USE_TIMER
        
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif
                if(zero_in){
                    

                    if ((widths.array() <= tol.array()).all()) {
                        final=interval1;
                            return true;
                    }
                    istack.emplace(interval1,split_i);
                }
                
            }
            // if(!f2&&f1){// vertices NOT calculated
            //     zero_in= Origin_in_function_bounding_box_double_vector(interval1,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,checked);
            //     if(zero_in){
            //         widths = width(interval1);

            //         if ((widths.array() <= tol.array()).all()) {
            //             final=interval1;
            //                 return true;
            //         }
            //         istack.emplace(interval1,split_i);
            //     }
                
            // }
            
        }
        else{
            current[split_i] = halves.second;
            widths(split_i) /=2;
#ifdef USE_TIMER
        
        timer.start();
        
#endif
            zero_in= Origin_in_function_bounding_box_double_vector(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,checked);
#ifdef USE_TIMER
        
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif                
                
                if(zero_in){
                    

                    if ((widths.array() <= tol.array()).all()) {
                        final=current;
                            return true;
                    }
                    
                    istack.emplace(current,split_i);
                }
            current[split_i] = halves.first;
#ifdef USE_TIMER
        timer.start();
#endif
            zero_in= Origin_in_function_bounding_box_double_vector_half(current,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e,check_vf,err_and_ms,prism,split_i,checked);
#ifdef USE_TIMER
        
        timer.stop();
        time20+=timer.getElapsedTimeInMicroSec();
#endif                
                if(zero_in){
                    
                    if ((widths.array() <= tol.array()).all()) {
                        final=current;
                            return true;
                    }
                    
                    istack.emplace(current,split_i);
                }

        }

    }
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
    std::cout<<"evaluate one dimension, "<<time23<<std::endl;
    std::cout<<"convert uvt, "<<time_convert<<std::endl;
    std::cout<<"time of call the vertex solving function, "<<time25<<std::endl;
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
