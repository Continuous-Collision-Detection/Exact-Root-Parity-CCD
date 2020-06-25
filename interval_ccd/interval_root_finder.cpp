// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"

#include <stack>

namespace intervalccd {

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
    double r=double(n.first)/pow(2,n.second);
    return r;
}
Eigen::VectorX3d width(const Interval3&x){
    Eigen::VectorX3d w;
    for(int i=0;i<3;i++){
        w[i]=Numccd2double(x[i].second)-Numccd2double(x[i].first);
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
bool Origin_in_function_bounding_box(
    const Interval3& paras,
    const std::function<Eigen::VectorX3I(const Numccd&, const Numccd&, const Numccd&)>& f){
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
                result=f(t[i],u[j],v[k]);
                if(interval_bounding_box_check(result,flag))
                    return true;
            }
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

    // interval is [k1/pow(2,n1), k2/pow(2,n2)], k1,k2,n1,n2 are all not negative
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
        k=k1*pow(2,n2-n1)+k2; assert(k%2==1);
        n=n2+1;
    }
    if(n2<n1){
        k=k1+k2*pow(2,n1-n2); assert(k%2==1);
        n=n1+1;
    }
    Numccd newnum(k,n);
    Singleinterval i1(low,newnum),i2(newnum,up);
    result.first=i1;result.second=i2;
    return result;
}
bool interval_root_finder_opt(
    const std::function<Eigen::VectorX3I(const Numccd&, const Numccd&, const Numccd&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,// initial interval, must be [0,1]x[0,1]x[0,1]
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
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
    while(!istack.empty()){
        current=istack.top().first;
        int last_split=istack.top().second;
        istack.pop();

        bool zero_in = Origin_in_function_bounding_box(current,f);
        if(!zero_in) continue;

        Eigen::VectorX3d widths = width(current);
        if ((widths.array() <= tol.array()).all()) {
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
        current[split_i] = halves.second;
        istack.emplace(current, split_i);
        current[split_i] = halves.first;
        istack.emplace(current, split_i);

    }
    return false;
    // std::stack<std::pair<Eigen::VectorX3d, int>> xs;
    // xs.emplace(x0, -1);
    // while (!xs.empty()) {
    //     x = xs.top().first;
    //     int last_split = xs.top().second;
    //     xs.pop();

    //     Eigen::VectorX3I y = f(x);

    //     if (!zero_in(y)) {
    //         continue;
    //     }

    //     Eigen::VectorX3d widths = width(x);
    //     if ((widths.array() <= tol.array()).all()) {
    //         if (constraint_predicate(x)) {
    //             return true;
    //         }
    //         continue;
    //     }

    //     // Bisect the next dimension that is greater than its tolerance
    //     int split_i;
    //     for (int i = 1; i <= x.size(); i++) {
    //         split_i = (last_split + i) % x.size();
    //         if (widths(split_i) > tol(split_i)) {
    //             break;
    //         }
    //     }
    //     std::pair<Interval, Interval> halves = bisect(x(split_i));
    //     Eigen::VectorX3I x1 = x;
    //     // Push the second half on first so it is examined after the first half
    //     if(check_vf){
    //         if(split_i==1){
    //             if(interval_satisfy_constrain(halves.second,x(2))){
    //                 x(split_i) = halves.second;
    //                 xs.emplace(x, split_i);
    //             }
    //             if(interval_satisfy_constrain(halves.first,x(2))){
    //                 x(split_i) = halves.first;
    //                 xs.emplace(x, split_i);
    //             }
    //         }
    //         if(split_i==2){
    //             if(interval_satisfy_constrain(halves.second,x(1))){
    //                 x(split_i) = halves.second;
    //                 xs.emplace(x, split_i);
    //             }
    //             if(interval_satisfy_constrain(halves.first,x(1))){
    //                 x(split_i) = halves.first;
    //                 xs.emplace(x, split_i);
    //             }
    //         }
    //         if(split_i==0){
    //             x(split_i) = halves.second;
    //             xs.emplace(x, split_i);
    //             x(split_i) = halves.first;
    //             xs.emplace(x, split_i);
    //         }
    //     }
    //     else{
    //         x(split_i) = halves.second;
    //         xs.emplace(x, split_i);
    //         x(split_i) = halves.first;
    //         xs.emplace(x, split_i);
    //     }
        
        
    // }
    // return false;
    
}
} // namespace ccd
