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
} // namespace ccd
