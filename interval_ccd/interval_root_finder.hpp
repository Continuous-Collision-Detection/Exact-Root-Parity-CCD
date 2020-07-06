// A root finder using interval arithmetic.
#pragma once

#include <functional>

#include "interval.hpp"
#include<vector>

namespace intervalccd {

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x);

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x,const bool check_vf);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x,const bool check_vf);
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x);
bool interval_root_finder_opt(
    const std::function<Eigen::VectorX3I(const Numccd&, const Numccd&, const Numccd&)>& f,
    const Eigen::VectorX3d& tol,
    Interval3& final,
    const bool check_vf);

bool interval_root_finder_double(
    const Eigen::VectorX3d& tol,
    //Eigen::VectorX3I& x,// result interval
    Interval3& final,
    const bool check_vf,
    const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e);
// return power t. n=result*2^t
int reduction(const int n, int& result);
std::pair<Singleinterval, Singleinterval> bisect(const Singleinterval& inter);
// calculate the sign of f. dim is the dimension we are evaluating.
template<typename T>
int function_f (
const Numccd&tpara, const Numccd&upara, const Numccd&vpara,const T& type, const int dim,
const Eigen::Vector3d& a0s,
    const Eigen::Vector3d& a1s,
    const Eigen::Vector3d& b0s,
    const Eigen::Vector3d& b1s,
    const Eigen::Vector3d& a0e,
    const Eigen::Vector3d& a1e,
    const Eigen::Vector3d& b0e,
    const Eigen::Vector3d& b1e );
void print_time_2();
double print_time_rational();
std::array<double, 3> get_numerical_error(const std::vector<Eigen::Vector3d> &vertices,const bool& check_vf);
} // namespace ccd
