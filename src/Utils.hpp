#pragma once

#include "Rational.hpp"

#include <Eigen/Core>

namespace ccd {

typedef Eigen::Matrix<Rational, 3, 1> Vector3r;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
static const int COPLANAR = -1;
static const int INTERSECTED = 1;
static const int NOT_INTERSECTED1 = 2;
static const int NOT_INTERSECTED2 = 3;
static const Vector3r ORIGIN = Vector3r(0, 0, 0);

static const int BI_DEGE_PLANE = 1;
static const int BI_DEGE_XOR_02 = 2;
static const int BI_DEGE_XOR_13 = 3;

template <typename V1, typename V2> Vector3r cross(const V1& v1, const V2& v2)
{
    Vector3r res;
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return res;
}
bool XOR(const bool a, const bool b)
{
    if (a && b)
        return false;
    if (!a && !b)
        return false;
    return true;
}
int int_XOR(const int a, const int b)
{
    if (a == -1 || b == -1)
        return -1;
    if (a == b)
        return 0;
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    std::cout << "impossible XOR cases" << std::endl;
}

template <typename V> void print(const V& v)
{
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}

void write(const Vector3d& v, std::ostream& out);
Vector3d read(std::istream& in);

int orient3d(
    const Vector3r& a, const Vector3r& b, const Vector3r& c, const Vector3r& d);
int orient2d(
    const Vector3r& a, const Vector3r& b, const Vector3r& c, const int axis);
int origin_ray_triangle_inter(
    const Vector3d& dir,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3);
bool segment_segment_intersection(
	const Vector3r& s0,
	const Vector3r& e0,
	const Vector3r& s1,
	const Vector3r& e1);
// 0 not intersected, 1 intersected, 2 s0 on segment
// can deal with degenerated cases
int ray_segment_intersection(
	const Vector3r& s0,
	const Vector3r& dir0,
	const Vector3r& s1,
	const Vector3r& e1);
bool segment_segment_inter(
    const Vector3r& s0,
    const Vector3r& e0,
    const Vector3r& s1,
    const Vector3r& e1,
    Vector3r& res,
    int axis);
// this function can also tell us if they are parallel and overlapped
// and also tell us if the parallel case has seg-seg overlapping:


int segment_triangle_inter(
    const Vector3d& e0,
    const Vector3d& e1,
    const Vector3d& t1,
    const Vector3d& t2,
    const Vector3d& t3);
int segment_triangle_inter(
    const Vector3r& e0,
    const Vector3r& e1,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3);
int line_triangle_inter(
    const Vector3r& e0,
    const Vector3r& e1,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3);
int ray_halfopen_triangle_inter(
    const Vector3r& p0,
    const Vector3r& dir,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3);
int ray_open_triangle_inter(
    const Vector3r& p0,
    const Vector3r& dir,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3);
int line_triangle_inter_return_t(
    const Vector3r& e0,
    const Vector3r& e1,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3,
    Rational& t);
bool same_point(const Vector3r& p1, const Vector3r& p2);
Vector3r tri_norm(const Vector3r& t0, const Vector3r& t1, const Vector3r& t2)
{
	Vector3r s1, s2;
	s1 = t1 - t0;
	s2 = t2 - t1;
	return cross(s1, s2);
}
} // namespace ccd
