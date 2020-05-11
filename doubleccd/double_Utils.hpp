#pragma once

#include <doubleCCD/Rational.hpp>
#include <vector>
#include <array>
#include <Eigen/Core>
#include <igl/Timer.h>
namespace doubleccd {

typedef Eigen::Matrix<Rational, 3, 1> Vector3r;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
static const int COPLANAR = -1;
static const int INTERSECTED = 1;
static const int NOT_INTERSECTED1 = 2;
static const int NOT_INTERSECTED2 = 3;
static const Vector3d ORIGIN = Vector3d(0, 0, 0);

static const int BI_DEGE_PLANE = 1;
static const int BI_DEGE_XOR_02 = 2;
static const int BI_DEGE_XOR_13 = 3;

class bilinear {
public:
	// v0, v1 are vertices of one triangle, v2, v3 are the vertices of another
	// one.
	bilinear(
		const Vector3d& v0,
		const Vector3d& v1,
		const Vector3d& v2,
		const Vector3d& v3);
	bool is_degenerated;
	std::vector<std::array<int, 3>> facets;
	std::array<int, 2> phi_f = { {2,2} };
	std::array<Vector3d, 4> v;

};
template <typename V1, typename V2> Vector3r cross(const V1& v1, const V2& v2)
{
    Vector3r res;
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return res;
}
int orient_3d(const Vector3d&p, const Vector3d&q, const Vector3d&r, const Vector3d& s);
int orient_2d(const Vector2d&p, const Vector2d&q, const Vector2d&r);

Rational func_g(
	const Vector3r& x,
	const std::array<Vector3d, 4>& corners,
	const std::array<int, 3>& indices);

Rational phi(const Vector3d x, const std::array<Vector3d, 4>& corners);
Rational phi(const Vector3r x, const std::array<Vector3d, 4>& corners);

void get_tet_phi(bilinear& bl);
//bool XOR(const bool a, const bool b)
//{
//    if (a && b)
//        return false;
//    if (!a && !b)
//        return false;
//    return true;
//}

// accept 0,1,2,3 as inputs
bool int_seg_XOR(const int a, const int b);

// accept -1,0,1,2,3 as inputs
int int_ray_XOR(const int a, const int b);
template <typename V> void print(const V& v)
{
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}

void write(const Vector3d& v, std::ostream& out);
Vector3d read(std::istream& in);


bool segment_segment_intersection_2d(
	const Vector2d& s0,
	const Vector2d& e0,
	const Vector2d& s1,
	const Vector2d& e1);

// 0 not intersected, 1 intersected, 2 s0 on segment
// can deal with degenerated cases
int ray_segment_intersection(
	const Vector3d& s0,
	const Vector3d& e0,
	const Vector3d& dir0,
	const Vector3d& s1,
	const Vector3d& e1);

// this function can also tell us if they are parallel and overlapped
// and also tell us if the parallel case has seg-seg overlapping:


bool is_triangle_degenerated(const Vector3d& t1, const Vector3d& t2, const Vector3d&t3);


bool same_point(const Vector3d& p1, const Vector3d& p2);


template<typename T>
static bool orient3D_LPI_prefilter_multiprecision(
	const T& px, const T& py, const T& pz, const T& qx, const T& qy, const T& qz,
	const T& rx, const T& ry, const T& rz, const T& sx, const T& sy, const T& sz, const T& tx, const T& ty, const T& tz,
	T& a11, T& a12, T& a13, T& d,T& n, const std::function<int(T)> &checker) {

	a11 = (px - qx);
	a12 = (py - qy);
	a13 = (pz - qz);
	T a21(sx - rx);
	T a22(sy - ry);
	T a23(sz - rz);
	T a31(tx - rx);
	T a32(ty - ry);
	T a33(tz - rz);
	T a2233((a22 * a33) - (a23 * a32));
	T a2133((a21 * a33) - (a23 * a31));
	T a2132((a21 * a32) - (a22 * a31));
	d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));
	int flag1 = checker(d);
	if (flag1 == -2 || flag1 == 0) {
		return false;// not enough precision
	}
	T px_rx(px - rx);
	T py_ry(py - ry);
	T pz_rz(pz - rz);

	n = ((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));// caution, this is -n actually

	a11 = a11 * n;
	a12 = a12 * n;
	a13 = a13 * n;
	return true;
}
template <typename T>
static int orient3D_LPI_postfilter_multiprecision(
    const T& a11,
    const T& a12,
    const T& a13,
    const T& d,
    const T& px,
    const T& py,
    const T& pz,
    const T& ax,
    const T& ay,
    const T& az,
    const T& bx,
    const T& by,
    const T& bz,
    const T& cx,
    const T& cy,
    const T& cz,
    const std::function<int(T)>& checker)
{

    T px_cx(px - cx);
    T py_cy(py - cy);
    T pz_cz(pz - cz);

    T d11((d * px_cx) + (a11));
    T d21(ax - cx);
    T d31(bx - cx);
    T d12((d * py_cy) + (a12));
    T d22(ay - cy);
    T d32(by - cy);
    T d13((d * pz_cz) + (a13));
    T d23(az - cz);
    T d33(bz - cz);

    T d2233(d22 * d33);
    T d2332(d23 * d32);
    T d2133(d21 * d33);
    T d2331(d23 * d31);
    T d2132(d21 * d32);
    T d2231(d22 * d31);

    T det(
        d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231));

    int flag2 = checker(det);
    if (flag2 == -2) {
        return 100; // not enough precision, only happens when using floating
                    // points
    }
    if (flag2 == 1) {
        if (d > 0) {
            return 1;
        }
        if (d < 0) {
            return -1;
        }
    }
    if (flag2 == -1) {
        if (d > 0) {
            return -1;
        }
        if (d < 0) {
            return 1;
        }
    }
    return 0;
}


static const   std::function<int(Rational)> check_rational = [](Rational v) {

	if (v.get_sign() > 0)
		return 1;
	if (v.get_sign() < 0)
		return -1;
	return 0;

};

int point_inter_triangle(
	const Vector3d&pt,
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	const bool& dege, const bool halfopen);

//already know lpi exist;
// 0 not intersected, 1 intersect open triangle, 2 shoot on edge, 3 shoot on edge t2-t3
int is_line_cut_triangle(
	const Vector3d& e0,
	const Vector3d& e1,
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	const bool halfopen);
int seg_triangle_inter_return_t(
	const Vector3d& e0,
	const Vector3d& e1,
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	Rational& t);
// if a line (going across pt, pt+dir) intersects triangle
// triangle is not degenerated
// 0 not intersected, 1 intersected, 3 intersected t2-t3 edge

// we check if triangle intersect segment,
// this function is used in cube edge--prism tri and cube edge--bilinear tri
// if halfopen= true, can tell us if intersect the edge t2-t3
// 0 not intersected, 1 intersected, 2 intersect edge, 3 intersect t2-t3 edge
int segment_triangle_intersection(
	const Vector3d& e0,
	const Vector3d& e1,
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	const bool halfopen);

int segment_triangle_intersection(
    const Vector3r& e0,
    const Vector3r& e1,
    const Vector3r& t1,
    const Vector3r& t2,
    const Vector3r& t3,
    const bool halfopen);
// 0 no intersection, 1 intersect, 2 point on triangle, 3 point or ray go to on t2-t3 edge, -1 shoot on border
int ray_triangle_intersection(
	const Vector3d& pt,
	const Vector3d& pt1,
	const Vector3d& dir,
	const Vector3d& t1,
	const Vector3d& t2,
	const Vector3d& t3,
	const bool halfopen);
void tri_bilinear(const bilinear &bl, int n, std::vector<std::array<Vector3r, 3>>& patch);
void save_obj(const std::string &name, const std::vector<std::array<Vector3r, 3>>& tris);
bool seg_discrete_bilinear_intersection(const bilinear &bl, int n,const Vector3d&s0,const Vector3d&s1);
int lpi_in_triangle(const Vector3d& p,const Vector3d& q,const Vector3d& r,const Vector3d& s,const Vector3d& t,bool halfopen);
bool lpi_rational(const Vector3d& p,const Vector3d& q,const Vector3d& r,const Vector3d& s,const Vector3d& t,
Rational&a11,Rational&a12,Rational&a13,Rational&d);
} // namespace ccd
