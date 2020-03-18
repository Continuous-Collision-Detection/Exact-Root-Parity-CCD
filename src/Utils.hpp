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
int int_XOR(const int a, const int b)// TODO
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

// this function can also tell us if they are parallel and overlapped
// and also tell us if the parallel case has seg-seg overlapping:




bool same_point(const Vector3r& p1, const Vector3r& p2);
Vector3r tri_norm(const Vector3r& t0, const Vector3r& t1, const Vector3r& t2)
{
	Vector3r s1, s2;
	s1 = t1 - t0;
	s2 = t2 - t1;
	return cross(s1, s2);
}
template<typename T>
static bool orient3D_LPI_prefilter_multiprecision(
	const T& px, const T& py, const T& pz, const T& qx, const T& qy, const T& qz,
	const T& rx, const T& ry, const T& rz, const T& sx, const T& sy, const T& sz, const T& tx, const T& ty, const T& tz,
	T& a11, T& a12, T& a13, T& d, const std::function<int(T)> &checker) {

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
	d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));//TODO maybe not safe
	int flag1 = checker(d);
	if (flag1 == -2 || flag1 == 0) {
		return false;// not enough precision
	}
	T px_rx(px - rx);
	T py_ry(py - ry);
	T pz_rz(pz - rz);

	T n((((py_ry)* a2133) - ((px_rx)* a2233)) - ((pz_rz)* a2132));

	a11 = a11 * n;
	a12 = a12 * n;
	a13 = a13 * n;
	return true;
}

template<typename T>
static bool orient3D_TPI_prefilter_multiprecision(
	const T& ov1x, const T& ov1y, const T& ov1z, const T& ov2x, const T& ov2y, const T& ov2z, const T& ov3x, const T& ov3y, const T& ov3z,
	const T& ow1x, const T& ow1y, const T& ow1z, const T& ow2x, const T& ow2y, const T& ow2z, const T& ow3x, const T& ow3y, const T& ow3z,
	const T& ou1x, const T& ou1y, const T& ou1z, const T& ou2x, const T& ou2y, const T& ou2z, const T& ou3x, const T& ou3y, const T& ou3z,
	T& d, T& n1, T& n2, T& n3, const std::function<int(T)> &checker
)
{
	::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

	T v3x(ov3x - ov2x);
	T v3y(ov3y - ov2y);
	T v3z(ov3z - ov2z);
	T v2x(ov2x - ov1x);
	T v2y(ov2y - ov1y);
	T v2z(ov2z - ov1z);
	T w3x(ow3x - ow2x);
	T w3y(ow3y - ow2y);
	T w3z(ow3z - ow2z);
	T w2x(ow2x - ow1x);
	T w2y(ow2y - ow1y);
	T w2z(ow2z - ow1z);
	T u3x(ou3x - ou2x);
	T u3y(ou3y - ou2y);
	T u3z(ou3z - ou2z);
	T u2x(ou2x - ou1x);
	T u2y(ou2y - ou1y);
	T u2z(ou2z - ou1z);

	T nvx(v2y * v3z - v2z * v3y);
	T nvy(v3x * v2z - v3z * v2x);
	T nvz(v2x * v3y - v2y * v3x);

	T nwx(w2y * w3z - w2z * w3y);
	T nwy(w3x * w2z - w3z * w2x);
	T nwz(w2x * w3y - w2y * w3x);

	T nux(u2y * u3z - u2z * u3y);
	T nuy(u3x * u2z - u3z * u2x);
	T nuz(u2x * u3y - u2y * u3x);

	T nwyuz(nwy * nuz - nwz * nuy);
	T nwxuz(nwx * nuz - nwz * nux);
	T nwxuy(nwx * nuy - nwy * nux);

	T nvyuz(nvy * nuz - nvz * nuy);
	T nvxuz(nvx * nuz - nvz * nux);
	T nvxuy(nvx * nuy - nvy * nux);

	T nvywz(nvy * nwz - nvz * nwy);
	T nvxwz(nvx * nwz - nvz * nwx);
	T nvxwy(nvx * nwy - nvy * nwx);

	d = (nvx * nwyuz - nvy * nwxuz + nvz * nwxuy);



	int flag1 = checker(d);
	if (flag1 == -2 || flag1 == 0) {
		return false;// not enough precision
	}

	T p1(nvx * ov1x + nvy * ov1y + nvz * ov1z);
	T p2(nwx * ow1x + nwy * ow1y + nwz * ow1z);
	T p3(nux * ou1x + nuy * ou1y + nuz * ou1z);

	n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
	n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
	n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;
	return true;
}

template<typename T>
static int orient3D_LPI_postfilter_multiprecision(
	const T& a11, const T& a12, const T& a13, const T& d,
	const T& px, const T& py, const T& pz,
	const T& ax, const T& ay, const T& az,
	const T& bx, const T& by, const T& bz,
	const T& cx, const T& cy, const T& cz, const std::function<int(T)> &checker) {

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

	T det(d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231));

	int flag2 = checker(det);
	if (flag2 == -2) {
		return 100;// not enough precision, only happens when using floating points
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

template<typename T>
static int orient3D_TPI_postfilter_multiprecision(
	const T& d, const T& n1, const T& n2, const T& n3,
	const T& q1x, const T& q1y, const T& q1z, const T& q2x, const T& q2y, const T& q2z, const T& q3x, const T& q3y, const T& q3z, const std::function<int(T)> &checker
)
{
	::feclearexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

	T dq3x(d * q3x);
	T dq3y(d * q3y);
	T dq3z(d * q3z);

	T a11(n1 - dq3x);
	T a12(n2 - dq3y);
	T a13(n3 - dq3z);
	T a21(q1x - q3x);
	T a22(q1y - q3y);
	T a23(q1z - q3z);
	T a31(q2x - q3x);
	T a32(q2y - q3y);
	T a33(q2z - q3z);

	T det(a11 * (a22*a33 - a23 * a32) - a12 * (a21*a33 - a23 * a31) + a13 * (a21*a32 - a22 * a31));

	int flag2 = checker(det);
	if (flag2 == -2) {
		return 100;// not enough precision
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

//already know lpi exist;
// 0 not intersected, 1 intersect open triangle, 2 shoot on edge, 3 shoot on edge t2-t3
int is_line_cut_triangle(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3, 
	const bool halfopen,  const Vector3r &norm) {
	
	////if (orient3d(n, t1, t2, t3) == 0) {
	//	//std::cout << "Degeneration happens" << std::endl;
	//	n = Vector3r(rand(), rand(), rand());
	//}

	Vector3r n = norm + t1;
	Rational a11, a12, a13, d;

	
	bool premulti = orient3D_LPI_prefilter_multiprecision(
		e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], 
		t1[0], t1[1], t1[2], t2[0], t2[1], t2[2], t3[0], t3[1], t3[2], a11, a12, a13, d, check_rational);
	if (premulti == false) return false;// point not exist, cannot happen

	int o1 = orient3D_LPI_postfilter_multiprecision(
		a11, a12, a13, d,
		e0[0], e0[1], e0[2],
		n[0], n[1], n[2],
		t1[0], t1[1], t1[2],
		t2[0], t2[1], t2[2],
		check_rational);
	int o2 = orient3D_LPI_postfilter_multiprecision(
		a11, a12, a13, d,
		e0[0], e0[1], e0[2],
		n[0], n[1], n[2],
		t2[0], t2[1], t2[2],
		t3[0], t3[1], t3[2],
		check_rational);// this edge
	int o3 = orient3D_LPI_postfilter_multiprecision(
		a11, a12, a13, d,
		e0[0], e0[1], e0[2],
		n[0], n[1], n[2],
		t3[0], t3[1], t3[2],
		t1[0], t1[1], t1[2],
		check_rational);// this edge

	if (halfopen) {
		if (o2 == 0 && o1 == o3)
			return 3;// on open edge t2-t3
	}

	if (o1 == o2 && o1 == o3)
		return 1;
	if (o1 == 0 && o2 * o3 >= 0)
		return 2;
	if (o2 == 0 && o1 * o3 >= 0)
		return 2;
	if (o3 == 0 && o2* o1 >= 0)
		return 2;

	return 0;
}
int line_triangle_inter_return_t(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	Rational& t)// TODO can not deal with degenerated triangle

{
	const Rational d = e0[0] * t1[1] * t2[2] - e0[0] * t1[1] * t3[2]
		- e0[0] * t1[2] * t2[1] + e0[0] * t1[2] * t3[1] + e0[0] * t2[1] * t3[2]
		- e0[0] * t2[2] * t3[1] - e0[1] * t1[0] * t2[2] + e0[1] * t1[0] * t3[2]
		+ e0[1] * t1[2] * t2[0] - e0[1] * t1[2] * t3[0] - e0[1] * t2[0] * t3[2]
		+ e0[1] * t2[2] * t3[0] + e0[2] * t1[0] * t2[1] - e0[2] * t1[0] * t3[1]
		- e0[2] * t1[1] * t2[0] + e0[2] * t1[1] * t3[0] + e0[2] * t2[0] * t3[1]
		- e0[2] * t2[1] * t3[0] - e1[0] * t1[1] * t2[2] + e1[0] * t1[1] * t3[2] 
		+ e1[0] * t1[2] * t2[1] - e1[0] * t1[2] * t3[1] - e1[0] * t2[1] * t3[2]
		+ e1[0] * t2[2] * t3[1] + e1[1] * t1[0] * t2[2] - e1[1] * t1[0] * t3[2]
		- e1[1] * t1[2] * t2[0] + e1[1] * t1[2] * t3[0] + e1[1] * t2[0] * t3[2]
		- e1[1] * t2[2] * t3[0] - e1[2] * t1[0] * t2[1] + e1[2] * t1[0] * t3[1]
		+ e1[2] * t1[1] * t2[0] - e1[2] * t1[1] * t3[0] - e1[2] * t2[0] * t3[1]
		+ e1[2] * t2[1] * t3[0];

	if (d.get_sign() == 0) // coplanar
		return -1;
	t = (e0[0] * t1[1] * t2[2] - e0[0] * t1[1] * t3[2] - e0[0] * t1[2] * t2[1]
		+ e0[0] * t1[2] * t3[1] + e0[0] * t2[1] * t3[2] - e0[0] * t2[2] * t3[1]
		- e0[1] * t1[0] * t2[2] + e0[1] * t1[0] * t3[2] + e0[1] * t1[2] * t2[0]
		- e0[1] * t1[2] * t3[0] - e0[1] * t2[0] * t3[2] + e0[1] * t2[2] * t3[0]
		+ e0[2] * t1[0] * t2[1] - e0[2] * t1[0] * t3[1] - e0[2] * t1[1] * t2[0]
		+ e0[2] * t1[1] * t3[0] + e0[2] * t2[0] * t3[1] - e0[2] * t2[1] * t3[0]
		- t1[0] * t2[1] * t3[2] + t1[0] * t2[2] * t3[1] + t1[1] * t2[0] * t3[2]
		- t1[1] * t2[2] * t3[0] - t1[2] * t2[0] * t3[1]
		+ t1[2] * t2[1] * t3[0])
		/ d;

	const Rational u = (-e0[0] * e1[1] * t1[2] + e0[0] * e1[1] * t3[2]
		+ e0[0] * e1[2] * t1[1] - e0[0] * e1[2] * t3[1]
		- e0[0] * t1[1] * t3[2] + e0[0] * t1[2] * t3[1]
		+ e0[1] * e1[0] * t1[2] - e0[1] * e1[0] * t3[2]
		- e0[1] * e1[2] * t1[0] + e0[1] * e1[2] * t3[0]
		+ e0[1] * t1[0] * t3[2] - e0[1] * t1[2] * t3[0]
		- e0[2] * e1[0] * t1[1] + e0[2] * e1[0] * t3[1]
		+ e0[2] * e1[1] * t1[0] - e0[2] * e1[1] * t3[0]
		- e0[2] * t1[0] * t3[1] + e0[2] * t1[1] * t3[0]
		+ e1[0] * t1[1] * t3[2] - e1[0] * t1[2] * t3[1]
		- e1[1] * t1[0] * t3[2] + e1[1] * t1[2] * t3[0]
		+ e1[2] * t1[0] * t3[1] - e1[2] * t1[1] * t3[0])
		/ d;
	const Rational v = (e0[0] * e1[1] * t1[2] - e0[0] * e1[1] * t2[2]
		- e0[0] * e1[2] * t1[1] + e0[0] * e1[2] * t2[1]
		+ e0[0] * t1[1] * t2[2] - e0[0] * t1[2] * t2[1]
		- e0[1] * e1[0] * t1[2] + e0[1] * e1[0] * t2[2]
		+ e0[1] * e1[2] * t1[0] - e0[1] * e1[2] * t2[0]
		- e0[1] * t1[0] * t2[2] + e0[1] * t1[2] * t2[0]
		+ e0[2] * e1[0] * t1[1] - e0[2] * e1[0] * t2[1]
		- e0[2] * e1[1] * t1[0] + e0[2] * e1[1] * t2[0]
		+ e0[2] * t1[0] * t2[1] - e0[2] * t1[1] * t2[0]
		- e1[0] * t1[1] * t2[2] + e1[0] * t1[2] * t2[1]
		+ e1[1] * t1[0] * t2[2] - e1[1] * t1[2] * t2[0]
		- e1[2] * t1[0] * t2[1] + e1[2] * t1[1] * t2[0])
		/ d;

	// std::cout << t << std::endl;

	// std::cout << u << std::endl;

	// std::cout << v << std::endl;

	if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u + v <= 1) {
		if (u == 0 || u == 1 || v == 0 || v == 1 || u + v == 1)
			return 2; // on the border
		return 1;
	}

	return 0;

}
// if a line (going across pt, pt+dir) intersects triangle
// triangle is not degenerated
// use ray_triangle_intersection twice, TODO modify this
// 0 not intersected, 1 intersected, 3 intersected t2-t3 edge 
int line_triangle_intersection(
	const Vector3r& pt,
	const Vector3r& dir,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen);
// we check if triangle intersect segment,
// this function is used in cube edge--prism tri and cube edge--bilinear tri
// if halfopen= true, can tell us if intersect the edge t2-t3
// 0 not intersected, 1 intersected, 2 one of pt is on plane, 3 intersect t2-t3 edge
int segment_triangle_intersection(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen);
// 0 no intersection, 1 intersect, 2 point on triangle, 3 point or ray go to on t2-t3 edge, -1 shoot on border
int ray_triangle_intersection(
	const Vector3r& pt,
	const Vector3r& dir,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen);
} // namespace ccd
