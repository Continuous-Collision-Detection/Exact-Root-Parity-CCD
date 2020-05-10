#include <Eigen/Dense>
#include <doubleCCD/double_Utils.hpp>
#include <doubleCCD/rayPlaneIntersection.h>
#include <fstream>
#include <predicates/indirect_predicates.h>

namespace doubleccd {
	bilinear::bilinear(
		const Vector3d& v0,
		const Vector3d& v1,
		const Vector3d& v2,
		const Vector3d& v3)
	{
		v = { { v0, v1, v2, v3 } };
		int ori = orient_3d(
			v0, v1, v2, v3);
		if (ori == 0) {
			is_degenerated = true;
		}
		else {
			is_degenerated = false;
		}
		if (ori >= 0) {
			facets.resize(4);// right hand outside order
			facets[0] = { { 1,2,0 } }; // 0,1 are one pair
			facets[1] = { { 3,0,2 } };

			facets[2] = { { 0,3,1 } }; // 2,3 are one pair
			facets[3] = { { 2,1,3 } };
		}
		if (ori == -1) {
			facets.resize(4);
			facets[0] = { { 1,0,2 } }; // 0,1 are one pair
			facets[1] = { { 3,2,0 } };

			facets[2] = { { 0,1,3 } }; // 2,3 are one pair
			facets[3] = { { 2,3,1 } };
		}
	}

	bool same_point(const Vector3d& p1, const Vector3d& p2)
	{
		if (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) {
			return true;
		}
		return false;
	}

	//int orient2d(
	//    const Vector3r& a, const Vector3r& b, const Vector3r& c, const int axis)
	//{
	//    Vector3r v1, v2;
	//    v1 = a - b;
	//    v2 = c - b;
	//    int i1, i2;
	//    if (axis == 0) {
	//        i1 = 1;
	//        i2 = 2;
	//    }
	//    if (axis == 1) {
	//        i1 = 0;
	//        i2 = 2;
	//    }
	//    if (axis == 2) {
	//        i1 = 0;
	//        i2 = 1;
	//    }
	//    const Rational det = v1[i1] * v2[i2] - v1[i2] * v2[i1];
	//    return det.get_sign();
	//}
	Rational phi(const Vector3d xd, const std::array<Vector3d, 4>& corners)
	{
		static const std::array<int, 4> vv = { { 0, 1, 2, 3 } };
		Vector3r x(xd[0], xd[1], xd[2]);
		const Rational g012 = func_g(x, corners, { { vv[0], vv[1], vv[2] } });
		const Rational g132 = func_g(x, corners, { { vv[1], vv[3], vv[2] } });
		const Rational g013 = func_g(x, corners, { { vv[0], vv[1], vv[3] } });
		const Rational g032 = func_g(x, corners, { { vv[0], vv[3], vv[2] } });

		const Rational h12 = g012 * g032;
		const Rational h03 = g132 * g013;

		const Rational phi = h12 - h03;

		return phi;
	}
	Rational phi(const Vector3r x, const std::array<Vector3d, 4>& corners)
	{
		static const std::array<int, 4> vv = { { 0, 1, 2, 3 } };
		const Rational g012 = func_g(x, corners, { { vv[0], vv[1], vv[2] } });
		const Rational g132 = func_g(x, corners, { { vv[1], vv[3], vv[2] } });
		const Rational g013 = func_g(x, corners, { { vv[0], vv[1], vv[3] } });
		const Rational g032 = func_g(x, corners, { { vv[0], vv[3], vv[2] } });

		const Rational h12 = g012 * g032;
		const Rational h03 = g132 * g013;

		const Rational phi = h12 - h03;

		return phi;
	}


	Rational func_g(
		const Vector3r& xr,
		const std::array<Vector3d, 4>& corners,
		const std::array<int, 3>& indices)
	{
		const int p = indices[0];
		const int q = indices[1];
		const int r = indices[2];
		Vector3r  pr(corners[p][0], corners[p][1], corners[p][2]),
			qr(corners[q][0], corners[q][1], corners[q][2]), rr(corners[r][0], corners[r][1], corners[r][2]);
		return (xr - pr)
			.dot(cross(qr - pr, rr - pr));// TODO maybe minus is not allowed
	}
	bool int_seg_XOR(const int a, const int b)
	{
		if (a == 2 || b == 2) return true;
		if (a == 0 && b == 1) return true;
		if (a == 1 && b == 0) return true;
		if (a == 3 || b == 3) return false;
		if (a == b) return false;
		return false;
		/*if (a == -1 || b == -1)
			return -1;
		if (a == b)
			return 0;
		if (a == 0)
			return b;
		if (b == 0)
			return a;
		std::cout << "impossible XOR cases" << std::endl;*/
	}
	int int_ray_XOR(const int a, const int b) {
		if (a == -1 || b == -1) return -1;
		if (a == 0) return b;
		if (b == 0) return a;
		if (a == b) return 0;
		if (a == 2 || b == 2) return 2;//this is case 2-3
		std::cout << "impossible to go here " << std::endl;
		return -1;
	}

	// parallel means not intersected
	int is_line_cut_triangle(
		const Vector3d& e0,
		const Vector3d& e1,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3,
		const bool halfopen) {


		explicitPoint3D p(e0[0], e0[1], e0[2]);
		explicitPoint3D q(e1[0], e1[1], e1[2]);
		explicitPoint3D a(t1[0], t1[1], t1[2]);
		explicitPoint3D b(t2[0], t2[1], t2[2]);
		explicitPoint3D c(t3[0], t3[1], t3[2]);
		implicitPoint3D_LPI l(p, q, a, b, c);
		explicitPoint3D ppp;
		///std::cout<<"in lpi, check if point exist"<<std::endl;

		if (!l.approxExplicit(ppp)) return 0;// this is important
		//std::cout<<"in lpi, already know point exist"<<std::endl;
		//int o1=orient_3d(e0,t1,t2,t3);
		//int o2=orient_3d(e1,t1,t2,t3);
		//std::cout<<"*orient "<<o1<<" "<<o2<<std::endl;
		//std::cout<<"orient "<<orient_3d(e0,e1,t1,t2)<<std::endl;
		//std::cout<<"orient "<<orient_3d(e0,e1,t2,t3)<<std::endl;
		//std::cout<<"orient "<<orient_3d(e0,e1,t1,t3)<<std::endl;
		//std::cout<<"inter bc? "<<genericPoint::pointInSegment(l, b, c)<<std::endl;
		{//this piece is problematic
		// if (genericPoint::pointInInnerTriangle(l, a, b, c))
			// 	return 1;
		}
		
			
		Vector3d ap=Vector3d::Random();
		while(orient_3d(ap,t1,t2,t3)==0){
			ap=Vector3d::Random();
		}
		explicitPoint3D ape(ap[0],ap[1],ap[2]);
		//std::cout<<"random point ori "<<orient_3d(ap,t1,t2,t3)<<std::endl;
		int op1=genericPoint::orient3D(l,ape,a,b);
		int op2=genericPoint::orient3D(l,ape,b,c);
		int op3=genericPoint::orient3D(l,ape,c,a);
		if(op1==op2&&op1==op3)//lpi in open triangle
            return 1;
                //std::cout<<"the orientations in double "<<op1<<" "<<op2<<" "<<op3<<std::endl;
			//std::cout<<"l\n"<<e0<<"\n\n"<<e1<<std::endl;
			//std::cout<<"a and b and c\n"<<t1<<"\n\n"<<t2<<"\n\n"<<t3<<std::endl;
		// if (genericPoint::pointInSegment(l, a, c)){
		// 	std::cout<<"l on ac"<<std::endl;
		// 	return 2;
		// }
		// if (genericPoint::pointInSegment(l, b, a)){
		// 	std::cout<<"seg seg check "<<genericPoint::segmentsCross(p,q,a,b)<<std::endl;
		// 	std::cout<<"l on ab"<<std::endl;
		// 	std::cout<<"l\n"<<e0<<"\n\n"<<e1<<std::endl;
		// 	std::cout<<"a and b and c\n"<<t1<<"\n\n"<<t2<<"\n\n"<<t3<<std::endl;
		// 	std::cout<<"orientation "<<orient_3d(e0,e1,t1,t2)<<std::endl;
		// 	return 2;
		// }
			
		
		if(orient_3d(e0,e1,t1,t3)==0&&genericPoint::segmentsCross(p,q,a,c)){
			//std::cout<<"l on ac"<<std::endl;
			return 2;
		}
		if (orient_3d(e0,e1,t1,t2)==0&&genericPoint::segmentsCross(p,q,a,b)){
			//std::cout<<"l on ab"<<std::endl;
			//std::cout<<"l\n"<<e0<<"\n\n"<<e1<<std::endl;
			//std::cout<<"a and b and c\n"<<t1<<"\n\n"<<t2<<"\n\n"<<t3<<std::endl;
			//std::cout<<"orientation "<<orient_3d(e0,e1,t1,t2)<<std::endl;
			return 2;
		}
		
			
		//if (genericPoint::pointInSegment(l, b, c))
		if (orient_3d(e0,e1,t3,t2)==0&&genericPoint::segmentsCross(p,q,c,b)){
			if (halfopen)
				return 3;// on open edge t2-t3
			else{
				//std::cout<<"return the last 2"<<std::endl;
				return 2;
			}
		}
			
				

		return 0;
	}

	//  when lpi does not exist, return 0
	int seg_triangle_inter_return_t(
		const Vector3d& e0,
		const Vector3d& e1,
		const Vector3d& t1d,
		const Vector3d& t2d,
		const Vector3d& t3d,
		Rational& t)

	{
		int o1 = orient_3d(e0, t1d, t2d, t3d);
		int o2 = orient_3d(e1, t1d, t2d, t3d);
		if (o1 == o2) return 0;//already know lpi exist, check if seg go cross the plane

		int inter = segment_triangle_intersection(e0, e1, t1d, t2d, t3d, false);
		if (inter == 0) return 0;// not intersected
		if (inter == 2) std::cout << "wrong here seg_triangle_inter_return_t" << std::endl;
		Rational a11, a12, a13, d, n;
		bool result = orient3D_LPI_prefilter_multiprecision(
			Rational(e0[0]), Rational(e0[1]), Rational(e0[2]), Rational(e1[0]), Rational(e1[1]), Rational(e1[2]),
			Rational(t1d[0]), Rational(t1d[1]), Rational(t1d[2]), Rational(t2d[0]), Rational(t2d[1]), Rational(t2d[2]),
			Rational(t3d[0]), Rational(t3d[1]), Rational(t3d[2]), a11, a12, a13, d, n, check_rational);

		t = -n / d;// in the function, n is actually -n
		return 1;

	}
	void get_tet_phi(bilinear& bl)
	{
		Rational v0, v1, v2;
		v0 = Rational(bl.v[0][0]) + Rational(bl.v[2][0]);
		v1 = Rational(bl.v[0][1]) + Rational(bl.v[2][1]);
		v2 = Rational(bl.v[0][2]) + Rational(bl.v[2][2]);
		Vector3r p02;
		p02[0] = v0 / Rational(2); p02[1] = v1 / Rational(2); p02[2] = v2 / Rational(2);


		Rational phi02 = phi(p02, bl.v);
		if (phi02.get_sign() > 0) {
			bl.phi_f[0] = 1;
			bl.phi_f[1] = -1;
			return;
		}
		else {
			bl.phi_f[0] = -1;
			bl.phi_f[1] = 1;
			return;
		}
		std::cout << "!!can not happen, get tet phi" << std::endl;
	}


	// point and seg are colinear, now check if point is on the segment(can deal with segment degeneration)
	bool colinear_point_on_segment(
		const Vector3d& pt, const Vector3d& s0, const Vector3d& s1)
	{
		if (std::min(s0[0], s1[0]) <= pt[0] <= std::max(s0[0], s1[0]))
			if (std::min(s0[1], s1[1]) <= pt[1] <= std::max(s0[1], s1[1]))
				if (std::min(s0[2], s1[2]) <= pt[2] <= std::max(s0[2], s1[2]))
					return true;
		return false;
	}
	bool point_on_segment(const Vector3d& pt, const Vector3d& s0, const Vector3d& s1) {
		if (!is_triangle_degenerated(pt, s0, s1))
			return false;
		return colinear_point_on_segment(pt, s0, s1);
	}


	int orient_2d(const Vector2d&p, const Vector2d&q, const Vector2d&r) {
		return orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
	}
	int orient_3d(const Vector3d&p, const Vector3d&q, const Vector3d&r, const Vector3d& s) {
		return -orient3d(p[0], p[1], p[2], q[0], q[1], q[2], r[0], r[1], r[2], s[0], s[1], s[2]);
	}
	//TODO orient 3d we flip the sign
	bool point_on_segment_2d(const Vector2d&p, const Vector2d&s0, const Vector2d&s1, bool need_check_colinear) {
		if (need_check_colinear) {
			if (orient_2d(p, s0, s1) != 0)
				return false;
		}
		if (std::min(s0[0], s1[0]) <= p[0] <= std::max(s0[0], s1[0]))
			if (std::min(s0[1], s1[1]) <= p[1] <= std::max(s0[1], s1[1]))
				return true;
		return false;
	}
	bool segment_segment_intersection_2d(
		const Vector2d& s0,
		const Vector2d& e0,
		const Vector2d& s1,
		const Vector2d& e1) {
		int o0 = orient_2d(s0, s1, e1);
		int o1 = orient_2d(e0, s1, e1);
		if (o0 == o1 && o0 + o1 != 0) return false;
		int o2 = orient_2d(s1, s0, e0);
		int o3 = orient_2d(e1, s0, e0);
		if (o2 == o3 && o2 + o3 != 0) return false;
		if (o0 != o1 && o2 != o3) return true;
		// the rest of cases: one seg degenerate, or they are colinear
		if (point_on_segment_2d(s0, s1, e1, false) || point_on_segment_2d(e0, s1, e1, false)
			|| point_on_segment_2d(s1, s0, e0, false) || point_on_segment_2d(e1, s0, e0, false))
			return true;
		return false;
	}

	// 0 not intersected; 1 intersected; 2 pt on s0
	int point_on_ray(const Vector3d& s0, const Vector3d& e0,
		const Vector3d& dir0, const Vector3d& pt) {
		if (same_point(s0, pt)|| same_point(e0, pt))
			return 2;
		if (!is_triangle_degenerated(s0, e0, pt)) return 0;

		// now the pt is on the line
		if (dir0[0] > 0) {
			if (pt[0] <= s0[0]) return 0;
		}
		if (dir0[0] < 0) {
			if (pt[0] >= s0[0]) return 0;
		}
		if (dir0[0] == 0) {
			if (pt[0] != s0[0]) return 0;
		}
		//
		if (dir0[1] > 0) {
			if (pt[1] <= s0[1]) return 0;
		}
		if (dir0[1] <= 0) {
			if (pt[1] >= s0[1]) return 0;
		}
		if (dir0[1] == 0) {
			if (pt[1] != s0[1]) return 0;
		}
		//
		if (dir0[2] > 0) {
			if (pt[2] <= s0[2]) return 0;
		}
		if (dir0[2] <= 0) {
			if (pt[2] >= s0[2]) return 0;
		}
		if (dir0[2] == 0) {
			if (pt[2] != s0[2]) return 0;
		}
		return 1;

	}
	// 0 not intersected, 1 intersected, 2 s0 on segment
	// can deal with degenerated cases
	int ray_segment_intersection(
		const Vector3d& s0,
		const Vector3d& e0,
		const Vector3d& dir0,
		const Vector3d& s1,
		const Vector3d& e1) {

		if (same_point(e1, s1))//degenerated case
			return point_on_ray(s0, e0, dir0, s1);

		/////////////////////////////////////
		if (orient_3d(s0, e0, s1, e1) != 0) return 0;

		if (point_on_segment(s0, s1, e1))
			return 2;
		if (!is_triangle_degenerated(s1, s0, e0)) {
			// we can get a point out of the plane
			Vector3d np = Vector3d::Random();
			while (orient_3d(np, s1, s0, e0) == 0) {// if coplanar, random
				np = Vector3d::Random();
			}
			int o1 = orient_3d(np, e1, s0, e0);
			if (o1 == 0) return point_on_ray(s0, e0, dir0, e1);
			int o2 = -1 * orient_3d(np, s1, s0, e0);// already know this is not 0
			int oo= orient_3d(np, e1, s0, s1);
			if (oo == 0) {// actually can directly return 0 because s0-s1-e0 is not degenerated
				if (point_on_ray(s0, e0, dir0, s1) > 0 || point_on_ray(s0, e0, dir0, e1) > 0)return 1;
				return 0;
			}
			if (o1 == oo && o2 == oo) return 1;
			return 0;
		}
		else {
			// s1 is on the line, if s1 is on ray, return 1; else return 0.(s0 is not on seg)
			if (point_on_ray(s0, e0, dir0, s1) == 1) return 1;
			return 0;
		}
		return 0;

	}

	void write(const Vector3d& v, std::ostream& out)
	{
		out.write(reinterpret_cast<const char*>(&v[0]), sizeof(v[0]));
		out.write(reinterpret_cast<const char*>(&v[1]), sizeof(v[1]));
		out.write(reinterpret_cast<const char*>(&v[2]), sizeof(v[2]));
	}

	Vector3d read(std::istream& in)
	{
		Vector3d res;
		double tmp;
		in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
		res[0] = tmp;

		in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
		res[1] = tmp;

		in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
		res[2] = tmp;

		return res;
	}
	// 2 on edge, 1 interior, 0 not intersect, 3 intersect OPEN edge t2-t3
	// norm follows right hand law
	int point_inter_triangle(
		const Vector3d&pt,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3,
		const bool& dege, const bool halfopen) {
		if (dege) {// check 2 edges are enough
			if (point_on_segment(pt, t1, t2))
				return 2;
			if (point_on_segment(pt, t1, t3))
				return 2;
			return 0;
		}
		else {
			if (orient_3d(pt, t1, t2, t3) != 0)
				return false;
			/*if (point_on_segment(pt, t1, t2))
				return 2;
			if (point_on_segment(pt, t1, t3))
				return 2;
			if (point_on_segment(pt, t2, t3))
				return 2;*///no need to do above
			Vector3d np = Vector3d::Random();
			while (orient_3d(np, t1, t2, t3) == 0) {// if coplanar, random
				np = Vector3d::Random();
			}
			int o1 = orient_3d(pt, np, t1, t2);
			int o2 = orient_3d(pt, np, t2, t3);// this edge
			int o3 = orient_3d(pt, np, t3, t1);
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
	}

	bool is_triangle_degenerated(const Vector3d& t1, const Vector3d& t2, const Vector3d&t3)
	{

		const auto to_2d = [](const Vector3d &p, int t)
		{
			return Vector2d(p[(t + 1) % 3], p[(t + 2) % 3]);
		};
		double r = ((t1 - t2).cross(t1 - t3)).norm();
		if (fabs(r) > 1e-8) return false;
		int ori;
		std::array<Vector2d, 3> p;
		for (int j = 0; j < 3; j++)
		{

			p[0] = to_2d(t1, j);
			p[1] = to_2d(t2, j);
			p[2] = to_2d(t3, j);

			ori = orient2d(p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1]);
			if (ori != 0)
			{
				return false;
			}
		}
		return true;
	}

	// we check if triangle intersect segment,
	// this function is used in cube edge--prism tri and cube edge--bilinear tri
	// if halfopen= true, can tell us if intersect the edge t2-t3
	// return 0, 1, 2, 3
	// this function is only used to check if seg intersect no degenerated bilinear
	//and if seg intersect opposite facets of tet.
	int segment_triangle_intersection(
		const Vector3d& e0,
		const Vector3d& e1,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3,
		const bool halfopen) {

		if (is_triangle_degenerated(t1, t2, t3))
			return 0;// since we already checked triangle edge - box segment

		int o1 = orient_3d(e0, t1, t2, t3);
		int o2 = orient_3d(e1, t1, t2, t3);
		if (o1 > 0 && o2 > 0)
			return 0;
		if (o1 < 0 && o2 < 0)
			return 0;

		if (o1 == 0 && o2 != 0)// e0 on the plane
			return point_inter_triangle(e0, t1, t2, t3, false, halfopen);
		if (o2 == 0 && o1 != 0)
			return point_inter_triangle(e1, t1, t2, t3, false, halfopen);
		if (o1 == 0 && o2 == 0)// two points are all on the plane. we already checked seg-two edges before
		{
			if (same_point(e0, e1))
				return point_inter_triangle(e0, t1, t2, t3, false, halfopen);
			else {
				int pinter0 = point_inter_triangle(e0, t1, t2, t3, false, halfopen);//2 is impossible
				int pinter1 = point_inter_triangle(e1, t1, t2, t3, false, halfopen);
				// two points all inside, inside; two outside, outside; others, intersect t2-t3
				if (pinter0 == 1 && pinter1 == 1) return 1; // interior
				if (pinter0 == 0 && pinter1 == 0) return 0;//  out
				if (halfopen)
					return 3;// intersect t2-t3
				else
					return 2;
			}
		}
		//std::cout<<"in xor2 go before line_trianlge test"<<std::endl;
		return is_line_cut_triangle(e0, e1, t1, t2, t3, halfopen);
	}
	// pt is not on the plane
	bool is_ray_intersect_plane(
		const Vector3d& pt,
		const Vector3d& dir,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3) {
		explicitPoint3D p(pt[0], pt[1], pt[2]), d(dir[0], dir[1], dir[2]);
		explicitPoint3D T1(t1[0], t1[1], t1[2]), T2(t2[0], t2[1], t2[2]), T3(t3[0], t3[1], t3[2]);
		return rayIntersectsPlane(p, d, T1, T2, T3);
	}
	Vector3r double2r(const Vector3d& pt){
		return Vector3r(pt[0],pt[1],pt[2]);
	}
	bool is_ray_intersect_plane_rational(const Vector3d& pt,
		const Vector3d& dir,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3){
			Vector3r ptr=double2r(pt),dirr=double2r(dir),t1r=double2r(t1),
			t2r=double2r(t2),t3r=double2r(t3),np;
			Vector3r l=dirr,l0=ptr,p0=t1r,n=cross(t1r-t2r,t1r-t3r);
			auto d=(p0-l0).dot(n)/(l.dot(n));
			return d>=0;
		}

	// 0 no intersection, 1 intersect, 2 point on triangle(including two edges), 3 point or ray shoot t2-t3 edge, -1 shoot on border
	int ray_triangle_intersection(
		const Vector3d& pt,
		const Vector3d& pt1,
		const Vector3d& dir,
		const Vector3d& t1,
		const Vector3d& t2,
		const Vector3d& t3,
		const bool halfopen) {


		if (is_triangle_degenerated(t1, t2, t3))// triangle degeneration
		{
			//std::cout<<"tri degenerated "<<std::endl;
			int inter1 = ray_segment_intersection(pt, pt1, dir, t1, t2);
			if (inter1 == 1) return -1;
			if (inter1 == 2) return 2;

			int inter2 = ray_segment_intersection(pt, pt1, dir, t1, t3);// check two segs is enough
			if (inter2 == 1) return -1;
			if (inter2 == 2) return 2;

			return 0;
		}

		int o1 = orient_3d(pt, t1, t2, t3);
		if (o1 == 0) {//point is on the plane
			int inter = point_inter_triangle(pt, t1, t2, t3, false, halfopen);
			if (inter == 1 || inter == 2) return 2;
			if (inter == 3) return 3;
			//if (inter == 0)
			else {// pt on the plane but not intersect triangle.
				if (orient_3d(pt1, t1, t2, t3) == 0) {// if ray is on the plane
					int inter1 = ray_segment_intersection(pt, pt1, dir, t1, t2);
					if (inter1 == 1) return -1;
					if (inter1 == 2) return 2;// acutally, cannot be 2 because already checked by pt_inter_tri

					int inter2 = ray_segment_intersection(pt, pt1, dir, t1, t3);
					if (inter2 == 1) return -1;
					if (inter2 == 2) return 2;
					//actually since point do not intersect triangle, check two segs are enough
					//int inter3 = ray_segment_intersection(pt, dir, t2, t3);
					//// ray overlaps t2-t3, shoot another ray, ray intersect it, shoot another one
					//if (inter3 == 1) return -1;
					//if (inter3 == 2) return 2;

					return 0;
				}


			}
			return 0;
		}

		// if point not on plane, and not point to plane, return 0
		//explicitPoint3D pte()
		bool in1=is_ray_intersect_plane_rational(pt, dir, t1, t2, t3);
		//bool in2=is_ray_intersect_plane(pt, dir, t1, t2, t3);
		//if(in1!=in2)std::cout<<"ray intersection don't match "<<in1<<" "<<in2<<std::endl;
		//TODO we are using Rational ray_tri_intersection
		if (!in1) return 0;
		//std::cout<<"ray inter plane"<<std::endl;
		// if ray go across the plane, then get lpi and 3 orientations
		int inter = is_line_cut_triangle(pt, pt1, t1, t2, t3, halfopen);
		//std::cout<<"is line cut triangle "<<inter<<std::endl;
		//std::cout<<"line cut tri? "<<inter<<std::endl;
		if (inter == 0)return 0;
		if (inter == 1)return 1;
		if (inter == 2)return -1;//shoot on edge
		if (inter == 3) return 3;

		return 0;
	}
	Vector3r sum(const Vector3r& a, const Vector3r& b) {
		Vector3r c;
		c[0] = a[0] + b[0];
		c[1] = a[1] + b[1];
		c[2] = a[2] + b[2];
		return c;
	}
	Vector3r para_p(const Vector3r&p1, const Vector3r&p2, int n, int i) {
		return sum(p1, (p2 - p1)*double(i) / double(n));
	};
	void tri_bilinear(const bilinear &bl, int n, std::vector<std::array<Vector3r, 3>>& patch) {
		patch.clear();

		const auto get_p = [](const Vector3r&v0, const Vector3r&v1, const Vector3r&v2, const Vector3r&v3,
			int n, int i, int j) {
			Vector3r p1, p2, p;
			p1 = para_p(v0, v1, n, i);
			p2 = para_p(v3, v2, n, i);
			p = para_p(p1, p2, n, j);
			return p;
		};
		std::array<Vector3r, 3> tri;
		Vector3r v0(bl.v[0][0], bl.v[0][1], bl.v[0][2]), v1(bl.v[1][0], bl.v[1][1], bl.v[1][2]), 
			v2(bl.v[2][0], bl.v[2][1], bl.v[2][2]), v3(bl.v[3][0], bl.v[3][1], bl.v[3][2]);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				tri[0] = get_p(v0, v1, v2, v3, n, i, j);
				tri[1] = get_p(v0, v1, v2, v3, n, i, j + 1);
				tri[2] = get_p(v0, v1, v2, v3, n, i + 1, j);
				patch.emplace_back(tri);
				tri[0] = get_p(v0, v1, v2, v3, n, i, j + 1);
				tri[1] = get_p(v0, v1, v2, v3, n, i + 1, j + 1);
				tri[2] = get_p(v0, v1, v2, v3, n, i + 1, j);
				patch.emplace_back(tri);
			}
		}
	}
	void save_obj(const std::string &name, const std::vector<std::array<Vector3r, 3>>& tris) {
		std::ofstream out(name);
		for (int i = 0; i < tris.size(); i++) {
			out << "v " << tris[i][0][0] << " " << tris[i][0][1] << " " << tris[i][0][2] << std::endl;
			out << "v " << tris[i][1][0] << " " << tris[i][1][1] << " " << tris[i][1][2] << std::endl;
			out << "v " << tris[i][2][0] << " " << tris[i][2][1] << " " << tris[i][2][2] << std::endl;
		}
		for (int i = 0; i < tris.size(); i++) {
			out << "f " << i * 3 + 1 << " " << i * 3 + 2 << " " << i * 3 + 3 << std::endl;
		}
		out.close();
	}
	int orient3d(
    const Vector3r& a, const Vector3r& b, const Vector3r& c, const Vector3r& d){
		const Rational det = (a - d).dot(cross(b - d, c - d));
    return det.get_sign();
	}
	int is_line_cut_triangle(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen, const Vector3r &norm) {

	////if (orient3d(n, t1, t2, t3) == 0) {
	//	//std::cout << "Degeneration happens" << std::endl;
	//	n = Vector3r(rand(), rand(), rand());
	//}

	Vector3r n = norm + t1;
	Rational a11, a12, a13, d,nr;


	bool premulti = orient3D_LPI_prefilter_multiprecision(
		e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
		t1[0], t1[1], t1[2], t2[0], t2[1], t2[2], t3[0], t3[1], t3[2], a11, a12, a13, d, nr,check_rational);
	if (premulti == false) return 0;

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
	int segment_triangle_intersection(
	const Vector3r& e0,
	const Vector3r& e1,
	const Vector3r& t1,
	const Vector3r& t2,
	const Vector3r& t3,
	const bool halfopen) {

	Vector3r norm = cross(t2 - t1, t3 - t2);
	if (norm[0]==0&&norm[1]==0&&norm[2]==0)// triangle degeneration
	{
		return 0;//we already checked triangle (at least two edges)edge against cube
		//TODO need to make it a general one?
	}

	int o1 = orient3d(e0, t1, t2, t3);
	int o2 = orient3d(e1, t1, t2, t3);
	if (o1 > 0 && o2 > 0)
		return 0;
	if (o1 < 0 && o2 < 0)
		return 0;

	if (o1 == 0 && o2 != 0)// e0 on the plane
		return is_line_cut_triangle(e0, e1, t1, t2, t3, halfopen, norm);
	if (o2 == 0 && o1 != 0)
		return is_line_cut_triangle(e0, e1, t1, t2, t3, halfopen, norm);
	if (o1 == 0 && o2 == 0)// two points are all on the plane. we already checked seg-two edges before
	{
		std::cout<<"seg and face parallel, we don't consider"<<std::endl;
		return 0;
		// if (same_point(e0,e1))
		// 	return point_inter_triangle(e0, t1, t2, t3, false, norm, halfopen);
		// else {
		// 	int pinter0 = point_inter_triangle(e0, t1, t2, t3, false, norm, halfopen);//2 is impossible
		// 	int pinter1 = point_inter_triangle(e1, t1, t2, t3, false, norm, halfopen);
		// 	// two points all inside, inside; two outside, outside; others, intersect t2-t3
		// 	if (pinter0 == 1 && pinter1 == 1) return 1; // interior
		// 	if (pinter0 == 0 && pinter1 == 0) return 0;//  out
		// 	if (halfopen)
		// 		return 3;// intersect t2-t3
		// 	else
		// 		return 2;
		// }
	}
	return is_line_cut_triangle(e0, e1, t1, t2, t3, halfopen, norm);
}

bool seg_discrete_bilinear_intersection(const bilinear &bl, int n,const Vector3d&s0,const Vector3d&s1){
std::vector<std::array<Vector3r, 3>> patch;
tri_bilinear(bl, n,patch );
Vector3r s0r(s0[0],s0[1],s0[2]),s1r(s1[0],s1[1],s1[2]);
for(int i=0;i<patch.size();i++){
	if(segment_triangle_intersection(s0r,s1r,patch[i][0],patch[i][1],patch[i][2],false)>0){
		return true;
	}

}
return false;
}


} // namespace ccd
