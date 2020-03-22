#pragma once
#include <Utils.hpp>
#include <subfunctions.h>
namespace ccd {

	// before going here, we already know the point can not be on the shape
	int ray_degenerated_bilinear_parity(
		const bilinear& bl,
		const Vector3r& pt,
		const Vector3r& dir,
		const int dege
	)
	{
		
		int r1, r2;
		if (dege == BI_DEGE_PLANE) {
			r1 = ray_triangle_intersection(// -1, 0, 1, 2
				pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
				bl.v[bl.facets[0][2]], true);
			if (r1 == 2)
				return 2;
			if (r1 == -1)
				return -1;
			if (r1 == 1)
				return 1;
			if (r1 == 3) {
				r1 = ray_triangle_intersection(// -1, 0, 1, 2
					pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
					bl.v[bl.facets[0][2]], false);//  this should have 3, if 3, see it as 1
				if (r1 == 2) return 2;// point on t2-t3 edge
				return 1;// ray go through t2-t3 edge
			}
			r2 = ray_triangle_intersection(
				pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
				bl.v[bl.facets[1][2]], false);
			if (r2 == 2)
				return 2;
			if (r2 == -1)
				return -1;
			if (r2 == 1)
				return 1;
			return 0;

		}
		else {

			if (dege == BI_DEGE_XOR_02) { // triangle 0-1-2 and 0-2-3
				r1 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
					bl.v[bl.facets[0][2]], true);//0: not hit, 1: hit on open triangle, 2: pt on halfopen T, since already checked, accept it 
				r2 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
					bl.v[bl.facets[1][2]], true);
				return int_ray_XOR(r1, r2);
			}

			if (dege == BI_DEGE_XOR_13) { // triangle 0-1-3 and 3-1-2
				r1 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[2][0]], bl.v[bl.facets[2][1]],
					bl.v[bl.facets[2][2]], true);//0: not hit, 1: hit on open triangle, 2: pt on halfopen T, since already checked, accept it 
				r2 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[3][0]], bl.v[bl.facets[3][1]],
					bl.v[bl.facets[3][2]], true);
				return int_ray_XOR(r1, r2);
			}
		}
		std::cout << "!! THIS CANNOT HAPPEN" << std::endl;
		return false;
	}
	// the point is inside of tet, tet is not degenerated
	//return: -1,1,0
	// phi_p is the phi of end point of ray
	int ray_correct_bilinear_face_pair_inter(
		const Vector3r& p,
		const Rational& phi_p,
		const Vector3r& dir,
		const bilinear& bl)
	{
		int r1, r2;
		/*if (phi_p.get_sign() == 0)
			return 2;*/
		if (bl.phi_f[0] * phi_p.get_sign() < 0) {
			r1 = ray_triangle_intersection(// -1,0,1,2,3
				p, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
				bl.v[bl.facets[0][2]], true);
			r2 = ray_triangle_intersection(
				p, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
				bl.v[bl.facets[1][2]], true);
			// when there is -1, -1(shoot on one of two edges); impossible to have 2; when there is 1, 1
			// when there is 3, pt is on t2-t3 edge(impossible), or ray go across that edge, parity +1, return 1
			// since pt is inside of tet, then impossible on any plane of the two triangles
			if (r1 == -1 || r2 == -1)
				return -1;
			if (r1 > 0 || r2 > 0)
				return 1; // cannot be degenerated, so this can work
			return 0;
		}
		else {
			r1 = ray_triangle_intersection(
				p, dir, bl.v[bl.facets[2][0]], bl.v[bl.facets[2][1]],
				bl.v[bl.facets[2][2]],true);
			r2 = ray_triangle_intersection(
				p, dir, bl.v[bl.facets[3][0]], bl.v[bl.facets[3][1]],
				bl.v[bl.facets[3][2]],true);
			if (r1 == -1 || r2 == -1)
				return -1;
			if (r1 > 0 || r2 > 0)
				return 1; // cannot be degenerated, so this can work
			return 0;
		}
	}

	int ray_bilinear_parity(
		bilinear& bl,
		const Vector3r& pt,
		const Vector3r& dir,
		const bool is_degenerated,
		const bool is_point_in_tet)// out of tet means no touch tet
	{
		
		if (!is_degenerated) {
			bool check = false;
			if (!is_point_in_tet) { // p out of tet,or on the border
				int r1, r2;
				//we need to know: if shoot on any edge?(two edges return -1, one edge return 1)
				
				r1 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
					bl.v[bl.facets[0][2]], true);
				r2 = ray_triangle_intersection(
					pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
					bl.v[bl.facets[1][2]], true);
				// idea is: if -1
				if (r1 == -1 || r2 == -1)
					return -1;
				if (r1 == 3 || r2 == 3) check = true;
				if (r1 == 2 && r2 == 0) check = true;
				if (r1 == 0 && r2 == 2) check = true;
				if (r1 == 1 && r2 == 1) return 0;
				if (r1 + r2 == 1) return 1;// 1-0 case
				if (r1 == 1 || r2 == 1) return 0;// 1-2 case
				if (r1 == 0 && r2 == 0) return 0;
				
				if (check == false) return 0;
				else {
					Vector3r norm0 = tri_norm(bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
						bl.v[bl.facets[0][2]]);
					Vector3r norm1 = tri_norm(bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
						bl.v[bl.facets[1][2]]);
					if (r1 == 3 || r2 == 3) {
						if (norm0.dot(dir) < 0 && norm1.dot(dir) < 0) return 1;
						else return 0;
					}
					if (r1 == 2) {
						if (norm0.dot(dir) < 0) return 1;
						else return 0;
					}
					if (r2 == 2) {
						if (norm1.dot(dir) < 0) return 1;
						else return 0;
					}
					std::cout << "impossible to go here " << std::endl;
					return 0;
				}
			}
			else { // p inside open tet 

				if (bl.phi_f[0] == 2) { // phi never calculated, need calculated
					get_tet_phi(bl);
				}
				Rational phip = phi(pt, bl.v);
				if (phip == 0)
					return 2;// point on bilinear
				return ray_correct_bilinear_face_pair_inter(pt, phip, dir, bl);
				
			}
		}
		else {// degenerated bilinear
			int degetype = bilinear_degeneration(bl);
			return ray_degenerated_bilinear_parity(bl, pt, dir, degetype);
		}
	}
	// -1 shoot another, 1 intersect, 2 point on triangle, 0 not intersect
	int ray_triangle_parity(
		const Vector3r& pt,
		const Vector3r& dir,
		const Vector3r& t0,
		const Vector3r& t1,
		const Vector3r& t2,
		const bool is_triangle_degenerated)
	{
		if (!is_triangle_degenerated) {
			return ray_triangle_intersection(pt, dir, t0, t1, t2, false);
			// 0 not hit, 1 hit on open triangle, -1 parallel or hit on edge, need
			// another shoot. 
		}
		else {
			// if pt on it (2), return 2; if 1(including overlap) return -1
			int i1 = ray_segment_intersection(pt, dir, t0, t1);// 2 means pt on the segment
			if (i1 == 2) return 2;
			if (i1 == 1) return -1;
			int i2 = ray_segment_intersection(pt, dir, t1, t2);
			if (i2 == 2) return 2;
			if (i2 == 1) return -1;
			/*int i3 = ray_segment_intersection(pt, dir, t2, t0);// if degenerated, check two edges is enough
			if (i3 == 2) return 2;
			if (i3 == 1) return -1;*/
			
			return 0; 
		}
	}




	// check if point has intersection with prism by counting parity
	// 1 -1 0
	int point_inside_prism(prism& psm, std::array<bilinear, 3> &bls, 
		const Vector3r& pt, const Vector3r& dir, const std::vector<bool>& is_pt_in_tet)
	{
		int S = 0;

		for (int patch = 0; patch < 3; ++patch) {

			int is_ray_patch = ray_bilinear_parity(
				bls[patch], pt, dir, bls[patch].is_degenerated, is_pt_in_tet[patch]);
			std::cout << "is_ray_patch " << is_ray_patch << std::endl;

			if (is_ray_patch == 2)
				return 1;

			if (is_ray_patch == -1)
				return -1;

			if (is_ray_patch == 1)
				S++;
		}

		int res;
		res = ray_triangle_parity(
			pt, dir, psm.p_vertices[0], psm.p_vertices[1], psm.p_vertices[2],
			psm.is_triangle_degenerated(0));

		if (res == 2)
			return 1;// it should be impossible 
		if (res == -1)
			return -1;

		if (res > 0)
			S++;
		res = ray_triangle_parity(
			pt, dir, psm.p_vertices[3], psm.p_vertices[4], psm.p_vertices[5],
			psm.is_triangle_degenerated(1));

		if (res == 2)
			return 1; // it should be impossible 
		if (res == -1)
			return -1;

		if (res > 0)
			S++;

		std::cout << "intersection nbrs " << S << std::endl;
		return ((S % 2) == 1) ? 1 : 0;
	}

	bool retrial_ccd(
		prism& psm, std::array<bilinear, 3>& bls,
		const Vector3r& pt,
		const std::vector<bool>& is_pt_in_tet)
	{
		static const int max_trials = 8;// TODO maybe dont need to set this
		Vector3r dir(1, 0, 0);

		int res = -1;
		int trials;
		for (trials = 0; trials < max_trials; ++trials) {
			res = point_inside_prism(psm, bls, pt, dir, is_pt_in_tet);

			if (res >= 0)
				break;

			Vector3d dir1 = Vector3d::Random();
			dir[0] = dir1[0];
			dir[1] = dir1[1];
			dir[2] = dir1[2];
		}

		if (trials == max_trials) {

			std::cout << "All rays are on edges, increase trials" << std::endl;
			throw "All rays are on edges, increase trials";
			return false;
		}

		return res >= 1; // >=1 means point inside of prism
	}

} // namespace ccd
