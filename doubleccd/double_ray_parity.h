#pragma once
//#include <Utils.hpp>
#include "double_subfunctions.h"
namespace ccd {

	// before going here, we already know the point can not be on the shape
	int ray_degenerated_bilinear_parity(
		const bilinear& bl,
		const Vector3r& pt,
		const Vector3r& dir,
		const int dege
	);
	// the point is inside of tet, tet is not degenerated
	//return: -1,1,0
	// phi_p is the phi of end point of ray
	int ray_correct_bilinear_face_pair_inter(
		const Vector3r& p,
		const Rational& phi_p,
		const Vector3r& dir,
		const bilinear& bl);

	int ray_bilinear_parity(
		bilinear& bl,
		const Vector3r& pt,
		const Vector3r& dir,
		const bool is_degenerated,
		const bool is_point_in_tet);// out of tet means no touch tet
		
	// -1 shoot another, 1 intersect, 2 point on triangle, 0 not intersect
	int ray_triangle_parity(
		const Vector3r& pt,
		const Vector3r& dir,
		const Vector3r& t0,
		const Vector3r& t1,
		const Vector3r& t2,
		const bool is_triangle_degenerated);




	// check if point has intersection with prism by counting parity
	// 1 -1 0
	int point_inside_prism(prism& psm, std::array<bilinear, 3> &bls,
		const Vector3r& pt, const Vector3r& dir, const std::vector<bool>& is_pt_in_tet);

	bool retrial_ccd(
		prism& psm, std::array<bilinear, 3>& bls,
		const Vector3r& pt,
		const std::vector<bool>& is_pt_in_tet);

} // namespace ccd
