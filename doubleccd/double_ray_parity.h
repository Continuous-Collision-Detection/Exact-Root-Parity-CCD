#pragma once
//#include <Utils.hpp>
#include <doubleCCD/double_subfunctions.h>

namespace doubleccd {

// before going here, we already know the point can not be on the shape
int ray_degenerated_bilinear_parity(
    const bilinear& bl,
    const Vector3d& pt,
    const Vector3d& pt1,
    const Vector3d& dir,
    const int dege);
// the point is inside of tet, tet is not degenerated
// return: -1,1,0
// phi_p is the phi of end point of ray
int ray_correct_bilinear_face_pair_inter(
    const Vector3d& p,
    const Vector3d& p1,
    const Rational& phi_p,
    const Vector3d& dir,
    const bilinear& bl);

int ray_bilinear_parity(
    bilinear& bl,
    const Vector3d& pt,
    const Vector3d& pt1,
    const Vector3d& dir,
    const bool is_degenerated,
    const bool is_point_in_tet); // out of tet means no touch tet

// -1 shoot another, 1 intersect, 2 point on triangle, 0 not intersect
int ray_triangle_parity(
    const Vector3d& pt,
    const Vector3d& pt1,
    const Vector3d& dir,
    const Vector3d& t0,
    const Vector3d& t1,
    const Vector3d& t2,
    const bool is_triangle_degenerated);

// check if point has intersection with prism by counting parity
// 1 -1 0
int point_inside_prism(
    prism& psm,
    std::array<bilinear, 3>& bls,
    const Vector3d& pt,
    const Vector3d& pt1,
    const Vector3d& dir,
    const std::vector<bool>& is_pt_in_tet);

bool retrial_ccd(
    prism& psm,
    std::array<bilinear, 3>& bls,
    const Vector3d& pt,
    const std::vector<bool>& is_pt_in_tet);
bool retrial_ccd_hex(
	std::array<bilinear, 6>& bls,
	const Vector3d& pt,
	const std::vector<bool>& is_pt_in_tet);
bool shoot_origin_ray_hex( std::array<bilinear, 6>& bls);
bool shoot_origin_ray_prism(prism& psm, std::array<bilinear, 3>& bls);
void ray_time();

} // namespace doubleccd
