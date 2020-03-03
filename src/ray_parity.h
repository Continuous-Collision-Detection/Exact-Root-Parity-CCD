#pragma once
#include <Utils.hpp>
#include <subfunctions.h>
namespace ccd {
Rational func_g(
    const Vector3r& x,
    const std::array<Vector3r, 4>& corners,
    const std::array<int, 3>& indices)
{
    const int p = indices[0];
    const int q = indices[1];
    const int r = indices[2];
    return (x - corners[p])
        .dot(cross(corners[q] - corners[p], corners[r] - corners[p]));
}
Rational phi(const Vector3r x, const std::array<Vector3r, 4>& corners)
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

void get_tet_phi(bilinear& bl)
{
    Vector3r p02 = (bl.v[0] + bl.v[2]) / 2;
    Rational phi02 = phi(p02, bl.v);
    if (phi02.get_sign() > 0) {
        bl.phi_f[0] = 1;
        bl.phi_f[1] = -1;
        return;
    }
    if (phi02.get_sign() < 0) {
        bl.phi_f[0] = -1;
        bl.phi_f[1] = 1;
        return;
    }
    std::cout << "!!can not happen, get tet phi" << std::endl;
}

int ray_bilinear_parity(
    bilinear& bl,
    const Vector3r& pt,
    const Vector3r& dir,
    const bool is_degenerated,
    const bool is_point_in_tet)
{
    bool inter1, inter2;
    if (!is_degenerated) {
        if (!is_point_in_tet) { // p out of tet
            if (ray_triangle_inter(
                    pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
                    bl.v[bl.facets[0][2]])
                    > 0
                || ray_triangle_inter(
                       pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
                       bl.v[bl.facets[1][2]])
                    > 0) { // choose one pair of facets
                return 1;
            }

            // TODO check half closed triangle shapes
        } else { // p inside tet

			if (bl.phi_f[0] == 2) {//phi never calculated, need calculated
                get_tet_phi(bl);
			}
                        Rational phip = phi(pt,bl.v);
                        if (phip==0)
            
        }
    } else {
    }
}
} // namespace ccd
