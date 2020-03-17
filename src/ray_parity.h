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
// before going here, we already know the point can not be on the shape
bool ray_degenerated_bilinear_parity(
    const bilinear& bl,
    const Vector3r& pt,
    const Vector3r& dir,
    const int dege
    )//TODO consider if it is correct when one of triangle is degenerated as a segment
{
    bool res;
    int r1, r2;
    if (dege == BI_DEGE_PLANE) {
        r1 = ray_halfopen_triangle_inter(// TODO need to check
            pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
            bl.v[bl.facets[0][2]]);
        if (r1 == 2)
            return 2;
        if (r1 == -1)
            return -1;
        if (r1 == 1)
            return 1;
        r2 = ray_halfopen_triangle_inter(
            pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
            bl.v[bl.facets[1][2]]);
        if (r2 == 2)
            return 2;
        if (r2 == -1)
            return -1;
        if (r2 == 1)
            return 1;
        return 0;
        
    } else {

        if (dege == BI_DEGE_XOR_02) { // triangle 0-1-2 and 0-2-3
            r1 = ray_halfopen_triangle_inter( 
                pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
                bl.v[bl.facets[0][2]]);//0: not hit, 1: hit on open triangle, 2: pt on halfopen T, since already checked, accept it 
            r2 = ray_halfopen_triangle_inter(
                pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
                bl.v[bl.facets[1][2]]);
            return int_XOR(r1, r2);
        }

        if (dege == BI_DEGE_XOR_13) { // triangle 0-1-3 and 3-1-2
            r1 = ray_halfopen_triangle_inter( 
                pt, dir, bl.v[bl.facets[2][0]], bl.v[bl.facets[2][1]],
                bl.v[bl.facets[2][2]]);//0: not hit, 1: hit on open triangle, 2: pt on halfopen T, since already checked, accept it 
            r2 = ray_halfopen_triangle_inter(
                pt, dir, bl.v[bl.facets[3][0]], bl.v[bl.facets[3][1]],
                bl.v[bl.facets[3][2]]);
            return int_XOR(r1, r2);
        }
    }
    std::cout << "!! THIS CANNOT HAPPEN" << std::endl;
    return false;
}
// the point is inside of tet, tet is not degenerated
//return: -1,1,0
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
        r1 = ray_halfopen_triangle_inter(
            p, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
            bl.v[bl.facets[0][2]]);
        r2 = ray_halfopen_triangle_inter(
            p, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
            bl.v[bl.facets[1][2]]);
        if (r1 == -1 || r2 == -1)
            return -1;
        if (r1 > 0 || r2 > 0)
            return 1; // cannot be degenerated, so this can work
        return 0;
    } else {
        r1 = ray_halfopen_triangle_inter(
            p, dir, bl.v[bl.facets[2][0]], bl.v[bl.facets[2][1]],
            bl.v[bl.facets[2][2]]);
        r2 = ray_halfopen_triangle_inter(
            p, dir, bl.v[bl.facets[3][0]], bl.v[bl.facets[3][1]],
            bl.v[bl.facets[3][2]]);
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
    bool inter1, inter2;
    if (!is_degenerated) {
        if (!is_point_in_tet) { // p out of tet
            int r1, r2;
            r1 = ray_halfopen_triangle_inter(
                pt, dir, bl.v[bl.facets[0][0]], bl.v[bl.facets[0][1]],
                bl.v[bl.facets[0][2]]);
            r2 = ray_halfopen_triangle_inter(
                pt, dir, bl.v[bl.facets[1][0]], bl.v[bl.facets[1][1]],
                bl.v[bl.facets[1][2]]);
            if (r1 == -1 || r2 == -1)
                return -1;
            if (r1 > 0 || r2 > 0)// since no touch tet, if ri cannot be 2 
                return 1;
            return 0;

            // TODO check half closed triangle shapes
        } else { // p inside tet

            if (bl.phi_f[0] == 2) { // phi never calculated, need calculated
                get_tet_phi(bl);
            }
            Rational phip = phi(pt, bl.v);
            if (phip == 0)
                return 2;// point on bilinear
            int res = ray_correct_bilinear_face_pair_inter(pt, phip, dir, bl);
            if (res == 1)
                return 1;
            if (res == -1)
                return -1;
            if (res == 0)
                return 0;
        }
    } else {// degenerated bilinear
        int degetype = bilinear_degeneration(bl);
        return ray_degenerated_bilinear_parity(bl, pt, dir, degetype);// TODO this can not fix the case one triangle totally segment
    }
}
int ray_triangle_parity(
    const Vector3r& pt,
    const Vector3r& dir,
    const Vector3r& t0,
    const Vector3r& t1,
    const Vector3r& t2,
    const bool is_triangle_degenerated)
{
    if (!is_triangle_degenerated) {
        return ray_open_triangle_inter(pt, dir, t0, t1, t2);
        // 0 not hit, 1 hit on open triangle, -1 parallel or hit on edge, need
        // another shoot. TODO with -1 need to cleaned
    } else {
        int i1 = ray_segment_inter(pt, dir, t0, t1);
        if (i1 != 0)
            return i1;
        int i2 = ray_segment_inter(pt, dir, t1, t2);
        if (i2 != 0)
            return i2;
        int i3 = ray_segment_inter(pt, dir, t2, t0);

        return i3; // TODO maybe we dont need return 1 because either ignor or
                   // shoot another ray
    }
}




    // check if point has intersection with prism by counting parity
int point_inside_prism( prism& psm, const Vector3r& pt, const Vector3r& dir, const std::vector<bool>& is_pt_in_tet)
{
    int S = 0;
    bool point_in_tet;
    bool is_degenerate;


    for (int patch = 0; patch < psm.bilinears.size(); ++patch) {
        bilinear tempbl(
            psm.bilinears[patch][0], psm.bilinears[patch][1],
            psm.bilinears[patch][2], psm.bilinears[patch][3]);

        int is_ray_patch = ray_bilinear_parity(
            tempbl, pt, dir, tempbl.is_degenerated, is_pt_in_tet[patch]);
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
        return 1;// it should be impossible TODO check if it is checked before
    if (res == -1)
        return -1;

    if (res > 0)
        S++;
    res = ray_triangle_parity(
        pt, dir, psm.p_vertices[3], psm.p_vertices[4], psm.p_vertices[5],
        psm.is_triangle_degenerated(1));

    if (res == 2)
        return 1; // it should be impossible TODO check if it is checked before
    if (res == -1)
        return -1;

    if (res > 0)
        S++;
    
    std::cout << "intersection nbrs " << S << std::endl;
    return ((S % 2) == 1) ? 1 : 0;
}
 
bool retrial_ccd(
    prism& psm,
    const Vector3r& pt,
    const std::vector<bool>& is_pt_in_tet)
{
    static const int max_trials = 8;// TODO maybe dont need to set this
    Vector3r dir(1, 0, 0);

    int res = -1;
    int trials;
    for (trials = 0; trials < max_trials; ++trials) {
        res = point_inside_prism(psm, pt, dir, is_pt_in_tet);

        if (res >= 0)
            break;

        dir = Vector3d::Random();
    }

    if (trials == max_trials) {

        std::cout << "All rays are on edges, increase trials" << std::endl;
        throw "All rays are on edges, increase trials";
        return false;
    }

    return res >= 1; // what is this mean?
}

} // namespace ccd