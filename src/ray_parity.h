#pragma once
#include <Utils.hpp>
#include <subfunctions.h>
namespace ccd {

int ray_bilinear_parity(
    bilinear& bl,
    const Vector3r& pt,
    const Vector3r& dir,
    const bool is_degenerated,
    const bool is_point_in_tet)
{
    if (!is_degenerated) {
        if (!is_point_in_tet) {
            ray_triangle_inter// TODO check half closed triangle shapes
		} else {
        
		}
	} else {
    }
}
}
