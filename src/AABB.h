#pragma once
#include "Utils.hpp"

#include <vector>
#include <array>


namespace ccd {

		//bool is_triangle_cut_bounding_box(const Vector3d &tri0, const Vector3d &tri1, const Vector3d &tri2, int index);
		//bool is_point_cut_bounding_box(const Vector3d &p,  int index);
		//bool is_segment_cut_bounding_box(const Vector3d &seg0, const Vector3d &seg1, int index);
		bool is_bbd_cut_bounding_box(const Vector3d &bbd0, const Vector3d &bbd1, int index);
	
}