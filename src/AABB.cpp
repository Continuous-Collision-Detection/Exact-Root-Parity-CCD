#include "AABB.h"

#include <cassert>

namespace ccd {

	/*
	bool is_triangle_cut_bounding_box(
		const Vector3d &tri0, const Vector3d &tri1, const Vector3d &tri2, int index)
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];
		Vector3 tmin, tmax;

		algorithms::get_tri_corners(tri0, tri1, tri2, tmin, tmax);
		bool cut = algorithms::box_box_intersection(tmin, tmax, bmin, bmax);
		if (cut == false) return false;

		if (cut) {

			std::array<Vector2, 3> tri;
			std::array<Vector2, 4> mp;
			int o0, o1, o2, o3, ori;
			for (int i = 0; i < 3; i++) {
				tri[0] = algorithms::to_2d(tri0, i);
				tri[1] = algorithms::to_2d(tri1, i);
				tri[2] = algorithms::to_2d(tri2, i);

				mp[0] = algorithms::to_2d(bmin, i);
				mp[1] = algorithms::to_2d(bmax, i);
				mp[2][0] = mp[0][0]; mp[2][1] = mp[1][1];
				mp[3][0] = mp[1][0]; mp[3][1] = mp[0][1];

				for (int j = 0; j < 3; j++) {
					o0 = fastEnvelope::Predicates::orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
					o1 = fastEnvelope::Predicates::orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
					o2 = fastEnvelope::Predicates::orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
					o3 = fastEnvelope::Predicates::orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
					ori = fastEnvelope::Predicates::orient_2d(tri[(j + 2) % 3], tri[j % 3], tri[(j + 1) % 3]);
					if (ori == 0) continue;
					if (ori*o0 <= 0 && ori*o1 <= 0 && ori*o2 <= 0 && ori*o3 <= 0) return false;
				}
			}
		}

		return cut;
	}
	bool is_point_cut_bounding_box(
		const Vector3 &p, int index)
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];
		if (p[0] < bmin[0] || p[1] < bmin[1] || p[2] < bmin[2]) return false;
		if (p[0] > bmax[0] || p[1] > bmax[1] || p[2] > bmax[2]) return false;
		return true;
	}
bool is_segment_cut_bounding_box(const Vector3 &seg0, const Vector3 &seg1, int
index)
{
        const auto &bmin = boxlist[index][0];
        const auto &bmax = boxlist[index][1];
        Scalar min[3], max[3];
        min[0] = std::min(seg0[0], seg1[0]);
        min[1] = std::min(seg0[1], seg1[1]);
        min[2] = std::min(seg0[2], seg1[2]);
        max[0] = std::max(seg0[0], seg1[0]);
        max[1] = std::max(seg0[1], seg1[1]);
        max[2] = std::max(seg0[2], seg1[2]);
        if (max[0] < bmin[0] || max[1] < bmin[1] || max[2] < bmin[2]) return
false; if (min[0] > bmax[0] || min[1] > bmax[1] || min[2] > bmax[2]) return
false; return true;
}
	*/
	

	
	bool is_bbd_cut_bounding_box(
		const Vector3 &bbd0, const Vector3 &bbd1, int index)
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];


		return algorithms::box_box_intersection(bbd0, bbd1, bmin, bmax);
	}
}
