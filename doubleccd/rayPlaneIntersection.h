#pragma once

#include <implicit_point.h>

bool rayIntersectsPlane(const explicitPoint3D& src, const explicitPoint3D& dir,
	const explicitPoint3D& v1, const explicitPoint3D& v2, const explicitPoint3D& v3);