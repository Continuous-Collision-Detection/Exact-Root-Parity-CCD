#include <catch2/catch.hpp>

#include <ccd.hpp>
#include <doubleccd.hpp>

static const double EPSILON = std::numeric_limits<float>::epsilon();
static const double MIN_DIST = 0.1 * EPSILON;
TEST_CASE(
    "Test Point-Triangle Continuous Collision Detection",
    "[ccd][point-triangle]")
{
    // point
    double v0z = GENERATE(0.0, -1.0);
    Eigen::Vector3d v0(0, 1, v0z);
    // triangle = (v1, v2, v3)
    Eigen::Vector3d v1(-1, 0, 1);
    Eigen::Vector3d v2(1, 0, 1);
    Eigen::Vector3d v3(0, 0, -1);

    // displacements
    double u0y
        = -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
    Eigen::Vector3d u0(0, u0y, u0z);
    double u1y
        = GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    Eigen::Vector3d u1(0, u1y, 0);

    SECTION("Clockwise triangle")
    {
        // already in clockwise order
    }
    SECTION("Counter-clockwise triangle") { std::swap(v1, v2); }

    bool hit = ccd::vertexFaceCCD(
        v0, v1, v2, v3, v0 + u0, v1 + u1, v2 + u1, v3 + u1, MIN_DIST);

    CAPTURE(v0z, u0y, u1y, u0z, EPSILON);
    CHECK(hit == ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z())));

    hit = doubleccd::vertexFaceCCD(
        v0, v1, v2, v3, v0 + u0, v1 + u1, v2 + u1, v3 + u1, MIN_DIST);

    CAPTURE(v0z, u0y, u1y, u0z, EPSILON);
    CHECK(hit == ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z())));
}

TEST_CASE("Test Edge-Edge Continuous Collision Detection", "[ccd][edge-edge]")
{
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    double e1x = GENERATE(
        -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
        1 + EPSILON);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement
        = GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
    Eigen::Vector3d u0(0, y_displacement, 0);
    Eigen::Vector3d u1(0, -y_displacement, 0);

    try {
        bool hit = ccd::edgeEdgeCCD(
            v0, v1, v2, v3, v0 + u0, v1 + u0, v2 + u1, v3 + u1, MIN_DIST);
        CAPTURE(y_displacement, e1x);
        CHECK(hit == (y_displacement >= 1.0 && e1x >= -1 && e1x <= 1));
    } catch (...) {
    }

    bool hit = doubleccd::edgeEdgeCCD(
        v0, v1, v2, v3, v0 + u0, v1 + u0, v2 + u1, v3 + u1, MIN_DIST);
    CAPTURE(y_displacement, e1x);
    CHECK(hit == (y_displacement >= 1.0 && e1x >= -1 && e1x <= 1));
}

TEST_CASE("Zhongshi test case", "[ccd][point-triangle]")
{
    double qy = GENERATE(-EPSILON, 0, EPSILON);

    Eigen::Vector3d q;
    q << 0, qy, 0;

    Eigen::Vector3d b0;
    b0 << 0, 0, 0;
    Eigen::Vector3d b1;
    b1 << 0, 1, 0;
    Eigen::Vector3d b2;
    b2 << 1, 0, 0;

    Eigen::Vector3d t0;
    t0 << 0, 0, 1;
    Eigen::Vector3d t1;
    t1 << 0, 1, 1;
    Eigen::Vector3d t2;
    t2 << 1, 0, 1;

    Eigen::Vector3d q1;
    q1 << 0, qy, 0;

    bool hit = ccd::vertexFaceCCD(q, b0, b1, b2, q1, t0, t1, t2, MIN_DIST);
    CAPTURE(qy);
    CHECK(hit == q.y() >= 0);

    hit = ccd::vertexFaceCCD(q, b0, b1, b2, q1, t0, t1, t2, MIN_DIST);
    CAPTURE(qy);
    CHECK(hit == q.y() >= 0);
}
