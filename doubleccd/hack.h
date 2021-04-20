//#include <doubleCCD/double_ray_parity.h>
#pragma once
#include <Eigen/Dense>
namespace doubleccd {
class hack {
public:
    Eigen::Vector3d dir;
    static hack& getInstance()
    {
        static hack instance; // Guaranteed to be destroyed.
                              // Instantiated on first use.
        return instance;
    }
};
} // namespace doubleccd
