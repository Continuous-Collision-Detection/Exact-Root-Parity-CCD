// Read collision data from a HDF5 file

#include <array>
#include <string>
#include <vector>

#include <Eigen/Core>
//#include <boost/filesystem.hpp>
#include <highfive/H5Easy.hpp>


/// Example: read_vertex_face_data("path/to/vf-data/", vertex_face_data);
void read_vertex_face_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& vertex_face_data);

void read_edge_edge_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& edge_edge_data);
