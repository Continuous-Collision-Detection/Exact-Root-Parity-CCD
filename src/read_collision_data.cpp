// Read collision data from a HDF5 file

#include <array>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <boost/filesystem.hpp>
#include <highfive/H5Easy.hpp>

#include <ccd.hpp>

void read_vertex_face_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& vertex_face_data)
{
    vertex_face_data.clear();

    for (auto& entry : boost::filesystem::directory_iterator(data_dir)) {
        H5Easy::File file(entry.path().string());

        HighFive::Group root = file.getGroup("/");

        const auto query_names = root.listObjectNames();
        for (const auto query_name : query_names) {
            vertex_face_data.push_back(
                H5Easy::load<Eigen::Matrix<double, 8, 3>>(file, query_name));
        }
    }
}

void read_edge_edge_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& edge_edge_data)
{
    edge_edge_data.clear();

    for (auto& entry : boost::filesystem::directory_iterator(data_dir)) {
        H5Easy::File file(entry.path().string());

        HighFive::Group root = file.getGroup("/");

        const auto query_names = root.listObjectNames();
        for (const auto query_name : query_names) {
            edge_edge_data.push_back(
                H5Easy::load<Eigen::Matrix<double, 8, 3>>(file, query_name));
        }
    }
}
