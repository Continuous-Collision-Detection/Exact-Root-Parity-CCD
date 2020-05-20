// Read collision data from a HDF5 file
#include "read_collision_data.hpp"

#include <array>
#include <string>
#include <vector>

#include <Eigen/Core>
//#include <boost/filesystem.hpp>
#include <highfive/H5Easy.hpp>

void read_vertex_face_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& vertex_face_data)
{
    vertex_face_data.clear();

    //for (auto& entry : boost::filesystem::directory_iterator(data_dir)) {
        H5Easy::File file(data_dir);

        HighFive::Group root = file.getGroup("/");

        const auto query_names = root.listObjectNames();
        for (const auto query_name : query_names) {
            vertex_face_data.push_back(
                H5Easy::load<Eigen::Matrix<double, 8, 3>>(file, query_name));
        }
    //}
}

void read_edge_edge_data(
    std::string data_dir,
    std::vector<Eigen::Matrix<double, 8, 3>>& edge_edge_data)
{
    edge_edge_data.clear();

    //for (auto& entry : boost::filesystem::directory_iterator(data_dir)) {
        //H5Easy::File file(entry.path().string());
    H5Easy::File file(data_dir);
        HighFive::Group root = file.getGroup("/");

        const auto query_names = root.listObjectNames();
        for (const auto query_name : query_names) {
            edge_edge_data.push_back(
                H5Easy::load<Eigen::Matrix<double, 8, 3>>(file, query_name));
        }
        std::cout<<"size of data "<<edge_edge_data.size()<<std::endl;
    //}
}
