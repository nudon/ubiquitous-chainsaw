#ifndef file_eigen_to_image_seen
#define file_eigen_to_image_seen

#include <Eigen/Dense>

using Eigen::MatrixXf;
void write_eigen_to_file(char* fn, MatrixXf mat);
void write_eigen_to_file(const char* fn, MatrixXf mat);
#endif
