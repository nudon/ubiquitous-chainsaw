#ifndef file_eigen_to_image_seen
#define file_eigen_to_image_seen

#include <Eigen/Dense>
#include <list>
#include "Chunk.h"
void write_eigen_to_file(std::string fn, Eigen::MatrixXf &mat);
void write_eigen_to_file(char* fn, Eigen::MatrixXf &mat);
void write_eigen_to_file(const char* fn, Eigen::MatrixXf &mat);

//for analysis of chunking process
void write_chunks_to_file(std::string fn, Eigen::MatrixXi &mat);
void write_chunks_to_file(const char* fn, Eigen::MatrixXi &mat);

void write_chunks_to_file(std::string fn, std::list<Chunk> &chunks, int rows, int cols);
void write_chunks_to_file(const char* fn, std::list<Chunk> &chunks, int rows, int cols);

#endif
