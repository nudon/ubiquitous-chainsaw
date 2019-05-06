#include <iostream>
#include <string.h>
#include <list>
#include <tuple>
#include <functional>
#include "eigen_to_image.h"
#include "Chunk.h"

#ifdef PICTURE
#include "pngwriter.h"
#endif

using std::tuple;
using Eigen::MatrixXi;
using Eigen::MatrixXf;

//assigns chunk id a random color
tuple<double, double ,double> hash_id(int id) {
  double r, g, b;
  std::hash<double> ihash;
  int min = 50;
  int max = 255;
  r = ihash(id) % max;
  g = ihash(id + 1) % max;
  b = ihash(id + 2) % max;
  int max_loop = 6;
  int count = 0;
  //don't want chunks to be too close to background color
  if (id <= 0) {
    return std::make_tuple(0.0,0.0,0.0);  
  }
  while(r < min && g < min &&  b < min) {
    if (count < max_loop) {
      r = ihash(r) % max;
      g = ihash(g) % max;
      b = ihash(b) % max;
      
    }
    else {
      //maybe the hash reached a fixed point, just set color to white
      r = max;
      g = max;
      b = max;
    }
    count++;
  }
  double d_max = max;
  return std::make_tuple(r / d_max ,g / d_max ,b / d_max);  
}


void write_chunks_to_file(std::string fn, Eigen::MatrixXi &mat) {
  write_chunks_to_file(fn.c_str(), mat);
}

void write_chunks_to_file(const char* fn, MatrixXi &mat) {
  #ifdef PICTURE
  ///std::cout << fn << "yoyoyo\n";
  int rows = mat.rows();
  int cols = mat.cols();
  int v = 0;
  pngwriter img(cols, rows, 0, fn);
  for (int row_i = 0; row_i < rows; row_i++) {
    for (int col_i = 0; col_i < cols; col_i++) {
      v = mat(row_i, col_i);
      auto rgb = hash_id(v);
      double r = std::get<0>(rgb);
      double g = std::get<1>(rgb);
      double b = std::get<2>(rgb);
      //std::cout << r << g << b << " \n";
      img.plot(col_i + 1, row_i + 1, r,g,b);
    }
  }
  img.close(); 
  #endif
}

void write_chunks_to_file(std::string fn, std::list<Chunk> &chunks, int rows, int cols) {
  write_chunks_to_file(fn.c_str(), chunks, rows, cols);
}

void write_chunks_to_file(const char* fn, std::list<Chunk> &chunks, int rows, int cols) {
  //convert chunks to a matrixXi
  MatrixXi mat(rows, cols);
  mat.setZero();
  int val;
  int x, y, xl, yl;
  for (Chunk a_chunk : chunks) {
    val = a_chunk.get_chunk_id();
    x = a_chunk.get_time_center() - a_chunk.get_time_margin();
    y = a_chunk.get_freq_center() - a_chunk.get_freq_margin();
    xl = a_chunk.get_time_margin() * 2;
    yl = a_chunk.get_freq_margin() * 2;
    mat.block(y,x, yl, xl) = MatrixXi::Constant(yl, xl, val);
  }
  write_chunks_to_file(fn, mat);
}


void write_eigen_to_file(std::string fn, Eigen::MatrixXf &mat) {
  write_eigen_to_file(fn.c_str(), mat);
}

//using Eigen::MatrixXf;
void write_eigen_to_file(const char* fn, MatrixXf &mat) {
  char* temp = strdup(fn);
  write_eigen_to_file(temp, mat);
  free(temp);
}


void write_eigen_to_file(char* fn, MatrixXf &mat) {
  #ifdef PICTURE
  int rows = mat.rows();
  int cols = mat.cols();
  double d = 0;
  float v = 0;
  int p = 0;
  pngwriter img(cols, rows, 0, fn);
  for (int row_i = 0; row_i < rows; row_i++) {
    for (int col_i = 0; col_i < cols; col_i++) {
      d = mat(row_i, col_i);
      v = mat(row_i, col_i);
      if (v < d || v < 0) {
	std::cout << v << " " << d << "\n";
      }
      img.plot(col_i + 1, row_i + 1, v,v,v);
    }
  }
  img.close(); 
  #endif
}



