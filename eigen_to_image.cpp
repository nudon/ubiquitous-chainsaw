#include <iostream>
#include <string.h>
#include "eigen_to_image.h"

#ifdef PICTURE
#include "pngwriter.h"
#endif




//using Eigen::MatrixXf;
void write_eigen_to_file(const char* fn, MatrixXf mat) {
  char* temp = strdup(fn);
  write_eigen_to_file(temp, mat);
  free(temp);
}


void write_eigen_to_file(char* fn, MatrixXf mat) {
  #ifdef PICTURE
  int rows = mat.rows();
  int cols = mat.cols();
  float v = 0;
  int p = 0;
  pngwriter img(cols, rows, 0, fn);
  for (int row_i = 0; row_i < rows; row_i++) {
    for (int col_i = 0; col_i < cols; col_i++) {
      v = mat(row_i, col_i);
      //maybe do scaling, or other things if you want spectogram to be a color
      p = v* 80;
      img.plot(col_i + 1, row_i + 1, p,p,p);
    }
  }
  img.close(); 
  #endif
}



