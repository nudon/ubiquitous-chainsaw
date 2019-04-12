#ifndef FILE_UTIL_SEEN
#define FILE_UTIL_SEEN

#include <stdio.h>
#include <FileLoop.h>
#include <FileWvOut.h>
#include <Fir.h>
#include <string.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>


//thing with fourier transform
Eigen::MatrixXf TFD_extract(stk::FileWvIn &input, int total_slices, int samples_per_slice, int channels);
Eigen::MatrixXf noread_TFD_extract(int total_slices, int samples_per_slice, int channels);
Eigen::MatrixXf filename_TFD_extract(std::string fn, int total_slices, int samples_per_slice, int channels);

//filtering methods
//using stk
void make_fir_bandpass_filter(double* coefs, int taps, int freq_center, int freq_margin,  int max_freq);
stk::Fir create_fir_from_coefs(double* coefs, int len);
stk::Fir create_1d_gaussian_filter(int length, double amplitude, double center, double stddev);
//using eigen

Eigen::MatrixXf create_1d_gaussian_filter_col(int length, double amplitude, double stddev);
Eigen::MatrixXf create_1d_gaussian_filter_row(int length, double amplitude, double stddev);
Eigen::MatrixXf one_d_convolve(Eigen::MatrixXf mat, Eigen::MatrixXf kern);
Eigen::MatrixXf dog(Eigen::MatrixXf full, Eigen::MatrixXf gauss);

//chunking functions
double snaz(Eigen::MatrixXf filt, int snazr);

#endif
