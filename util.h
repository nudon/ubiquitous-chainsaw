#ifndef FILE_UTIL_SEEN
#define FILE_UTIL_SEEN

#include <stdio.h>
#include <list>
#include <FileLoop.h>
#include <FileWvOut.h>
#include <Fir.h>
#include <string.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "Chunk.h"
#include "ChunkMatch.h"

//fourier transform
Eigen::MatrixXf TFD_extract(stk::FileWvIn &input, int total_slices, int samples_per_slice, int channels);
void kiss_to_stk(kiss_fft_scalar* in, stk::StkFrames &out);
void stk_to_kiss(stk::StkFrames &in, kiss_fft_scalar* out);

//filtering methods
//using stk
void make_fir_bandpass_filter(double* coefs, int taps, int freq_center, int freq_margin,  int max_freq);
stk::Fir create_fir_from_coefs(double* coefs, int len);
stk::Fir create_1d_gaussian_filter(int length, double amplitude, double center, double stddev);
//using eigen

Eigen::MatrixXf create_1d_gaussian_filter_col(int length, double amplitude, double stddev);
Eigen::MatrixXf create_1d_gaussian_filter_row(int length, double amplitude, double stddev);
Eigen::MatrixXf one_d_convolve(Eigen::MatrixXf &mat, Eigen::MatrixXf &kern);
Eigen::MatrixXf dog(Eigen::MatrixXf &full, Eigen::MatrixXf &gauss);

//chunking functions
double snaz(Eigen::MatrixXf &filt, int snazr);
//double snaz(std::list<Chunk> &many_chunks, int snazr);
//generic sub-non-zero average zeroing
//though it's not really zeroing because you are removing things from a list
template<class T>
double gsnaz(std::list<T> &many_t, double (*func)(T), int snazr) {
  double nz_tot = 0;
  double nz_avg = 0;
  double val = 0;
  int nz_count = 0;
  std::list<T> temp1 (many_t);
  std::list<T> temp2;
  for (int i = 0; i <= snazr; i++) {
    nz_tot = 0;
    nz_count = 0;
    for (T a_t : temp1) {
      val = func(a_t);
      if (val > nz_avg) {
	nz_tot += val;
	nz_count++;
	temp2.push_front(a_t);
      }
    }
    if (i < snazr) {
      nz_avg = nz_tot / nz_count;
      temp1.assign(temp2.begin(), temp2.end());
      temp2.erase(temp2.begin(), temp2.end());
    }
  }
  many_t.assign(temp2.begin(), temp2.end());
  return nz_avg;
}

//okay so I have a case where good things are lower(match scores)
//similar idea to sub-non-zero average zeroing
//but instead it's anything above the non-zero average gets zerod/removed
//super non-zero average zeroing? same acronym as sub variant
//just add a reverse in the name
template<class T>
double reverse_gsnaz(std::list<T> &many_t, double (*func)(T), int snazr) {
  double nz_tot = 0;
  double nz_avg = std::numeric_limits<double>::infinity();
  double val = 0;
  int nz_count = 0;
  std::list<T> temp1 (many_t);
  std::list<T> temp2;
  for (int i = 0; i <= snazr; i++) {
    nz_tot = 0;
    nz_count = 0;
    for (T a_t : temp1) {
      val = func(a_t);
      if (val <= nz_avg) {
	nz_tot += val;
	nz_count++;
	temp2.push_front(a_t);
      }
    }
    if (i < snazr) {
      nz_avg = nz_tot / nz_count;
      temp1.assign(temp2.begin(), temp2.end());
      temp2.erase(temp2.begin(), temp2.end());
    }
  }
  many_t.assign(temp2.begin(), temp2.end());
  return nz_avg;
}
#endif
