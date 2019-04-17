#include <FileRead.h>
#include <FileWvIn.h>
#include "Song.h"
#include "util.h"

using Eigen::MatrixXf;
using Eigen::MatrixXi;
using namespace stk;
Song::Song(std::string fn) {
  file_name = fn;
  FileRead song_file;
  FileWvIn song_wv;
  song_wv.openFile(fn);
  channels = song_wv.channelsOut();
  song_file.open(fn, channels);
  samples = song_file.fileSize();
  file_rate = song_file.fileRate();
  song_file.close();
  song_wv.closeFile();
}

MatrixXf Song::spectrogram(int samples_per_slice) {
  //just setup and apply fft to song
  FileWvIn input;
  input.openFile(get_file_name());
  //int nyquist = get_file_rate() / 2;
  int mono_sample_length = get_tot_samples() / get_channels();
  int total_slices = ceil(mono_sample_length / samples_per_slice);
  //MatrixXf zero (samples_per_slice, total_slices);
  //zero.setZero();
  //MatrixXf rand = MatrixXf::Random(samples_per_slice, total_slices);
  MatrixXf spec = TFD_extract(input,total_slices, samples_per_slice, get_channels());
  input.closeFile();
  return spec;    
}
