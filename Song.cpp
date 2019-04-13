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
  //MatrixXf noread = noread_TFD_extract(total_slices, samples_per_slice, get_channels());
  //okay figured out the error. could have guessed it from valgrind logs but wan't smart enough
  //so, valgrind had logs about errors in closing the input
  //somewhere it looked like it had previously freed it, and was doing it twice
  //error was handing in FileWvIn to tfd_extract. because I handed in whole structure
  //it made a shallow copy, probably also copying over some file handle pointer or something
  //when it reached the end of tfd_extract, it freed the local copy
  //then also went to close/free some shared resource
  //then in this function, it double freed.
  //which resulted in lots of weird stuff like list iterators suddenly getting clobbered and calls to stk addtime and tick to blow up
  //so handing in a filename and having function open it up
  //as well as passing by reference, made the error go away
  //well along the way I read that handing thing by reference(&arg) in functions was a good practice anyway
  //MatrixXf fn = filename_TFD_extract(get_file_name(), total_slices, samples_per_slice, get_channels());
  MatrixXf ret = spec;
  input.closeFile();
  return ret;    
}
