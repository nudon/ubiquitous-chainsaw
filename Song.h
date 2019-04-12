#ifndef FILE_SONG_SEEN
#define FILE_SONG_SEEN

//#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <list>
#include "Chunk.h"


class Song {
 public:
  Song(){}
  Song(std::string fn);
  ~Song(){}
  
  int get_channels() { return channels; }
  int get_tot_samples() { return samples; }
  int get_samples_per_channel() { return get_tot_samples() / get_channels(); }
  int get_file_rate() { return file_rate; }
  const std::string get_file_name() { return (const std::string)file_name;}
  Eigen::MatrixXf spectrogram(int samples_per_slice);
  
  
 private:
  std::string file_name;
  int channels;
  int samples;
  int file_rate;
};

#endif
