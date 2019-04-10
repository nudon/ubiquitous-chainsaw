#ifndef FILE_SONGEMBEDDER_SEEN
#define FILE_SONGEMBEDDER_SEEN
#include <list>
#include "Song.h"

class SongEmbedder{
 public:
  SongEmbedder(){}
  SongEmbedder(Song a, Song b);
  ~SongEmbedder(){}
  
  void funk();
  
 private:
  Eigen::MatrixXf spectrogram(int samples_per_slice);
  Song reciptor;
  Song replacer;
  
};

#endif
