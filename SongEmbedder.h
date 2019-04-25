#ifndef FILE_SONGEMBEDDER_SEEN
#define FILE_SONGEMBEDDER_SEEN
#include <list>
#include "Song.h"
#include "ChunkMatch.h"

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

  std::list<ChunkMatch> get_matches(std::list<Chunk> &reciptor_chunks, std::list<Chunk> &replacer_chunks);

  
  void output_remix(int sps,std::list<ChunkMatch> &matches, std::string fn);
  
};

#endif
