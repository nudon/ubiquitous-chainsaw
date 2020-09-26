#ifndef FILE_SONGEMBEDDER_SEEN
#define FILE_SONGEMBEDDER_SEEN
#include <list>
#include <Stk.h>
#include <FileWvIn.h>
#include <LentPitShift.h>
#include "Song.h"
#include "ChunkMatch.h"
#include "SoundPatch.h"


class SongEmbedder{
 public:
  SongEmbedder(){}
  SongEmbedder(Song a, Song b, int sample_size, int edge_snazr, int chunk_snazr, int match_snazr);
  ~SongEmbedder(){}
  

  void funk();
  void extract(Song &s);
  void output_remix(std::list<ChunkMatch> &matches, std::string fn);
  
 private:
  Song reciptor;
  Song replacer;
  int samples_per_slice;
  int edge_snazr;
  int chunk_snazr;
  int match_snazr;

  std::list<Chunk> make_chunks(Song input);
  //Eigen::MatrixXf spectrogram(int samples_per_slice);
  
  std::list<SoundPatch> individual_sound_patches(std::list<Chunk>& res, std::list<Chunk>& src);
  std::list<SoundPatch> group_sound_patches(std::list<Chunk>& res, std::list<Chunk>& src);

  std::list<ChunkMatch> get_matches(std::list<Chunk> &reciptor_chunks, std::list<Chunk> &replacer_chunks);
  std::list<ChunkGroupMatch> get_matches(std::list<ChunkGroup> &res, std::list<ChunkGroup>& src);

  void output_remix(std::list<SoundPatch>& patches, std::string fn);

  void normal_filt(ChunkMatch& chunk_match, stk::FileWvIn &sample_src, stk::StkFrames &frame_in, stk::StkFrames &out, int i);
  void alt_filt(ChunkMatch& chunk_match, stk::FileWvIn &sample_src, stk::StkFrames &out, int i, stk::LentPitShift &lent);


  
  
};

#endif
