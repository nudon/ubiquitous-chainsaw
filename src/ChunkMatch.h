#ifndef FILE_CHUNKMATCH_SEEN
#define FILE_CHUNKMATCH_SEEN

#include <list>
#include "Chunk.h"
#include "ChunkCompare.h"

class ChunkMatch {
 public :
  ChunkMatch() {}
  ChunkMatch(Chunk to_match);
  ~ChunkMatch() {}
  
  //template< template <class T> class container>
  void best_match_chunk(std::list<Chunk> &many_chunks, ChunkCompare &comp);
  int get_active_start();
  int get_active_end();
  double get_score() { return score; }

  Chunk get_orig_chunk() { return orig; }
  Chunk get_match_chunk() { return match; }

  static std::list<ChunkMatch> self_match(std::list<Chunk>&rec);
  static std::list<ChunkMatch> best_match(std::list<Chunk>&rec, std::list<Chunk> &rep, ChunkCompare& comp);
  static std::list<ChunkMatch> quick_match(std::list<Chunk>&rec, std::list<Chunk> &rep, ChunkCompare& comp);
 
  static bool comp_active_start(ChunkMatch &a, ChunkMatch &b);
  static bool comp_active_end(ChunkMatch &a, ChunkMatch &b);

  static bool comp_orig_start(ChunkMatch &a, ChunkMatch &b);
  static bool comp_orig_end(ChunkMatch &a, ChunkMatch &b);

  static int match_hash( ChunkMatch &arg);

  

  bool operator == (const ChunkMatch other) const;
  
 private:
  Chunk orig;
  Chunk match;
  double score;
};

#endif
