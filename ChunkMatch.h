#ifndef FILE_CHUNKMATCH_SEEN
#define FILE_CHUNKMATCH_SEEN

#include <list>
#include "Chunk.h"
#include "ChunkCompare.h"
//for chunk matching, have 2 ideas
//iterate over entire list, find absolute best match
//or, create multiple linear collections, sorted by some property
//main ones I'm thinking are frequency center, frequency margin, and time_margin
//find closest match, and return a list of match and n adjacent chunks
//join all subsets, find maximum match over those
//should still implement full matching for comparison, but other method might also work
//and be quicker
//yeah full match takes forever even on smaller (1:30 min) songs

class ChunkMatch {
 public :
  ChunkMatch() {}
  ChunkMatch(Chunk to_match);
  ~ChunkMatch() {}
  
  //template< template <class T> class container>
  void best_match_chunk(std::list<Chunk> many_chunks, ChunkCompare comp);
  int get_active_start();
  int get_active_end();

  Chunk get_orig_chunk() { return orig; }
  Chunk get_match_chunk() { return match; }

  static bool comp_active_start(ChunkMatch a, ChunkMatch b);
  static bool comp_active_end(ChunkMatch a, ChunkMatch b);

  bool operator == (ChunkMatch other);
  
 private:
  Chunk orig;
  Chunk match;
  double score;

};


#endif
