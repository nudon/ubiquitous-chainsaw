#include <limits>
#include <iostream>
#include "ChunkMatch.h"

using stk::StkFrames;
ChunkMatch::ChunkMatch(Chunk to_match) {
  orig = to_match;
}
//template< template <class T> class container>
void ChunkMatch::best_match_chunk(std::list<Chunk> &many_chunks, ChunkCompare comp) {
  double lowest = std::numeric_limits<double>::infinity();
  double temp = 0;
  Chunk best_match;
  for (Chunk b : many_chunks) {
    temp = comp.compare(orig, b);
    if (temp < lowest) {
      best_match = b;
      lowest = temp;
      if (lowest < 0.00001) {
	break;
      }
    }
  }
  match = best_match;
  score = lowest;
  match.make_chunk_filter();
}
int ChunkMatch::get_active_start() {
  return orig.get_time_center() - match.get_time_margin();
}

int ChunkMatch::get_active_end() {
  //std::cout << orig.get_time_center() + match.get_time_margin() << "\n";
  return orig.get_time_center() + match.get_time_margin();
}



bool ChunkMatch::comp_active_start(ChunkMatch &a, ChunkMatch &b) {
  return a.get_active_start() < b.get_active_start();
}

bool ChunkMatch::comp_active_end(ChunkMatch &a, ChunkMatch &b) {
  return a.get_active_end() < b.get_active_end();
}

bool ChunkMatch::comp_orig_start(ChunkMatch &a, ChunkMatch &b) {
  return a.orig.get_time_start() < b.orig.get_time_start();
}
bool ChunkMatch::comp_orig_end(ChunkMatch &a, ChunkMatch &b) {
  return a.orig.get_time_end() < b.orig.get_time_end();
}


bool  ChunkMatch:: operator == (ChunkMatch other) {
  return (get_orig_chunk() == other.get_orig_chunk() &&
	  get_match_chunk() == other.get_match_chunk());
}
