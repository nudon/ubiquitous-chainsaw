#include <limits>
#include <iostream>
#include "ChunkMatch.h"

using stk::StkFrames;
ChunkMatch::ChunkMatch(Chunk to_match) {
  orig = to_match;
  score = std::numeric_limits<double>::infinity();
}
//template< template <class T> class container>
void ChunkMatch::best_match_chunk(std::list<Chunk> &many_chunks, ChunkCompare &comp) {
  double lowest = score;
  double temp = 0;
  Chunk best_match = match;
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
  match.make_chunk_filter();
  score = lowest;
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




int ChunkMatch::match_hash( ChunkMatch &arg) {
  std::hash<double> hashy;
  double n;
  n = (arg.orig.get_time_start() - arg.match.get_time_end()) * (arg.orig.get_time_end() - arg.match.get_time_start());
  //d = arg.orig.get_rel_freq_center() / arg.match.get_rel_freq_center();
  return hashy(n);
}

bool  ChunkMatch:: operator == (const ChunkMatch other) const  {
  return (orig == other.orig && match == other.match);
}

  
