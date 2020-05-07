#include <limits>
#include "ChunkGroupMatch.h"

using std::list;

ChunkGroupMatch::ChunkGroupMatch(ChunkGroup& to_match) {
  orig = to_match;
  score = std::numeric_limits<double>::infinity();
}

void ChunkGroupMatch::best_match(list<ChunkGroup>& pool, ChunkGroupCompare& comp) {
  double lowest = score;
  double temp = 0;
  ChunkGroup best_match = match;
  for (ChunkGroup b : pool) {
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
}

int ChunkGroupMatch::get_active_start() {
  return orig.get_super_chunk().get_time_center() - match.get_super_chunk().get_time_margin();
}

int ChunkGroupMatch::get_active_end() {
  return orig.get_super_chunk().get_time_center() + match.get_super_chunk().get_time_margin();
}
