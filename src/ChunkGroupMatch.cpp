#include <limits>
#include "ChunkGroupMatch.h"

using std::list;

ChunkGroupMatch::ChunkGroupMatch(ChunkGroup& to_match) {
  orig = to_match;
  score = std::numeric_limits<double>::infinity();
}

list<ChunkGroupMatch> ChunkGroupMatch::best_match(list<ChunkGroup> &res ,list<ChunkGroup> &src, ChunkGroupCompare &comp) {
  list<ChunkGroupMatch> matches;
  ChunkGroupMatch a_match;
  for (ChunkGroup& g : res) {
    a_match = ChunkGroupMatch(g);
    a_match.best_group_match(src, comp);
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}


void ChunkGroupMatch::best_group_match(list<ChunkGroup>&pool, ChunkGroupCompare& comp) {
  double lowest = score;
  double temp = 0;
  ChunkGroup best_match = match;
  for (ChunkGroup& b : pool) {
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
