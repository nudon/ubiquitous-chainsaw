#include <limits>
#include "ChunkMatch.h"

ChunkMatch::ChunkMatch(Chunk to_match) {
  orig = to_match;
}
//template< template <class T> class container>
void ChunkMatch::best_match_chunk(std::list<Chunk> many_chunks, ChunkCompare comp) {
  double lowest = std::numeric_limits<double>::infinity();
  double temp = 0;
  Chunk best_match;
  for (Chunk b : many_chunks) {
    temp = comp.compare(orig, b);
    if (temp < lowest) {
      best_match = b;
      lowest = temp;
    }
  }
  match = best_match;
  score = lowest;
}
