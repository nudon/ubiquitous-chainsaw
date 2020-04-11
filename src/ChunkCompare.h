#ifndef FILE_CHUNKCOMPARE_SEEN
#define FILE_CHUNKCOMPARE_SEEN

#include "Chunk.h"

class ChunkCompare {
 public:
  ChunkCompare() {}
  ChunkCompare(double w1, double w2, double w3);
  ~ChunkCompare() {}
  double compare(Chunk a, Chunk b);

 private:
  double freq_center_weight;
  double freq_margin_weight;
  double time_margin_weight;  
};

#endif
