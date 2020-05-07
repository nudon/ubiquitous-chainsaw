#ifndef FILE_CHUNKGROUPCOMPARE_SEEN
#define FILE_CHUNKGROUPCOMPARE_SEEN

#include "ChunkGroup.h"

class ChunkGroupCompare {
 public:
  ChunkGroupCompare() {}
  ChunkGroupCompare(double r, double b, double s);
  ~ChunkGroupCompare() {}
  double compare(ChunkGroup& a, ChunkGroup& b);

 private:
  double ratio_w;
  double base_freq_w;
  double size_w;
};


#endif
