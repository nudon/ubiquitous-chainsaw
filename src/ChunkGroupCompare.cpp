#include "ChunkGroupCompare.h"
#include "util.h"

ChunkGroupCompare::ChunkGroupCompare(double r, double b, double s) {
  ratio_w = r;
  base_freq_w = b;
  size_w = s;
}

double ChunkGroupCompare::compare(ChunkGroup& a, ChunkGroup& b) {
  double diff = 0;
  diff += ratio_w * square_diff(a.get_ratio(), b.get_ratio());
  diff += base_freq_w * square_diff(a.get_rel_base_freq(), b.get_rel_base_freq());
  diff += size_w * square_diff(a.get_size() ,  b.get_size());
  return diff;
}
