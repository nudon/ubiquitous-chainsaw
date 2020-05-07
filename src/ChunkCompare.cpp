#include "ChunkCompare.h"
#include "util.h"

ChunkCompare::ChunkCompare(double w1, double w2, double w3) {
  freq_center_weight = w1;
  freq_margin_weight = w2;
  time_margin_weight = w3;  
}

double ChunkCompare::compare(Chunk a, Chunk b) {
  double diff = 0;
  diff += freq_center_weight * abs_log_diff(a.get_rel_freq_center(), b.get_rel_freq_center());
  diff += freq_margin_weight * square_diff(a.get_freq_margin(), b.get_freq_margin());
  diff += time_margin_weight * square_diff(a.get_time_margin(), b.get_time_margin());
  
  return diff;
}
