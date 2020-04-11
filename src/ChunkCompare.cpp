#include <math.h>
#include <tgmath.h>
#include "ChunkCompare.h"

double abs_log_diff(double a, double b);
double abs_diff(double a, double b);
double square_diff(double a, double b);

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

double abs_log_diff(double a, double b) {
  if (a == 0 || b == 0) {
    a += 0.1;
    b += 0.1;
  }
  return abs(log2(a) - log2(b));
}

double abs_diff(double a, double b) {
  return abs(a - b);
}

double square_diff(double a, double b) {
  return pow(a - b, 2);
}