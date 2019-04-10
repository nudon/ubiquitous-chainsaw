#include <math.h>
#include "ChunkCompare.h"

double abs_diff(double a, double b);
double squared_diff(double a, double b);

ChunkCompare::ChunkCompare(double w1, double w2, double w3) {
  freq_center_weight = w1;
  freq_margin_weight = w2;
  time_margin_weight = w3;  
}

double ChunkCompare::compare(Chunk a, Chunk b) {
  double diff = 0;
  double (*func)(double , double ) = abs_diff;
  diff += freq_center_weight * func(a.get_freq_center(), b.get_freq_center());
  diff += freq_margin_weight * func(a.get_freq_margin(), b.get_freq_margin());
  diff += time_margin_weight * func(a.get_time_margin(), b.get_time_margin());
  
  return diff;
}

double abs_diff(double a, double b) {
  return abs(a - b);
}

double squared_diff(double a, double b) {
  return pow(a - b, 2);
}
