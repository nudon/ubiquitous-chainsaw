#ifndef FILE_CHUNKSTATS_NOTSEEN
#define FILE_CHUNKSTATS_NOTSEEN

//#include <stdio.h>
#include <Eigen/Dense>

class ChunkStats {
public :
  ChunkStats();
  ChunkStats(Eigen::MatrixXi chunk_ids);
  ~ChunkStats();
  
private:
  static const int stat_fields;
  static const int chunk_size_i;
  static const int min_freq_i;
  static const int max_freq_i;
  static const int min_time_i;
  static const int max_time_i;

  Eigen::MatrixXi stats;
  Eigen::MatrixXi make_stats(Eigen::MatrixXi chunk_ids);
  Eigen::MatrixXi get_size();
  Eigen::MatrixXi get_min_freq();
  Eigen::MatrixXi get_max_freq();
  Eigen::MatrixXi get_min_time();
  Eigen::MatrixXi get_max_time();
  Eigen::MatrixXi get_delta_freq();
  Eigen::MatrixXi get_delta_time();

};

#endif
