#ifndef FILE_CHUNKSTATS_SEEN
#define FILE_CHUNKSTATS_SEEN

//#include <stdio.h>
#include <list>
#include <Eigen/Dense>
#include "Chunk.h"


class ChunkStats {
public :
  ChunkStats();
  ChunkStats(Eigen::MatrixXi chunk_ids);
  ~ChunkStats();

  //would be nice to make these read-only
  Eigen::MatrixXi get_size() {  return stats.col(chunk_size_i); }
  Eigen::MatrixXi get_min_freq() {  return stats.col(min_freq_i); }
  Eigen::MatrixXi get_max_freq() {  return stats.col(max_freq_i); }
  Eigen::MatrixXi get_min_time() {  return stats.col(min_time_i); }
  Eigen::MatrixXi get_max_time() {  return stats.col(max_time_i); }
  Eigen::MatrixXi get_delta_freq() {  return get_max_freq() - get_min_freq();}
  Eigen::MatrixXi get_delta_time() {  return get_max_time() - get_min_time();}
  std::list<Chunk> cull_chunks();
  
private:
  static const int stat_fields;
  static const int chunk_size_i;
  static const int min_freq_i;
  static const int max_freq_i;
  static const int min_time_i;
  static const int max_time_i;

  Eigen::MatrixXi stats;
  Eigen::MatrixXi make_stats(Eigen::MatrixXi chunk_ids);



};

#endif
