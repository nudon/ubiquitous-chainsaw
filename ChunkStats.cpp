#include <iostream>
#include "ChunkStats.h"
using Eigen::MatrixXi;

ChunkStats::ChunkStats(MatrixXi chunk_ids) {
  stats = make_stats(chunk_ids);
}

ChunkStats::~ChunkStats() {

}

MatrixXi ChunkStats::make_stats(MatrixXi chunk_ids) {
  int ci = 0, ri = 0;
  chunk_ids.maxCoeff(&ri, &ci);
  int num_chunks = chunk_ids(ri, ci);
  std::cout << num_chunks << " many chunks! \n";
  //potentially construct this to be row_major
  MatrixXi stats = MatrixXi::Constant(num_chunks + 1, stat_fields, -1);
  stats.col(chunk_size_i) = MatrixXi::Constant(num_chunks + 1, 1, 0);

  int freq_range = chunk_ids.rows();
  int time_range = chunk_ids.cols();
  int chunk = 0;
  for (int time_i = 0; time_i < time_range; time_i++) {
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      chunk = chunk_ids(freq_i, time_i);
      stats(chunk, chunk_size_i) += 1;
      //count for no-group/silence might be interesting
      //that min/max freq/time will not be interesting
      if (chunk > 0) {
	stats(chunk, min_freq_i) = std::min(freq_i, stats(chunk, min_freq_i));
	stats(chunk, max_freq_i) = std::max(freq_i, stats(chunk, max_freq_i));
	
	stats(chunk, min_time_i) = std::min(time_i, stats(chunk, min_time_i));
	stats(chunk, max_time_i) = std::max(time_i, stats(chunk, max_time_i));
      }
    }
  }
  stats.col(chunk_size_i).minCoeff(&ri);
  std::cout << "smallest size is " << stats(ri, chunk_size_i) << "\n";
  return stats;
}

MatrixXi ChunkStats::get_size() {
  return stats.col(chunk_size_i);
}

MatrixXi ChunkStats::get_min_freq() {
  return stats.col(min_freq_i);
}

MatrixXi ChunkStats::get_max_freq() {
  return stats.col(max_freq_i);
}

MatrixXi ChunkStats::get_min_time() {
  return stats.col(min_time_i);
}

MatrixXi ChunkStats::get_max_time(){
  return stats.col(max_time_i);
}

MatrixXi ChunkStats::get_delta_freq() {
  return get_max_freq() - get_min_freq();
}

MatrixXi ChunkStats::get_delta_time() {
  return get_max_time() - get_min_time();
}


const int ChunkStats::stat_fields = 5;
const int ChunkStats::chunk_size_i = 0;
const int ChunkStats::min_freq_i = 1;
const int ChunkStats::max_freq_i = 2;
const int ChunkStats::min_time_i = 3;
const int ChunkStats::max_time_i = 4;

