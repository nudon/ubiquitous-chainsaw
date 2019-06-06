#include <iostream>
#include <limits>
#include "ChunkStats.h"
#include "util.h"
using Eigen::MatrixXi;

double get_size_wrap(Chunk a);
double get_freq_wrap(Chunk a);
double get_time_wrap(Chunk a);

ChunkStats::ChunkStats(MatrixXi chunk_ids) {
  bin_size = chunk_ids.rows();
  bin_len = chunk_ids.cols();
  stats = make_stats(chunk_ids);
}

ChunkStats::~ChunkStats() {

}

MatrixXi ChunkStats::make_stats(MatrixXi chunk_ids) {
  int ci = 0, ri = 0;
  chunk_ids.maxCoeff(&ri, &ci);
  int num_chunks = chunk_ids(ri, ci);
  //potentially construct this to be row_major
  MatrixXi stats = MatrixXi::Constant(num_chunks + 1, stat_fields, -1);
  double max_val = std::numeric_limits<int>::max();
  double min_val = -1 * std::numeric_limits<int>::max() ; 
  stats.col(chunk_size_i) = MatrixXi::Constant(num_chunks + 1, 1, 0);
  stats.col(min_freq_i) = MatrixXi::Constant(num_chunks + 1, 1, max_val);
  stats.col(min_time_i) = MatrixXi::Constant(num_chunks + 1, 1, max_val);
  stats.col(max_freq_i) = MatrixXi::Constant(num_chunks + 1, 1, min_val);
  stats.col(max_time_i) = MatrixXi::Constant(num_chunks + 1, 1, min_val);

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
  return stats;
}

std::list<Chunk> ChunkStats::cull_chunks(int snazr, int opt) {
  std::list<Chunk> chunk_list;
  MatrixXi sizes = get_size();
  MatrixXi minf = get_min_freq();
  MatrixXi maxf = get_max_freq();
  MatrixXi mint = get_min_time();
  MatrixXi maxt = get_max_time();
  int chunk_count = 0;
  int chunks = sizes.rows();
  double average_size = 0;
  
  //chunk id zero is not a chunk , just stats for everything that didn't get assigned one
  //ignore it when taking the average and when looping
  //generalize this into another snaz
  for (int i = 1; i < chunks; i++) {
    if (sizes(i) > average_size) {
      chunk_count++;
      Chunk temp = Chunk(minf(i), maxf(i), mint(i), maxt(i), bin_size, bin_len,i);
      // std::cout << minf(i) << " " << maxf(i) << " " << mint(i) << " " << maxt(i) << "\n";
      chunk_list.insert(chunk_list.begin(), temp);
    }
  }
  //average_size = snaz(chunk_list, snazr);
  double (*the_func)(Chunk) = NULL;
  if (opt == 0) {
    the_func = get_size_wrap;
  }
  else if (opt == CHUNK_TIME_OPT) {
    the_func = get_time_wrap;
  }
  else if (opt == CHUNK_FREQ_OPT) {
    the_func = get_freq_wrap;
  }
  average_size = gsnaz(chunk_list, the_func, snazr);
  //std::cout << "original size was " << chunks << " culled size is " << chunk_list.size() << "\n";
  //std::cout << "average size was " << average_size << "\n";
  return chunk_list; 
}

double get_size_wrap(Chunk a) {
  return a.get_chunk_size();
}

double get_freq_wrap(Chunk a) {
  return a.get_freq_margin();
}

double get_time_wrap(Chunk a) {
  return a.get_time_margin();
}


const int ChunkStats::stat_fields = 5;
const int ChunkStats::chunk_size_i = 0;
const int ChunkStats::min_freq_i = 1;
const int ChunkStats::max_freq_i = 2;
const int ChunkStats::min_time_i = 3;
const int ChunkStats::max_time_i = 4;

