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
  return stats;
}

std::list<Chunk> ChunkStats::cull_chunks() {
  //wait I don't even need the chunkids anymore 
  
  //thinking of returning a linked list of chunk objects
  //because working with matrixes will probably be annoying past this point


  std::list<Chunk> chunk_list;
  MatrixXi sizes = get_size();
  int chunks = sizes.rows();
  //ignoring chunk_id 0
  sizes = get_size().block(1,0, chunks - 1,1);
  chunks = sizes.rows();

  double average_size = (float)sizes.sum() / sizes.sum();
  MatrixXi minf = get_min_freq();
  MatrixXi maxf = get_max_freq();
  MatrixXi mint = get_min_time();
  MatrixXi maxt = get_max_time();
  int chunk_count = 0;
  for (int i = 0; i < chunks; i++) {
    if (sizes(i) > average_size) {
      //gather min, max, delta of time and frequency, construct a chunk
      //find a better name for chunks to many things have chunk in the name
      chunk_count++;
      Chunk temp = Chunk(minf(i), maxf(i), mint(i), maxt(i));
      chunk_list.insert(chunk_list.begin(), temp);
    }
  }
  std::cout << "original size was " << chunks << " culled size is " << chunk_count << "\n";
  return chunk_list; 
}



const int ChunkStats::stat_fields = 5;
const int ChunkStats::chunk_size_i = 0;
const int ChunkStats::min_freq_i = 1;
const int ChunkStats::max_freq_i = 2;
const int ChunkStats::min_time_i = 3;
const int ChunkStats::max_time_i = 4;

