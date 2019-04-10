#include <iostream>
#include <Eigen/Eigen>
#include <list>
#include "SongEmbedder.h"
#include "Chunk.h"
#include "ChunkStats.h"
#include "ChunkMatch.h"
#include "util.h"


using std::list;
using Eigen::MatrixXf;
using Eigen::MatrixXi;

list<Chunk> make_chunks(int samples_per_slice, Song input);
std::pair<MatrixXf, MatrixXf> filter(MatrixXf input);
MatrixXi chunkify(MatrixXf input, int vertical_range, int horizontal_range);
std::list<Chunk> get_important_chunks(MatrixXi chunks);

SongEmbedder::SongEmbedder(Song src, Song dst) {
  reciptor = src;
  replacer = dst;
}

//placeholder heavy lifting function
void SongEmbedder::funk() {
  int samples_per_slice = 1024;
  if (reciptor.get_file_rate() == replacer.get_file_rate()) {
    //generate tfd's of both songs
    //take x/y sobel, apply difference of gaussian
    //then call chunkify, ChunkStats, then cull_chunks, join x/y chunks
    list<Chunk> reciptor_chunks = make_chunks(samples_per_slice, reciptor);
    list<Chunk> replacer_chunks = make_chunks(samples_per_slice, replacer);

    //construct a comparer
    //weights for frequency center, frequency margin, and time maring
    ChunkCompare comp (1,1,1);
    
    //iterate over 1 chunk list, create and fill a bunch of match structs
    list<ChunkMatch> matches;
    ChunkMatch a_match;
    for (Chunk a_chunk : reciptor_chunks) {
      a_match = ChunkMatch(a_chunk);
      a_match.best_match_chunk(replacer_chunks, comp);
      matches.insert(matches.begin(), a_match);
    }
    
    //then, need to mix in matches into an output file
    
    //get start/end times for when chunk needs to be mixed into output
    
    //have 2 sorted lists of chunks based on start/end times
    
    //have some container of active chunks
    
    //if current time slice == first in start time, move chunk from starttimelist into active chunks, repeat until bool condition evals to false
    
    //then, for each thing in active chunks, take input tick, apply a filter to extract a part of each active chunk, mix/sum into some output tick
    
    //then check if time slice == first in end time. if true, remove chunk from active_chunks and endtimelist, repeat until bool evals to false
    std::cout << "made it to the end!";
  }
  else {
    //ERROR
    //would have to do weird stuff to work with this
    //because the different sampling rates would manifest in fft
    //would have the frequency bins of both spectrograms not line up
  }
}



//template< template <class T> class container>
list<Chunk> make_chunks(int samples_per_slice, Song input) {
  MatrixXf spec = input.spectrogram(samples_per_slice);
  std::pair<MatrixXf, MatrixXf> xy_chunks = filter(spec);
  int x_vert = 1;
  int x_horz = 2;
  int y_vert = x_vert;
  int y_horz = x_horz;
  MatrixXi x_chunks = chunkify(std::get<0>(xy_chunks), x_vert, x_horz);
  MatrixXi y_chunks = chunkify(std::get<1>(xy_chunks), y_vert, y_horz);
  list<Chunk> part1 = get_important_chunks(x_chunks);
  list<Chunk> part2 = get_important_chunks(y_chunks);
  list<Chunk> ret;
  ret.insert(ret.begin(), part1.begin(), part1.end());
  ret.insert(ret.begin(), part2.begin(), part2.end());
  return ret;
}



std::pair<MatrixXf, MatrixXf> filter(MatrixXf input) {
  MatrixXf gauss_col = create_1d_gaussian_filter_col(16, 1, 1);
  MatrixXf gauss_row = create_1d_gaussian_filter_row(16, 1, 1);
  MatrixXf gaussian_blur = one_d_convolve(input, gauss_col);
  gaussian_blur = one_d_convolve(gaussian_blur, gauss_row);
  MatrixXf dogg = dog(input, gaussian_blur);
  MatrixXf sobel_x (1,3);
  MatrixXf sobel_y (3,1);
  sobel_x << -1, 0, 1;
  sobel_y << -1, 0, 1;

  MatrixXf x_sobeled_dog = one_d_convolve(dogg, sobel_x);
  MatrixXf y_sobeled_dog = one_d_convolve(dogg, sobel_y);
  std::pair<MatrixXf, MatrixXf> ret(x_sobeled_dog, y_sobeled_dog);
  return ret;
}


MatrixXi chunkify(MatrixXf input, int vertical_range, int horizontal_range) {
  int chunk_id = 0;
  int freq_range = input.rows();
  int time_range = input.cols();
  int temp_i = -1;
  int temp = 0;
  float nz_avg = 0;
  MatrixXf filt (freq_range, time_range);
  MatrixXi chunk_ids (freq_range, time_range);
  filt << input;
  chunk_ids.setZero( );
  //threshold values < average to zero
  //do this to not care about very minor sounds
  nz_avg = snaz(filt, 1);
  
  //isn't a perfect grouping algorithim, since the local grouping is kind of naive
  //hoping further culling methods down the line will get rid of bad groups
  //but loop could be improved by having a more robust/expensive joining loop(s)
  for (int time_i = 0; time_i < time_range; time_i++) {
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      if (filt(freq_i, time_i) > nz_avg) {
	//check if you want to join it vertically
	for (int freq_off = 1; freq_off <= vertical_range; freq_off++) {
	  //group into chunks within freq_off indexes below current cell 
	  temp_i = freq_i - freq_off;
	  if (temp_i >= 0) {
	    temp = chunk_ids(temp_i, time_i);
	    if (temp != 0) {
	      chunk_ids(freq_i, time_i) = temp;
	      break;
	    }
	  }
	  else {
	    break;
	  }
	}
	//check if you want to join it horizontally, potentially overwriting vertical joins
	for (int time_off = 1; time_off <= horizontal_range; time_off++) {
	  temp_i = time_i - time_off;
	  if (temp_i >= 0) {
	    temp = chunk_ids(freq_i, temp_i);
	    if (temp != 0) {
	      chunk_ids(freq_i, time_i) = temp;
	      break;
	    }
	  }
	  else {
	    break;
	  }
	}
	//else assigning to a new group
	if (chunk_ids(freq_i, time_i) == 0) {
	  chunk_ids(freq_i, time_i) = ++chunk_id;
	}
      }
    }
  }
  //basic error checking, making sure sufficiently low values aren't being assigned chunk_ids
  for (int time_i = 0; time_i < time_range; time_i++) {
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      if (chunk_ids(freq_i, time_i) != 00 && filt(freq_i, time_i) < nz_avg) {
	std::cout << filt(freq_i, time_i) <<  "bad! ";
      }
    }
  }
  return chunk_ids;
}



//template< template <class T> class container>
std::list<Chunk> get_important_chunks(MatrixXi chunks) {
  ChunkStats stats (chunks);
  return stats.cull_chunks();
}
