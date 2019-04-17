#include <iostream>
#include <Eigen/Eigen>
#include <list>
#include <math.h>
#include "SongEmbedder.h"
#include "Chunk.h"
#include "ChunkStats.h"
#include "ChunkMatch.h"
#include "util.h"
#include "eigen_to_image.h"


using std::list;
using std::string;
using std::max;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
using stk::FileWvIn;
using stk::FileWvOut;
using stk::StkFrames;

list<Chunk> make_chunks(int samples_per_slice, Song input);
std::pair<MatrixXf, MatrixXf> filter(MatrixXf input);
MatrixXi chunkify(MatrixXf input, int vertical_range, int horizontal_range);
list<Chunk> get_important_chunks(MatrixXi chunks);

MatrixXi fill_chunkify(MatrixXf &input, int thresh);
void recurse_fill(MatrixXf &input, MatrixXi &ids, int x, int y, int cid);


SongEmbedder::SongEmbedder(Song src, Song dst) {
  reciptor = src;
  replacer = dst;
}

double chunkmatch_get_score_wrapper(ChunkMatch a) {
  return a.get_score();
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
    list<ChunkMatch>matches = get_matches(reciptor_chunks, replacer_chunks);
    //don't actually want gsnaz
    //good match scores are zero, bad are high
    //straight snaz results in only making bad matches
    //could fix by subtracting out highest/worst match score then multiplying by -1
    //gsnaz(matches, chunkmatch_get_score_wrapper, 1);
    std::cout << "orig len is " << matches.size();
    reverse_gsnaz(matches, chunkmatch_get_score_wrapper, 1);
    std::cout << "rsnazed len is " << matches.size() << "\n";

    list<ChunkMatch> active_start;
    list<ChunkMatch> active_end;
    for (ChunkMatch a_match : matches) {
      active_start.push_front(a_match);
      active_end.push_front(a_match);
    }

    active_start.sort(ChunkMatch::comp_active_start);
    active_end.sort(ChunkMatch::comp_active_end);
    
    std::string fn = "muxed";
    output_remix(samples_per_slice, active_start, active_end, fn);
    std::cout << "made it to the end!\n";
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
  std::pair<MatrixXf, MatrixXf> xy_filt = filter(spec);
  MatrixXf x_filt = std::get<0>(xy_filt);
  MatrixXf y_filt = std::get<1>(xy_filt);
  int x_vert = 1;
  int x_horz = 2;
  int y_vert = x_vert;
  int y_horz = x_horz;
  int snazr = 1;
  snaz(x_filt, snazr);
  snaz(y_filt, snazr);
  MatrixXi x_chunks = fill_chunkify(x_filt, 0.0);
  MatrixXi y_chunks = fill_chunkify(y_filt, 0.0);
  list<Chunk> part1 = get_important_chunks(x_chunks);
  list<Chunk> part2 = get_important_chunks(y_chunks);
  list<Chunk> ret;
  ret.insert(ret.begin(), part1.begin(), part1.end());
  ret.insert(ret.begin(), part2.begin(), part2.end());

  //std::cout << "making pictures\n";
  string bn = input.get_file_name();
  string post = ".png";
  write_eigen_to_file(bn + ".spec" + post, spec);
  write_eigen_to_file(bn + ".filX" + post, x_filt);
  write_eigen_to_file(bn + ".filY" + post, y_filt);

  write_chunks_to_file(bn + ".chunkX" + post, x_chunks);
  write_chunks_to_file(bn + ".chunkY" + post, y_chunks);
  int xcy = x_chunks.rows();
  int xcx = x_chunks.cols();
  int ycy = y_chunks.rows();
  int ycx = y_chunks.cols();
  write_chunks_to_file(bn + ".cull_chunkX" + post, part1, xcy, xcx);
  write_chunks_to_file(bn + ".cull_chunkY" + post, part2, ycy, ycx);
  write_chunks_to_file(bn + ".cull_chunkTOT" + post, ret, max(xcy, ycy), max(xcx, ycx));

  /*
  int maxsnaz = 3;
  for (int i = 0; i < maxsnaz; i++) {
    std::string number = std::to_string(i);
    snaz(x_filt , 1);
    snaz(y_filt , 1);
    write_eigen_to_file(bn + ".filX" + number + post, x_filt);
    write_eigen_to_file(bn + ".filY" + number + post, y_filt);
  }
  */
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
  int snazr = 1;
  nz_avg = snaz(filt, snazr);
  
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

//the original chunkify had some problems
//this will just assign any non-zero input values the same chunk id as it's adjacent non-zero input values
MatrixXi fill_chunkify(MatrixXf &input, int thresh) {
  int chunk_id = 1;
  int freq_range = input.rows();
  int time_range = input.cols();
  MatrixXi ids(freq_range, time_range);
  ids.setZero();
  for (int time_i = 0; time_i < time_range; time_i++) {
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      if (input(freq_i, time_i) > thresh) {
	recurse_fill(input, ids, time_i, freq_i, chunk_id++);
      }
    }
  }
  return ids;
}

void recurse_fill(MatrixXf &input, MatrixXi &ids, int x, int y, int cid) {
  int rows = input.rows();
  int cols = input.cols();
  if (0 <= x && x < cols && 0 <= y && y < rows) {
    if (input(y,x) > 0.0 && ids(y,x) == 0) {
      ids(y,x) = cid;
      recurse_fill(input, ids, x + 1, y, cid);
      recurse_fill(input, ids, x - 1, y, cid);
      recurse_fill(input, ids, x, y + 1, cid);
      recurse_fill(input, ids, x, y - 1, cid);
    }
  } 
}



//template< template <class T> class container>
std::list<Chunk> get_important_chunks(MatrixXi chunks) {
  ChunkStats stats (chunks);
  int snazr = 1;
  list<Chunk> ret = stats.cull_chunks(snazr);
  return ret;
}


list<ChunkMatch> SongEmbedder::get_matches(list<Chunk> reciptor_chunks, list<Chunk> replacer_chunks) {
  //construct a comparer
  //weights for frequency center, frequency margin, and time maring
  ChunkCompare comp (1,1,1);
  list<ChunkMatch> matches;
  ChunkMatch a_match;
  std::cout << "matching sound chunks! This may take a while.\n";
  int count = 0;
  int tot = reciptor_chunks.size();
  for (Chunk a_chunk : reciptor_chunks) {
    a_match = ChunkMatch(a_chunk);
    a_match.best_match_chunk(replacer_chunks, comp);
    matches.insert(matches.begin(), a_match);
    //std::cout << "out of " << tot << " matches completed " << ++count << "\n";
  }
  std::cout << "Done with matching!\n";
  return matches;
}


void SongEmbedder::output_remix(int samples_per_slice, list<ChunkMatch> active_start, list<ChunkMatch> active_end, std::string fn) {
  int new_slices = reciptor.get_samples_per_channel() / samples_per_slice;
  int channels = reciptor.get_channels();
  float offset = -1;
  FileWvIn sample_src;
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames filt_frames (samples_per_slice,1);
  StkFrames temp_frames (samples_per_slice,1);
  StkFrames zero (0, samples_per_slice, channels);
    
  sample_src.openFile(replacer.get_file_name());
  output.openFile("muxed.wav", channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  bool print_some_samples = false;
  list<ChunkMatch> active_chunks;
  for (int i = 0; i < new_slices; i++) {
    //add newly activated chunks
    while(active_start.size() > 0 && active_start.front().get_active_start() <= i) {
      active_chunks.insert(active_chunks.begin(), active_start.front());
      active_start.pop_front();
    }   
    //take active chunks, apply filter to replacer at specified time
    //sum into some output
    frame_out *= zero;
    for (ChunkMatch chunk_match : active_chunks) {
      Chunk match = chunk_match.get_match_chunk();
      ChunkFilter filt = match.get_filter();
      sample_src.reset();
      //need to find time offset into other file
      //offset into chunk + start of chunk in sample_src
      offset = (i - chunk_match.get_active_start()) + (match.get_time_center() - match.get_time_margin());
      offset *= samples_per_slice;
      //not sure which of these or any are even correct
      //does skipping ahead 4 samples skip 4 samples in all channels?
      //does it split them evenly? what about skipping odd amounts?
      //offset /= replacer.get_tot_samples();
      //offset *= replacer.get_tot_samples();
      sample_src.addTime(offset);
      
      sample_src.tick(frame_in);
      for (int i = 0; i < channels; i++) {
	frame_in.getChannel(i, temp_frames,0);
	filt_frames = filt.fir_filter_frame(temp_frames);
	frame_out.getChannel(i, temp_frames, 0);
	temp_frames += filt_frames;
	frame_out.setChannel(i, temp_frames, 0);
      }
      if (print_some_samples && i % 777 == 0) {
	for (int i = 0; i < samples_per_slice; i++) {
	  double out = frame_out[i];
	  double in = frame_in[i];
	  std::cout << "in :" << in << " out :" << out << "\n ";
	}

      }
    }

    output.tick(frame_out);
      
    //remove chunks that just ended
    while(active_end.size() > 0 && active_end.front().get_active_end() <= i) {
      active_chunks.remove(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}
