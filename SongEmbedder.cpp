#include <iostream>
#include <Eigen/Eigen>
#include <list>
#include <map>
#include <unordered_set>
#include <math.h>
#include <LentPitShift.h>
#include "SongEmbedder.h"
#include "Chunk.h"
#include "ChunkStats.h"
#include "ChunkMatch.h"
#include "util.h"
#include "eigen_to_image.h"


using std::list;
using std::map;
using std::unordered_set;
using std::pair;
using std::string;
using std::max;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
using stk::FileWvIn;
using stk::FileWvOut;
using stk::StkFrames;
using stk::LentPitShift;

list<Chunk> make_chunks(int samples_per_slice, Song input);
pair<MatrixXf, MatrixXf> filter(MatrixXf &input);
MatrixXi chunkify(MatrixXf input, int vertical_range, int horizontal_range);
list<Chunk> get_important_chunks(MatrixXi chunks, int snazr);
list<ChunkMatch> best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
list<Chunk> get_closest_bin(map<int, list<Chunk>> bins, int start);
MatrixXi fill_chunkify(MatrixXf &input, int thresh);
void recurse_fill(MatrixXf &input, MatrixXi &ids, int x, int y, int cid);




void test_spec_filter();
void test_sound_filter(int samples_per_slice, std::string fn);


void print_frame(StkFrames in);

void test_filter();

SongEmbedder::SongEmbedder(Song src, Song dst) {
  reciptor = src;
  replacer = dst;
}

double chunkmatch_get_score_wrapper(ChunkMatch a) {
  return a.get_score();
}



void SongEmbedder::funk() {
  int samples_per_slice = 1024 * 2;
  //FilterTest(samples_per_slice, reciptor.get_file_name());
  if (reciptor.get_file_rate() == replacer.get_file_rate()) {
    list<Chunk> reciptor_chunks = make_chunks(samples_per_slice, reciptor);
    list<Chunk> replacer_chunks = make_chunks(samples_per_slice, replacer);
    list<ChunkMatch>matches = get_matches(reciptor_chunks, replacer_chunks);
    
    std::string fn = "muxed";
    std::cout << "outputting song\n";
    output_remix(samples_per_slice, matches, fn);
    //output_remix_alt(samples_per_slice, matches, "bombs");
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
  int filt_snazr = 1;
  int chunk_snazr = 1;
  std::cout << "Building spectrogram for " << input.get_file_name() << "\n";
  MatrixXf spec = input.spectrogram(samples_per_slice);
  foobar_spec(spec);
  std::cout << "filtering spectrogram\n";
  pair<MatrixXf, MatrixXf> xy_filt = filter(spec);
  MatrixXf time_edges = std::get<0>(xy_filt);
  MatrixXf freq_edges = std::get<1>(xy_filt);


  snaz(time_edges, filt_snazr);
  snaz(freq_edges, filt_snazr);
  std::cout << "extracting chunks from spectrogram\n";
  MatrixXi time_chunks = fill_chunkify(time_edges, 0.0);
  MatrixXi freq_chunks = fill_chunkify(freq_edges, 0.0);
  //carry on
  std::cout << "culling chunks\n";
  ChunkStats time_stats (time_chunks);
  ChunkStats freq_stats (freq_chunks);
  list<Chunk> part1 = time_stats.cull_chunks(chunk_snazr, CHUNK_TIME_OPT);
  list<Chunk> part2 = freq_stats.cull_chunks(chunk_snazr, CHUNK_FREQ_OPT);
  list<Chunk> ret;
  ret.insert(ret.begin(), part1.begin(), part1.end());
  ret.insert(ret.begin(), part2.begin(), part2.end());
  int tot_count = time_stats.get_chunk_count() + freq_stats.get_chunk_count();
  std::cout << tot_count << " many chunk(s) extracted from " << input.get_file_name() << "\n";
  std::cout << ret.size() << " many chunk(s) left after culling " << input.get_file_name() << "\n";

  //std::cout << "making pictures\n";
  string bn = input.get_file_name();
  string post = ".png";
  write_eigen_to_file(bn + ".time_edges" + post, time_edges);
  write_eigen_to_file(bn + ".freq_edges" + post, freq_edges);
  write_eigen_to_file(bn + ".spec" + post, spec);
  write_chunks_to_file(bn + ".time_chunks" + post, time_chunks);
  write_chunks_to_file(bn + ".freq_chunks" + post, freq_chunks);
  write_chunks_to_file(bn + ".time_chunks_culled" + post, part1, time_chunks.rows(), time_chunks.cols());
  write_chunks_to_file(bn + ".freq_chunks_culled" + post, part2, freq_chunks.rows(), freq_chunks.cols());
  return ret;
}

pair<MatrixXf, MatrixXf> filter(MatrixXf &input) {
  //gauss col is a collunm vec, row is a row vec
  int gauss_col_dim = 33;
  int gauss_row_dim = 17;
  double div = 9.0;
  //the stddev is roughly set with div
  //only tested with gauss dims around ~3-100
  //going with higher divs seems to have the gradient fall of to zero sooner
  // div = 9 gives an approximate 1->0 gradient from center to edge
  // div = 4 gives an approximate 1->0.5 gradient from center to edge
  MatrixXf gauss_col = create_1d_gaussian_filter_col(gauss_col_dim, 1, ceil(gauss_col_dim * gauss_col_dim / div));
  MatrixXf gauss_row = create_1d_gaussian_filter_row(gauss_row_dim, 1, ceil(gauss_row_dim * gauss_row_dim / div));
  write_eigen_to_file("gauss_col.png", gauss_col);
  write_eigen_to_file("gauss_row.png", gauss_row);
  //x gauss is image blurred horizontally, y is image blurred vertically
  MatrixXf x_gauss = one_d_convolve(input, gauss_row);
  MatrixXf y_gauss = one_d_convolve(input, gauss_col);

  //extract difference of gaussian
  MatrixXf x_dog = input - x_gauss;
  MatrixXf y_dog = input - y_gauss;
  
  //smoothes edges so barely seperated chunks are joined
  //x_dog picked out vertical edges, so smooth vertically to join. vice versa for y and horizontal
  MatrixXf x_dog_smooth = one_d_convolve(x_dog, gauss_col);
  MatrixXf y_dog_smooth = one_d_convolve(y_dog, gauss_row);
  pair<MatrixXf, MatrixXf> gret(x_dog_smooth, y_dog_smooth);
  return gret;
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
std::list<Chunk> get_important_chunks(MatrixXi chunks, int snazr) {
  ChunkStats stats (chunks);
  list<Chunk> ret = stats.cull_chunks(snazr);
  return ret;
}


list<ChunkMatch> SongEmbedder::get_matches(list<Chunk> &reciptor_chunks, list<Chunk> &replacer_chunks) {
  //construct a comparer
  //weights for frequency center, frequency margin, and time maring
  ChunkCompare comp (50,1,10);

  int snazr = 1;
  std::cout << "matching sound chunks!\n";
  //list<ChunkMatch> matches = best_match(reciptor_chunks, replacer_chunks, comp);
  list<ChunkMatch> matches = quick_match(reciptor_chunks, replacer_chunks, comp);
  std::cout << "Done with matching!\n";
  if (snazr > 0) {
    std::cout << "culling bad matches\n";
    std::cout << "orig len is " << matches.size() << "\n";
    double avg = reverse_gsnaz(matches, chunkmatch_get_score_wrapper, snazr);
    std::cout << "only accepted scores <= " << avg << "\n";
    std::cout << "culled len is " << matches.size() << "\n";
  }

  return matches;
}

list<ChunkMatch> best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
  list<ChunkMatch> matches;
  ChunkMatch a_match;
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    a_match.best_match_chunk(rep, comp);
    a_match.get_match_chunk().make_chunk_filter();
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
  map<int, list<Chunk>> freq_cent_map;
  map<int, list<Chunk>> freq_marg_map;
  map<int, list<Chunk>> time_marg_map;
  list<ChunkMatch> matches;
  //map<int, int> 
  int fc_ind, fm_ind, tm_ind;
  //populate maps
  list<Chunk> sub_list;
  for (Chunk a_chunk : rep) {
    fm_ind = a_chunk.get_freq_length();
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();
    
    freq_marg_map[fm_ind].push_front(a_chunk);
    time_marg_map[tm_ind].push_front(a_chunk);
    freq_cent_map[fc_ind].push_front(a_chunk);
  }

  ChunkMatch a_match;
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    fm_ind = a_chunk.get_freq_length();
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();

    sub_list = freq_marg_map[fm_ind];
    //sub_list = get_closest_bin(freq_marg_map, fm_ind);
    a_match.best_match_chunk(sub_list, comp);
    
    sub_list = time_marg_map[tm_ind];
    //sub_list = get_closest_bin(time_marg_map, tm_ind);
    a_match.best_match_chunk(sub_list, comp);
    
    sub_list = freq_cent_map[fc_ind];
    //sub_list = get_closest_bin(freq_cent_map, fc_ind);
    a_match.best_match_chunk(sub_list, comp);
    
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

list<Chunk> get_closest_bin(map<int, list<Chunk>> bins, int start) {
  //various bins may have zero items in them
  //find the closest non-empty one
  bool done = false;
  int offset = 0;
  int ind = -1;
  while(!done) {
    if (bins[start + offset].size() > 0) {
      done = true;
      ind = start + offset;
    }
    else if (start - offset >= 0) {
      if (bins[start - offset].size() > 0) {
	done = true;
	ind = start - offset;
      }
    }
    offset++;
  }
  return bins[ind];
}


struct cm_hash {
  size_t operator ()(const ChunkMatch & cm) const {
    return ChunkMatch::match_hash((ChunkMatch&)cm);
    //return 55;
  }
};

void SongEmbedder::output_remix(int samples_per_slice, list<ChunkMatch> &matches, std::string fn) {
  list<ChunkMatch> active_start;
  list<ChunkMatch> active_end;
  for (ChunkMatch a_match : matches) {
    active_start.push_front(a_match);
    active_end.push_front(a_match);
  }
  
  active_start.sort(ChunkMatch::comp_active_start);
  active_end.sort(ChunkMatch::comp_active_end);
  int new_slices = reciptor.get_samples_per_channel() / samples_per_slice;
  int channels = reciptor.get_channels();
  float offset = -1;
  FileWvIn sample_src;
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  
  StkFrames temp (0, samples_per_slice, channels);
  StkFrames rev_feather (0, samples_per_slice, channels);
  StkFrames reg_feather (0, samples_per_slice, channels);
  int feather_margin = samples_per_slice;
  float feather_val = 0;
  for (int i = 0; i < samples_per_slice; i++) {
    feather_val = std::min(1.0f, (float)i / feather_margin);
    for (int chan = 0; chan < channels; chan++) {
      reg_feather(i, chan) = feather_val;
      rev_feather(samples_per_slice - 1 - i, chan) = feather_val;
    }
  }
    
  sample_src.openFile(replacer.get_file_name());
  output.openFile("muxed.wav", channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  //list<ChunkMatch> active_chunks;
  unordered_set<ChunkMatch,cm_hash> active_chunks; 
  for (int i = 0; i < new_slices; i++) {
    //add newly activated chunks
    while(active_start.size() > 0 && active_start.front().get_active_start() == i) {
      active_chunks.insert(active_start.front());
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
      offset = (i - chunk_match.get_active_start()) + match.get_time_start();
      offset *= samples_per_slice;
      sample_src.addTime(offset);
      
      sample_src.tick(frame_in);

      temp *= zero;
      filt.fir_filter_frame(frame_in, temp);
      if (chunk_match.get_active_start() == i) {
	temp *= reg_feather;
      }
      if (chunk_match.get_active_end() == i) {
	temp *= rev_feather;
      }
      frame_out += temp;
    }

    output.tick(frame_out);
    
    //remove chunks that just ended
    while(active_end.size() > 0 && active_end.front().get_active_end() == i) {
      active_chunks.erase(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}




//unfinished time compress version
void SongEmbedder::output_remix_alt(int samples_per_slice, list<ChunkMatch> &matches, std::string fn) {
  list<ChunkMatch> active_start;
  list<ChunkMatch> active_end;
  list<ChunkMatch> active_chunks;
  for (ChunkMatch a_match : matches) {
    active_start.push_front(a_match);
    active_end.push_front(a_match);
  }

  //new sort methods, just original reciptor time
  active_start.sort(ChunkMatch::comp_orig_start);
  active_end.sort(ChunkMatch::comp_orig_end);
  
  int new_slices = reciptor.get_samples_per_channel() / samples_per_slice;
  int channels = reciptor.get_channels();
  float offset = -1;
  FileWvIn sample_src;
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  sample_src.openFile(replacer.get_file_name());
  output.openFile(fn, channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);

  LentPitShift lent(1, samples_per_slice);
  
  for (int i = 0; i < new_slices; i++) {
    while(active_start.size() > 0 && active_start.front().get_active_start() <= i) {
      active_chunks.insert(active_chunks.begin(), active_start.front());
      active_start.pop_front();
    }   
    frame_out *= zero;
    for (ChunkMatch chunk_match : active_chunks) {
      Chunk match = chunk_match.get_match_chunk();
      ChunkFilter filt = match.get_filter();
      sample_src.reset();
      Chunk orig = chunk_match.get_orig_chunk();
      double offset_ratio =  (double)(i - orig.get_time_start()) / (orig.get_time_length());
      offset = (offset_ratio * match.get_time_length() + match.get_time_start()) * samples_per_slice;
      sample_src.addTime(offset);

      //then, extract a variable amount of samples
      double next_offset_ratio =  (double)(i + 1  - orig.get_time_start()) / (orig.get_time_length());
      double next_offset = (next_offset_ratio * match.get_time_length() + match.get_time_start()) * samples_per_slice;
      int size1 = next_offset - offset;
      double compression_ratio = (double)orig.get_time_length() / match.get_time_length();
      int size2 = samples_per_slice / compression_ratio;
      //tick size samples from sample_src at curr offset, then time stretch it
      //result should be same size as frame_out if the correct amount of samples got correctly stretched
      int size = size2;
      StkFrames temp_in (size, channels);
      StkFrames filt_in (size, channels);
      sample_src.tick(temp_in);
      
      filt.fir_filter_frame(temp_in, filt_in);
      reshape_chunk(filt_in, frame_out, orig, match, lent);
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


void test_spec_filter() {
  string post = "_test.png";
  int samp_rows = 7;
  int samp_cols = 21;
  MatrixXf sample(samp_rows, samp_cols);
  sample.setZero();
  sample.row((samp_rows - 1) / 2) = MatrixXf::Constant(1, samp_cols, 1);
  sample.col((samp_cols - 1) / 2) = MatrixXf::Constant(samp_rows, 1, 1);

  write_eigen_to_file("samples" + post, sample);
  //gauss col is a collunm vec, row is a row vec
  MatrixXf gauss_col = create_1d_gaussian_filter_col(5, 1, 1);
  MatrixXf gauss_row = create_1d_gaussian_filter_row(5, 1, 1);
  write_eigen_to_file("gausscol" + post, gauss_col);
  write_eigen_to_file("gaussrow" + post, gauss_row);

  //x gauss is image blurred horizontally, y is image blurred vertically
  MatrixXf x_gauss = one_d_convolve(sample, gauss_row);
  MatrixXf y_gauss = one_d_convolve(sample, gauss_col);
  write_eigen_to_file("x_gauss" + post, x_gauss);
  write_eigen_to_file("y_gauss" + post, y_gauss);
  

  MatrixXf x_dog = (sample - x_gauss);
  MatrixXf y_dog = (sample - y_gauss);
  write_eigen_to_file("x_dog" + post, x_dog);
  write_eigen_to_file("y_dog" + post, y_dog);
  
  //smoothes edges so barely seperated chunks are joined
  MatrixXf x_dog_smooth = one_d_convolve(x_dog, gauss_col);
  MatrixXf y_dog_smooth = one_d_convolve(y_dog, gauss_row);
  write_eigen_to_file("x_dog_smooth" + post, x_dog_smooth);
  write_eigen_to_file("y_dog_smooth" + post, y_dog_smooth);
}


void test_sound_filter(int samples_per_slice, std::string fn) {
  FileWvIn sample_src;
  sample_src.openFile(fn);
  int channels = sample_src.channelsOut();
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames filt (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  StkFrames temp (samples_per_slice, 1);
  StkFrames temp_out (samples_per_slice, 1);
  
  output.openFile("filt_test.wav", channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  float cent = samples_per_slice / 4;
  float marg = samples_per_slice / 100;
  ChunkFilter filter1 (1 * cent, marg, samples_per_slice);
  ChunkFilter filter2 (2 * cent, marg, samples_per_slice);
  ChunkFilter filter3 (3 * cent, marg, samples_per_slice);
  list<ChunkFilter> filters;
  filters.push_front(filter1);
  filters.push_front(filter2);
  filters.push_front(filter3);
   
  while(!sample_src.isFinished()) {
    frame_out *= zero;
    filt *= zero;
    sample_src.tick(frame_in);
    for (ChunkFilter a_filter : filters) {
      a_filter.fir_filter_frame(frame_in, filt);
      frame_out += filt;
    }
    output.tick(frame_out);
  }
  output.closeFile();
}


