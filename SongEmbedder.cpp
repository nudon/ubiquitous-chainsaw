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


pair<MatrixXf, MatrixXf> extract_edges(MatrixXf &input);
MatrixXi chunkify(MatrixXf input, int vertical_range, int horizontal_range);
list<ChunkMatch> best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
int get_closest_bin_index(map<int, list<Chunk>> bins, int start);
MatrixXi fill_chunkify(MatrixXf &input, int thresh);
void recurse_fill(MatrixXf &input, MatrixXi &ids, int x, int y, int cid);


void test_filter();
void test_spec_filter();
void test_sound_filter(int samples_per_slice, std::string fn);

void test_center_grouping(MatrixXi &song_chunks);
list<Chunk> group_chunks(list<Chunk> &chunks);

void test_new_time_freq_extract(MatrixXf &raw_spec, MatrixXf &time_edges);
void test_new_time_recurse(MatrixXf& raw_spec, MatrixXf &time_edges, MatrixXi &groups, double raw_thresh, int assign_id, int freq_i, int time_i);



SongEmbedder::SongEmbedder(Song src, Song dst, int sps, int edge_r, int chunk_r, int match_r) {
  reciptor = src;
  replacer = dst;
  samples_per_slice = sps;
  edge_snazr = edge_r;
  chunk_snazr = chunk_r;
  match_snazr = match_r;  
}

double chunkmatch_get_score_wrapper(ChunkMatch a) {
  return a.get_score();
}



void SongEmbedder::funk() {
  if (reciptor.get_file_rate() == replacer.get_file_rate()) {
    list<Chunk> reciptor_chunks = make_chunks(reciptor);
    list<Chunk> replacer_chunks = make_chunks(replacer);
    list<ChunkMatch>matches = get_matches(reciptor_chunks, replacer_chunks);
    
    std::string fn = "muxed.wav";
    std::cout << "outputting song\n";
    output_remix(matches, fn);
    //output_remix_alt(samples_per_slice, matches, "muxed_alt.wav");
    std::cout << "completed output\n";
  }
  else {
    //ERROR
    //would have to do weird stuff to work with this
    //because the different sampling rates would manifest in fft
    //would have the frequency bins of both spectrograms not line up
    std::cerr << "unable to mix files, sample rates are different, please use ffmpeg or alternative software to convert songs to same sample rate\n"; 
  }
}


list<Chunk> SongEmbedder::make_chunks(Song input) {
  std::cout << "Building spectrogram for " << input.get_file_name() << "\n";
  MatrixXf spec = input.spectrogram(samples_per_slice);
  foobar_spec(spec);
  
  std::cout << "filtering spectrogram\n";
  pair<MatrixXf, MatrixXf> xy_filt = extract_edges(spec);
  MatrixXf freq_edges = std::get<1>(xy_filt);
  
  snaz(freq_edges, edge_snazr);
 
  std::cout << "extracting chunks from spectrogram\n";
  MatrixXi freq_chunks = fill_chunkify(freq_edges, 0.0);

  std::cout << "culling chunks\n";
  ChunkStats freq_stats (freq_chunks);
  list<Chunk> part2 = freq_stats.cull_chunks(chunk_snazr, CHUNK_FREQ_OPT);

  list<Chunk> ret;
  ret.insert(ret.begin(), part2.begin(), part2.end());
  int tot_count = 0;
  tot_count += freq_stats.get_chunk_count();
  std::cout << tot_count << " many chunk(s) extracted from " << input.get_file_name() << "\n";
  std::cout << ret.size() << " many chunk(s) left after culling " << input.get_file_name() << "\n";

   string bn = input.get_file_name();
  string post = ".png";

  write_eigen_to_file(bn + ".freq_edges" + post, freq_edges);
  write_eigen_to_file(bn + ".spec" + post, spec);
  write_chunks_to_file(bn + ".freq_chunks" + post, freq_chunks);
  write_chunks_to_file(bn + ".freq_chunks_culled" + post, part2, freq_chunks.rows(), freq_chunks.cols());

  /* time chunks mostly pick out garbage sound, though some of it may be useful in other ways
     MatrixXf time_edges = std::get<0>(xy_filt);
     snaz(time_edges, edge_snazr);
     MatrixXi time_chunks = fill_chunkify(time_edges, 0.0);
     ChunkStats time_stats (time_chunks);
     list<Chunk> part1 = time_stats.cull_chunks(chunk_snazr, CHUNK_TIME_OPT); 
     write_eigen_to_file(bn + ".time_edges" + post, time_edges);
     write_chunks_to_file(bn + ".time_chunks" + post, time_chunks);
     write_chunks_to_file(bn + ".time_chunks_culled" + post, part1, time_chunks.rows(), time_chunks.cols());
  */
  return ret;
}

pair<MatrixXf, MatrixXf> extract_edges(MatrixXf &input) {
  //gauss col is a collunm vec, row is a row vec
  int gauss_col_dim = 5;
  int gauss_row_dim = 11;
  int dog_smooth_col_dim = 11;
  int dog_smooth_row_dim = 3;

  double div = 9.0;
  //the stddev is roughly set with div
  //only tested with gauss dims around ~3-100
  //going with higher divs seems to have the gradient fall of to zero sooner
  // div = 9 gives an approximate 1->0 gradient from center to edge
  // div = 4 gives an approximate 1->0.5 gradient from center to edge
  MatrixXf gauss_col = create_1d_gaussian_filter_col(gauss_col_dim, 1, ceil(gauss_col_dim * gauss_col_dim / div));
  MatrixXf gauss_row = create_1d_gaussian_filter_row(gauss_row_dim, 1, ceil(gauss_row_dim * gauss_row_dim / div));
  MatrixXf dog_smooth_row = create_1d_gaussian_filter_row(dog_smooth_row_dim, 1, ceil(dog_smooth_row_dim * dog_smooth_row_dim / div));
  MatrixXf dog_smooth_col = create_1d_gaussian_filter_col(dog_smooth_col_dim, 1, ceil(dog_smooth_col_dim * dog_smooth_col_dim / div));
  
  //x gauss is image blurred horizontally, y is image blurred vertically
  MatrixXf x_gauss = one_d_convolve(input, gauss_row);
  MatrixXf y_gauss = one_d_convolve(input, gauss_col);

  //extract difference of gaussian
  MatrixXf x_dog = input - x_gauss;
  MatrixXf y_dog = input - y_gauss;
  
  //smoothes edges so barely seperated chunks are joined
  //x_dog picked out vertical edges, so smooth vertically to join. vice versa for y and horizontal
  MatrixXf x_dog_smooth = one_d_convolve(x_dog, dog_smooth_col);
  MatrixXf y_dog_smooth = one_d_convolve(y_dog, dog_smooth_row);
  pair<MatrixXf, MatrixXf> ret(x_dog_smooth, y_dog_smooth);
  return ret;
}

//will assign adjacent non-zero cells the same chunk ID
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

list<ChunkMatch> SongEmbedder::get_matches(list<Chunk> &reciptor_chunks, list<Chunk> &replacer_chunks) {
  //weights for frequency center, frequency margin, and time marin
  ChunkCompare comp (50,1,10);

  std::cout << "matching sound chunks!\n";
  list<ChunkMatch> matches = quick_match(reciptor_chunks, replacer_chunks, comp);
  std::cout << "Done with matching!\n";
  if (match_snazr > 0) {
    std::cout << "culling bad matches\n";
    std::cout << "orig len is " << matches.size() << "\n";
    double avg = reverse_gsnaz(matches, chunkmatch_get_score_wrapper, match_snazr);
    std::cout << "only accepted scores <= " << avg << "\n";
    std::cout << "culled len is " << matches.size() << "\n";
  }

  return matches;
}

//iterates over all chunks for matches
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

//iterates over chunks which has at least one identical stat for matches
list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
  map<int, list<Chunk>> freq_cent_map;
  map<int, list<Chunk>> time_marg_map;
  list<ChunkMatch> matches;
  int fc_ind, tm_ind, temp;
  
  //populate maps
  list<Chunk> sub_list;
  for (Chunk a_chunk : rep) {
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();
    
    time_marg_map[tm_ind].push_front(a_chunk);
    freq_cent_map[fc_ind].push_front(a_chunk);
  }

  ChunkMatch a_match;
  map<int, int> closest_freq_cent;
  map<int, int> closest_time_marg;
  
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();

    if (closest_freq_cent.count(fc_ind) == 0) {
      temp = get_closest_bin_index(freq_cent_map, fc_ind);
      closest_freq_cent.insert(std::pair<int,int>(fc_ind, temp));
    }
    fc_ind = closest_freq_cent[fc_ind];

    if (closest_time_marg.count(tm_ind) == 0) {
      temp = get_closest_bin_index(time_marg_map, tm_ind);
      closest_time_marg.insert(std::pair<int,int>(tm_ind, temp));
    }
    tm_ind = closest_time_marg[fc_ind];
    
    sub_list = time_marg_map[tm_ind];
    //a_match.best_match_chunk(sub_list, comp);
    
    sub_list = freq_cent_map[fc_ind];
    a_match.best_match_chunk(sub_list, comp);
    
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

int get_closest_bin_index(map<int, list<Chunk>> bins, int start) {
  //various bins may have zero items in them
  //find the closest non-empty one
  bool done = false;
  int offset = 0;
  int ind = -1;
  while(!done) {
    if (bins.count(start + offset) == 1) {
      done = true;
      ind = start + offset;
    }
    else if (start - offset >= 0) {
      if (bins.count(start - offset) == 1) {
	done = true;
	ind = start - offset;
      }
    }
    offset++;
  }
  return ind;
}


struct cm_hash {
  size_t operator ()(const ChunkMatch & cm) const {
    return ChunkMatch::match_hash((ChunkMatch&)cm);
  }
};

void SongEmbedder::output_remix(list<ChunkMatch> &matches, std::string fn) {
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
  FileWvIn sample_src;
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  StkFrames temp (0, samples_per_slice, channels);

  LentPitShift lent(1, samples_per_slice);
  
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
  output.openFile(fn, channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  
  unordered_set<ChunkMatch,cm_hash> active_chunks; 
  for (int i = 0; i < new_slices; i++) {
    while(active_start.size() > 0 && active_start.front().get_active_start() <= i) {
      if (active_start.front().get_active_start() == i) {
	active_chunks.insert(active_start.front());
	active_start.pop_front();
      }
      else if (active_start.front().get_active_start() < 0) {
	//odd case, happens when song matches a chunk which is long enough and close enough to start
	//that it has a negative start time.
	//just tossing it away for now
	active_start.pop_front();
      }
    }   

    frame_out *= zero;
    for (ChunkMatch chunk_match : active_chunks) {
      temp *= zero;
      normal_filt(chunk_match, sample_src, frame_in, temp, i);
      //alt_filt(chunk_match, sample_src, temp, i, lent);

      if (chunk_match.get_active_start() == i) {
	temp *= reg_feather;
      }
      if (chunk_match.get_active_end() == i) {
	temp *= rev_feather;
      }
      
      frame_out += temp;
    }

    output.tick(frame_out);
    
    while(active_end.size() > 0 && active_end.front().get_active_end() == i) {
      active_chunks.erase(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}

void SongEmbedder::normal_filt(ChunkMatch& chunk_match, FileWvIn &sample_src, StkFrames &frame_in, StkFrames &out, int i) {
  Chunk match = chunk_match.get_match_chunk();
  ChunkFilter filt = match.get_filter();
  int offset = (i - chunk_match.get_active_start()) + match.get_time_start();
  offset *= samples_per_slice;
  sample_src.reset();
  sample_src.addTime(offset);
  
  sample_src.tick(frame_in);
  
  //filt.fir_filter_frame(frame_in, out);
  filt.fir_filter_frame(frame_in);
  out += frame_in;
}

void SongEmbedder::alt_filt(ChunkMatch& chunk_match, FileWvIn &sample_src, StkFrames &out, int i, LentPitShift &lent)  {
  int channels = sample_src.channelsOut();;
  Chunk match = chunk_match.get_match_chunk();
  Chunk orig = chunk_match.get_orig_chunk();
  ChunkFilter filt = match.get_filter();
  double offset_ratio =  (double)(i - orig.get_time_start()) / (orig.get_time_length());
  int offset = (offset_ratio * match.get_time_length() + match.get_time_start()) * samples_per_slice;
  sample_src.reset();
  sample_src.addTime(offset);

  double compression_ratio = (double)orig.get_time_length() / match.get_time_length();
  int size = samples_per_slice / compression_ratio;
  StkFrames temp_in (size, channels);
  sample_src.tick(temp_in);
     
  filt.fir_filter_frame(temp_in);
  reshape_chunk(temp_in, out, orig, match, lent);
}


void test_sound_filter(int samples_per_slice, std::string fn) {
  //setup 3 filters which filter out low, middle, and high frequencies. 
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

  std::string filt_test_fn = "filt_test" + fn;
  output.openFile(filt_test_fn, channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
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

  Song filt_song(filt_test_fn);
  Song raw_song(fn);
  MatrixXf filt_spec = filt_song.spectrogram(samples_per_slice);
  MatrixXf raw_spec = raw_song.spectrogram(samples_per_slice);
  write_eigen_to_file("raw_spec.png", raw_spec);
  write_eigen_to_file("filt_spec.png", filt_spec);
}


void test_center_grouping(MatrixXi &song_chunks) {
  printf("starting group test\n");
  ChunkStats stats (song_chunks);
  list<Chunk> clist = stats.cull_chunks(0, CHUNK_FREQ_OPT);
  list<Chunk> groups = group_chunks(clist);

  write_chunks_to_file("center_groups_in.png", clist, song_chunks.rows(), song_chunks.cols());
  write_chunks_to_file("center_groups_out.png", groups, song_chunks.rows(), song_chunks.cols());
  printf("done with grouping test\n");
}

list<Chunk> group_chunks(list<Chunk> &chunks) {
  list<Chunk> temp (chunks);
  temp.sort(Chunk::comp_time_center);
  double center = std::numeric_limits<double>::infinity();
  double center_sum = 0;
  double center_count = 0;
  //instead of a constant, could be the average distance between adjacent elements in center_sorted list
  double tolerance = 5;
  int group_id = 0;
  for (Chunk &aChunk : temp) {
    if ( std::abs(aChunk.get_time_center() - center) > tolerance) {
      //new group
      group_id++;
      center_sum = 0;
      center_count = 0;
    }
    //assign to current group
    center_sum += aChunk.get_time_center();
    center_count++;
    aChunk.set_chunk_id(group_id);
    center = center_sum / center_count;
  }
  return temp;
  
}
