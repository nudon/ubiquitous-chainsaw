#include <cmath>
#include <FileRead.h>
#include <FileWvIn.h>
#include "Song.h"
#include "ChunkStats.h"
#include "util.h"
#include "eigen_to_image.h"
#include "SongEmbedder.h"
#include "ChunkGroup.h"

using Eigen::MatrixXf;
using Eigen::MatrixXi;
using std::list;
using std::pair;
using std::string;
using namespace stk;

pair<MatrixXf, MatrixXf> extract_edges(MatrixXf &input);
MatrixXi fill_chunkify(MatrixXf &input, int thresh);
void recurse_fill(MatrixXf &input, MatrixXi &ids, int x, int y, int cid);

Song::Song(string fn) {
  file_name = fn;
  FileRead song_file;
  FileWvIn song_wv;
  song_wv.openFile(fn);
  channels = song_wv.channelsOut();
  song_file.open(fn, channels);
  samples = song_file.fileSize();
  file_rate = song_file.fileRate();
  song_file.close();
  song_wv.closeFile();
}

MatrixXf Song::spectrogram(int samples_per_slice) {
  //just setup and apply fft to song
  FileWvIn input;
  input.openFile(get_file_name());
  int mono_sample_length = get_tot_samples() / get_channels();
  int total_slices = ceil(mono_sample_length / samples_per_slice);
  MatrixXf spec = TFD_extract(input,total_slices, samples_per_slice, get_channels());
  input.closeFile();
  return spec;    
}

list<Chunk> Song::make_chunks(int samples_per_slice, int edge_snazr, int chunk_snazr) {
  std::cout << "Building spectrogram for " << get_file_name() << "\n";
  MatrixXf spec = spectrogram(samples_per_slice);
  foobar_spec(spec);
  
  std::cout << "filtering spectrogram\n";
  pair<MatrixXf, MatrixXf> xy_filt = extract_edges(spec);
  MatrixXf freq_edges = std::get<1>(xy_filt);
  
  snaz(freq_edges, edge_snazr);
 
  std::cout << "extracting chunks from spectrogram\n";
  MatrixXi freq_chunks = fill_chunkify(freq_edges, 0.0);

  std::cout << "culling chunks\n";
  ChunkStats freq_stats (freq_chunks);
  list<Chunk> freq_list = freq_stats.cull_chunks(chunk_snazr, CHUNK_FREQ_OPT);

  int tot_count = freq_stats.get_chunk_count();
  std::cout << tot_count << " many chunk(s) extracted from " << get_file_name() << "\n";
  std::cout << freq_list.size() << " many chunk(s) left after culling " << get_file_name() << "\n";

  string bn = get_file_name();
  string post = ".png";
  
  write_eigen_to_file(bn + ".freq_edges" + post, freq_edges);
  write_eigen_to_file(bn + ".spec" + post, spec);
  write_chunks_to_file(bn + ".freq_chunks" + post, freq_chunks);
  write_chunks_to_file(bn + ".freq_chunks_culled" + post, freq_list, freq_chunks.rows(), freq_chunks.cols());
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
  return freq_list;
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

void test_sound_filter(int samples_per_slice, string fn) {
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

  string filt_test_fn = "filt_test" + fn;
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

void Song::test_chunk_grouping(int samples_per_slice, int edge_snazr, int chunk_snazr) {
  list<Chunk> chunks = make_chunks(samples_per_slice, edge_snazr, chunk_snazr);

  printf("starting group test\n");
  //groups based on time
  //further seperates groups based on frequency
  list<ChunkGroup> cel = group_chunks(chunks);
  
  printf("yeah %lu\n", cel.size());
  printf("done with grouping test\n");
}
