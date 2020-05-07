#include <math.h>
#include <vector>
#include <string>
#include <LentPitShift.h>
#include "util.h"
#include "kiss_fftr.h"
#include "FIRFilterCode.h"
#include "eigen_to_image.h"
#include "Song.h"
#include "SongEmbedder.h"

using Eigen::MatrixXf;
using Eigen::MatrixXi;
using std::max;
using std::list;
using std::string;
using namespace stk;

static list<Chunk> time_grouper(list<Chunk>& chunks);
static list<ChunkGroup> freq_grouper(list<Chunk>& chunks);
static double chunk_ratio(Chunk& a, Chunk& b);
static list<ChunkGroup> build_freq_groups(list<Chunk>& group);
static void chain_build(list<Chunk>& group, list<Chunk*>& chain, int pos, double ratio, double tolerance);


MatrixXf TFD_extract(FileWvIn &input, int total_slices, int samples_per_slice, int channels) {
  //int nfft = samples_per_slice;
  bool is_inverse = false;
  //number of bins fft puts things into. apparently it's the same as nfft.
  //since nothing is placed in the upper half, divide by 2
  //though if you want to play around with inverse fft, set to nfft
  //int freq_chunks = nfft / 2;
  //freq_chunks = samples_per_slice;
  kiss_fftr_cfg cfg = kiss_fftr_alloc( samples_per_slice ,is_inverse,0,0 );
  size_t input_size = sizeof(kiss_fft_scalar) * (samples_per_slice + 1);
  size_t output_size = sizeof(kiss_fft_cpx) * (samples_per_slice + 1);
  kiss_fft_scalar* fft_input =(kiss_fft_scalar*)KISS_FFT_MALLOC(input_size);
  kiss_fft_cpx* fft_output =(kiss_fft_cpx*)KISS_FFT_MALLOC(output_size);
  StkFrames slice(samples_per_slice, channels);
  StkFrames single_chan(samples_per_slice, 1);

  for (int i = 0 ; i < samples_per_slice; i++) {
    fft_input[i] = 0;
    fft_output[i] = (kiss_fft_cpx){.r = 0, .i = 0};
  }
  int rows = samples_per_slice / 2;
  int cols = total_slices;
  MatrixXf tfd = MatrixXf(rows, cols); 
  tfd.setZero();
  
  for (int slice_i = 0; slice_i < cols; slice_i++) {
    input.tick(slice);
    for (int i = 0; i < channels; i++) {
      slice.getChannel(i, single_chan, 0);
      stk_to_kiss(single_chan, fft_input);
      kiss_fftr(cfg, fft_input, fft_output);
      for (int freq_i = 0; freq_i < rows ; freq_i++) {
	kiss_fft_cpx a = fft_output[freq_i];
	double modulus =  pow(pow(a.r, 2) + pow(a.i, 2), 2);
	tfd(freq_i, slice_i) += modulus;
      }
    }
  }
  KISS_FFT_FREE(fft_input);
  KISS_FFT_FREE(fft_output);
  kiss_fftr_free(cfg);
  return tfd;
}


//data transform because my base spec output looked weird
//used some math that foobar uses for it's visualization.
/* 
   original code found at g_preprocess_chunk in 
   https://github.com/stengerh/foo_vis_spectrogram/blob/master/foo_vis_spectrogram/CSpectrogramView.cpp
 */
void foobar_spec(MatrixXf &tfd) {
  float temp = -1;
  float pow = -1;
  float minf = tfd.minCoeff();
  float maxf = tfd.maxCoeff();
  float low = std::max(20 * log10(minf), -800.f);
  float high = 20 * log10(maxf);
  float div = 1.0 / (high - low);
  for (int col_i = 0; col_i < tfd.cols(); col_i++) {
    for (int row_i = 0; row_i < tfd.rows(); row_i++) {
      temp = tfd(row_i, col_i);
      pow = 20.0 * log10(temp);
      temp = (pow - low) * div;
      temp = std::min(1.0f, std::max(0.0f, temp));
      tfd(row_i, col_i) = temp;
    }
  }
}

void kiss_to_stk(kiss_fft_scalar* in, StkFrames &out) {
  int size = out.frames();
  for(int i = 0; i < size; i++) {
    out[i] = in[i];
  }
}

void stk_to_kiss(StkFrames &in, kiss_fft_scalar* out) {
  int size = in.frames();
  for(int i = 0; i < size; i++) {
    out[i] = in[i];
  }
}

void fetch_frame(int frame_n, int samples_per_slice, FileWvIn& wav, StkFrames& out) {
  int offset = frame_n * samples_per_slice;
  wav.reset();
  wav.addTime(offset);
  wav.tick(out);
}

void filter_frame(Chunk& c, StkFrames& in, StkFrames& out) {
  ChunkFilter filt = c.get_filter();
  //filt.fir_filter_frame(in, out);
  filt.fir_filter_frame(in);
  out += in;
}


void resample_frame(StkFrames &cur_frames,StkFrames &new_frames) {
  int cur_len = cur_frames.frames();
  int new_len = new_frames.frames();
  float cur_i = -1;
  for (int new_i = 0; new_i < new_len; new_i++) {
    cur_i = cur_len * (double)new_i / new_len;
    new_frames[new_i] = cur_frames.interpolate(cur_i);
  }
}

void reverse_frame(StkFrames &cur_frames, StkFrames &rev_frames) {
  int rev_i = -1;
  int channels = cur_frames.channels();
  int len = cur_frames.frames();
  for (int i = 0; i < len; i++) {
    rev_i = len - 1 - i;
    for (int chan = 0; chan < channels; chan++) {
      rev_frames(rev_i, chan) = cur_frames(i,chan);
    }
  }
}

void reshape_chunk(StkFrames &frame_in, StkFrames &frame_out, Chunk &src, Chunk &shape, LentPitShift &lent) {
  //reshape in time/frequency space to match another chunks shape
  //change the sampling rate first to fit the time
  //calculate how resampling changed pitch, then make/use a pitch shifter to correct to desired pitch center
  
  int channels = frame_in.channels();
  int orig_size = frame_in.frames();
  int new_size = frame_out.frames();

  StkFrames temp_in (orig_size, 1);
  StkFrames temp_out (new_size, 1);
  StkFrames garb (new_size, 1);
  StkFrames reshape (new_size, channels);

  //forgot about my experiments with lentPitShifts
  //there might be a delay before output starts getting output
  //could either try repeating input, handing in a previouse frame, or see what ticking in zeros does

  double resample_freq_scale = (double)orig_size / new_size;
  double resampled_freq = src.get_freq_center() * resample_freq_scale;
  double pitch_correct = shape.get_freq_center() / resampled_freq;

  lent.setShift(pitch_correct);
  //std::cout << "shift is " << pitch_correct << "\n";
  //first step, resample at a different rate
  for (int chan_i = 0; chan_i < channels; chan_i++) {
    frame_in.getChannel(chan_i, temp_in, 0);
    resample_frame(temp_in, temp_out);
    lent.tick(temp_out);
    //print_frame(temp_out);
    reshape.setChannel(chan_i, temp_out, 0);
  }
  frame_out += reshape;
}

//creates a 1-dimensional gaussian filter, as an stk finite impulse response filter
//  int filter_size = 20;
//  Fir gaussian_filt = create_1d_gaussian_filter(20, 1, filter_size / 2, 1);
stk::Fir create_1d_gaussian_filter(int length, double amplitude, double center, double stddev) {
  std::vector<StkFloat> kernel(0);
  StkFloat val = 0;
  double exp = 0;
  for (int i = 0; i < length; i++) {
    exp = -1 * pow(i - center, 2) / stddev;
    val = amplitude * pow(M_E, exp);
    kernel.insert(kernel.begin(), val);
  }
  Fir ret (kernel);
  return ret;
}

MatrixXf create_1d_gaussian_filter_col(int length, double amplitude, double stddev) {
  MatrixXf kern (length, 1);
  float val = 0;
  double exp = 0;
  double center = (length - 1) / 2;
  for (int i = 0; i < length; i++) {
    exp = -1 * pow(i - center, 2) / stddev;
    val = amplitude * pow(M_E, exp);
    kern(i,0) = val;
  }
  return kern;;
}

MatrixXf create_1d_gaussian_filter_row(int length, double amplitude, double stddev) {
  return create_1d_gaussian_filter_col(length, amplitude, stddev).transpose();
}


void make_fir_bandpass_filter(double* coefs, int taps, int freq_center, int freq_margin,  int max_freq) {
  double rel_center = (double)freq_center / max_freq;
  double rel_margin = (double)freq_margin / max_freq;
  TFIRPassTypes band = firBPF;
  TWindowType wt = wtKAISER;
  double beta = 5;
  RectWinFIR(coefs, taps, band, rel_center, rel_margin);
  FIRFilterWindow(coefs, taps, wt, beta);
}

stk::Fir create_fir_from_coefs(double* coefs, int len) {
  std::vector<StkFloat> kernel(0);
  StkFloat val = 0;
  for (int i = 0; i < len; i++) {
    val = coefs[i];
    kernel.insert(kernel.begin(), val);
  }
  stk::Fir ret(kernel);
  return ret;
}

MatrixXf one_d_convolve_valid(MatrixXf &mat, MatrixXf &kern) {
  //perform valid 1d convolution on matrix
  //or cross-corelation, planning on using symetric kernels so doesn't matter
  int res_dim;
  int kern_dim = -1;
  MatrixXf temp, result;
  int d_start, d_end, d_width;
  
  if (kern.rows() == 1) {
    result = MatrixXf(mat.rows(), mat.cols() - (kern.cols() - 1));
    result.setZero();
    res_dim = result.cols();
    d_width = mat.rows();
    kern_dim = kern.cols();
    for (int i = 0; i < kern_dim; i++) {
      d_start = i;
      d_end = res_dim;
      temp = mat.block(0,d_start, d_width, d_end);
      result += temp * kern(0,i);
    }
  }
  else if (kern.cols() == 1) {
    result = MatrixXf(mat.rows() - (kern.rows() - 1), mat.cols());
    result.setZero();
    res_dim = result.rows();
    d_width = mat.cols();
    kern_dim = kern.rows();
    for (int i = 0; i < kern_dim; i++) {
      d_start = i;
      d_end = res_dim;
      temp = mat.block(d_start, 0, d_end, d_width);
      result += temp * kern(i, 0);
    }
  }
  else {
    //error, kernel is not 1d
  }
  if (kern_dim % 2 == 0) {
    //kernel has even dim and doesn't have a perfect center index, funky stuff may happen
    std::cerr << "Warning, one-d-convolving with an even-length kernel\n";
  }
  return result;
}

MatrixXf one_d_convolve(MatrixXf &mat, MatrixXf &kern) {
  //perform same 1d convolution on matrix
  int kern_dim = -1;
  MatrixXf temp, ref, result;
  //reflecting edges
  int pad = -1;
  if (kern.rows() == 1) {
    pad  = kern.cols() / 2;
  }
  else if (kern.cols() == 1) {
    pad  = kern.rows() / 2;
  }
  int row_i = 0;
  int col_i = 0;
  int low_row_bound = pad;
  int high_row_bound = mat.rows() + pad - 1;
  int low_col_bound = pad;
  int high_col_bound = mat.cols() + pad - 1;
  int ref_row = -1;
  int ref_col = -1;
  float sum = kern.sum();
  result = MatrixXf(mat.rows(), mat.cols());
  result.setZero();
  ref = MatrixXf(mat.rows() + 2 * pad, mat.cols() + 2 * pad);
  ref.setZero();
  bool col_region = false;
  bool row_region = false;
  //reflect edges
  row_i = 0;
  while(row_i < ref.rows()) {
    if (row_i < low_row_bound || row_i > high_row_bound) {
      //inside pad boundary
      //pad things
      ref_row = pad - row_i;
      if (ref_row < 0) {
	ref_row = 2 * mat.rows() - (row_i - pad) - 1;
      }
      row_region = true;
    }
    else {
      ref_row = row_i - pad;
      row_region = false;
    }
    col_i = 0;
    while(col_i < ref.cols()) {
      if (col_i < low_col_bound || col_i > high_col_bound) {
	//inside pad boundary
	//pad things
	ref_col = pad - col_i;
	if (ref_col < 0) {
	  ref_col = 2 * mat.cols() - (col_i - pad) - 1;
	}
	col_region = true;
      }
      else {
	ref_col = col_i - pad;
	col_region = false;
      }
      if (col_region || row_region) {
	ref(row_i, col_i) = mat(ref_row, ref_col);
	col_i++;
      }
      else {
	//not reflecting anything, skip some columns
	col_i = high_col_bound + 1;
      }
    }
     row_i++;
  }
  //fill center
  ref.block(pad, pad, mat.rows(), mat.cols()) = mat;

  if (kern.rows() == 1) {
    kern_dim = kern.cols();
    for (int i = 0; i < kern_dim; i++) {
      temp = ref.block(pad,i,mat.rows(), mat.cols());
      result += temp * (kern(0,i) / sum);
    }
  }
  else if (kern.cols() == 1) {
    kern_dim = kern.rows();
    for (int i = 0; i < kern_dim; i++) {
      temp = ref.block(i,pad,mat.rows(), mat.cols());
      result += temp * (kern(i, 0) / sum);
    }
  }
  else {
    std::cerr << "Warning, one-d-convolving with a 2d kernel\n";
  }
  if (kern_dim % 2 == 0) {
    //kernel has even dim and doesn't have a perfect center index, funky stuff may happen
    std::cerr << "Warning, one-d-convolving with an even-length kernel\n";
  }
  return result;
}

//difference of gaussian
//only really a function because gauss can be smaller
//if using regular/same convolution, no need to use this
MatrixXf dog(MatrixXf &full, MatrixXf &gauss) {
  int f_rows = full.rows();
  int f_cols = full.cols();
  int g_rows = gauss.rows();
  int g_cols = gauss.cols();
  MatrixXf ret;
  if (f_rows == g_rows && f_cols == g_cols) {
    ret = full - gauss;
  }
  else {
    MatrixXf temp = full.block((f_rows - g_rows) / 2, (f_cols - g_cols) / 2, g_rows, g_cols);
    ret = temp - gauss;
  }
  return ret.cwiseAbs();
}

//stands for sub nonzero-average zeroing
//take average of all non-zero elements, zero anything below average
//takes matrix and rounds to apply snaz, zeros values in matrix
//returns average of penultimate round
//having snazr as anything above 1 is typically excessive, consider at most using 2
double snaz(MatrixXf &filt, int snazr) {
  int freq_range = filt.rows();
  int time_range = filt.cols();
  double nz_tot = 0;
  double nz_avg = 0.0;
  int nz_count = 0;
  float val = 0;
  for (int i = 0; i <= snazr; i++) {
    nz_tot = 0;
    nz_count = 0;
    //don't want to include zeros in average, can't use builtin sum methods, have to traverse :(
    //combined two tasks in one, do the sub_average zeroing then compute new average
    //need to to 1 final loop after outermost for loop finishes to do final sub_average zeroing
    //accomplished by using i <= snazr and only recomputing nz_avg if i < snazr
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      for (int time_i = 0; time_i < time_range; time_i++) {
	val = filt(freq_i, time_i);
	if (val <= nz_avg) {
	  filt(freq_i, time_i) = 0.0;
	}
	else {
	  if (i < snazr) {
	    nz_tot += val;
	    nz_count++;
	  }
	}
      }
    }
    if (i < snazr) {
      nz_avg = nz_tot / nz_count;
      //std::cout << "non-zero average for round " << i << " is " << nz_avg << "\n";
    }
  }
  return nz_avg;
}

double snaz(list<Chunk> &many_chunks, int snazr) {
  double nz_tot = 0;
  double nz_avg = 0.0;
  double val = 0;
  int nz_count = 0;
  list<Chunk> temp1 (many_chunks);
  list<Chunk> temp2;
  for (int i = 0; i <= snazr; i++) {
    nz_tot = 0;
    nz_count = 0;
    for (Chunk a_chunk : temp1) {
      val = a_chunk.get_chunk_size();
      if (val > nz_avg) {
	nz_tot += val;
	nz_count++;
	temp2.push_front(a_chunk);
      }
    }
    if (i < snazr) {
      nz_avg = nz_tot / nz_count;
      temp1.assign(temp2.begin(), temp2.end());
      temp2.erase(temp2.begin(), temp2.end());
    }
  }
  many_chunks.assign(temp2.begin(), temp2.end());
  return nz_avg;
}

/*
grouping
2 stages.
first, group by time by having tome tolerance between chunks (avg dist beween chunks)
if a chunk is within tolerance to anything in a group then add chunk to group
could also use average distance between chunk and everything in group

second, take groups  from stage 1 and split each group in more based groups on frequency
have some comparision between chunks ( either a scalar or an offset)
then search through group making other subgroups where items in group have roughly the same comparision with their next element

*/
void output_groups(list<Chunk>& chunks, int samples_per_slice, string fn) {
  int curr_id = 0;
  list<Chunk> group;
  list<ChunkMatch> matches;
  Song test_song(fn);
  chunks.sort(Chunk::comp_chunk_id);
  SongEmbedder output (test_song, test_song, samples_per_slice, -1,-1,-1);
  for (Chunk& c : chunks) {
    if (c.get_chunk_id() != curr_id) {
      if (group.size() > 0) {
	matches = ChunkMatch::self_match(group);
	output.output_remix(matches, "group" + std::to_string(curr_id));
	matches.clear();	
      }
      //end
      group.clear();
      curr_id = c.get_chunk_id();
    }
    else {
      group.push_front(c);
    }
  }
  int sample_length =test_song. get_samples_per_channel();
  int output_rows = samples_per_slice / 2;
  int output_cols = ceil((double)sample_length / samples_per_slice);
  write_chunks_to_file("center_groups.png", chunks, output_rows, output_cols);
}

list<ChunkGroup> group_chunks(list<Chunk>& pool) {
  list<Chunk> temp = time_grouper(pool);
  list<ChunkGroup> ret = freq_grouper(temp);
  return ret;
}

list<Chunk> time_grouper(list<Chunk>& chunks) {
  list<Chunk> temp (chunks);
  temp.sort(Chunk::comp_time_center);
  double center = std::numeric_limits<double>::infinity();
  double center_sum = 0;
  double center_count = 0;
  //instead of a constant, could be the average distance between adjacent elements in center_sorted list
  double total_dist = 0;
  Chunk* prev = NULL;
  for (Chunk& c: temp) {
    if (prev != NULL) {
      total_dist += abs(c.get_time_center() - prev->get_time_center());
    }
    prev = &c;
  }
  double tolerance = 5;
  tolerance = total_dist / temp.size();
  int group_id = 0;
  for (Chunk& aChunk : temp) {
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
  printf("Created %d groups out of %ld chunks\n", group_id, chunks.size());
  return temp;
}

//testing for further frequency grouping
//for chunks within the same group
//split chunks into groups where the fn * group_freq = fn+1
list<ChunkGroup> freq_grouper(list<Chunk>& chunks) {
  int max_id = -1;
  for (Chunk& c : chunks) {
    max_id = max(c.get_chunk_id(), max_id);
  }
  
  chunks.sort(Chunk::comp_chunk_id);

  int curr_id = 0;
  list<Chunk> group;
  list<ChunkGroup> ret;
  for (Chunk& c : chunks) {
    if (c.get_chunk_id() != curr_id) {
      //start of new time chunk
      //do something with previouse, completed chunks
      if (group.size() > 1) {
	ret.splice(ret.end(), build_freq_groups(group));
      }
      //end
      group.clear();
      curr_id = c.get_chunk_id();
    }
    else {
      group.push_front(c);
    }
  }
  return ret;
}

double chunk_ratio(Chunk& a, Chunk& b) {
  return (b.get_freq_center() + 0.5) / (a.get_freq_center() + 0.5);
}

static int chain_id = 1;
list<ChunkGroup> build_freq_groups(list<Chunk>& group) {
  printf("Group is %ld big \n", group.size());
  group.sort(Chunk::comp_chunk_freq);

  int base_count = 0;
  int memb_count = 0;
  double rat;
  double tolerance = 0.05;
  list<Chunk> sub_group;
  list<Chunk*> chain;
  list<ChunkGroup> ret;
  for (Chunk& base : group) {
    for (Chunk& memb: group) {
      if (memb_count > base_count) {
	chain.push_back(&base);
	chain.push_back(&memb);
	rat = chunk_ratio(base, memb);
	chain_build(group, chain, memb_count, rat, tolerance);
	if (chain.size() > 2) {
	  printf("Built a chain of size %ld \n", chain.size());
	  printf("Ratio for chain is %f\n", rat);
	  for (Chunk *link : chain) {
	    Chunk temp = *link;
	    //need a way of making chian id's
	    temp.set_chunk_id(chain_id);
	    sub_group.push_back(temp);
	  }
	  ret.push_back(ChunkGroup(sub_group, rat));
	  chain_id++;
	  sub_group.clear();
	}
	chain.clear();
      }
      memb_count++;
    }
    base_count++;
  }
  return ret;
}

//continues building a chain  where f(n+1)/f(n) = ratio
//tolerance is a percentile tolerance for the ratio
//builds up multiple continuations for the original chain than only returns the largest chain
void chain_build(list<Chunk>& group, list<Chunk*>& chain, int pos, double ratio, double tolerance) {
  int idx = 0;
  Chunk prev_chunk;
  double cur_ratio = 0;
  unsigned long max = chain.size();
  list<Chunk*> ret = chain;
  list<list<Chunk*>> collec;
  for (Chunk& c: group) {
    if (idx >= pos) {
      if (idx == pos) {
	prev_chunk = c;
      }
      else {
	cur_ratio = chunk_ratio(prev_chunk, c);
	if ( cur_ratio > ratio * (1 + tolerance)) {
	  break;
	}
	if ( abs(cur_ratio - ratio) < ratio * tolerance ) {
	  list<Chunk*> aChain = chain;
	  aChain.push_back(&c);
	  chain_build(group, aChain, idx, ratio, tolerance);
	  collec.push_back(aChain);
	}
      }
    }
    idx++;
  }
  for (list<Chunk*> aChain : collec) {
    if (aChain.size() > max) {
      ret = aChain;
      max = aChain.size();
    }
  }
  chain = ret;
}


void print_frame(StkFrames in) {
  int len = in.frames();
  for (int i = 0; i < len; i++) {
    std::cout << in[i] << " ";
  }
  std::cout << "\n";
}


double abs_log_diff(double a, double b) {
  if (a == 0 || b == 0) {
    a += 0.1;
    b += 0.1;
  }
  return abs(log2(a) - log2(b));
}

double abs_diff(double a, double b) {
  return abs(a - b);
}

double square_diff(double a, double b) {
  return pow(a - b, 2);
}

