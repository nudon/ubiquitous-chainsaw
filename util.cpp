#include <math.h>
#include <vector>
#include <LentPitShift.h>
#include "util.h"
#include "kiss_fftr.h"
#include "FIRFilterCode.h"

using Eigen::MatrixXf;
using Eigen::MatrixXi;
using std::list;
using namespace stk;

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

void resample_frame(StkFrames &cur_frames,StkFrames &new_frames) {
  int cur_len = cur_frames.frames();
  int new_len = new_frames.frames();
  float cur_i = -1;
  for (int new_i = 0; new_i < new_len; new_i++) {
    cur_i = cur_len * (double)new_i / new_len;
    new_frames[new_i] = cur_frames.interpolate(cur_i);
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
  //also hand this in as an arguement since construction will be same across all func calls
  //int tMax = new_size;
  //LentPitShift lent(1, tMax);
  //LentPitShift lent = src.getPitShift();

  double resample_freq_scale = (double)orig_size / new_size;
  double resampled_freq = src.get_freq_center() * resample_freq_scale;
  double pitch_correct = shape.get_freq_center() / resampled_freq;

  lent.setShift(pitch_correct);
  //std::cout << "shift is " << pitch_correct << "\n";
  //first step, resample at a different rate
  for (int chan_i = 0; chan_i < channels; chan_i++) {
    frame_in.getChannel(chan_i, temp_in, 0);
    resample_frame(temp_in, temp_out);
    //lent.tick(temp_out, garb);
    //lent.tick(temp_out, garb);
    lent.tick(temp_out);
    //std::cout << "size is " << temp_out.frames() << "\n";
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
  MatrixXf kern (1, length);
  float val = 0;
  double exp = 0;
  double center = (double)length / 2;
  for (int i = 0; i < length; i++) {
    exp = -1 * pow(i - center, 2) / stddev;
    val = amplitude * pow(M_E, exp);
    kern(0,i) = val;
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
  return Fir(kernel);
}

MatrixXf one_d_convolve(MatrixXf &mat, MatrixXf &kern) {
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
  }
  return result;
}

//difference of gaussian
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
  //may not need this but worried about breaking iterators
  //std::list<Chunk> temp;
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



void print_frame(StkFrames in) {
  int len = in.frames();
  for (int i = 0; i < len; i++) {
    std::cout << in[i] << " ";
  }
  std::cout << "\n";
}
