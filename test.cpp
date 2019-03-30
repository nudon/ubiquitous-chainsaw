#include <stdio.h>
#include <FileLoop.h>
#include <FileWvOut.h>
#include <Fir.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

#include "kiss_fftr.h"
#include "FIRFilterCode.h"
using namespace stk;
using Eigen::MatrixXf;


void process_song(std::string fn_in, std::string fn_out);
MatrixXf TFD_extract(FileWvIn input, int nfft, int total_slices, int samples_per_slice, int channels);
void TFDI_extract(kiss_fft_cpx* fft_input, int nfft, int samples_per_slice, StkFrames output);

stk::Fir create_1d_gaussian_filter(int length, double amplitude, double center, double stddev);

void make_fir_bandpass_filter(double* coefs, int taps, int freq_center, int freq_margin,  int max_freq);
stk::Fir create_fir_from_coefs(double* coefs, int len);


MatrixXf create_1d_gaussian_filter_col(int length, double amplitude, double stddev);

MatrixXf create_1d_gaussian_filter_row(int length, double amplitude, double stddev);
MatrixXf one_d_convolve(MatrixXf mat, MatrixXf kern);


int main() {
  std::string file_in = "./sound/tsu.wav";
  std::string file_out = "./pet.wav";

  process_song(file_in, file_out);

  return 0;
}


void process_song(std::string fn_in, std::string fn_out) {
  FileRead song;
  FileWvIn input;
  FileWvOut output;
  fn_in = "./sound/tsu.wav";
  fn_out = "pet.wav";
  song.open(fn_in, 2);
  input.openFile(fn_in);
  output.openFile(fn_out, 2, FileWrite::FILE_WAV, Stk::STK_SINT16);

  int channels = input.channelsOut();
  long samples = song.fileSize();
  long samples_length = samples / channels;
  double file_rate = input.getFileRate();
  double len_in_sec = samples_length / file_rate;
  
  
  
  int samples_per_slice = 1024;
  double slices_per_second = file_rate / samples_per_slice;
  int total_slices = ceil(samples_length / samples_per_slice);

  int max_freq = file_rate / 2;
  //kiss_fft mentions that the fftr only populates nfft/2 + 1 things 
  //guessing that means I initialize with nfft, then only look at nfft/2 + 1 results?
  //actualyl nfft should just be length of sample. is also the resulting frequency resolution.


  bool filter_test = false;
  if (filter_test) {
    FileWvOut filteredOut;
    filteredOut.openFile("hi.wav", 2, FileWrite::FILE_WAV, Stk::STK_SINT16);
  
    StkFrames slice(samples_per_slice, channels);
    int taps = 60;
    double coefs[taps];
    make_fir_bandpass_filter(coefs, taps, max_freq / 2, max_freq / 64, max_freq);
    stk::Fir bp = create_fir_from_coefs(coefs, taps);
    input.reset();

    for (int slice_i = 0; slice_i < total_slices; slice_i++) {
      input.tick(slice);
      for (int chan_i = 0; chan_i < channels; chan_i++) {
	bp.tick(slice, chan_i);
      }
      filteredOut.tick( slice );
    }
    input.reset();
  }


  int nfft = samples_per_slice;

  //generate time-frequency distribution
  MatrixXf tfd = TFD_extract(input, nfft, total_slices, samples_per_slice, channels);

  //test out convolution code
  MatrixXf guass_col = create_1d_gaussian_filter_col(5, 1, 1);
  MatrixXf guass_row = create_1d_gaussian_filter_row(5, 1, 1);
  MatrixXf result = one_d_convolve(tfd, guass_col);
  result = one_d_convolve(result, guass_row);
  
  /// close files
  song.close();
  input.closeFile();
  output.closeFile();
}


//extract the time-frequency distribution of a song
//going to have it return eigen 
MatrixXf TFD_extract(FileWvIn input, int nfft, int total_slices, int samples_per_slice, int channels) {
  bool is_inverse = false;
  //numbe of binds fft puts things into. apparently it's the same as nfft. 
  int freq_chunks = nfft;
  kiss_fftr_cfg cfg = kiss_fftr_alloc( nfft ,is_inverse,0,0 );
  size_t input_size = sizeof(kiss_fft_scalar) * samples_per_slice;
  size_t output_size = sizeof(kiss_fft_cpx) * freq_chunks;
  kiss_fft_scalar* fft_input =(kiss_fft_scalar*)KISS_FFT_MALLOC(input_size);
  kiss_fft_cpx* fft_output =(kiss_fft_cpx*)KISS_FFT_MALLOC(output_size);
  bool print = false;
  StkFrames slice(samples_per_slice, channels);
  StkFrames single_chan(samples_per_slice, 1);

  MatrixXf tfd = MatrixXf(freq_chunks, total_slices); 
  tfd.Zero(freq_chunks, total_slices);
  for (int slice_i = 0; slice_i < total_slices; slice_i++) {
    input.tick(slice);
    for (int i = 0; i < channels; i++) {
      slice.getChannel(i, single_chan, 0);
      for (int i = 0; i < samples_per_slice; i++) {
	fft_input[i] = single_chan[i];
      }
      kiss_fftr(cfg, fft_input, fft_output);
      for (int freq_i = 0; freq_i < freq_chunks; freq_i++) {
	kiss_fft_cpx a = fft_output[freq_i];
	double modulus =  pow(a.r, 2) + pow(a.i, 2);
	if (print) {
	  std::cout << modulus << " ";
	  }
	//plus equals to get both channels
	tfd(freq_i, slice_i) += modulus;
      }
      if (print) {
	std::cout << "\n";
      }
    }
  }
  kiss_fftr_free(cfg);
  return tfd;
}

//do the inverse of the TFD, which I guess just returns the original input amplitudes?
//call this inside the regular TFD extract, hand in cpx, output some scalar value into some passed in stk-frame output, then in tfd, tick that to some output file
void TFDI_extract(kiss_fft_cpx* fft_input, int nfft, int samples_per_slice, StkFrames output) {
  bool print = true;
  bool is_inverse = true;
  int freq_chunks = nfft / 2 + 1;
  kiss_fftr_cfg cfg = kiss_fftr_alloc( nfft ,is_inverse,0,0 );
  //just swapped the input/output sizes 
  //actually don't even need the input size since data is handed in
  //size_t input_size = sizeof(kiss_fft_cpx) * freq_chunks;
  size_t output_size = sizeof(kiss_fft_scalar) * samples_per_slice;

  kiss_fft_scalar* fft_output =(kiss_fft_scalar*)KISS_FFT_MALLOC(output_size);

  
  //kiss_fftri(cfg, fft_input, fft_output);
  //output is a bunch of scalars, put into some stack frame
  for (int freq_i = 0; freq_i < samples_per_slice; freq_i++) {
    kiss_fft_scalar a = fft_output[freq_i];
    if (print) {
      std::cout << a << " ";
    }
    output[freq_i] = a;
  }
  if (print) {
    std::cout << "\n";
  }
  kiss_fftr_free(cfg);
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
  TWindowType wt = wtSINC;
  double some_window_param = 10;
  RectWinFIR(coefs, taps, band, rel_center, rel_margin);
  FIRFilterWindow(coefs, taps, wt, some_window_param);
}

stk::Fir create_fir_from_coefs(double* coefs, int len) {
  std::vector<StkFloat> kernel(0);
  StkFloat val = 0;
  for (int i = 0; i < len; i++) {
    val = coefs[i];
    kernel.insert(kernel.begin(), val);
  }
  Fir ret (kernel);
  return ret;
}

MatrixXf one_d_convolve(MatrixXf mat, MatrixXf kern) {
  //perform valid 1d convolution on matrix
  //or cross-corelation, planning on using symetric kernels so doesn't matter
  int res_dim;
  int kern_dim = -1;
  MatrixXf temp, result;
  int d_start, d_end, d_width;
  
  if (kern.rows() == 1) {
    result = MatrixXf(mat.rows(), mat.cols() - (kern.cols() - 1));
    result.Zero(result.rows(), result.cols());
    res_dim = result.cols();
    d_width = mat.rows();
    kern_dim = kern.cols();
    for (int i = 0; i < kern_dim; i++) {
      d_start = i;
      d_end = res_dim;
      temp = mat.block(0,d_start, d_width, d_end);
      temp = temp * kern(0,i);
      result += temp;
    }
  }
  else if (kern.cols() == 1) {
    result = MatrixXf(mat.rows() - (kern.rows() - 1), mat.cols());
    result.Zero(result.rows(), result.cols());
    res_dim = result.rows();
    d_width = mat.cols();
    kern_dim = kern.rows();
    for (int i = 0; i < kern_dim; i++) {
      d_start = i;
      d_end = res_dim;
      temp = mat.block(d_start, 0, d_end, d_width);
      temp = temp * kern(i, 0);
      result += temp;
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
