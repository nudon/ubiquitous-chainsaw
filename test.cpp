#include <stdio.h>
#include <FileLoop.h>
#include <FileWvOut.h>
#include <Fir.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>

#include "kiss_fftr.h"
#include "FIRFilterCode.h"
#include "eigen_to_image.h"

#include "ChunkStats.h"
#include "Chunk.h"

using namespace stk;

using Eigen::MatrixXf;
using Eigen::MatrixXi;
using std::list;

//main test
void process_song(std::string fn_in, std::string fn_out);

//thing with fourier transform
MatrixXf TFD_extract(FileWvIn input, int nfft, int total_slices, int samples_per_slice, int channels);
void TFDI_extract(kiss_fft_cpx* fft_input, int nfft, int samples_per_slice, StkFrames output);


//filtering methods
//using stk
void make_fir_bandpass_filter(double* coefs, int taps, int freq_center, int freq_margin,  int max_freq);
stk::Fir create_fir_from_coefs(double* coefs, int len);
stk::Fir create_1d_gaussian_filter(int length, double amplitude, double center, double stddev);
//using eigen
MatrixXf create_1d_gaussian_filter_col(int length, double amplitude, double stddev);
MatrixXf create_1d_gaussian_filter_row(int length, double amplitude, double stddev);
MatrixXf one_d_convolve(MatrixXf mat, MatrixXf kern);
MatrixXf dog(MatrixXf full, MatrixXf gauss);

//chunking functions
double snaz(MatrixXf filt, int snazr);
MatrixXi chunkify(MatrixXf tfd, int vert_range, int horz_range);
//list(Chunk) cull_chunks(MatrixXi chunk_ids, ChunkStats stats)
list<Chunk> cull_chunks(ChunkStats stats);

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
  
  write_eigen_to_file("some_spectogram.png", tfd);
  //test out convolution code
  MatrixXf guass_col = create_1d_gaussian_filter_col(16, 1, 1);
  MatrixXf guass_row = create_1d_gaussian_filter_row(16, 1, 1);
  MatrixXf sobel_x (1,3);
  MatrixXf sobel_y (3,1);
  //3x3 sobels but I only implemented 1d convolution :(
  //sobel_x << -1, 0, 1, -3, 0, 3, -1, 0, 1;
  //sobel_y << -1, -3, -1, 0, 0, 0, 1, 3, 1;
  sobel_x << -1, 0, 1;
  sobel_y << -1, 0, 1;

  MatrixXf tfd_y = one_d_convolve(tfd, sobel_y);
  MatrixXf tfd_x = one_d_convolve(tfd, sobel_x);

  write_eigen_to_file("y_sobel_spectogram.png", tfd_y);
  write_eigen_to_file("x_sobel_spectogram.png", tfd_x);
  
  
  MatrixXf sobel = one_d_convolve(tfd, sobel_y);
  sobel = one_d_convolve(sobel, sobel_x);
  write_eigen_to_file("abs_sobel_spectogram.png", sobel);
  
  //chunkify(sobel, 1,2);

  MatrixXf result = one_d_convolve(tfd, guass_col);
  result = one_d_convolve(result, guass_row);
  write_eigen_to_file("gauss_spectogram.png", result);
  MatrixXf dogg = dog(tfd, result);
  write_eigen_to_file("dog_spectogram.png", dogg);

  MatrixXf x_sobeled_dog = one_d_convolve(dogg, sobel_x);
  MatrixXf y_sobeled_dog = one_d_convolve(dogg, sobel_y);
  
  write_eigen_to_file("x_sobel_dog_spectogram.png", x_sobeled_dog * 100);
  write_eigen_to_file("y_sobel_dog_spectogram.png", y_sobeled_dog * 100);
  int vert = 1;
  int horz = 2;
  MatrixXi chunk_ids = chunkify(y_sobeled_dog, vert,horz);
  //MatrixXi stats = chunk_stats(chunk_ids);
  ChunkStats stats = ChunkStats(chunk_ids);
  list<Chunk> y_imps = cull_chunks(stats);
  chunk_ids = chunkify(x_sobeled_dog, vert,horz);
  stats = ChunkStats(chunk_ids);
  list<Chunk> x_imps = cull_chunks(stats);
  list<Chunk> all_imps;
  all_imps.insert(all_imps.begin(), x_imps.begin(), x_imps.end());
  all_imps.insert(all_imps.begin(), y_imps.begin(), y_imps.end());
  std::cout << "length of joined list is " << all_imps.size() << "\n";
						    
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
  //since nothing is placed in the upper half, divide by 2
  //though if you want to play around with inverse fft, just set to nfft I think?
  int freq_chunks = nfft / 2;
  kiss_fftr_cfg cfg = kiss_fftr_alloc( nfft ,is_inverse,0,0 );
  size_t input_size = sizeof(kiss_fft_scalar) * samples_per_slice;
  size_t output_size = sizeof(kiss_fft_cpx) * freq_chunks;
  kiss_fft_scalar* fft_input =(kiss_fft_scalar*)KISS_FFT_MALLOC(input_size);
  kiss_fft_cpx* fft_output =(kiss_fft_cpx*)KISS_FFT_MALLOC(output_size);
  bool print = false;
  StkFrames slice(samples_per_slice, channels);
  StkFrames single_chan(samples_per_slice, 1);

  MatrixXf tfd = MatrixXf(freq_chunks, total_slices); 
  tfd.setZero();
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
  std::cout << tfd(tfd.rows() / 4, tfd.cols() - 1) << "\n";
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


MatrixXf dog(MatrixXf full, MatrixXf gauss) {
  int f_rows = full.rows();
  int f_cols = full.cols();
  int g_rows = gauss.rows();
  int g_cols = gauss.cols();
  if (f_rows == g_rows && f_cols == g_cols) {
    return full - gauss;
  }
  else {
    MatrixXf temp = full.block((f_rows - g_rows) / 2, (f_cols - g_cols) / 2, g_rows, g_cols);
    temp = temp - gauss;
    return temp.cwiseAbs();
  }
}

//stands for sub nonzero-average zeroing
//take average of all non-zero elements, zero anything below average
//takes matrix and rounds to apply snaz, zeros values in matrix
//returns average of penultimate round
//having snazr as anything above 1 is typically excessive, consider at most using 2
double snaz(MatrixXf filt, int snazr) {
  int freq_range = filt.rows();
  int time_range = filt.cols();
  double nz_tot = 0;
  double nz_avg = 0;
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
	if (val < nz_avg) {
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
      std::cout << "non-zero average for round " << i << " is " << nz_avg << "\n";
    }
  }
  return nz_avg;
}

MatrixXi chunkify(MatrixXf tfd, int vert_range, int horz_range) {
  //chunkify the input matrix
  //defining chunkify as assign groups of sound in tfd that are signifigant together
  //do this by locally grouping cells within proximity in matrix
  //output is a new matrix, same dimensions, values are positive intergers
  //zero indicated cell wasn't assigned a chunk
  //positive integer represents some group id to which the cell got assigned to
  int chunk_id = 0;
  int freq_range = tfd.rows();
  int time_range = tfd.cols();
  int temp_i = -1;
  int temp = 0;
  float nz_avg = 0;
  MatrixXf filt (freq_range, time_range);
  MatrixXi chunk_ids (freq_range, time_range);
  filt << tfd;
  chunk_ids.setZero( );
  //threshold values < average to zero
  //do this to not care about very minor sounds
  //sub-nonzero-average-zeroing rounds
  nz_avg = snaz(filt, 1);
  
  //isn't a perfect grouping algorithim, since the local grouping is kind of naive
  //hoping further culling methods down the line will get rid of bad groups
  //but loop could be improved by having a more robust/expensive joining loop(s)
  for (int time_i = 0; time_i < time_range; time_i++) {
    for (int freq_i = 0; freq_i < freq_range; freq_i++) {
      if (filt(freq_i, time_i) > nz_avg) {
	//check if you want to join it vertically
	for (int freq_off = 1; freq_off <= vert_range; freq_off++) {
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
	for (int time_off = 1; time_off <= horz_range; time_off++) {
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

list<Chunk> cull_chunks(ChunkStats stats) {
  //wait I don't even need the chunkids anymore 
  
  //thinking of returning a linked list of chunk objects
  //because working with matrixes will probably be annoying past this point


  list<Chunk> chunk_list;
  MatrixXi sizes = stats.get_size();
  int chunks = sizes.rows();
  //ignoring chunk_id 0
  sizes = stats.get_size().block(1,0, chunks - 1,1);
  chunks = sizes.rows();

  double average_size = (float)sizes.sum() / sizes.sum();
  MatrixXi minf = stats.get_min_freq();
  MatrixXi maxf = stats.get_max_freq();
  MatrixXi mint = stats.get_min_time();
  MatrixXi maxt = stats.get_max_time();
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
