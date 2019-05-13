#include "util.h"
#include "ChunkFilter.h"
#include <iostream>
using stk::StkFrames;
using stk::Fir;

//underlying library for setting coeffs allows taps up to 256
#define TAPS 60


ChunkFilter::ChunkFilter(float freq_center, int freq_margin, int size) {
  //just create an FIR filter based on chunk dimensions
  center = freq_center;
  margin = freq_margin;
  tot_size = size;
  double coefs[TAPS];

  make_fir_bandpass_filter(coefs, TAPS, freq_center, freq_margin, tot_size);
  fir_filter = create_fir_from_coefs(coefs, TAPS);
  //std::cout << fir_filter << "\n";
}


StkFrames ChunkFilter::fir_filter_frame(StkFrames &input) {
  //apply fir_filter to input

  int channels = input.channels();
  int samples = input.frames();
  StkFrames filtered (0,samples, channels);
  //copy input to filtered
  filtered += input;
  for (int i = 0; i < channels; i++) {
    fir_filter.tick(filtered, i);
  }
  return filtered;
}

int ChunkFilter::fir_filter_frame(StkFrames &input, StkFrames &output) {
  int ret = 0;
  if (input.channels() != output.channels()) {
    std::cerr << "filter input and output have a different amount of channels\n";
    ret = 1;
  }
  if (input.frames() != output.frames()) {
    std::cerr << "filter input and output have different sizes\n";
    ret = 2;
  }
  //apply fir_filter to input

  if (ret == 0) {
    ret = fir_filter_multi(input,output);
  }
  return ret;
}

int ChunkFilter::fir_filter_multi(StkFrames &input, StkFrames &output) {
  //apply fir_filter to input
  int ret = 0;
  int channels = input.channels();
  int samples = input.frames();
  StkFrames temp_frames (samples,1);
  //StkFrames frame_out (samples, 1);
  StkFrames filt(samples, channels);
  for (int i = 0; i < channels; i++) {
    fir_filter.tick(input, temp_frames, i, 0);
    filt.setChannel(i, temp_frames, 0);
  }
  output += filt;
  return ret;
}


int ChunkFilter::fft_filter_frame(StkFrames &input, StkFrames &output) {
  int ret = 0;
  if (input.channels() != output.channels()) {
    std::cerr << "filter input and output have a different amount of channels\n";
    ret = 1;
  }
  if (input.frames() != output.frames()) {
    std::cerr << "filter input and output have different sizes\n";
    ret = 2;
  }
  //apply fir_filter to input

  if (ret == 0) {
    kiss_fftr_cfg fft_d = kiss_fftr_alloc(2, false, 0,0);
    ret = fft_filter_multi(input,output, fft_d, fft_d);
  }
  return ret;
}


int ChunkFilter::fft_filter_multi(stk::StkFrames &input, stk::StkFrames &output, kiss_fftr_cfg &fft, kiss_fftr_cfg &ifft) {
  //so, take a frame of input, get a tfd out of it,
  //try manipulating entries of tfd then calling inverse, store in output

  //code for initing cfgs, put outside of function body
  int samples_per_slice = input.frames();
  int channels = input.channels();
  kiss_fftr_cfg fft_normal = kiss_fftr_alloc(samples_per_slice, false, 0,0);
  kiss_fftr_cfg fft_inverse = kiss_fftr_alloc(samples_per_slice, true, 0,0);
  
  size_t input_size = sizeof(kiss_fft_scalar) * (samples_per_slice + 1);
  size_t output_size = sizeof(kiss_fft_cpx) * (samples_per_slice + 1);
  kiss_fft_scalar* fft_input =(kiss_fft_scalar*)KISS_FFT_MALLOC(input_size);
  kiss_fft_cpx* fft_output =(kiss_fft_cpx*)KISS_FFT_MALLOC(output_size);
  StkFrames single_chan(samples_per_slice, 1);
  StkFrames filt(samples_per_slice, channels);
  //generate tfd for input
  for (int i = 0; i < channels; i++) {
    input.getChannel(i, single_chan, 0);
    stk_to_kiss(single_chan, fft_input);
    kiss_fftr(fft_normal, fft_input, fft_output);
    //do things to fft_output
    filter_kiss_cpx(fft_output, samples_per_slice);
    kiss_fftri(fft_inverse, fft_output, fft_input);
    kiss_to_stk(fft_input, single_chan);
    filt.setChannel(i, single_chan, 0);
  }
  output += filt;
  return 0;
}

void ChunkFilter::filter_kiss_cpx(kiss_fft_cpx* data, int size) {
  kiss_fft_cpx zero = (kiss_fft_cpx){.r = 0, .i = 0};
  float conv_margin = margin;
  float conv_center = center;
  int nyquist = size / 2;
  if (size != tot_size) {
    conv_margin = size * (margin / tot_size);
    conv_center = size * (center / tot_size);
  }
  for (int i = 0; i < nyquist; i++) {
    if (abs(i - conv_center) > conv_margin) {
      data[i] = zero;
    }
  }
}
