#include "util.h"
#include "ChunkFilter.h"
#include <iostream>
using stk::StkFrames;
using stk::Fir;
#define TAPS 120


ChunkFilter::ChunkFilter(float freq_center, int freq_margin, int size) {
  //just create an FIR filter based on chunk dimensions
  center = freq_center;
  margin = freq_margin;
  tot_size = size;
  int taps = 60;
  //was getting heap-corruption like behavior
  //seems to have been through coefs
  //give it a constant size so if 
  double coefs[taps];

  make_fir_bandpass_filter(coefs, taps, freq_center, freq_margin, tot_size);
  fir_filter = create_fir_from_coefs(coefs, taps);
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
  StkFrames frame_out (samples, 1);
  for (int i = 0; i < channels; i++) {
    fir_filter.tick(input, temp_frames, i, 0);
    output.getChannel(i, frame_out, 0);
    frame_out += temp_frames;
    output.setChannel(i, frame_out, 0);
  }
  return ret;
}



StkFrames fft_filter_frame(StkFrames &input) {
  //do an fft
  //zero out non-chunk parts
  //compute ifft of that
  //shove into stkFrames
  int channels = input.channels();
  int samples = input.frames();
  StkFrames filtered (0,samples, channels);
  return filtered;
}
