#include "util.h"
#include "ChunkFilter.h"
using stk::StkFrames;
#define TAPS 60
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
