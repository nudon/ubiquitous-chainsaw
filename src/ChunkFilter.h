#ifndef FILE_FILTER_SEEN
#define FILE_FILTER_SEEN
#include <Stk.h>
#include <Fir.h>
#include "kiss_fftr.h"
#include "FIRFilterCode.h"

class ChunkFilter {
 public:
  ChunkFilter() { center = -1; margin = -1; tot_size = -1;}
  ChunkFilter(float freq_center, int freq_margin, int tot_size);
  ~ChunkFilter() {}

  int fir_filter_frame(stk::StkFrames &input, stk::StkFrames &output);
  int fft_filter_frame(stk::StkFrames &input, stk::StkFrames &output);

  float get_center() { return center; }
  int get_margin() {return margin; }
  int get_size() { return tot_size; }

 private:
  float center;
  int margin;
  int tot_size;

  int fir_filter_multi(stk::StkFrames &input, stk::StkFrames &output);
  int fft_filter_multi(stk::StkFrames &input, stk::StkFrames &output, kiss_fftr_cfg &fft, kiss_fftr_cfg &ifft);
  void filter_kiss_cpx(kiss_fft_cpx* data, int size);
  
  stk::Fir fir_filter;
  
  
};

#endif
