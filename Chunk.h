#ifndef FILE_CHUNK_SEEN
#define FILE_CHUNK_SEEN

class Chunk {
 public:
  Chunk() {
    min_freq = -1;
    max_freq = -1;
    min_time = -1;
    max_time = -1;
  }
  Chunk(int minf, int maxf, int mint, int maxt) {
    min_freq = minf;
    max_freq = maxf;
    min_time = mint;
    max_time = maxt;
  }
  
  ~Chunk() {}

  float get_freq_center() { return (max_freq + min_freq) / 2.0; }
  float get_time_center() { return (max_time + min_time) / 2.0; }

  int get_freq_margin() { return max_freq - min_freq; }
  int get_time_margin() { return max_time - min_time; }
  
 private:
  int min_freq;
  int max_freq;
  int min_time;
  int max_time;
};



#endif
