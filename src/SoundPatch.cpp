#include <FileWvIn.h>
#include <list>
#include "SoundPatch.h"
#include "util.h"

using stk::FileWvIn;
using stk::StkFrames;
using std::list;
//for constructors, iterate over chunnk(s) and save things



SoundPatch::SoundPatch(ChunkMatch& match, Song& src) {
  FileWvIn wav;
  wav.openFile(src.get_file_name());
  Chunk c = match.get_match_chunk();
  c.make_chunk_filter();
  int samples_per_slice = c.get_bin_size() * 2;
  StkFrames slice(samples_per_slice, src.get_channels());
  StkFrames zero(0, samples_per_slice, src.get_channels());
  int start = c.get_time_start();
  int end = c.get_time_end();
  for (int i = start; i <= end; i++) {
    fetch_frame(i, samples_per_slice, wav, slice);
    //same in/out might do weird things
    filter_frame(c, slice, slice);
    //might need some manual memory allocation here
    sav.push_back(slice);
    slice *= zero;
  }
  wav.closeFile();
  
  this->start_time = match.get_active_start();
  this->end_time = match.get_active_end();
}

SoundPatch::SoundPatch(ChunkGroupMatch& match, Song& src) {
  FileWvIn wav;
  wav.openFile(src.get_file_name());
  int samples_per_slice = match.get_orig().get_chunks().front().get_bin_size() * 2;
  StkFrames total(samples_per_slice, src.get_channels());
  StkFrames slice(samples_per_slice, src.get_channels());
   StkFrames zero(0, samples_per_slice, src.get_channels());
  int start = -1;
  int end = -1;
  start = match.get_match().get_super_chunk().get_time_start();
  end = match.get_match().get_super_chunk().get_time_end();
  list<Chunk> pool = match.get_match().get_chunks();
  double val = 1.0 / pool.size();
  StkFrames norm(val, samples_per_slice, src.get_channels());
  for (int i = start; i <= end; i++) {
    for (Chunk& c : pool) {
      c.make_chunk_filter();
      
      slice *= zero;
      fetch_frame(i, samples_per_slice, wav, slice);
      filter_frame(c, slice, total);
      //total += slice;
    }
    //save it

    total *= norm;
    sav.push_back(total);
    total *= zero;
  }
  wav.closeFile();

  this->start_time = match.get_active_start();
  this->end_time = match.get_active_end();
}

bool SoundPatch::comp_active_start(SoundPatch &a, SoundPatch &b) {
  return a.get_active_start() < b.get_active_start();
}

bool SoundPatch::comp_active_end(SoundPatch &a, SoundPatch &b) {
  return a.get_active_end() < b.get_active_end();
}

//time_i is time index in song output, 
void SoundPatch::tick_out(int time_i, stk::StkFrames& out) {
  int rel_i = time_i - get_active_start();
  if (time_i < get_active_start() || time_i > get_active_end()) {
    fprintf(stderr, "Error, time index fell outside active range of sound patch\n");
  }
  else {
    get_frame(rel_i, out);
  }
}

//i is a direct index to the stkframes
void SoundPatch::get_frame(int i, stk::StkFrames& out) {
  if (i < 0 || i >= (int)sav.size()) {
    fprintf(stderr, "Error, sound patch index out of range\n");
  }
  else {
    out += sav[i];
  }
}

int SoundPatch::patch_hash(const SoundPatch &arg) {
  std::hash<double> hashy;
  double n;
  n = arg.start_time * arg.end_time;
  return hashy(n);
}
