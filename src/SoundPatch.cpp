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
  for (int i = start; i <= end; i++) {
    for (Chunk c : pool) {
    //need to be able to iterate over ChunkGroup members
      slice *= zero;
      fetch_frame(i, samples_per_slice, wav, slice);
      filter_frame(c, slice, total);
    }
    //save it
    sav.push_back(total);
    slice *= zero;
  }
  wav.closeFile();

  this->start_time = match.get_active_start();
  this->end_time = match.get_active_end();
}

//time_i is time index in song output, 
void tick_out(int time_i, stk::StkFrames& out);

bool comp_active_start(ChunkMatch &a, ChunkMatch &b);
bool comp_active_end(ChunkMatch &a, ChunkMatch &b);

//i is a direct index to the stkframes
void get_frame(int i, stk::StkFrames& out);
