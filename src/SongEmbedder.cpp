#include <iostream>
#include <Eigen/Eigen>
#include <list>
#include <map>
#include <unordered_set>
#include <math.h>
#include <LentPitShift.h>

#include "SongEmbedder.h"
#include "Chunk.h"
#include "util.h"


using std::list;
using std::map;
using std::unordered_set;
using std::pair;
using std::string;
using std::max;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
using stk::FileWvIn;
using stk::FileWvOut;
using stk::StkFrames;
using stk::LentPitShift;


int get_closest_bin_index(map<int, list<Chunk>> bins, int start);


SongEmbedder::SongEmbedder(Song src, Song dst, int sps, int edge_r, int chunk_r, int match_r) {
  reciptor = src;
  replacer = dst;
  samples_per_slice = sps;
  edge_snazr = edge_r;
  chunk_snazr = chunk_r;
  match_snazr = match_r;  
}

double chunkmatch_get_score_wrapper(ChunkMatch& a) {
  return a.get_score();
}

double chunkgroupmatch_get_score_wrapper(ChunkGroupMatch& a) {
  return a.get_score();
}

list<SoundPatch> SongEmbedder::individual_sound_patches(list<Chunk>& res, list<Chunk>& src) {
  list<ChunkMatch> matches = get_matches(res, src);
  list<SoundPatch> ret;
  for (ChunkMatch& m : matches) {
    ret.push_back(SoundPatch(m, replacer));
  }
  return ret;
}


list<SoundPatch> SongEmbedder::group_sound_patches(list<Chunk>& res, list<Chunk>& src) {
  
  list<ChunkGroup> res_pool = group_chunks(res);
  printf("For result chunks, put %lu chunks into %lu groups\n", res.size(), res_pool.size());
  list<ChunkGroup> src_pool = group_chunks(src);
  printf("For source chunks, put %lu chunks into %lu groups\n", src.size(), src_pool.size());
  list<ChunkGroupMatch> matches = get_matches(res_pool, src_pool);
  list<SoundPatch> ret;
  for (ChunkGroupMatch& m : matches) {
    ret.push_back(SoundPatch(m, replacer));
  }
  return ret;
}

list<ChunkMatch> SongEmbedder::get_matches(list<Chunk> &reciptor_chunks, list<Chunk> &replacer_chunks) {
  //weights for frequency center, frequency margin, and time marin
  ChunkCompare comp (50,1,10);

  std::cout << "matching sound chunks!\n";
  list<ChunkMatch> matches = ChunkMatch::quick_match(reciptor_chunks, replacer_chunks, comp);
  std::cout << "Done with matching!\n";
  if (match_snazr > 0) {
    std::cout << "culling bad matches\n";
    std::cout << "orig len is " << matches.size() << "\n";
    double avg = reverse_gsnaz(matches, chunkmatch_get_score_wrapper, match_snazr);
    std::cout << "only accepted scores <= " << avg << "\n";
    std::cout << "culled len is " << matches.size() << "\n";
  }

  return matches;
}

list<ChunkGroupMatch> SongEmbedder::get_matches(list<ChunkGroup> &res, list<ChunkGroup>& src) {
  //weights for frequency center, frequency margin, and time marin
  ChunkGroupCompare comp (3,9,1);

  std::cout << "matching chunk groups!\n";
  list<ChunkGroupMatch> matches = ChunkGroupMatch::best_match(res, src, comp);
  std::cout << "Done with matching!\n";
  if (match_snazr > 0) {
    std::cout << "culling bad matches\n";
    std::cout << "orig len is " << matches.size() << "\n";
    double avg = reverse_gsnaz(matches, chunkgroupmatch_get_score_wrapper, match_snazr);
    std::cout << "only accepted scores <= " << avg << "\n";
    std::cout << "culled len is " << matches.size() << "\n";
  }

  return matches;
}



void SongEmbedder::funk() {
  if (reciptor.get_file_rate() == replacer.get_file_rate()) {
    list<Chunk> reciptor_chunks = reciptor.make_chunks(samples_per_slice, edge_snazr, chunk_snazr);
    list<Chunk> replacer_chunks = replacer.make_chunks(samples_per_slice, edge_snazr, chunk_snazr);
    
    std::string fn = "muxed.wav";
    bool useSoundPatch = true;
    bool useGroupMatch = true;
    if (useSoundPatch) {
      list<SoundPatch> patches;
      if (useGroupMatch) {
	patches = group_sound_patches(reciptor_chunks, replacer_chunks);
      }
      else {
	patches = individual_sound_patches(reciptor_chunks, replacer_chunks);
      } 
      output_remix(patches, fn);
    }
    else {
      list<ChunkMatch>matches = get_matches(reciptor_chunks, replacer_chunks);
      output_remix(matches, fn);
    }
    //output_remix_alt(samples_per_slice, matches, "muxed_alt.wav");
    std::cout << "completed output\n";
  }
  else {
    //ERROR
    //would have to do weird stuff to work with this
    //because the different sampling rates would manifest in fft
    //would have the frequency bins of both spectrograms not line up
    std::cerr << "unable to mix files, sample rates are different, please use ffmpeg or alternative software to convert songs to same sample rate\n"; 
  }
}

struct cm_hash {
  size_t operator ()(const ChunkMatch & cm) const {
    return ChunkMatch::match_hash((ChunkMatch&)cm);
  }
};

void SongEmbedder::output_remix(list<ChunkMatch> &matches, std::string fn) {
  list<ChunkMatch> active_start;
  list<ChunkMatch> active_end;
  for (ChunkMatch a_match : matches) {
    active_start.push_front(a_match);
    active_end.push_front(a_match);
  }
  
  active_start.sort(ChunkMatch::comp_active_start);
  active_end.sort(ChunkMatch::comp_active_end);
  int new_slices = reciptor.get_samples_per_channel() / samples_per_slice;
  int channels = reciptor.get_channels();
  FileWvIn sample_src;
  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  StkFrames temp (0, samples_per_slice, channels);

  LentPitShift lent(1, samples_per_slice);
  
  StkFrames rev_feather (0, samples_per_slice, channels);
  StkFrames reg_feather (0, samples_per_slice, channels);
  int feather_margin = samples_per_slice;
  float feather_val = 0;
  for (int i = 0; i < samples_per_slice; i++) {
    feather_val = std::min(1.0f, (float)i / feather_margin);
    for (int chan = 0; chan < channels; chan++) {
      reg_feather(i, chan) = feather_val;
      rev_feather(samples_per_slice - 1 - i, chan) = feather_val;
    }
  }
    
  sample_src.openFile(replacer.get_file_name());
  output.openFile(fn, channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  
  unordered_set<ChunkMatch,cm_hash> active_chunks; 
  for (int i = 0; i < new_slices; i++) {
    while(active_start.size() > 0 && active_start.front().get_active_start() <= i) {
      if (active_start.front().get_active_start() == i) {
	active_chunks.insert(active_start.front());
	active_start.pop_front();
      }
      else if (active_start.front().get_active_start() < 0) {
	//odd case, happens when song matches a chunk which is long enough and close enough to start
	//that it has a negative start time.
	//just tossing it away for now
	active_start.pop_front();
      }
    }   

    frame_out *= zero;
    for (ChunkMatch chunk_match : active_chunks) {
      temp *= zero;
      normal_filt(chunk_match, sample_src, frame_in, temp, i);
      //alt_filt(chunk_match, sample_src, temp, i, lent);

      if (chunk_match.get_active_start() == i) {
	temp *= reg_feather;
      }
      if (chunk_match.get_active_end() == i) {
	temp *= rev_feather;
      }
      
      frame_out += temp;
    }
    //skips sections of silence from start and end of song
    if (active_start.size() != matches.size() && active_end.size() != 0) {
      output.tick(frame_out);
    }
    
    while(active_end.size() > 0 && active_end.front().get_active_end() == i) {
      active_chunks.erase(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}

void SongEmbedder::normal_filt(ChunkMatch& chunk_match, FileWvIn &sample_src, StkFrames &frame_in, StkFrames &out, int i) {
  Chunk match = chunk_match.get_match_chunk();
  int frame_n = (i - chunk_match.get_active_start()) + match.get_time_start();
  fetch_frame(frame_n, samples_per_slice, sample_src, frame_in);
  filter_frame(match, frame_in, out);
}


void SongEmbedder::alt_filt(ChunkMatch& chunk_match, FileWvIn &sample_src, StkFrames &out, int i, LentPitShift &lent)  {
  int channels = sample_src.channelsOut();;
  Chunk match = chunk_match.get_match_chunk();
  Chunk orig = chunk_match.get_orig_chunk();
  ChunkFilter filt = match.get_filter();
  double offset_ratio =  (double)(i - orig.get_time_start()) / (orig.get_time_length());
  int offset = (offset_ratio * match.get_time_length() + match.get_time_start()) * samples_per_slice;
  sample_src.reset();
  sample_src.addTime(offset);

  double compression_ratio = (double)orig.get_time_length() / match.get_time_length();
  int size = samples_per_slice / compression_ratio;
  StkFrames temp_in (size, channels);
  sample_src.tick(temp_in);
     
  filt.fir_filter_frame(temp_in, temp_in);
  reshape_chunk(temp_in, out, orig, match, lent);
}

static bool sound_p_act_start_wrap(SoundPatch* a, SoundPatch* b) {
  return SoundPatch::comp_active_start(*a,*b);
}

static bool sound_p_act_end_wrap(SoundPatch* a, SoundPatch* b) {
  return SoundPatch::comp_active_end(*a,*b);
}

struct sp_hash {
  size_t operator ()(const SoundPatch* sp) const {
    return SoundPatch::patch_hash(*sp);
  }
};

void SongEmbedder::output_remix(list<SoundPatch>& patches, std::string fn) {
  std::cout << "outputting song\n";
  list<SoundPatch*> active_start;
  list<SoundPatch*> active_end;
  unordered_set<SoundPatch*,sp_hash> active_patches; 
  for (SoundPatch& p : patches) {
    active_start.push_front(&p);
    active_end.push_front(&p);
  }
  
  active_start.sort(sound_p_act_start_wrap);
  active_end.sort(sound_p_act_end_wrap);
  
  int new_slices = reciptor.get_samples_per_channel() / samples_per_slice;
  int channels = reciptor.get_channels();

  FileWvOut output;
  StkFrames frame_in (samples_per_slice, channels);
  StkFrames frame_out (samples_per_slice, channels);
  StkFrames zero (0, samples_per_slice, channels);
  StkFrames temp (0, samples_per_slice, channels);

  LentPitShift lent(1, samples_per_slice);
  
  StkFrames rev_feather (0, samples_per_slice, channels);
  StkFrames reg_feather (0, samples_per_slice, channels);
  int feather_margin = samples_per_slice;
  float feather_val = 0;
  for (int i = 0; i < samples_per_slice; i++) {
    feather_val = std::min(1.0f, (float)i / feather_margin);
    for (int chan = 0; chan < channels; chan++) {
      reg_feather(i, chan) = feather_val;
      rev_feather(samples_per_slice - 1 - i, chan) = feather_val;
    }
  }
    
  output.openFile(fn, channels, stk::FileWrite::FILE_WAV, stk::Stk::STK_SINT16);
  

  for (int i = 0; i < new_slices; i++) {
    while(active_start.size() > 0 && active_start.front()->get_active_start() <= i) {
      if (active_start.front()->get_active_start() == i) {
	active_patches.insert(active_start.front());
	active_start.pop_front();
      }
      else if (active_start.front()->get_active_start() < 0) {
	//odd case, happens when song matches a chunk which is long enough and close enough to start
	//that it has a negative start time.
	//just tossing it away for now
	active_start.pop_front();
      }
    }   

    frame_out *= zero;
    for (SoundPatch* p : active_patches) {
      temp *= zero;

      p->tick_out(i, temp);

      if (p->get_active_start() == i) {
	temp *= reg_feather;
      }
      if (p->get_active_end() == i) {
	temp *= rev_feather;
      }
      frame_out += temp;
    }
    //skips sections of silence from start and end of song
    if (active_start.size() != patches.size() && active_end.size() != 0) {
      output.tick(frame_out);
    }
    
    while(active_end.size() > 0 && active_end.front()->get_active_end() == i) {
      active_patches.erase(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}
