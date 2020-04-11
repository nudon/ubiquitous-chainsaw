#include <iostream>
#include <Eigen/Eigen>
#include <list>
#include <map>
#include <unordered_set>
#include <math.h>
#include <LentPitShift.h>

#include "SongEmbedder.h"
#include "Chunk.h"
#include "ChunkMatch.h"
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


list<ChunkMatch> best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp);
int get_closest_bin_index(map<int, list<Chunk>> bins, int start);


SongEmbedder::SongEmbedder(Song src, Song dst, int sps, int edge_r, int chunk_r, int match_r) {
  reciptor = src;
  replacer = dst;
  samples_per_slice = sps;
  edge_snazr = edge_r;
  chunk_snazr = chunk_r;
  match_snazr = match_r;  
}

double chunkmatch_get_score_wrapper(ChunkMatch a) {
  return a.get_score();
}



void SongEmbedder::funk() {
  if (reciptor.get_file_rate() == replacer.get_file_rate()) {
    list<Chunk> reciptor_chunks = reciptor.make_chunks(samples_per_slice, edge_snazr, chunk_snazr);
    list<Chunk> replacer_chunks = replacer.make_chunks(samples_per_slice, edge_snazr, chunk_snazr);
    list<ChunkMatch>matches = get_matches(reciptor_chunks, replacer_chunks);
    
    std::string fn = "muxed.wav";
    std::cout << "outputting song\n";
    output_remix(matches, fn);
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

list<ChunkMatch> SongEmbedder::get_matches(list<Chunk> &reciptor_chunks, list<Chunk> &replacer_chunks) {
  //weights for frequency center, frequency margin, and time marin
  ChunkCompare comp (50,1,10);

  std::cout << "matching sound chunks!\n";
  list<ChunkMatch> matches = quick_match(reciptor_chunks, replacer_chunks, comp);
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

//iterates over all chunks for matches
list<ChunkMatch> best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
  list<ChunkMatch> matches;
  ChunkMatch a_match;
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    a_match.best_match_chunk(rep, comp);
    a_match.get_match_chunk().make_chunk_filter();
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

//iterates over chunks which has at least one identical stat for matches
list<ChunkMatch> quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
  map<int, list<Chunk>> freq_cent_map;
  map<int, list<Chunk>> time_marg_map;
  list<ChunkMatch> matches;
  int fc_ind, tm_ind, temp;
  
  //populate maps
  list<Chunk> sub_list;
  for (Chunk a_chunk : rep) {
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();
    
    time_marg_map[tm_ind].push_front(a_chunk);
    freq_cent_map[fc_ind].push_front(a_chunk);
  }

  ChunkMatch a_match;
  map<int, int> closest_freq_cent;
  map<int, int> closest_time_marg;
  
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    
    fc_ind = int(round(2 * a_chunk.get_freq_center()));
    tm_ind = a_chunk.get_time_length();

    if (closest_freq_cent.count(fc_ind) == 0) {
      temp = get_closest_bin_index(freq_cent_map, fc_ind);
      closest_freq_cent.insert(std::pair<int,int>(fc_ind, temp));
    }
    fc_ind = closest_freq_cent[fc_ind];

    if (closest_time_marg.count(tm_ind) == 0) {
      temp = get_closest_bin_index(time_marg_map, tm_ind);
      closest_time_marg.insert(std::pair<int,int>(tm_ind, temp));
    }
    tm_ind = closest_time_marg[fc_ind];
    
    sub_list = time_marg_map[tm_ind];
    //a_match.best_match_chunk(sub_list, comp);
    
    sub_list = freq_cent_map[fc_ind];
    a_match.best_match_chunk(sub_list, comp);
    
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

int get_closest_bin_index(map<int, list<Chunk>> bins, int start) {
  //various bins may have zero items in them
  //find the closest non-empty one
  bool done = false;
  int offset = 0;
  int ind = -1;
  while(!done) {
    if (bins.count(start + offset) == 1) {
      done = true;
      ind = start + offset;
    }
    else if (start - offset >= 0) {
      if (bins.count(start - offset) == 1) {
	done = true;
	ind = start - offset;
      }
    }
    offset++;
  }
  return ind;
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

    output.tick(frame_out);
    
    while(active_end.size() > 0 && active_end.front().get_active_end() == i) {
      active_chunks.erase(active_end.front());
      active_end.pop_front();
    }
  }
  output.closeFile();
}

void SongEmbedder::normal_filt(ChunkMatch& chunk_match, FileWvIn &sample_src, StkFrames &frame_in, StkFrames &out, int i) {
  Chunk match = chunk_match.get_match_chunk();
  ChunkFilter filt = match.get_filter();
  int offset = (i - chunk_match.get_active_start()) + match.get_time_start();
  offset *= samples_per_slice;
  sample_src.reset();
  sample_src.addTime(offset);
  
  sample_src.tick(frame_in);
  
  //filt.fir_filter_frame(frame_in, out);
  filt.fir_filter_frame(frame_in);
  out += frame_in;
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
     
  filt.fir_filter_frame(temp_in);
  reshape_chunk(temp_in, out, orig, match, lent);
}
