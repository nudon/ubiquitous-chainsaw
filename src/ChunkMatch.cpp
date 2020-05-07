#include <limits>
#include <iostream>
#include <list>
#include <map>
#include "ChunkMatch.h"

using stk::StkFrames;
using std::map;
using std::list;

static int get_closest_bin_index(map<int, list<Chunk>> bins, int start);
 
ChunkMatch::ChunkMatch(Chunk to_match) {
  orig = to_match;
  score = std::numeric_limits<double>::infinity();
}

void ChunkMatch::best_match_chunk(std::list<Chunk> &many_chunks, ChunkCompare &comp) {
  double lowest = score;
  double temp = 0;
  Chunk best_match = match;
  for (Chunk b : many_chunks) {
    temp = comp.compare(orig, b);
    if (temp < lowest) {
      best_match = b;
      lowest = temp;
      if (lowest < 0.00001) {
	break;
      }
    }
  }
  match = best_match;
  match.make_chunk_filter();
  score = lowest;
}

list<ChunkMatch> ChunkMatch::self_match(list<Chunk>&rec) {
  printf("self matching\n");
  list<ChunkMatch> matches;
  ChunkMatch a_match;
  for (Chunk a_chunk : rec) {
    a_match = ChunkMatch(a_chunk);
    a_match.match = a_chunk;
    a_match.get_match_chunk().make_chunk_filter();
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}
 
//iterates over all chunks for matches
list<ChunkMatch> ChunkMatch::best_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
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
list<ChunkMatch> ChunkMatch::quick_match(list<Chunk>&rec, list<Chunk> &rep, ChunkCompare& comp) {
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
    a_match.best_match_chunk(sub_list, comp);
    
    sub_list = freq_cent_map[fc_ind];
    a_match.best_match_chunk(sub_list, comp);
    
    matches.insert(matches.begin(), a_match);
  }
  return matches;
}

static int get_closest_bin_index(map<int, list<Chunk>> bins, int start) {
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

int ChunkMatch::get_active_start() {
  return orig.get_time_center() - match.get_time_margin();
}

int ChunkMatch::get_active_end() {
  return orig.get_time_center() + match.get_time_margin();
}

bool ChunkMatch::comp_active_start(ChunkMatch &a, ChunkMatch &b) {
  return a.get_active_start() < b.get_active_start();
}

bool ChunkMatch::comp_active_end(ChunkMatch &a, ChunkMatch &b) {
  return a.get_active_end() < b.get_active_end();
}

bool ChunkMatch::comp_orig_start(ChunkMatch &a, ChunkMatch &b) {
  return a.orig.get_time_start() < b.orig.get_time_start();
}
bool ChunkMatch::comp_orig_end(ChunkMatch &a, ChunkMatch &b) {
  return a.orig.get_time_end() < b.orig.get_time_end();
}


int ChunkMatch::match_hash( ChunkMatch &arg) {
  std::hash<double> hashy;
  double n;
  n = (arg.orig.get_time_start() - arg.match.get_time_end()) * (arg.orig.get_time_end() - arg.match.get_time_start());
  //d = arg.orig.get_rel_freq_center() / arg.match.get_rel_freq_center();
  return hashy(n);
}

bool  ChunkMatch:: operator == (const ChunkMatch other) const  {
  return (orig == other.orig && match == other.match);
}

  
