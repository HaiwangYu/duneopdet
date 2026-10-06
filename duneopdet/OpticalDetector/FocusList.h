//=========================================================
// FocusList.h
// The sample ranges ("focus ranges") of a simulated photon
// detector waveform where a digitizer adds line noise and
// evaluates its self-trigger: each photoelectron pulse plus
// padding on both sides.
//
// The ranges are kept disjoint. Overlapping ranges make
// every loop over them visit the same samples more than
// once: WaveformDigitizerSim then stores those samples in
// two readout windows (duplicate OpHits downstream) and adds
// the line noise to them twice. With pile-up (radiological
// backgrounds) the legacy AddRange left ranges overlapping
// whenever a new range enclosed an existing one or an
// extension reached a neighbouring range.
//=========================================================

#ifndef DUNEOPDET_OPTICALDETECTOR_FOCUSLIST_H
#define DUNEOPDET_OPTICALDETECTOR_FOCUSLIST_H

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace opdet {

  typedef std::vector< std::pair< size_t, size_t > > Ranges_t;

  class FocusList
  {
  public:
    // mergeOverlaps = false reproduces the legacy behaviour, which can leave
    // overlapping ranges; kept only to reproduce existing samples.
    FocusList(size_t nSamples, size_t padding, bool mergeOverlaps = true)
      : fNSamples(nSamples), fPadding(padding), fMergeOverlaps(mergeOverlaps)
    {}

    // Add the range [from, to] (inclusive, in ticks), padded on both sides
    // and clipped to the waveform.
    void AddRange(size_t from, size_t to)
    {
      from -= std::min(from, fPadding);
      to   =  std::min(to+fPadding, fNSamples-1);

      if (!fMergeOverlaps) {
        AddRangeLegacy(from, to);
        return;
      }

      // Merge the new range with every range it overlaps. The union goes into
      // the first overlapping range in list order (the one the legacy code
      // extends) and the others are removed. The existing ranges are disjoint,
      // so one pass finds all of them. When nothing is enclosed or bridged the
      // result, including the order of the ranges, is the legacy one.
      auto merged = ranges.end();
      for (auto it = ranges.begin(); it != ranges.end(); ) {
        if (from <= it->second && to >= it->first) {
          from = std::min(from, it->first);
          to   = std::max(to, it->second);
          if (merged == ranges.end()) {
            merged = it;
            ++it;
          }
          else {
            it = ranges.erase(it);
          }
        }
        else {
          ++it;
        }
      }

      if (merged == ranges.end()) ranges.emplace_back(from, to);
      else                        *merged = std::make_pair(from, to);
    }

    void everything()
    {
      ranges.clear();
      ranges.emplace_back(0, fNSamples-1);
    }

    Ranges_t ranges;

  protected:
    size_t fNSamples;
    size_t fPadding;
    bool   fMergeOverlaps;

  private:
    void AddRangeLegacy(size_t from, size_t to)
    {
      for(size_t i = 0; i < ranges.size(); ++i){
        std::pair<size_t, size_t>& r = ranges[i];
        // Completely nested, discard
        if(from >= r.first && to <= r.second) return;
        // Extend end
        if(from >= r.first && from <= r.second){
          r.second = to;
          return;
        }
        // Extend front
        if(to >= r.first && to <= r.second){
          r.first = from;
          return;
        }
      }
      // Discontiguous, add
      ranges.emplace_back(from, to);
    }
  };

}

#endif
