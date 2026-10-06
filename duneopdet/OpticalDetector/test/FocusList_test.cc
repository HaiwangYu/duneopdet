// Unit tests for duneopdet/OpticalDetector/FocusList.h: the focus ranges of
// WaveformDigitizerSim must stay disjoint, and must be the legacy ones
// whenever the legacy code kept them disjoint.

#define BOOST_TEST_MODULE FocusList_test
#include <boost/test/unit_test.hpp>

#include "duneopdet/OpticalDetector/FocusList.h"

#include <random>
#include <vector>

using opdet::FocusList;
using opdet::Ranges_t;

namespace {

  // Samples covered by the ranges, counting a sample once per range covering it
  std::vector<int> coverage(Ranges_t const& ranges, size_t nSamples)
  {
    std::vector<int> cov(nSamples, 0);
    for (auto const& r : ranges)
      for (size_t k = r.first; k <= r.second; ++k) ++cov[k];
    return cov;
  }

  bool disjoint(Ranges_t const& ranges, size_t nSamples)
  {
    for (int c : coverage(ranges, nSamples))
      if (c > 1) return false;
    return true;
  }

}

BOOST_AUTO_TEST_CASE(padding_and_clipping)
{
  FocusList fl(1000, 10);
  fl.AddRange(5, 20);     // front clipped at 0
  fl.AddRange(990, 995);  // back clipped at nSamples-1
  fl.AddRange(500, 510);
  Ranges_t const expected{{0, 30}, {980, 999}, {490, 520}};
  BOOST_TEST(fl.ranges == expected);
}

BOOST_AUTO_TEST_CASE(nested_extend_and_disjoint_as_legacy)
{
  // No range encloses or bridges another: the legacy result, order included
  FocusList fixed(1000, 0), legacy(1000, 0, false);
  for (auto* fl : {&fixed, &legacy}) {
    fl->AddRange(100, 200);
    fl->AddRange(500, 600);
    fl->AddRange(150, 180);  // nested
    fl->AddRange(190, 250);  // extends the end of the first
    fl->AddRange(450, 520);  // extends the front of the second
    fl->AddRange(300, 310);  // disjoint
  }
  Ranges_t const expected{{100, 250}, {450, 600}, {300, 310}};
  BOOST_TEST(fixed.ranges == expected);
  BOOST_TEST(legacy.ranges == expected);
}

BOOST_AUTO_TEST_CASE(enclosing_range_is_merged)
{
  FocusList fixed(1000, 0), legacy(1000, 0, false);
  for (auto* fl : {&fixed, &legacy}) {
    fl->AddRange(100, 200);
    fl->AddRange(50, 300);   // encloses [100, 200]
  }
  BOOST_TEST(fixed.ranges == (Ranges_t{{50, 300}}));
  // The legacy bug: appended as a second, overlapping range
  BOOST_TEST(legacy.ranges == (Ranges_t{{100, 200}, {50, 300}}));
}

BOOST_AUTO_TEST_CASE(bridging_range_is_merged)
{
  FocusList fixed(1000, 0), legacy(1000, 0, false);
  for (auto* fl : {&fixed, &legacy}) {
    fl->AddRange(100, 200);
    fl->AddRange(300, 400);
    fl->AddRange(500, 600);
    fl->AddRange(150, 550);  // reaches the second and third ranges
  }
  // Merged into the first overlapping range, the others removed
  BOOST_TEST(fixed.ranges == (Ranges_t{{100, 600}}));
  // The legacy bug: the first range is extended over the other two
  BOOST_TEST(legacy.ranges == (Ranges_t{{100, 550}, {300, 400}, {500, 600}}));
}

BOOST_AUTO_TEST_CASE(touching_ranges_stay_separate)
{
  // Ranges sharing no sample are not merged, as in the legacy code
  FocusList fixed(1000, 0);
  fixed.AddRange(100, 200);
  fixed.AddRange(201, 300);
  BOOST_TEST(fixed.ranges == (Ranges_t{{100, 200}, {201, 300}}));
  fixed.AddRange(200, 201);
  BOOST_TEST(fixed.ranges == (Ranges_t{{100, 300}}));
}

BOOST_AUTO_TEST_CASE(everything_covers_the_waveform)
{
  FocusList fl(1000, 10);
  fl.AddRange(100, 200);
  fl.everything();
  BOOST_TEST(fl.ranges == (Ranges_t{{0, 999}}));
}

BOOST_AUTO_TEST_CASE(random_pileup)
{
  // Pulses at random times, as in WaveformDigitizerSim: 325-tick pulses plus
  // 100 ticks of padding. The fixed ranges must be disjoint and cover exactly
  // the padded pulses; whenever the legacy ranges came out disjoint they must
  // be identical to the fixed ones.
  size_t const nSamples = 20000, padding = 100, pulse = 325;
  std::mt19937 gen(12345);
  int nLegacyDisjoint = 0, nLegacyOverlap = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    size_t const nPulses = 1 + trial % 60;
    std::uniform_int_distribution<size_t> start(0, nSamples - 1);
    FocusList fixed(nSamples, padding), legacy(nSamples, padding, false);
    std::vector<bool> truth(nSamples, false);
    for (size_t p = 0; p < nPulses; ++p) {
      size_t const from = start(gen);
      size_t const to = std::min(from + pulse - 1, nSamples - 1);
      fixed.AddRange(from, to);
      legacy.AddRange(from, to);
      for (size_t k = (from > padding ? from - padding : 0); k <= std::min(to + padding, nSamples - 1); ++k)
        truth[k] = true;
    }
    auto const cov = coverage(fixed.ranges, nSamples);
    bool sameCover = true, isDisjoint = true;
    for (size_t k = 0; k < nSamples; ++k) {
      if ((cov[k] > 0) != truth[k]) sameCover = false;
      if (cov[k] > 1) isDisjoint = false;
    }
    BOOST_TEST(isDisjoint);
    BOOST_TEST(sameCover);
    if (disjoint(legacy.ranges, nSamples)) {
      ++nLegacyDisjoint;
      BOOST_TEST(fixed.ranges == legacy.ranges);
    }
    else
      ++nLegacyOverlap;
  }
  // Both regimes must have been exercised
  BOOST_TEST(nLegacyDisjoint > 100);
  BOOST_TEST(nLegacyOverlap > 100);
}
