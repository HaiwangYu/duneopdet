//=========================================================
// SIPMOpSensorSim_module.cc
// This module produces detected photons (creating OpDetDivRec)
// from photon detectors taking OpDetBacktrackerRecords as input.
// Applies quantum efficiency and cross-talk
//
// Gleb Sinev, Duke, 2015
// Anne Christensen, CSU, 2019
// Alex Himmel, FNAL, 2021
//=========================================================

#ifndef SIPMOpSensorSim_h
#define SIPMOpSensorSim_h

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h" //vitor
#include "art_root_io/TFileDirectory.h"



// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes

#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "duneopdet/OpticalDetector/AlgoSSPLeadingEdge.h"
#include "dunecore/DuneObj/OpDetDivRec.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"


// CLHEP includes

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandBinomial.h"

// C++ includes

#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <memory>

// ROOT includes

#include "TTree.h"

namespace opdet {

  class SIPMOpSensorSim : public art::EDProducer{

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>  InputTag            { Name("InputTag"),            Comment("Input tag for OpDetBacktrackerRecords") };
      fhicl::Atom<double>         QuantumEfficiency   { Name("QuantumEfficiency"),   Comment("Probabilityof recording a photon") };
      fhicl::Atom<double>         DarkNoiseRate       { Name("DarkNoiseRate"),       Comment("Rate in Hz") };
      fhicl::Atom<double>         CrossTalk           { Name("CrossTalk"),           Comment("Cross talk (1->2 PE) probability") };
      fhicl::Atom<double>         Correction          { Name("Correction"),          Comment("Adjust the amount of total light. Kept seprate from QE for clarity."), 1};
      fhicl::OptionalAtom<double> LateLightCorrection { Name("LateLightCorrection"), Comment("Adjust the amount of late light")};
      fhicl::OptionalAtom<double> LateLightBoundary   { Name("LateLightBoundary"),   Comment("Boundary that defines late light (ns)") };
      fhicl::Atom<std::string>    LateLightReference  { Name("LateLightReference"),
                                                        Comment("Time that late light is counted from, on each OpDet: \"Event\" = first photon of the "
                                                                "whole event (legacy; with pile-up every later interaction's light is late), "
                                                                "\"Track\" = first photon of the same track"),
                                                        "Event" };
      fhicl::Atom<bool>           DarkNoiseTimeInNs   { Name("DarkNoiseTimeInNs"),
                                                        Comment("Store dark-count times in ns, like the photon times. false reproduces the legacy "
                                                                "behaviour, which stores the us value as ns and puts every dark count within "
                                                                "a few us of t = 0. Default true."),
                                                        true };
      fhicl::Atom<bool>           BinomialQE          { Name("BinomialQE"),
                                                        Comment("Detect each photon of an SDP with probability QE (binomial). false reproduces the "
                                                                "legacy Poisson(QE*n), which on top of PDFastSim's own Poisson draw inflates the PE "
                                                                "variance and allows more PE than photons. Default true."),
                                                        true };

    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit SIPMOpSensorSim(Parameters const & config);
    void produce(art::Event&) override;

  private:

    // The parameters read from the FHiCL file
    art::InputTag fInputTag;            // Input tag for OpDet collection
    art::ProductToken< std::vector<sim::OpDetBacktrackerRecord> > fInputToken;
    double        fQE;
    double        fDarkNoiseRate;	      // In Hz
    double        fCrossTalk;           // Probability of SiPM producing 2 PE signal

    double        fTimeBegin;           // Earliest and latest possible times for reading out data
    double        fTimeEnd;             // Used for defining the dark noise range

    // unused double        fCorrection;          // Correct light yield. Kept separate for clarity
    bool          fCorrectLateLight;    // Do we apply a late light correction?
    double        fLateLightCorrection; // How much to correct the late light
    double        fLateLightBoundary;   // What is the boundary which defines late light?
    bool          fLateLightPerTrack;   // Is late light counted from each track's first photon (else the event's)?
    bool          fDarkNoiseTimeInNs;   // Are dark-count times converted to ns (else legacy: us stored as ns)?
    bool          fBinomialQE;          // Is detection binomial per photon (else legacy: Poisson on the mean)?

    // Random number engines
    CLHEP::HepRandomEngine& fSIPMEngine;
    CLHEP::RandExponential  fRandExponential;
    CLHEP::RandFlat         fRandFlat;
    CLHEP::RandPoissonQ     fRandPoissPhot;
    CLHEP::RandBinomial     fRandBinomPhot;
  
    // Produce waveform on one of the optical detectors
    // These cannot be const due to the random number generators
    void PhotonsToPE(sim::OpDetBacktrackerRecord const& btr, sim::OpDetDivRec& dr_plusnoise);
    void AddDarkNoise(sim::OpDetDivRec &);
    unsigned short CrossTalk();

    void SetBeginEndTimes();

  };
}

#endif

namespace opdet {

  DEFINE_ART_MODULE(SIPMOpSensorSim)

}

namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  SIPMOpSensorSim::SIPMOpSensorSim(Parameters const & config)
    : art::EDProducer{config}
    , fInputTag{     config().InputTag()}
    , fInputToken{   consumes< std::vector<sim::OpDetBacktrackerRecord> >(fInputTag) }
    , fDarkNoiseRate{config().DarkNoiseRate()}
    , fCrossTalk{    config().CrossTalk()}
    , fSIPMEngine(
        art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
          createEngine(0, "HepJamesRandom", "sipm"),
          "HepJamesRandom",
          "sipm", 
          config.get_PSet(), 
          "SeedSiPM"))
    , fRandExponential(fSIPMEngine)
    , fRandFlat(fSIPMEngine)
    , fRandPoissPhot(fSIPMEngine)
    , fRandBinomPhot(fSIPMEngine)
  {

    // This module produces (infrastructure piece)
    produces< std::vector< sim::OpDetDivRec > >();


    if (fDarkNoiseRate < 0.0) {
      throw art::Exception(art::errors::Configuration) 
        << "fDarkNoiseRate: " << fDarkNoiseRate << '\n'
        << "Dark noise rate should be non-negative!\n";
    }


    double tempQE = config().QuantumEfficiency() * config().Correction();
    if (tempQE < 0 || tempQE > 1) {
      throw art::Exception(art::errors::Configuration) 
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "       leading to total efficiency: " << tempQE << "\n"
        << "It should be between 0 and 1!\n";
    }

    // Correct out the prescaling applied during simulation
    auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fQE = tempQE / LarProp->ScintPreScale();

    if (fQE > 1.0001 ) {
      throw art::Exception(art::errors::Configuration)
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "       leading to total efficiency: " << tempQE << "\n"
        << "is too large.\n"
        << "It is larger than the prescaling applied during simulation, " << LarProp->ScintPreScale() << ".\n"
        << "Final QE must be equal to or smaller than the QE applied at simulation time.\n";
    }

    //This will be handled on the WaveformDigitizerSim_module.cc 
    // Check for non-trivial channel mapping which is not supported
/*    art::ServiceHandle< geo::Geometry > geometry;
    for (unsigned int opDet = 0; opDet < geometry->NOpDets() ; ++opDet) {
      if (geometry->NOpHardwareChannels(opDet) > 1)
        throw art::Exception(art::errors::Configuration)
          << "OpDet #" << opDet << " has " << geometry->NOpHardwareChannels(opDet) 
          << " channels associated with it. \n"
          << "This kind of channel mapping is not supported by SIPMOpSensorSim.\n"
          << "You need to use the legacy OpDetDigitizerDUNE instead.\n";
    }*/

    // Set time ranges if needed for dark noise
    if (fDarkNoiseRate > 0) SetBeginEndTimes();

    // Check for optional parameters for adjusting late light.
    // Apply adjustments if both are given.
    // Throw an error if only one is specified.
    bool haveCorrection = config().LateLightCorrection(fLateLightCorrection);
    bool haveBoundary   = config().LateLightBoundary(fLateLightBoundary);
    fCorrectLateLight = haveBoundary && haveCorrection;

    if (haveBoundary != haveCorrection) {
      throw art::Exception(art::errors::Configuration)
        << "Must specify both LateLightCorrection and LateLightBoundary in order to apply correction.\n"
        << "Leave both out to not adjust late light.\n";
    }

    std::string const lateLightReference = config().LateLightReference();
    if (lateLightReference != "Event" && lateLightReference != "Track") {
      throw art::Exception(art::errors::Configuration)
        << "LateLightReference: \"" << lateLightReference << "\"\n"
        << "Must be \"Event\" or \"Track\".\n";
    }
    fLateLightPerTrack = (lateLightReference == "Track");
    fDarkNoiseTimeInNs = config().DarkNoiseTimeInNs();
    fBinomialQE        = config().BinomialQE();

    if (fCorrectLateLight && fLateLightCorrection*tempQE > LarProp->ScintPreScale())
    {
      throw art::Exception(art::errors::Configuration)
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "         and late light correction: " << fLateLightCorrection << "\n"
        << "       leading to total efficiency: " << tempQE * fLateLightCorrection << "\n"
        << "is too large.\n"
        << "It is larger than the prescaling applied during simulation, " << LarProp->ScintPreScale() << ".\n"
        << "Final QE must be equal to or smaller than the QE applied at simulation time.\n";
    }
  }


  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::produce(art::Event& event)
  {
    // A pointer that will store produced OpDetDivRec
    auto OpDetDivRecPtr = std::make_unique< std::vector< sim::OpDetDivRec > >();

    // Get OpDetBacktrackerRecord from the event
    auto const & btr_handle = event.getValidHandle(fInputToken);

    // For every optical detector:
    for (auto const& btr : *btr_handle) {
      int opDet = btr.OpDetNum();

      auto DivRecPlusNoise = sim::OpDetDivRec(opDet);

      PhotonsToPE(btr, DivRecPlusNoise);

      // Generate dark noise
      if (fDarkNoiseRate > 0.0) AddDarkNoise(DivRecPlusNoise);

      // AddAfterPulsing();

      OpDetDivRecPtr->emplace_back(DivRecPlusNoise);
    }

    event.put(std::move(OpDetDivRecPtr));
  }

  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::PhotonsToPE(sim::OpDetBacktrackerRecord const& btr,
                                    sim::OpDetDivRec& dr_plusnoise)
  {
    // Don't do anything without any records
    if (btr.timePDclockSDPsMap().size() == 0)
      return;

    // Late light is light arriving more than fLateLightBoundary after a
    // reference time on this OpDet. The BTR covers the whole event, so with
    // pile-up the first photon of the event usually comes from an earlier
    // interaction, and all light of the later ones, prompt included, counts
    // as late. Per track, each track's light is counted from that track's
    // first photon: tracks of one interaction start within ns of each other,
    // and delayed processes (captures, decays) get their own reference.
    double firstTime = 0.;
    std::unordered_map<int, double> firstTimeOfTrack;
    if (fCorrectLateLight) {
      firstTime = btr.timePDclockSDPsMap()[0].first;
      for (auto const& [time, sdps]: btr.timePDclockSDPsMap()) {
        if (time < firstTime) firstTime = time;
        if (!fLateLightPerTrack) continue;
        for (auto const& sdp : sdps) {
          auto const [it, inserted] = firstTimeOfTrack.try_emplace(sdp.trackID, time);
          if (!inserted && time < it->second) it->second = time;
        }
      }
    }

    // Loop through times in vector< pair< arrival time (in ns), vector< SDP > > >
    for (auto const& [time, sdps]: btr.timePDclockSDPsMap()) {

      double lateScale = 1.;
      if (fCorrectLateLight && !fLateLightPerTrack) {
        if (time > firstTime + fLateLightBoundary)
          lateScale = fLateLightCorrection;
      }

      // Loop through SDPs
      for(auto const& sdp : sdps) {

        if (fCorrectLateLight && fLateLightPerTrack)
          lateScale = (time > firstTimeOfTrack.at(sdp.trackID) + fLateLightBoundary) ? fLateLightCorrection : 1.;

        // Reduce true photons by QE. sdp.numPhotons is already a Poisson draw
        // (PDFastSim), so each photon is detected with probability QE; a second
        // Poisson draw on the mean would give Var/mean = 1+QE instead of 1-QE.
        int nphot = 0;
        if (fBinomialQE) {
          double const p = std::min(1., fQE * lateScale);
          if (sdp.numPhotons > 0)
            nphot = static_cast<int>(fRandBinomPhot.fire(sdp.numPhotons, p));
        }
        else
          nphot = fRandPoissPhot.fire(fQE * (double)sdp.numPhotons * lateScale);

        // For each true photon detected
        for(int truePh=0; truePh<nphot; ++truePh) {

          // Determine actual PE with Cross Talk
          unsigned int PE = 1+CrossTalk();
          for(unsigned int i = 0; i < PE; i++) {
            // Add to collection
            dr_plusnoise.AddPhoton(btr.OpDetNum(), // Channel
                                   sdp.trackID,    // TrackID
                                   time);          // Time
          }
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  unsigned short SIPMOpSensorSim::CrossTalk()
  {
    // Use recursion to allow cross talk to create cross talk
    if      (fCrossTalk <= 0.0)                 return 0;
    else if (fRandFlat.fire(1.0) > fCrossTalk)  return 0;
    else                                        return 1+CrossTalk();
  }

  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::AddDarkNoise(sim::OpDetDivRec & dr_plusnoise)
  {

    double tau = 1.0/fDarkNoiseRate*1.e6; // Typical time between Dark counts in us
                                          // Multiply by 10^6 since fDarkNoiseRate is in Hz

    double darkNoiseTime = fRandExponential.fire(tau) + fTimeBegin;
    while (darkNoiseTime < fTimeEnd) {
	    //Currently using all zeros as an indicator of Dark Noise for the trackID.
      int PE = 1+CrossTalk();
      for(int j = 0; j < PE; j++) {
        // darkNoiseTime is in us (fTimeBegin, fTimeEnd, tau); DivRec times are ns
        dr_plusnoise.AddPhoton(dr_plusnoise.OpDetNum(), 0,
                               fDarkNoiseTimeInNs ? darkNoiseTime*1000. : darkNoiseTime);
      }

      // Find next time to simulate a single PE pulse
      darkNoiseTime += fRandExponential.fire(tau);
    }
  }

  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::SetBeginEndTimes()
  {

    // Get the window start/end time from the detector properties services
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    auto const geometry  = art::ServiceHandle< geo::Geometry >();

    double maxDrift = 0.0;
    for (geo::TPCGeo const& tpc : geometry->Iterate<geo::TPCGeo>())
      if (maxDrift < tpc.DriftDistance()) maxDrift = tpc.DriftDistance();

    // Start at -1 drift window
    fTimeBegin = -1*maxDrift/detProp.DriftVelocity();

    // Take the TPC readout window size and convert
    // to us with the electronics clock frequency
    fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();
  }


} // end opdet namespace 
