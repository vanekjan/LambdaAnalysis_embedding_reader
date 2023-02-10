#include "TFile.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TSystem.h"

#include "StarClassLibrary/StParticleDefinition.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StParticleTypes.hh"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEventTypes.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StTrackTopologyMap.h"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcEvent.hh"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"

#include "StMcAnaCuts.h"
#include "StMcAnalysisMaker.h"

ClassImp(StMcAnalysisMaker);

StMcAnalysisMaker::StMcAnalysisMaker(const char *name, const char *title): StMaker(name),
   mMuDst(NULL), mField(-999), mFillTpcHitsNtuple(false), 
   mFile(NULL), mTracks(NULL), mEventCount(NULL), mMcEvent(NULL), mEvent(NULL), mAssoc(NULL)
{
   LOG_INFO << "StMcAnalysisMaker() - DONE" << endm;
}
//__________________________________
int StMcAnalysisMaker::Init()
{
   if (!mOutfileName.Length())
   {
      // StBFChain* bfChain = (StBFChain *) StMaker::GetChain();
      //
      // if (!bfChain) return kStFatal;
      //
      // mOutfileName = bfChain->GetFileIn();

      if (mOutfileName.Length())
      {
         LOG_INFO << mOutfileName << endm;
         mOutfileName = gSystem->BaseName(mOutfileName.Data());
         mOutfileName = mOutfileName.ReplaceAll(".event.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".geant.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".MuDst.root", "");
      }
      else
      {
         mOutfileName = "mcAnalysis";
      }
   }

   mOutfileName = mOutfileName.ReplaceAll(".root", "");
   mFile = new TFile(Form("%s.McAna.root", mOutfileName.Data()), "recreate");
   assert(mFile && !mFile->IsZombie());

   mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");
   if (!mAssoc)
   {
      LOG_ERROR << "Could not get StAssociationMaker" << endm;
      return kStErr;
   }
/*
   for (int ii = 0; ii < McAnaCuts::maxNumberOfTriggers; ++ii)
   {
      firedTriggersIndices.push_back(-1);
   };
*/
   mEventCount = new TNtuple("eventCount", "eventCount", "runId:eventId:mcVx:mcVy:mcVz:vx:vy:vz:vzVpd:"
                             "zdc:bbc:nMcTracks:nRTracks:magField");

   mTracks = new TNtuple("tracks", "", "pt:p:eta:y:phi:geantId:eventGenLabel:startVtxX:startVtxY:startVtxZ:stopVtxX:stopVtxY:stopVtxZ:" // MC
                         "gPt:gEta:gPhi:nFit:nMax:nCom:nDedx:dedx:nSigPi:nSigK:nSigP:dca:dcaXY:dcaZ"); // global

   if (mFillTpcHitsNtuple)
   {
      hTpcHitsDiffXVsPadrowVsSector = new TH3F("hTpcHitsDiffXVsPadrowVsSector", "hTpcHitsDiffXVsPadrowVsSector", 50, 0, 50, 26, 0, 26, 200, -2, 2);
      hTpcHitsDiffYVsPadrowVsSector = new TH3F("hTpcHitsDiffYVsPadrowVsSector", "hTpcHitsDiffYVsPadrowVsSector", 50, 0, 50, 26, 0, 26, 200, -2, 2);
      hTpcHitsDiffZVsPadrowVsSector = new TH3F("hTpcHitsDiffZVsPadrowVsSector", "hTpcHitsDiffZVsPadrowVsSector", 50, 0, 50, 26, 0, 26, 200, -2, 2);

      mTpcNtuple = new TNtuple("tpc", "tpc", "eta:phi:mcX:mcY:mcZ:mcTb:rcX:rcY:rcZ:rcTb:sector:padrow:pad:dE:adc:mcHitToRcTrackX:mcHitToRcTrackY:mcHitToRcTrackZ");
   }

   LOG_INFO << "Init() - DONE" << endm;

   return kStOk;
}

//__________________________________
int StMcAnalysisMaker::Make()
{
   StMuDstMaker* muDstMaker = (StMuDstMaker*)GetMaker("MuDst");

   if (!muDstMaker)
   {
      LOG_WARN << " No MuDstMaker, will try to take all event information from StEvent" << endm;
      mMuDst = NULL;
   }
   else
   {
     mMuDst = (StMuDst*)muDstMaker->muDst();
   }

   if(!mMuDst || !mMuDst->event())
   {
     LOG_WARN << "MuDst or mMuDst->event() is missing, will try to take all event information from StEvent" << endm;
     mMuDst = NULL;
   }

   mMcEvent = (StMcEvent*)GetDataSet("StMcEvent");

   if (!mMcEvent)
   {
      LOG_WARN << "No StMcEvent" << endm;
      return kStWarn;
   }

   mEvent = (StEvent*)GetDataSet("StEvent");
   if (!mEvent)
   {
      LOG_WARN << "No StEvent" << endm;
      return kStWarn;
   }

   mField = (float)mEvent->summary()->magneticField();

   
   // Fill
   int nRTracks = -1;
   int nMcTracks = -1;

   int fillTracksStatus = kStOk;
   
   //for pp trigger check not needed?
   fillTracksStatus = fillTracks(nRTracks, nMcTracks);
   
/*   
   if (passTrigger())
   {
      fillTracksStatus = fillTracks(nRTracks, nMcTracks);
   }
   else
   {
      LOG_INFO << "No interesting triggers. Counting event then skipping." << endm;
   }
*/
   int fillEventCountStatus = fillEventCounts((float)nRTracks, (float)nMcTracks);

   return fillTracksStatus && fillEventCountStatus;
}

int StMcAnalysisMaker::fillTracks(int& nRTracks, int& nMcTracks)
{
   nRTracks = 0;
   nMcTracks = 0;

   LOG_INFO << "Filling " << mMcEvent->tracks().size() << " mcTracks..." << "\n";

   for (unsigned int iTrk = 0;  iTrk < mMcEvent->tracks().size(); ++iTrk)
   {
      StMcTrack* const mcTrack = mMcEvent->tracks()[iTrk];

      if (!mcTrack)
      {
         LOG_WARN << "Empty mcTrack container" << endm;
         continue;
      }

      if (!isGoodMcTrack(mcTrack)) continue;
      ++nMcTracks;

      int nCommonHits = 0;
      StTrack const* const rcTrack = findPartner(mcTrack, nCommonHits);

      float array[220];
      for (int ii = 0; ii < 220; ++ii) array[ii] = -999;

      int idx = 0;
      

      fillMcTrack(array, idx, mcTrack);

      if (rcTrack)
      {
         ++nRTracks;
         fillRcTrack(array, idx, mcTrack, rcTrack, nCommonHits);
         if (mFillTpcHitsNtuple) fillTpcNtuple(mcTrack, rcTrack);
      }

      mTracks->Fill(array);
   }

   LOG_INFO << endm;

   return kStOk;
}

void StMcAnalysisMaker::fillMcTrack(float* array, int& idx, StMcTrack const* const mcTrk)
{
   array[idx++] = mcTrk->pt();
   array[idx++] = mcTrk->momentum().mag();
   array[idx++] = mcTrk->pseudoRapidity();
   array[idx++] = mcTrk->rapidity();
   array[idx++] = mcTrk->momentum().phi();
   array[idx++] = mcTrk->geantId();
   array[idx++] = mcTrk->eventGenLabel();
   array[idx++] = mcTrk->startVertex()->position().x();
   array[idx++] = mcTrk->startVertex()->position().y();
   array[idx++] = mcTrk->startVertex()->position().z();

   if (mcTrk->stopVertex())
   {
      array[idx++] = mcTrk->stopVertex()->position().x();
      array[idx++] = mcTrk->stopVertex()->position().y();
      array[idx++] = mcTrk->stopVertex()->position().z();
   }
   else
   {
      idx += 3;
   }
}

void StMcAnalysisMaker::fillRcTrack(float* array, int& idx, StMcTrack const* const mcTrack, StTrack const* const rcTrack, int const ncom)
{
   if (!rcTrack) return;

   array[idx++] = rcTrack->geometry()->momentum().perp();
   array[idx++] = rcTrack->geometry()->momentum().pseudoRapidity();
   array[idx++] = rcTrack->geometry()->momentum().phi();
   array[idx++] = rcTrack->fitTraits().numberOfFitPoints(kTpcId);
   array[idx++] = rcTrack->numberOfPossiblePoints(kTpcId);
   array[idx++] = ncom;

   // dedx info
   float nDedxPts = -9999;
   float dedx = -9999;
   float nSigPi = -9999;
   float nSigK = -9999;
   float nSigP = -9999;
   static StTpcDedxPidAlgorithm aplus(McAnaCuts::dedxMethod);
   static StPionPlus* Pion = StPionPlus::instance();
   static StKaonPlus* Kaon = StKaonPlus::instance();
   static StProton* Proton = StProton::instance();
   StParticleDefinition const* prtcl = rcTrack->pidTraits(aplus);

   if (prtcl)
   {
      nDedxPts = aplus.traits()->numberOfPoints();
      dedx = aplus.traits()->mean();
      nSigPi = aplus.numberOfSigma(Pion);
      nSigK = aplus.numberOfSigma(Kaon);
      nSigP = aplus.numberOfSigma(Proton);
   }

   array[idx++] = getNHitsDedx(rcTrack);
   array[idx++] = dedx;
   array[idx++] = nSigPi;
   array[idx++] = nSigK;
   array[idx++] = nSigP;


   float dca = -999.;
   float dcaXY = -999.;
   float dcaZ = -999.;

   getDca(rcTrack, dca, dcaXY, dcaZ);

   array[idx++] = dca;
   array[idx++] = dcaXY;
   array[idx++] = dcaZ;
   
}

void StMcAnalysisMaker::fillTpcNtuple(StMcTrack const* const mcTrack, StTrack const* const rcTrack)
{
   if (!rcTrack) return;

   StPtrVecHit rcTpcHits = rcTrack->detectorInfo()->hits(kTpcId);
   StPhysicalHelixD helix = rcTrack->geometry()->helix();

   for (size_t ih = 0; ih < rcTpcHits.size(); ++ih)
   {
      StTpcHit* rcHit = dynamic_cast<StTpcHit*>(rcTpcHits[ih]);
      if (!rcHit) continue;

      pair<rcTpcHitMapIter, rcTpcHitMapIter>  bounds = mAssoc->rcTpcHitMap()->equal_range(rcHit);

      // loop over all mcHits associated with this rcHit
      bool found = false;
      for (rcTpcHitMapIter iter = bounds.first; iter != bounds.second; iter++)
      {
         const StMcTpcHit* mcHit = (*iter).second;

         // fill histograms if this mcHit belongs to this mcTrack
         if (mcHit->parentTrack()->key() == mcTrack->key())
         {
            StThreeVectorF mcHitTotRcTrack = helix.at(helix.pathLength(mcHit->position().x(), mcHit->position().y())) - mcHit->position();

            float arr[50];
            int iArr = 0;
            arr[iArr++] = mcTrack->pseudoRapidity();
            arr[iArr++] = mcTrack->momentum().phi();
            arr[iArr++] = mcHit->position().x();
            arr[iArr++] = mcHit->position().y();
            arr[iArr++] = mcHit->position().z();
            arr[iArr++] = mcHit->timeBucket();
            arr[iArr++] = rcHit->position().x();
            arr[iArr++] = rcHit->position().y();
            arr[iArr++] = rcHit->position().z();
            arr[iArr++] = rcHit->timeBucket();
            arr[iArr++] = mcHit->sector();
            arr[iArr++] = mcHit->padrow();
            arr[iArr++] = mcHit->pad();
            arr[iArr++] = mcHit->dE();
            arr[iArr++] = rcHit->adc();
            arr[iArr++] = mcHitTotRcTrack.x();
            arr[iArr++] = mcHitTotRcTrack.y();
            arr[iArr++] = mcHitTotRcTrack.z();

            mTpcNtuple->Fill(arr);

            hTpcHitsDiffXVsPadrowVsSector->Fill(mcHit->padrow(), mcHit->sector(), mcHit->position().x() - rcHit->position().x());
            hTpcHitsDiffYVsPadrowVsSector->Fill(mcHit->padrow(), mcHit->sector(), mcHit->position().y() - rcHit->position().y());
            hTpcHitsDiffZVsPadrowVsSector->Fill(mcHit->padrow(), mcHit->sector(), mcHit->position().z() - rcHit->position().z());
            found = true;
         }
         else
         {
            cout << "Not a good candidate" << endl;
         }
      }

      if (!found)
      {
         LOG_WARN << "No mc hit was found for this rc Hit!!!!" << endm;
      }
   }
}

bool StMcAnalysisMaker::isGoodMcTrack(StMcTrack const* const mcTrack) const
{
   return ( mcTrack->parent()->geantId() == McAnaCuts::motherGeantId_1 || mcTrack->parent()->geantId() == McAnaCuts::motherGeantId_2 )
          && mcTrack->geantId() == McAnaCuts::geantId_1
          && mcTrack->startVertex()->position().perp() < McAnaCuts::mcTrackStartVtxR;
}

int StMcAnalysisMaker::fillEventCounts(float nRTracks, float nMcTracks)
{
   float vars[50];

   float vpdVz = -999;
   StBTofHeader* tofheader = 0;
   if (mEvent->btofCollection())  tofheader = mEvent->btofCollection()->tofHeader();
   if (tofheader) vpdVz = tofheader->vpdVz();

   int iVar = 0;
   vars[iVar++] = (float)mEvent->runId();
   vars[iVar++] = (float)mEvent->id();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().z();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().z();
   vars[iVar++] = vpdVz;
   vars[iVar++] = (float)mEvent->runInfo()->zdcCoincidenceRate();
   vars[iVar++] = (float)mEvent->runInfo()->bbcCoincidenceRate();
   vars[iVar++] = nMcTracks;
   vars[iVar++] = nRTracks;
   vars[iVar++] = (float)mEvent->summary()->magneticField() / 10;

   mEventCount->Fill(vars);

   return kStOk;
}
/*
bool StMcAnalysisMaker::passTrigger()
{
   LOG_INFO << "Checking triggers..." << endm;
   bool interesting_event = false;

   if (!mEvent)
   {
      LOG_FATAL << "mEvent doesn't exist" << endm;
   }

   if (McAnaCuts::interesting_triggers.size() == 0)
   {
      LOG_WARN << "No triggers in McAnaCuts::interesting_triggers ... accepting event anyway" << endm;
      return true;
   }

   const StTriggerId* st_trgid = mEvent->triggerIdCollection()->nominal();

   for (unsigned int ii = 0; ii < firedTriggersIndices.size(); ++ii)
   {
      firedTriggersIndices[ii] = -1;
   }

   // Fill interesting triggers
   LOG_INFO << "Interesting fired triggers: " << "\n";

   for (unsigned int ii = 0; ii < st_trgid->maxTriggerIds(); ++ii)
   {
      unsigned int id = st_trgid->triggerId(ii);

      int trgIndex = -1;

      for (unsigned int jj = 0; jj < McAnaCuts::interesting_triggers.size(); ++jj)
      {
         if (McAnaCuts::interesting_triggers[jj] == id)
         {
            trgIndex = jj;
            interesting_event = true;
            LOG_INFO << id << " ";
            break;
         }
      }

      if (trgIndex != -1) firedTriggersIndices[trgIndex] = 1.0;
   }

   LOG_INFO << endm;

   return interesting_event;
}
*/
StTrack const* StMcAnalysisMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits) const
{
   //..StMcTrack find partner from the StTracks
   pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);

   const StTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (mcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();
      const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

void StMcAnalysisMaker::getDca(StTrack const* const rcTrack, float& dca, float& dcaXY, float& dcaZ) const
{
   StPhysicalHelixD helix = rcTrack->geometry()->helix();
   dca = helix.distance(mEvent->primaryVertex()->position());
   dcaXY = helix.geometricSignedDistance(mEvent->primaryVertex()->position().x(), mEvent->primaryVertex()->position().y());

   StThreeVectorF dcaPoint = helix.at(helix.pathLength(mEvent->primaryVertex()->position()));
   dcaZ = dcaPoint.z() - mEvent->primaryVertex()->position().z();
}


StMcTrack const* StMcAnalysisMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits) const
{
   //.. StGlobalTracks find partner from StMcTracks.
   //.. See example from StRoot/StMcAnalysisMaker
   pair<rcTrackMapIter, rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);

   const StMcTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (rcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();

      const StMcTrack* track = k->second->partnerMcTrack();

      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

int StMcAnalysisMaker::Finish()
{
   mFile->cd();
   mFile->Write();
   mFile->Close();
   return kStOk;
}

StDedxPidTraits const* StMcAnalysisMaker::findDedxPidTraits(StTrack const* const rcTrack) const
{
   StDedxPidTraits* pid = 0;
   StPtrVecTrackPidTraits traits = rcTrack->pidTraits(kTpcId);

   for (unsigned int ii = 0; ii < traits.size(); ++ii)
   {
      pid = dynamic_cast<StDedxPidTraits*>(traits[ii]);
      if (pid && pid->method() == McAnaCuts::dedxMethod) break;
   }

   return pid;
}

int StMcAnalysisMaker::getNHitsDedx(StTrack const* const t) const
{
   int ndedx = -1;
   StPtrVecTrackPidTraits pidTraits = t->pidTraits(kTpcId);

   if (pidTraits.size())
   {
      StDedxPidTraits* pid;
      for (unsigned int ii = 0; ii < pidTraits.size(); ++ii)
      {
         pid = dynamic_cast<StDedxPidTraits*>(pidTraits[ii]);

         if (pid && (pid->method() == McAnaCuts::dedxMethod))
         {
            ndedx = pid->numberOfPoints();            //number of dedx hits
            break;
         }
      }
   }

   return ndedx;
}
