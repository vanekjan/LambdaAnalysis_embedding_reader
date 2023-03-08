#include "TFile.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

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

using namespace std;

ClassImp(StMcAnalysisMaker);

StMcAnalysisMaker::StMcAnalysisMaker(const char *name, const char *title): StMaker(name),
   mMuDst(NULL), mField(-999), mFillTpcHitsNtuple(false), 
   mFile(NULL), mTracks(NULL), mLambda(NULL),mEventCount(NULL), mMcEvent(NULL), mEvent(NULL), mAssoc(NULL)
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
                         
   mLambda = new TNtuple("Lambda", "Lambda", "px_pi:py_pi:pz_pi:E_pi:y_pi:" // MC - pi
                         "px_p:py_p:pz_p:E_p:y_p:" // MC - p
                         "px_L:py_L:pz_L:E_L:y_L:M_L:decL_L:Id:" // MC - Lambda
                         "gPx_pi:gPy_pi:gPz_pi:nFit_pi:nMax_pi:nCom_pi:nDedx_pi:dedx_pi:nSigPi_pi:nSigK_pi:nSigP_pi:dca_pi:dcaXY_pi:dcaZ_pi:" // RC - pi
                         "gPx_p:gPy_p:gPz_p:nFit_p:nMax_p:nCom_p:nDedx_p:dedx_p:nSigPi_p:nSigK_p:nSigP_p:dca_p:dcaXY_p:dcaZ_p:" // RC - p
                         "gPx_L:gPy_L:gPz_L:gDecL_L:gCosTheta_L:pairDCA_L:eventId"); // RC - L
                         
   mPDGdata = new TDatabasePDG();

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
   
   fillLambdas();
   
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
//    ( mcTrack->parent()->geantId() == McAnaCuts::motherGeantId_1 || mcTrack->parent()->geantId() == McAnaCuts::motherGeantId_2 )
   return  ( mcTrack->geantId() == McAnaCuts::geantId_1 || mcTrack->geantId() == McAnaCuts::geantId_2 || mcTrack->geantId() == McAnaCuts::geantId_3 || mcTrack->geantId() == McAnaCuts::geantId_4)
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

float StMcAnalysisMaker::getPairDca(StTrack const* trk1, StTrack const* trk2) const
{
  StPhysicalHelixD helix1 = trk1->geometry()->helix();
  StPhysicalHelixD helix2 = trk2->geometry()->helix();
  
  pair<double, double> const ss = helix1.pathLengths(helix2);
  
  StThreeVectorF dcaPoint1 = helix1.at(ss.first);
  StThreeVectorF dcaPoint2 = helix2.at(ss.second);

  return (dcaPoint1 - dcaPoint2).mag();

}

float StMcAnalysisMaker::getDecayLength(StTrack const* trk1, StTrack const* trk2) const
{
  StPhysicalHelixD helix1 = trk1->geometry()->helix();
  StPhysicalHelixD helix2 = trk2->geometry()->helix();
  
  pair<double, double> const ss = helix1.pathLengths(helix2);
  
  StThreeVectorF dcaPoint1 = helix1.at(ss.first);
  StThreeVectorF dcaPoint2 = helix2.at(ss.second);
  
  StThreeVectorF secondaryVertex = (dcaPoint1 + dcaPoint2) * 0.5;

  return secondaryVertex.mag();

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

int StMcAnalysisMaker::fillLambdas()
{
  //vectors of pions and p
  //save index of particle in MC event
  vector<int> pi_plus_index;
  vector<int> pi_minus_index;
  
  vector<int> p_index;
  vector<int> p_bar_index;

  //first fill pi and p vectors
  for (unsigned int iTrk = 0;  iTrk < mMcEvent->tracks().size(); ++iTrk)
  {
    StMcTrack* const mcTrack = mMcEvent->tracks()[iTrk];

    if (!mcTrack)
    {
       LOG_WARN << "Empty mcTrack container" << endm;
       continue;
    }

    if (!isGoodMcTrack(mcTrack)) continue;
    
    
    if( mcTrack->geantId() == 8 ) pi_plus_index.push_back(iTrk);
    if( mcTrack->geantId() == 9 ) pi_minus_index.push_back(iTrk);
    
    if( mcTrack->geantId() == 14 ) p_index.push_back(iTrk);
    if( mcTrack->geantId() == 15 ) p_bar_index.push_back(iTrk);
    
  
  }//end for loop over MC tracks


  //analyze vectors
  //find pi and p from decays of Lambda
  //save info olny for MC Lambdas
  //then find corresponding RC tracks
  
  float L_tree_vars[100];  
  for(unsigned int i = 0; i < 100; i++) L_tree_vars[i] = -999; //default values
  
  
  TParticlePDG *Lambda_PDG = mPDGdata->GetParticle(3122);
  
  TParticlePDG *pi_plus_PDG = mPDGdata->GetParticle(211);
  TParticlePDG *pi_minus_PDG = mPDGdata->GetParticle(-211);
  
  TParticlePDG *p_PDG = mPDGdata->GetParticle(2212);
  TParticlePDG *p_bar_PDG = mPDGdata->GetParticle(-2212);
  
  //Lambda
  
  for(unsigned int i_p = 0; i_p < p_index.size(); i_p++)
  {
    StMcTrack* const proton = mMcEvent->tracks()[p_index.at(i_p)];
    
    TLorentzVector p_4mom;
    //p_4mom.SetPxPyPzE(proton->momentum().x(), proton->momentum().y(), proton->momentum().z(), proton->energy());
    p_4mom.SetXYZM(proton->momentum().x(), proton->momentum().y(), proton->momentum().z(), p_PDG->Mass());
    
    int nCommonHits_p = 0;
    StTrack const* const rcProton = findPartner(proton, nCommonHits_p);
    
    
    for(unsigned int i_pi = 0; i_pi < pi_minus_index.size(); i_pi++)
    {
      StMcTrack* const pion = mMcEvent->tracks()[pi_minus_index.at(i_pi)];
      
      TLorentzVector pi_4mom;
      //pi_4mom.SetPxPyPzE(pion->momentum().x(), pion->momentum().y(), pion->momentum().z(), pion->energy());
      pi_4mom.SetXYZM(pion->momentum().x(), pion->momentum().y(), pion->momentum().z(), pi_minus_PDG->Mass());
      
      int nCommonHits_pi = 0;
      StTrack const* const rcPion = findPartner(pion, nCommonHits_pi);
      
      
      TLorentzVector L_4mom = p_4mom+pi_4mom;
      
      if( fabs(L_4mom.M() - Lambda_PDG->Mass()) < 1e-5  )
      //if( fabs(proton->startVertex()->position().mag() - pion->startVertex()->position().mag()) < 1e-5 && fabs(L_4mom.M() - Lambda_PDG->Mass()) < 1e-5 )
      {
        //cout<<"Lambda found!"<<endl;
        
        int iVar = 0;
        
        L_tree_vars[iVar++] = pion->momentum().x();
        L_tree_vars[iVar++] = pion->momentum().y();
        L_tree_vars[iVar++] = pion->momentum().z();
        L_tree_vars[iVar++] = pion->energy();
        L_tree_vars[iVar++] = pion->rapidity();
        
        L_tree_vars[iVar++] = proton->momentum().x();
        L_tree_vars[iVar++] = proton->momentum().y();
        L_tree_vars[iVar++] = proton->momentum().z();
        L_tree_vars[iVar++] = proton->energy();
        L_tree_vars[iVar++] = proton->rapidity();
        
        L_tree_vars[iVar++] = L_4mom.Px();
        L_tree_vars[iVar++] = L_4mom.Py();
        L_tree_vars[iVar++] = L_4mom.Pz();
        L_tree_vars[iVar++] = L_4mom.E();
        L_tree_vars[iVar++] = L_4mom.Rapidity();
        L_tree_vars[iVar++] = L_4mom.M();
        L_tree_vars[iVar++] = pion->startVertex()->position().mag(); //L MC decay vertex based on pion origin
        
        L_tree_vars[iVar++] = 3122; //L PDG id
        
               
        //get RC partner
               
        if(rcPion)
        {
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().x();
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().y();
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().z();
          L_tree_vars[iVar++] = rcPion->fitTraits().numberOfFitPoints(kTpcId);
          L_tree_vars[iVar++] = rcPion->numberOfPossiblePoints(kTpcId);
          L_tree_vars[iVar++] = nCommonHits_pi;
          
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
         StParticleDefinition const* prtcl = rcPion->pidTraits(aplus);

         if (prtcl)
         {
            nDedxPts = aplus.traits()->numberOfPoints();
            dedx = aplus.traits()->mean();
            nSigPi = aplus.numberOfSigma(Pion);
            nSigK = aplus.numberOfSigma(Kaon);
            nSigP = aplus.numberOfSigma(Proton);
         }

         L_tree_vars[iVar++] = getNHitsDedx(rcPion);
         L_tree_vars[iVar++] = dedx;
         L_tree_vars[iVar++] = nSigPi;
         L_tree_vars[iVar++] = nSigK;
         L_tree_vars[iVar++] = nSigP;


         float dca = -999.;
         float dcaXY = -999.;
         float dcaZ = -999.;

         getDca(rcPion, dca, dcaXY, dcaZ);

         L_tree_vars[iVar++] = dca;
         L_tree_vars[iVar++] = dcaXY;
         L_tree_vars[iVar++] = dcaZ;                             
          
        
        }
        else
        {
          iVar += 14;        
        }
        //________________________________
        
        if(rcProton)
        {
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().x();
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().y();
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().z();
          L_tree_vars[iVar++] = rcProton->fitTraits().numberOfFitPoints(kTpcId);
          L_tree_vars[iVar++] = rcProton->numberOfPossiblePoints(kTpcId);
          L_tree_vars[iVar++] = nCommonHits_p;
          
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
         StParticleDefinition const* prtcl = rcProton->pidTraits(aplus);

         if (prtcl)
         {
            nDedxPts = aplus.traits()->numberOfPoints();
            dedx = aplus.traits()->mean();
            nSigPi = aplus.numberOfSigma(Pion);
            nSigK = aplus.numberOfSigma(Kaon);
            nSigP = aplus.numberOfSigma(Proton);
         }

         L_tree_vars[iVar++] = getNHitsDedx(rcProton);
         L_tree_vars[iVar++] = dedx;
         L_tree_vars[iVar++] = nSigPi;
         L_tree_vars[iVar++] = nSigK;
         L_tree_vars[iVar++] = nSigP;


         float dca = -999.;
         float dcaXY = -999.;
         float dcaZ = -999.;

         getDca(rcProton, dca, dcaXY, dcaZ);

         L_tree_vars[iVar++] = dca;
         L_tree_vars[iVar++] = dcaXY;
         L_tree_vars[iVar++] = dcaZ;                             
          
        
        }
        else
        {
          iVar += 14;        
        }
        
        if(rcPion && rcProton)
        {
        
          TVector3 pi_rc_mom(rcPion->geometry()->momentum().x(), rcPion->geometry()->momentum().y(), rcPion->geometry()->momentum().z());
          TVector3 p_rc_mom(rcProton->geometry()->momentum().x(), rcProton->geometry()->momentum().y(), rcProton->geometry()->momentum().z());
          
          TVector3 L_rc_mom = pi_rc_mom+p_rc_mom;
          
          L_tree_vars[iVar++] = L_rc_mom.X();
          L_tree_vars[iVar++] = L_rc_mom.Y();
          L_tree_vars[iVar++] = L_rc_mom.Z();
          
          L_tree_vars[iVar++] = getDecayLength(rcPion, rcProton);
          L_tree_vars[iVar++] = getPairDca(rcPion, rcProton);
        
        }
        else
        {        
          iVar += 5;        
        }
        
        L_tree_vars[iVar++] = mMcEvent->eventNumber();
        
        mLambda->Fill(L_tree_vars);

      }
    
    } //end for pions
  
  } //end for protons
  
  
  //L-bar
  for(unsigned int i_p = 0; i_p < p_bar_index.size(); i_p++)
  {
    StMcTrack* const proton = mMcEvent->tracks()[p_bar_index.at(i_p)];
    
    TLorentzVector p_4mom;
    //p_4mom.SetPxPyPzE(proton->momentum().x(), proton->momentum().y(), proton->momentum().z(), proton->energy());
    p_4mom.SetXYZM(proton->momentum().x(), proton->momentum().y(), proton->momentum().z(), p_PDG->Mass());
    
    int nCommonHits_p = 0;
    StTrack const* const rcProton = findPartner(proton, nCommonHits_p);
    
    
    for(unsigned int i_pi = 0; i_pi < pi_plus_index.size(); i_pi++)
    {
      StMcTrack* const pion = mMcEvent->tracks()[pi_plus_index.at(i_pi)];
      
      TLorentzVector pi_4mom;
      //pi_4mom.SetPxPyPzE(pion->momentum().x(), pion->momentum().y(), pion->momentum().z(), pion->energy());
      pi_4mom.SetXYZM(pion->momentum().x(), pion->momentum().y(), pion->momentum().z(), pi_minus_PDG->Mass());
      
      int nCommonHits_pi = 0;
      StTrack const* const rcPion = findPartner(pion, nCommonHits_pi);
      
      
      TLorentzVector L_4mom = p_4mom+pi_4mom;
      
      if( fabs(L_4mom.M() - Lambda_PDG->Mass()) < 1e-5  )
      //if( fabs(proton->startVertex()->position().mag() - pion->startVertex()->position().mag()) < 1e-5 && fabs(L_4mom.M() - Lambda_PDG->Mass()) < 1e-5 )
      {
        //cout<<"Lambda found!"<<endl;
        
        int iVar = 0;
        
        L_tree_vars[iVar++] = pion->momentum().x();
        L_tree_vars[iVar++] = pion->momentum().y();
        L_tree_vars[iVar++] = pion->momentum().z();
        L_tree_vars[iVar++] = pion->energy();
        L_tree_vars[iVar++] = pion->rapidity();
        
        L_tree_vars[iVar++] = proton->momentum().x();
        L_tree_vars[iVar++] = proton->momentum().y();
        L_tree_vars[iVar++] = proton->momentum().z();
        L_tree_vars[iVar++] = proton->energy();
        L_tree_vars[iVar++] = proton->rapidity();
        
        L_tree_vars[iVar++] = L_4mom.Px();
        L_tree_vars[iVar++] = L_4mom.Py();
        L_tree_vars[iVar++] = L_4mom.Pz();
        L_tree_vars[iVar++] = L_4mom.E();
        L_tree_vars[iVar++] = L_4mom.Rapidity();
        L_tree_vars[iVar++] = L_4mom.M();
        L_tree_vars[iVar++] = pion->startVertex()->position().mag(); //L MC decay vertex based on pion origin             
                
        L_tree_vars[iVar++] = -3122; //L-bar PDG id
        
               
        //get RC partner
               
        if(rcPion)
        {
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().x();
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().y();
          L_tree_vars[iVar++] = rcPion->geometry()->momentum().z();
          L_tree_vars[iVar++] = rcPion->fitTraits().numberOfFitPoints(kTpcId);
          L_tree_vars[iVar++] = rcPion->numberOfPossiblePoints(kTpcId);
          L_tree_vars[iVar++] = nCommonHits_pi;
          
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
         StParticleDefinition const* prtcl = rcPion->pidTraits(aplus);

         if (prtcl)
         {
            nDedxPts = aplus.traits()->numberOfPoints();
            dedx = aplus.traits()->mean();
            nSigPi = aplus.numberOfSigma(Pion);
            nSigK = aplus.numberOfSigma(Kaon);
            nSigP = aplus.numberOfSigma(Proton);
         }

         L_tree_vars[iVar++] = getNHitsDedx(rcPion);
         L_tree_vars[iVar++] = dedx;
         L_tree_vars[iVar++] = nSigPi;
         L_tree_vars[iVar++] = nSigK;
         L_tree_vars[iVar++] = nSigP;


         float dca = -999.;
         float dcaXY = -999.;
         float dcaZ = -999.;

         getDca(rcPion, dca, dcaXY, dcaZ);

         L_tree_vars[iVar++] = dca;
         L_tree_vars[iVar++] = dcaXY;
         L_tree_vars[iVar++] = dcaZ;                             
          
        
        }
        else
        {
          iVar += 14;        
        }
        //________________________________
        
        if(rcProton)
        {
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().x();
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().y();
          L_tree_vars[iVar++] = rcProton->geometry()->momentum().z();
          L_tree_vars[iVar++] = rcProton->fitTraits().numberOfFitPoints(kTpcId);
          L_tree_vars[iVar++] = rcProton->numberOfPossiblePoints(kTpcId);
          L_tree_vars[iVar++] = nCommonHits_p;
          
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
         StParticleDefinition const* prtcl = rcProton->pidTraits(aplus);

         if (prtcl)
         {
            nDedxPts = aplus.traits()->numberOfPoints();
            dedx = aplus.traits()->mean();
            nSigPi = aplus.numberOfSigma(Pion);
            nSigK = aplus.numberOfSigma(Kaon);
            nSigP = aplus.numberOfSigma(Proton);
         }

         L_tree_vars[iVar++] = getNHitsDedx(rcProton);
         L_tree_vars[iVar++] = dedx;
         L_tree_vars[iVar++] = nSigPi;
         L_tree_vars[iVar++] = nSigK;
         L_tree_vars[iVar++] = nSigP;


         float dca = -999.;
         float dcaXY = -999.;
         float dcaZ = -999.;

         getDca(rcProton, dca, dcaXY, dcaZ);

         L_tree_vars[iVar++] = dca;
         L_tree_vars[iVar++] = dcaXY;
         L_tree_vars[iVar++] = dcaZ;     
         
                  
        
        }
        else
        {
          iVar += 14;        
        }
        
        if(rcPion && rcProton)
        {
        
          TVector3 pi_rc_mom(rcPion->geometry()->momentum().x(), rcPion->geometry()->momentum().y(), rcPion->geometry()->momentum().z());
          TVector3 p_rc_mom(rcProton->geometry()->momentum().x(), rcProton->geometry()->momentum().y(), rcProton->geometry()->momentum().z());
          
          TVector3 L_rc_mom = pi_rc_mom+p_rc_mom;
          
          L_tree_vars[iVar++] = L_rc_mom.X();
          L_tree_vars[iVar++] = L_rc_mom.Y();
          L_tree_vars[iVar++] = L_rc_mom.Z();
          
          L_tree_vars[iVar++] = getDecayLength(rcPion, rcProton);
          L_tree_vars[iVar++] = getPairDca(rcPion, rcProton);
        
        }
        else
        {        
          iVar += 5;
        }
        
        L_tree_vars[iVar++] = mMcEvent->eventNumber();
        
        mLambda->Fill(L_tree_vars);

      }
    
    } //end for pions
  
  } //end for protons
  
  
/*  
   "px_pi:py_pi:pz_pi:E_pi:y_pi:" // MC - pi
   "px_p:py_p:pz_p:E_p:y_p:" // MC - p
   "px_L:py_L:pz_L:E_L:y_L:Id:" // MC - Lambda
   "gPx_pi:gPy_pi:gPz_pi:nFit_pi:nMax_pi:nCom_pi:nDedx_pi:dedx_pi:nSigPi_pi:nSigK_pi:nSigP_pi:dca_pi:dcaXY_pi:dcaZ_pi:" // RC - pi
   "gPx_p:gPy_p:gPz_p:nFit_p:nMax_p:nCom_p:nDedx_p:dedx_p:nSigPi_p:nSigK_p:nSigP_p:dca_p:dcaXY_p:dcaZ_p:" // RC - p
   "gPx_L:gPy_L:gPz_L:eventId"
   
   array[idx++] = rcTrack->geometry()->momentum().perp();
   array[idx++] = rcTrack->geometry()->momentum().pseudoRapidity();
   array[idx++] = rcTrack->geometry()->momentum().phi();
   array[idx++] = rcTrack->fitTraits().numberOfFitPoints(kTpcId);
   array[idx++] = rcTrack->numberOfPossiblePoints(kTpcId);
   array[idx++] = ncom;
  
*/
}


