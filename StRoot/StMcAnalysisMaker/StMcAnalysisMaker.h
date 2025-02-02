#ifndef ST_MCANALYSISMAKER_H
#define ST_MCANALYSISMAKER_H

#include <vector>
#include "TString.h"
#include "TDatabasePDG.h"

#include"StThreeVectorF.hh"

#include "StChain/StMaker.h"

class TFile;
class TNtuple;
class TH3F;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StMcEvent;
class StEvent;
class StMuDst;
class StDedxPidTraits;

class StMcAnalysisMaker : public StMaker
{
private:
   TString mOutfileName;
   StMuDst*       mMuDst;
   //std::vector<float> firedTriggersIndices;
   double mField; //.. magnetic field
   int    mCentrality;
   bool mFillTpcHitsNtuple;

   TFile* mFile;
   TNtuple* mTracks;
   TNtuple* mEventCount; //.. For counting purposes
   TNtuple* mTpcNtuple;
   TNtuple* mLambda;
   TNtuple* mLambdaRC;
   
   TDatabasePDG *mPDGdata;

   TH3F* hTpcHitsDiffXVsPadrowVsSector;                                                                                                                                                                  
   TH3F* hTpcHitsDiffYVsPadrowVsSector;
   TH3F* hTpcHitsDiffZVsPadrowVsSector;
   
   TH1F *nLambdas_hist;
   TH1F *nLambdaBars_hist;
   
   TH1F *nLambdaMothers_hist;
   TH1F *nLambdaBarMothers_hist;

   StMcEvent* mMcEvent;
   StEvent* mEvent;
   StAssociationMaker* mAssoc;

   StTrack const* findPartner(StMcTrack*, int&) const;
   StMcTrack const* findPartner(StGlobalTrack*, int&) const;
   StDedxPidTraits const* findDedxPidTraits(StTrack const*) const;
   int getNHitsDedx(StTrack const*) const;

   //bool passTrigger();
   int  fillEventCounts(float nRTracks = -1, float nMcTracks = -1);
   int  fillTracks(int& nRTracks, int& nMcTracks);
   void fillMcTrack(float* array,int& idx,StMcTrack const*);
   void fillRcTrack(float* array,int& idx,StMcTrack const*,StTrack const*,int const ncom);
   
   int fillLambdas_direct();
   int  fillLambdas_pairs(); //find MC Lambdas and match them to RC
   int fillRCLambdas(); //fill RC p-pi pairs, like as is done in data. This includes combinatorial background.
   //void fillMCLambdas( float* array, StMcTrack const* );
   //void fillRCLambdas(float* array,int& idx,StMcTrack const*,StTrack const*,int const ncom);
   
   void getDca(StTrack const*,float& dca, float& dcaXY, float& dcaZ) const;
   
   float getPairDca(StTrack const* trk1, StTrack const* trk2) const;
   float getPairDcaStraightLine(StTrack const* trk1, StTrack const* trk2) const;
   
   StThreeVectorF getSecondaryVertex(StTrack const* trk1, StTrack const* trk2) const;
   StThreeVectorF getSecondaryVertexStraightLine(StTrack const* trk1, StTrack const* trk2) const;
   
   StThreeVectorF getMomentumAt(StTrack const* trk, StThreeVectorF vertex) const;

   bool isGoodMcTrack(StMcTrack const*) const;
   
   bool isGoodRcPion(StTrack const*, int nCommonHits) const;
   bool isGoodRcProton(StTrack const*, int nCommonHits) const;

   void fillTpcNtuple(StMcTrack const* const,StTrack const* const);
   
   int mJobIndex;

public:
   StMcAnalysisMaker (const char *name="StMcAnalysisMaker", const char *title="event/StMcAnalysisMaker");

   int Init();
   int Make();
   int Finish();

   void setOutFileName(std::string);
   void fillTpcHitsNtuple(bool t=true);
   
   void setJobIndex(int jobIndex); //to set uniqure eventId for each submitted job

   ClassDef(StMcAnalysisMaker, 0)
};

inline void StMcAnalysisMaker::fillTpcHitsNtuple(bool t){ mFillTpcHitsNtuple=t;}
inline void StMcAnalysisMaker::setOutFileName(std::string s){ mOutfileName = s.c_str();}
inline void StMcAnalysisMaker::setJobIndex( int jobIndex ){ mJobIndex = jobIndex; }

#endif
