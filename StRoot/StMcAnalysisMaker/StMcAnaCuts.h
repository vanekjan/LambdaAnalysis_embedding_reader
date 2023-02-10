#ifndef StMcAnaCuts_H
#define StMcAnaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include <vector>

#include "Rtypes.h"
#include "StEvent/StEnumerations.h"

namespace McAnaCuts
{
/*
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    
    //Run16 Au+Au MB triggerIds
    t.push_back(520001);
    t.push_back(520011); 
    t.push_back(520021); 
    t.push_back(520031); 
    t.push_back(520041);
    t.push_back(520051);
    
    

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();
*/

  float const mcTrackStartVtxR = 30.0; // maximum
  
  
  int const geantId_1 = 8; //pi+
  int const geantId_2 = 14; //p
  
  
  int const motherGeantId_1 = 18; //Lambda
  int const motherGeantId_2 = 26; //Lambda-bar

  StDedxMethod dedxMethod = kLikelihoodFitId;

  //int const maxNumberOfTriggers = 6;
}
#endif
