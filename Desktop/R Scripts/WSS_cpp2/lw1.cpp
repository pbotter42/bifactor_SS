#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]

List lw1(List allTlines,
            int nItems) {
  // part two the lw alg
  List lwIter(nItems); // global list containing all the lw iterations
  List lwFinalIter(1); // the final iteration 
  
  // startrp first iter info
  List rpList(nItems); // list containing the rp
  List rpListTemp1(2);
  std::vector<double> tempQ;
  std::vector<double> tempP;
  tempQ.push_back(as<List>(as<List>(as<List>(allTlines[0])[0]))[1]);
  tempP.push_back(as<List>(as<List>(as<List>(allTlines[0])[1]))[1]);
  rpListTemp1[0] = tempQ;
  rpListTemp1[1] = tempP;  
  rpList[0] = rpListTemp1;
  // end rp first iter info
  
  /*
  I will need to go back eventually an modify the lwIter and lwFinalIter relationship
  eventually to save memory. with something like a switch 
  */
  lwIter[0] = allTlines[0];
  for(int nLwIter = 1; nLwIter < nItems; ++nLwIter) {
    
    List prevIter = lwIter[nLwIter-1];
    int prevIterSize = prevIter.size();
    List lwIterTemp(prevIterSize*2);
    
    // start rp
    List rpPrev = rpList[nLwIter-1];
    List rpCurrent(prevIterSize*2);
    // end rp
    
    //lwIter[nLwIter] = lwIterTemp;
    int iterN = 0;
    for(int prev = 0; prev < prevIterSize; ++prev) {
      
      double prevWeight = as<List>(as<List>(as<List>(lwIter[nLwIter-1])[prev]))[1];      
      NumericVector prevTline = as<List>(as<List>(as<List>(lwIter[nLwIter-1])[prev]))[0];
      
      for(int cur = 0; cur < 2; ++cur) {
        
        double curWeight = as<List>(as<List>(as<List>(allTlines[nLwIter])[cur]))[1];
        NumericVector curTline = as<List>(as<List>(as<List>(allTlines[nLwIter])[cur]))[0];
        
        double iterWeight = prevWeight + curWeight;
        NumericVector iterTline = prevTline * curTline;
        List lwIterTemp2(2);
        lwIterTemp2[0] = iterTline;
        lwIterTemp2[1] = iterWeight;
        lwIterTemp[iterN] = lwIterTemp2;
        
        // start rp
        //rpTemp.insert(std::end(rpTemp), curWeight);
        std::vector<double> rpTemp = as<List>(as<List>(rpPrev))[prev];
        rpTemp.insert(std::end(rpTemp), curWeight);
        rpCurrent[iterN] = rpTemp;
        // end rp
        
        iterN = iterN + 1;
        
      }
    }
    lwIter[nLwIter] = lwIterTemp;
    lwFinalIter[0] = lwIterTemp; // again i need to change this to save memory
    
    // start rp
    rpList[nLwIter] = rpCurrent;
    // end rp
  }
  List out(2);
  out[0] = lwFinalIter[0];
  out[1] = rpList; 
  return out;
}

/*
 I can not calculate the sum of weights at each step 
 and instead just sum the RP. as n grows its hard to immagine
 how this would not be a little faster
 
*/