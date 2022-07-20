/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/
#ifndef ROO_NA60_SHAPE
#define ROO_NA60_SHAPE
  
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
  
class RooRealVar;
  
class NA60Shape : public RooAbsPdf {
 public:
  NA60Shape() {} ;
  NA60Shape(const char *name, const char *title, RooAbsReal& _m,
                RooAbsReal& _m0, RooAbsReal& _sigma,
                RooAbsReal& _alpha, RooAbsReal& _p1,
                RooAbsReal& _p2, RooAbsReal& _p3,
                RooAbsReal& _alpha2, RooAbsReal& _p12,
                RooAbsReal& _p22, RooAbsReal& _p32);
  
  NA60Shape(const NA60Shape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new NA60Shape(*this,newname); }
  
  inline virtual ~NA60Shape() { }
  
  virtual Int_t getAnalyticalIntegral( RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0 ) const;
  virtual Double_t analyticalIntegral( Int_t code, const char* rangeName=0 ) const;
  
  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;
  
 protected:
  
  Double_t ApproxErf(Double_t arg) const ;
  
  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alpha;
  RooRealProxy p1;
  RooRealProxy p2;
  RooRealProxy p3;
  RooRealProxy alpha2;
  RooRealProxy p12;
  RooRealProxy p22;
  RooRealProxy p32;
  
  Double_t evaluate() const;
  Double_t evaluateLoc(double x) const;
  
 private:
  
  ClassDef(NA60Shape,1) // Extended Crystal Ball lineshape PDF
};
  
#endif
