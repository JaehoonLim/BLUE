//-----------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
// 
// Copyright (C) 2012-2014, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.1.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//-----------------------------------------------------------------------------
#ifndef Blue_H
#define Blue_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDatime.h"
 
class Blue: public TObject {

public:

  //----------------------------------------------------------------------------
  // The Interface of the Class
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Constructor
  //----------------------------------------------------------------------------
  Blue(const Int_t NumEst, const Int_t NumUnc);
  Blue(const Int_t NumEst, const Int_t NumUnc, 
       const Int_t NumObs, const Int_t* IntObs);

  Blue(const Int_t NumEst, const Int_t NumUnc, const Int_t* IntFac);
  Blue(const Int_t NumEst, const Int_t NumUnc, const Int_t NumObs, 
       const Int_t* IntObs, const Int_t* IntFac);

  //----------------------------------------------------------------------------
  // Destructor
  //----------------------------------------------------------------------------
  ~Blue();
  
  //----------------------------------------------------------------------------
  // Solver
  //----------------------------------------------------------------------------
  void FixInp();
  void ReleaseInp();
  void ResetInp();
  void Solve();
  void SolveRelUnc(const Double_t Dx);
  void SolveAccImp(const Double_t Dx);
  void SolveAccImp(const Int_t ImpFla, const Double_t Dx);

  void SolveScaRho();
  void SolveScaRho(const Int_t RhoFla);
  void SolveScaRho(const Int_t RhoFla, const Double_t *MinRho, const Double_t *MaxRho);

  void SolveInfWei();
  void SolvePosWei();
  void SolveMaxVar(const Int_t IFuRho);
  
  //----------------------------------------------------------------------------
  // Filler
  //----------------------------------------------------------------------------
  void FillEst(const Int_t i, const Double_t *x);
  void FillCor(const Int_t k, const Double_t *x);
  void FillCor(const Int_t k, const Double_t rho);
  void FillNamEst(const TString* NamEst);
  void FillNamUnc(const TString* NamUnc);
  void FillNamObs(const TString* NamObs);
  
  //----------------------------------------------------------------------------
  // Geters
  //----------------------------------------------------------------------------
  // After FixInp()
  Int_t GetActEst() const;
  Int_t GetActEst(const Int_t n) const;
  Int_t GetActUnc() const;
  Int_t GetActObs() const;

  TString GetNamEst(const Int_t i) const;
  TString GetNamUnc(const Int_t k) const;
  TString GetNamObs(const Int_t n) const;

  Int_t GetIndEst(Int_t *IndEst) const;
  Int_t GetIndUnc(Int_t *IndUnc) const;
  Int_t GetIndObs(Int_t *IndObs) const;

  Int_t GetPreEst(const Int_t n) const;

  Int_t GetCov(TMatrixD *UseCov) const;
  Int_t GetCov(Double_t *RetCov) const;
  Int_t GetCovInvert(TMatrixD *UseCovI) const;
  Int_t GetCovInvert(Double_t *RetCovI) const;
  Int_t GetRho(TMatrixD *UseRho) const;
  Int_t GetRho(Double_t *RetRho) const;
  Int_t GetEst(TMatrixD *UseEst) const;
  Int_t GetEst(Double_t *RetEst) const;
  Int_t GetEstVal(TMatrixD *UseEstVal) const;
  Int_t GetEstVal(Double_t *RetEstVal) const;
  Int_t GetEstUnc(TMatrixD *UseEstUnc) const;
  Int_t GetEstUnc(Double_t *RetEstUnc) const;
  Int_t GetParams(const Int_t Ifl, TMatrixD *UseParams) const;
  Int_t GetParams(const Int_t Ifl, Double_t *RetParams) const;

  // After Solve()
  Double_t GetChiq() const;
  Int_t    GetNdof() const;
  Double_t GetProb() const;
  Double_t GetPull(const Int_t i) const;

  Int_t GetCovRes(TMatrixD *UseCovRes) const;
  Int_t GetCovRes(Double_t *RetCovRes) const;
  Int_t GetRhoRes(TMatrixD *UseRhoRes) const;
  Int_t GetRhoRes(Double_t *RetRhoRes) const;
  Int_t GetWeight(TMatrixD *UseWeight) const;
  Int_t GetWeight(Double_t *RetWeight) const;
  Int_t GetResult(TMatrixD *UseResult) const;
  Int_t GetResult(Double_t *RetResult) const;
  Int_t GetUncert(TMatrixD *UseUncert) const;
  Int_t GetUncert(Double_t *RetUncert) const;

  // For InspectLike
  Int_t GetInspectLike(TMatrixD *UseInsLik) const;
  Int_t GetInspectLike(Double_t *RetInsLik) const;

  // After SolveXXX()
  Int_t GetAccImpLasEst(const Int_t n) const;
  Int_t GetAccImpIndEst(const Int_t n, Int_t *IndEst) const;

  Int_t GetNumScaFac() const;
  Int_t GetNumScaRho() const;

  Int_t GetScaVal(const Int_t n, TMatrixD *UseScaVal) const;
  Int_t GetScaVal(const Int_t n, Double_t *RetScaVal) const;
  Int_t GetScaUnc(const Int_t n, TMatrixD *UseScaUnc) const;
  Int_t GetScaUnc(const Int_t n, Double_t *RetScaUnc) const;

  //----------------------------------------------------------------------------
  // Setters
  //----------------------------------------------------------------------------
  void SetActiveEst(const Int_t i);
  void SetActiveUnc(const Int_t k);

  void SetInActiveEst(const Int_t i);
  void SetInActiveUnc(const Int_t k);

  void SetRhoValUnc(const Double_t RhoVal);
  void SetRhoValUnc(const Int_t k, const Double_t RhoVal);
  void SetRhoValUnc(const Int_t k, const Int_t l, const Double_t RhoVal);
  void SetNotRhoValUnc();
  void SetNotRhoValUnc(const Int_t k);

  void SetRhoFacUnc(const Double_t RhoFac);
  void SetRhoFacUnc(const Int_t k, const Double_t RhoFac);
  void SetRhoFacUnc(const Int_t k, const Int_t l, const Double_t RhoFac);
  void SetNotRhoFacUnc();
  void SetNotRhoFacUnc(const Int_t k);

  void SetRhoRedUnc();
  void SetRhoRedUnc(const Int_t k);
  void SetNotRhoRedUnc();
  void SetNotRhoRedUnc(const Int_t k);

  void SetRelUnc();
  void SetRelUnc(const Int_t k);
  void SetRelUnc(const Int_t i, const Int_t k, const Double_t *ActCof);
  void SetNotRelUnc();
  void SetNotRelUnc(const Int_t k);

  //----------------------------------------------------------------------------
  // Control Printout
  //----------------------------------------------------------------------------
  void SetPrintLevel(const Int_t p);
  void SetQuiet();
  void SetNotQuiet();

  //----------------------------------------------------------------------------
  // Print Out
  //----------------------------------------------------------------------------
  void PrintMatrix(const TMatrixD *TryMat) const;
  void PrintMatrix(const TMatrixD *TryMat, const TString ForVal) const;
  void PrintMatrix(const TMatrixD *TryMat, const Int_t NumRow, const Int_t NumCol) const;
  void PrintMatrix(const TMatrixD *TryMat, const Int_t NumRow, const Int_t NumCol, 
		   const TString ForVal) const;

  void PrintDouble(const Double_t *TryDou, const Int_t NumRow, const Int_t NumCol) const;
  void PrintDouble(const Double_t *TryDou, const Int_t NumRow, const Int_t NumCol,
		   const TString ForVal) const;

  void PrintListEst() const;
  void PrintListUnc() const;

  void PrintNamEst() const;
  void PrintNamUnc() const;

  void PrintEst(const Int_t i) const;
  void PrintEst() const;

  void PrintCofRelUnc(const Int_t k) const;
  void PrintCofRelUnc() const;

  void PrintCor(const Int_t k) const;
  void PrintCor() const;
  void PrintCov(const Int_t k) const;
  void PrintCov() const;
  void PrintCovInvert() const;
  void PrintRho() const;
  void PrintCompatEst(const TString FilNam) const;
  void PrintCompatEst() const;

  void PrintParams(const Int_t Ifl) const;
  void PrintParams() const;

  void PrintPull(const Int_t i) const;
  void PrintPull() const;

  void PrintListObs() const;
  void PrintNamObs() const;

  void PrintCovRes() const;
  void PrintRhoRes() const;
  void PrintWeight() const;
  void PrintResult() const; 
  void PrintCompatObs() const;

  void PrintChiPro() const;

  void PrintInspectLike() const;

  // After SolveXXX()
  void PrintAccImp() const;

  void PrintScaRho(const TString FilNam) const;
  void PrintScaRho() const;

  void PrintInfWei() const;
  void PrintMaxVar() const;

  void PrintStatus() const;

  //----------------------------------------------------------------------------
  // Utilities
  //----------------------------------------------------------------------------
  Double_t GetPara(const Int_t ifl, const Double_t rho, const Double_t zva) const;
  Double_t FunPara(const Double_t* x, const Double_t* par) const;

  void InspectPair(const Int_t i, const Int_t j) const;
  void InspectPair(const Int_t i, const Int_t j, const TString FilNam) const;
  void InspectPair(const Int_t i, const Int_t j, const TString FilNam, const Int_t IndFig) const;

  void DrawSens(const Double_t xv1, const Double_t xv2, const Double_t sv1, 
		const Double_t sv2, const Double_t rho, const TString FilNam) const;
  void DrawSens(const Double_t xv1, const Double_t xv2, const Double_t sv1, 
		const Double_t sv2, const Double_t rho, const TString FilNam, const Int_t IndFig) const;

  void InspectLike(const Int_t n);
  void InspectLike(const Int_t n, TString FilNam);

  Int_t InspectResult() const;

  void LatexResult(const TString FilNam) const;
  void LatexResult(const TString FilNam, const TString ForVal, const TString ForUnc,
		   const TString ForWei, const TString ForRho, const TString ForPul) const;
  void DisplayResult(const Int_t n, const TString FilNam) const;
  void DisplayResult(const Int_t n, const TString FilNam, 
		     const TString ForVal, const TString ForUnc) const;

  // After SolveXXX()
  void DisplayAccImp(const Int_t n, const TString FilNam) const;
  void DisplayAccImp(const Int_t n, const TString FilNam, 
		     const TString ForVal, const TString ForUnc) const;

private:  

  //TString *_name;

  // The version
  TString Versio;

  // The date
  TDatime* Date;
  TString  Today;
  TString* Months;

  // Actual values of estimates, uncertainties and observables
  Int_t InpEst, InpUnc, InpObs;

  // Actual names of estimates, uncertainties and observables
  TString* EstNam;
  TString* UncNam;
  TString* ObsNam;

  // Arrays of informations 0/1 = yes/no for filled or active
  Int_t* EstAct;
  Int_t* UncAct;
  Int_t* ObsAct;

  Int_t* EstFil;
  Int_t* UncFil;

  // Variables for RhoScaleFactors
  Int_t     FlaFac;
  Int_t     InrFac;
  Int_t     InpFac;
  TMatrixD* MatFac;
  TMatrixD* ActFac;
  TMatrixD* MinFac;
  TMatrixD* MaxFac;

  TMatrixD** ValSca;
  TMatrixD** SigSca;
  TMatrixD** VtoSca;
  TMatrixD** StoSca;
  Int_t      FaiFac;

  // Arrays of indices of changed, scaled or reduced correlations
  Int_t* UncCha;
  Int_t* UncFac;
  Int_t* UncRed;

  // Array of indices of relative uncertainties 
  Int_t* UncRel;

  // Filling flags for names of Est, Unc and Obs
  Int_t IsFillEstNam;
  Int_t IsFillUncNam;
  Int_t IsFillObsNam;

  // Relation of estimates to observables
  // EstObs gives the observable an estimate determines
  Int_t* EstObs;

  // LisXXX is the list of active xxx = Est, Unc, Obs
  Int_t* LisEst;
  Int_t* LisUnc;
  Int_t* LisObs;

  // Control flag for all input filled and fixed
  Int_t InpFil, InpFix;

  // Control flag for parameters caculated
  Int_t IsCalcParams;

  // Control flags for combination performed
  Int_t IsSolve;
  Int_t IsSolveRelUnc;
  Int_t IsSolveAccImp;
  Int_t IsSolveScaRho;
  Int_t IsSolveInfWei;
  Int_t IsSolvePosWei;
  Int_t IsSolveMaxVar;

  // Control flag and matrix of coefficients for iterative Blue
  Int_t IsRelUncMode;
  static const Int_t MaxCof = 3;
  TMatrixD* Cof;

  // The list of importance of estimates for SolveAccImp == abs(weight) 
  Double_t  PreAcc;
  Int_t     IntAcc;
  Int_t*    LisImp;
  TMatrixD* LasImp;
  TMatrixD* IndImp;
  TMatrixD* ValImp;
  TMatrixD* UncImp;
  TMatrixD* StaImp;
  TMatrixD* SysImp;

  // The arrays for DisplayAny
  Int_t*    Indx;
  Int_t*    Colo;
  TString*  Name;
  Double_t* Valu;
  Double_t* Stat;
  Double_t* Syst;
  Double_t* Full;

  // The list of:
  // 1) initial/final variances/correlations for SolveMaxVar
  // 2) the rho values 
  // 3) the steering for the fijk futch 
  TVectorD* VarMax;
  TVectorD* VarMin;
  TMatrixD* RhoMax;
  TMatrixD* RhoMin;
  TMatrixD* RhoOoz;

  Double_t  RhoFco;
  TVectorD* RhoFpk;
  TMatrixD* RhoFij;
  Int_t     IFuCor;
  Int_t     IFuFla;
  Int_t*    IFuFai;

  // Control flag for print out level
  Int_t IPrint;

  // Control flag for quiet mode
  Int_t IQuiet;

  // Save harbor for all input
  Int_t     InpEstOrig, InpUncOrig, InpObsOrig;
  TVectorD* XvaOrig;
  TVectorD* SigOrig;
  TMatrixD* UncOrig;
  TMatrixD* CorOrig;
  TMatrixD* UmaOrig;
  TMatrixD* UtrOrig;
  
  TString* EstNamOrig;
  TString* UncNamOrig;
  TString* ObsNamOrig;

  // Vectors to do the job
  TVectorD* Xva;
  TVectorD* Sig;

  // Matrices to do the job
  TMatrixD* Unc;
  TMatrixD* Cor;
  TMatrixD* Cov;
  TMatrixD* CovI;
  TMatrixD* Rho;
  TMatrixD* Lam;
  TVectorD* Pul;
  TMatrixD* Uma;
  TMatrixD* Utr;

  // Matrices for the results
  TVectorD* XvaRes;
  TMatrixD* CorRes;
  TMatrixD* CovRes;
  TMatrixD* RhoRes;

  // Chiq information of the results
  Double_t ChiQua;
  Double_t ChiPro;
  Int_t    NumDof;

  // Vectors for the informatioon weights
  TVectorD* VarInd;
  TVectorD* IntWei;
  TVectorD* MarWei;
  TVectorD* BluWei;

  // Matrices for the parameters
  TMatrixD* SRat;
  TMatrixD* Beta;
  TMatrixD* Sigx;
  TMatrixD* DBdr;
  TMatrixD* DSdr;
  TMatrixD* DBdz;
  TMatrixD* DSdz;

  // Matrix for control of the relative uncertainties 
  TMatrixD*  IndRel;

  // Flag to disable not yet released parts
  Int_t    NotRel;

  // Matrices for the InspectLike results
  TMatrixD* ChiRes;
  static const Int_t LikDim = 7;
  TMatrixD* LikRes;
  Int_t     LikFla;

  // Control flags for InspectLike
  Int_t IsInspectLike;

private:
  //----------------------------------------------------------------------------
  // Private Functions do do the job
  // -- Warning -- NOT PRIVATE in ACLiC --- Refrain from using these functions
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Constructor
  //----------------------------------------------------------------------------
  void Construct(const Int_t NumEst, const Int_t NumUnc, const Int_t NumObs, 
		 const Int_t *IntObs, const Int_t *IntFac);
  void SetupRoot();

  //----------------------------------------------------------------------------
  // Solver
  //----------------------------------------------------------------------------
  void SolveScaRho(const Int_t RhoFla, const TMatrixD* MinFac, const TMatrixD* MaxFac);

  void SetIsSolved(const Int_t l);
  void SetIsSolvedRelUnc(const Int_t l);
  void SetIsSolvedAccImp(const Int_t l);
  void SetIsSolvedScaRho(const Int_t l);
  void SetIsSolvedInfWei(const Int_t l);
  void SetIsSolvedPosWei(const Int_t l);
  void SetIsSolvedMaxVar(const Int_t l);
  void SetIsInspectLike(const Int_t l);

  //----------------------------------------------------------------------------
  // Fillers
  //----------------------------------------------------------------------------
  void FillCov();
  void FillCovInvert();
  void FillSig();
  void FillRho();

  //----------------------------------------------------------------------------
  // Getters
  //----------------------------------------------------------------------------
  Int_t IsAllowedEst(const Int_t i) const;
  Int_t IsActiveEst(const Int_t i) const;
  Int_t IsFilledEst(const Int_t i) const;
  Int_t IsFilledEst() const;
  Int_t IsWhichEst(const Int_t ia) const;
  Int_t IsIndexEst(const Int_t i) const;
  Int_t EstWhichObs(const Int_t i) const;

  Int_t IsAllowedUnc(const Int_t k) const;
  Int_t IsActiveUnc(const Int_t k) const;
  Int_t IsFilledUnc(const Int_t k) const;
  Int_t IsFilledUnc() const;

  Int_t IsAllowedFac(const Int_t l) const;
  Int_t IsRhoValUnc(const Int_t k) const;
  Int_t IsRhoFacUnc(const Int_t k) const;
  Int_t IsRhoRedUnc(const Int_t k) const;
  Int_t IsRelValUnc(const Int_t k) const;

  Int_t IsWhichUnc(const Int_t ka) const;
  Int_t IsIndexUnc(const Int_t k) const;

  Int_t IsAllowedObs(const Int_t n) const;
  Int_t IsActiveObs(const Int_t n) const;
  Int_t IsWhichObs(const Int_t na) const;
  Int_t IsIndexObs(const Int_t n) const;

  Int_t IsFilledNamEst() const;
  Int_t IsFilledNamUnc() const;
  Int_t IsFilledNamObs() const;

  Int_t IsFilledInp() const;
  Int_t IsFixedInp() const;

  Int_t IsCalcedParams() const;

  Int_t IsSolved() const;
  Int_t IsSolvedRelUnc() const;
  Int_t IsSolvedAccImp() const;
  Int_t IsSolvedScaRho() const;
  Int_t IsSolvedInfWei() const;
  Int_t IsSolvedPosWei() const;
  Int_t IsSolvedMaxVar() const;
  Int_t IsSolvedAnyWay() const;

  Int_t IsInspectedLike() const;

  Int_t IsRelValUnc() const;

  Int_t IsPrintLevel() const;
  Int_t IsQuiet() const;

  Int_t IsFutchCor() const;

  Int_t GetMatEigen(const TMatrixD *TryMat) const;

  //----------------------------------------------------------------------------
  // Setters
  //----------------------------------------------------------------------------
  void SetActiveEst(const Int_t i, const Int_t l);
  void SetFilledEst(const Int_t i);
  void SetFilledEst();

  void SetActiveUnc(const Int_t i, const Int_t l);
  void SetFilledUnc(const Int_t k);
  void SetFilledUnc();

  void SetRhoValUnc(const Int_t k, const Int_t l, const Double_t RhoVal, const Int_t m);
  void SetRhoFacUnc(const Int_t k, const Int_t l, const Double_t RhoFac, const Int_t m);
  void SetRhoRedUnc(const Int_t i, const Int_t l);

  void SetRelValUnc(const Int_t l);

  void SetFilledInp();
  void SetFixedInp(const Int_t l);

  void SetQuiet(const Int_t l);

  void SetCalcedParams(const Int_t l);

  void SetFutchCor(const Int_t i, const Int_t j, const Int_t l, const Double_t RhoFut);

  //----------------------------------------------------------------------------
  // PrintOut
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Calculate
  //----------------------------------------------------------------------------
  Double_t CalcRelUnc(const Int_t i, const Int_t k, const Double_t ValCom) const;
  Double_t Likelihood(const Double_t* x, const Double_t* par) const;
  void     CalcChiRes(const Int_t n, const Int_t IsLike, const Double_t xt) const;

  void CalcParams();

  void     ResetScaRho(const Int_t IFlag);
  Double_t CalcActFac(const Int_t i, const Int_t j, const Int_t k);

  Double_t CalcFutchCor(const Int_t i, const Int_t j);

  Int_t IsWhichMatrix(const TMatrixD *TryMat) const;

  void MatrixtoDouble(const TMatrixD *InpMat, Double_t *OutDou) const;

  void WriteInput(const Int_t IWhat, const TString FilNam, 
		  const TString ForVal, const TString ForRho) const;

  //----------------------------------------------------------------------------
  // Display
  //----------------------------------------------------------------------------
  void PlotRes(const Int_t N, const TString FilNam, 
	       const TString ForVal, const TString ForUnc) const;

  // The Class definition
  ClassDef(Blue,0);

};

#endif // Blue_H
