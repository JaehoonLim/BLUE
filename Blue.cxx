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
#include "TObject.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TMatrixDEigen.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TDatime.h"
#include "Blue.h"

// ----> First the implementation of the public member functions

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

ClassImp(Blue);

Blue::Blue(const Int_t NumEst, const Int_t NumUnc) : TObject() {

  // The NumObs=1 default
  Int_t NumObs = 1;
  Int_t IWhichObs[NumEst];
  for(Int_t i = 0; i<NumEst; i++){IWhichObs[i] = 0;}; 

  // The NumFac=1 default
  Int_t NumFac = 1;
  Int_t IWhichFac[NumEst*NumEst];
  for(Int_t i = 0; i<NumEst*NumEst; i++){IWhichFac[i] = 0;}; 

  Construct(NumEst, NumUnc, NumObs, &IWhichObs[0], &IWhichFac[0]);
};

//------------------------------------------------------------------------------

Blue::Blue(const Int_t NumEst, const Int_t NumUnc, const Int_t NumObs, 
	   const Int_t* IntObs){

  // The NumFac=1 default
  Int_t NumFac = 1;
  Int_t IWhichFac[NumEst*NumEst];
  for(Int_t i = 0; i<NumEst*NumEst; i++){IWhichFac[i] = 0;}; 

  Construct(NumEst, NumUnc, NumObs, &IntObs[0], &IWhichFac[0]);
};

//------------------------------------------------------------------------------

Blue::Blue(const Int_t NumEst, const Int_t NumUnc, const Int_t* IntFac) : TObject() {

  // The NumObs=1 default
  Int_t NumObs = 1;
  Int_t IWhichObs[NumEst];
  for(Int_t i = 0; i<NumEst; i++){IWhichObs[i] = 0;}; 

  Construct(NumEst, NumUnc, NumObs, &IWhichObs[0], &IntFac[0]);
};

//------------------------------------------------------------------------------

Blue::Blue(const Int_t NumEst, const Int_t NumUnc, const Int_t NumObs, 
	   const Int_t* IntObs, const Int_t* IntFac){

  Construct(NumEst, NumUnc, NumObs, &IntObs[0], &IntFac[0]);
};

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------

Blue::~Blue(){
  //
  // Clean up all local objects
  //

  // The Date
  Date->~TDatime(); Date = NULL;
  Months->~TString(); Months = NULL;
  
  // The names of estimates, uncertainties and observables
  EstNam->~TString(); EstNam = NULL;
  UncNam->~TString(); UncNam = NULL;
  ObsNam->~TString(); ObsNam = NULL;
  EstNamOrig->~TString(); EstNamOrig = NULL;
  UncNamOrig->~TString(); UncNamOrig = NULL;
  ObsNamOrig->~TString(); ObsNamOrig = NULL;
  
  // The array of which estimates determines which observable
  delete EstObs; EstObs = NULL;

  // The matrix of scale factors
  MatFac->Delete(); MatFac = NULL;

  // The actual factors and ranges for SolveScaRho
  ActFac->Delete(); ActFac = NULL;
  MinFac->Delete(); MinFac = NULL;
  MaxFac->Delete(); MaxFac = NULL;

  // The number of scan values and the result matrices
  for(Int_t n = 0; n<InpObsOrig; n++){
    ValSca[n]->Delete(); ValSca[n] = NULL;
    SigSca[n]->Delete(); SigSca[n] = NULL;
    VtoSca[n]->Delete(); VtoSca[n] = NULL;
    StoSca[n]->Delete(); StoSca[n] = NULL;
  }

  // The lists of active estimates and uncertainties
  delete LisEst; LisEst = NULL;
  delete LisUnc; LisUnc = NULL;

  // The U-Matrices for the combination
  Uma->Delete(); Uma = NULL;
  Utr->Delete(); Utr = NULL;

  // Arrays for controlling the estimates
  delete EstAct; EstAct = NULL;
  delete EstFil; EstFil = NULL;
  delete UncAct; UncAct = NULL;
  delete UncFil; UncFil = NULL;
  delete UncCha; UncCha = NULL;
  delete UncFac; UncFac = NULL;
  delete UncRed; UncRed = NULL;
  delete UncRel; UncRel = NULL;
  
  // List of active observables and list of observables
  delete LisObs; LisObs = NULL;
  delete ObsAct; ObsAct = NULL;

  // Array of coefficients forBLUE with relative uncertainties
  Cof->Delete(); Cof = NULL;

  // The list of importance and the steering flag for SolveAccImp
  delete LisImp; LisImp = NULL;
  delete LasImp; LasImp = NULL;
  delete IndImp; IndImp = NULL;
  delete ValImp; ValImp = NULL;
  delete UncImp; UncImp = NULL;
  delete StaImp; StaImp = NULL;
  delete SysImp; SysImp = NULL;

  // The arrays for PlotRes
  delete Indx; Indx = NULL;
  delete Colo; Colo = NULL;
  Name->~TString(); Name = NULL;
  delete Valu; Valu = NULL;
  delete Stat; Stat = NULL;
  delete Syst; Syst = NULL;
  delete Full; Full = NULL;

  // The list of initial/final variances/correlations for SolveMaxVar
  VarMax->Delete(); VarMax = NULL;
  VarMin->Delete(); VarMin = NULL;
  RhoMax->Delete(); RhoMax = NULL;
  RhoMin->Delete(); RhoMin = NULL;
  RhoOoz->Delete(); RhoOoz = NULL;
  RhoFpk->Delete(); RhoFpk = NULL;
  RhoFij->Delete(); RhoFij = NULL;
  delete IFuFai; IFuFai = NULL;

  // Save harbor for all inputs
  XvaOrig->Delete(); XvaOrig = NULL;
  SigOrig->Delete(); SigOrig = NULL;
  UncOrig->Delete(); UncOrig = NULL;
  CorOrig->Delete(); CorOrig = NULL;
  UmaOrig->Delete(); UmaOrig = NULL;
  UtrOrig->Delete(); UtrOrig = NULL;

   // Vectors to do the job
  Xva->Delete(); Xva = NULL;
  Sig->Delete(); Sig = NULL;

  // Matrices to do the job
  Unc->Delete(); Unc = NULL;
  Cor->Delete(); Cor = NULL;
  Cov->Delete(); Cov = NULL;
  CovI->Delete(); CovI = NULL;
  Rho->Delete(); Rho = NULL;
  Lam->Delete(); Lam = NULL;
  Pul->Delete(); Pul = NULL;

  // Matrices for the results
  XvaRes->Delete(); XvaRes = NULL;
  CorRes->Delete(); CorRes = NULL;
  CovRes->Delete(); CovRes = NULL;
  RhoRes->Delete(); RhoRes = NULL;

  // Vectors for the information weights
  VarInd->Delete(); VarInd = NULL;
  IntWei->Delete(); IntWei = NULL;
  MarWei->Delete(); MarWei = NULL;
  BluWei->Delete(); BluWei = NULL;

   // Matrices for the parameters
  SRat->Delete(); SRat = NULL;
  Beta->Delete(); Beta = NULL;
  Sigx->Delete(); Sigx = NULL;
  DBdr->Delete(); DBdr = NULL;
  DSdr->Delete(); DSdr = NULL;
  DBdz->Delete(); DBdz = NULL;
  DSdz->Delete(); DSdz = NULL;

  // Matrix for control of the relative uncertainties 
  IndRel->Delete(); IndRel = NULL;
  
  // Matrices for the InspectLike results
  ChiRes->Delete(); ChiRes = NULL;
  LikRes->Delete(); LikRes = NULL;

  // if(_name){delete _name; _name = 0;};

  // Say that we are done
  if(IsQuiet() == 0){
    printf("... Blue~Blue: Delete object for InpEst=%3i, InpUnc=%3i, InpObs=%2i\n",
	   InpEstOrig, InpUncOrig, InpObsOrig);
  }

  return;
};

//------------------------------------------------------------------------------
// Solver
//------------------------------------------------------------------------------

void Blue::FixInp(){
  if(IsFixedInp() == 1){
    printf("... Blue->FixInp(): IGNORED Nothing changed since last combination.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  } 
  if(IsFilledInp() == 0){
    printf("... Blue->FixInp(): IGNORED input is not yet filled \n");
    return;
  }

  // Check whether relative uncertainties are properly filled
  Int_t IFail = 0;
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      if(IsRelValUnc(k) == -1){
	printf("... Blue->FixInp(): Relative Unc %2i has the follwing estimates not set", k);
	for(Int_t i = 0; i<InpEstOrig; i++){
	  if(IsActiveEst(i) == 1){
	    if(IndRel->operator()(i,k) == 0.){printf(" %2i",i);};
	  }
	}
	printf("\n");
	IFail = IFail + 1;
      }
    }
  }
  if(IFail >= 1){
    printf("... Blue->FixInp(): ---------------------------------------------------------- \n");
    printf("... Blue->FixInp(): !!! %2i relative uncertainties are not properly filled !!! \n", IFail);
    printf("... Blue->FixInp(): !!! NO COMBINATION POSSIBLE I will not fix the input   !!! \n");
    printf("... Blue->FixInp(): ---------------------------------------------------------- \n");
    return;
  }

  // Fill number and list of active estimates
  InpEst = 0;
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i) == 1){
      LisEst[InpEst] = i;
      InpEst = InpEst + 1;
    }
  }
  
  // Get number of active uncertainties
  InpUnc = 0;
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      LisUnc[InpUnc] = k;
      InpUnc = InpUnc + 1;
    }
  }

  // Fill number and list of active observables
  //printf("... Blue->FixInp(): UmaOrig \n"); UmaOrig->Print();
  Int_t mm = 0;
  InpObs = 0;
  for(Int_t n = 0; n<InpObsOrig; n++){
    ObsAct[n] = 0;
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i) == 1){
	if(UmaOrig->operator()(i,n) == 1){
	  ObsAct[n] = 1;
	}
      }
    }
    if(IsActiveObs(n) == 1){
      InpObs = InpObs + 1;
      LisObs[mm] = n;
      mm = mm + 1;
    }
  }

  // Fill U and UT matricees
  TMatrixD *F  = new TMatrixD(InpEst,InpObs);
  Int_t j = 0;
  Int_t m = 0;
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i) == 1){
      for(Int_t n = 0; n<InpObsOrig; n++){
	if(IsActiveObs(n) == 1){
	  F->operator()(j,m) = UmaOrig->operator()(i,n);
	  m = m + 1;
	}
      }
      j = j + 1;
      m = 0;
    }
  }
  Uma->Delete(); Uma = new TMatrixD(InpEst,InpObs);
  Utr->Delete(); Utr = new TMatrixD(InpObs,InpEst);
  Uma->operator=(*F);
  Utr->Transpose(*Uma);
  //printf("... Blue->FixInp(): Uma \n"); Uma->Print();
  //printf("... Blue->FixInp(): Utr \n"); Utr->Print();
  F->Delete(); F = NULL;

  // Fill vector of estimates
  Xva->Delete(); Xva = new TVectorD(InpEst);
  j = 0;
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i) == 1){
      Xva->operator()(j) = XvaOrig->operator()(i);
      j = j + 1;
    }
  }
  //printf("... Blue->FixInp(): Xva \n"); Xva->Print();

  // Fill matrix of uncertainties
  Unc->Delete(); Unc  = new TMatrixD(InpEst*InpUnc,InpEst*InpUnc);
  Int_t IRow, IRowOrig;
  j = 0;
  Int_t l = 0;
  Double_t ValCom = 0, UncCom = 0;
  for(Int_t k = 0; k<InpUncOrig; k++){ 
    if(IsActiveUnc(k) == 1){
      j = 0;
      for(Int_t i = 0; i<InpEstOrig; i++){
	if(IsActiveEst(i) == 1){
	  IRowOrig=i+k*InpEstOrig;
	  IRow=j+l*InpEst;
	  UncCom = UncOrig->operator()(IRowOrig,IRowOrig);
	  // Rescale relative uncertainties
	  if(IsRelValUnc() == 1 && IsRelValUnc(k) == 1){
	    ValCom = 0.;
	    for(Int_t n = 0; n<InpObs; n++){
	      //23.10.13 (i,n) -> (j,n);
	      ValCom = ValCom + Uma->operator()(j,n)*XvaRes->operator()(n);
	    }
	    UncCom = CalcRelUnc(i, k, ValCom);
	    if(UncCom < 0)IFail = IFail + 1;
	  }
	  Unc->operator()(IRow,IRow) = UncCom;
	  j = j + 1;
	}
      }
      l = l + 1;
    }
  }
  // printf("... Blue->FixInp(): Unc \n"); Unc->Print();
  if(IFail >= 1){
    printf("... Blue->FixInp(): ---------------------------------------------------------- \n");
    printf("... Blue->FixInp(): !!! %2i relative uncertainties are not properly filled !!! \n", IFail);
    printf("... Blue->FixInp(): !!! NO COMBINATION POSSIBLE I will not fix the input   !!! \n");
    printf("... Blue->FixInp(): ---------------------------------------------------------- \n");
    return;
  }
  
  // Fill matrix of correlations
  //printf("... Blue->FixInp(): CorOrig \n"); CorOrig->Print();
  TMatrixD *FN  = new TMatrixD(InpEst,InpEst);
  TMatrixD *FO  = new TMatrixD(InpEstOrig,InpEstOrig);
  TMatrixD *FS  = new TMatrixD(InpEstOrig,InpEstOrig);
  TMatrixD *CC  = new TMatrixD(InpEst*InpUnc,InpEst*InpUnc);
  
  Int_t ii = 0;
  Int_t jj = 0;
  Int_t kk = 0;
  Int_t IRowLow, IRowHig, IColLow, IColHig;
  for(Int_t k = 0; k<InpUncOrig; k++){ 
    if(IsActiveUnc(k) == 1){
      IRowLow = k*InpEstOrig; 
      IRowHig = IRowLow + InpEstOrig-1;
      IColLow = IRowLow;
      IColHig = IRowHig;
      CorOrig->GetSub(IRowLow, IRowHig, IColLow, IColHig, *FO, "S");
      //printf("... Blue->FixInp(): F0 \n"); FO->Print();

      // This is the place to do 1) changed, 2) scaled or 3) reduced
      // correlations by changing FO
      // 1) Do the change
      if(IsRhoValUnc(k) == 1){
	for(Int_t i1 = 0; i1<InpEstOrig; i1++){
	  for(Int_t j1 = i1+1; j1<InpEstOrig; j1++){
	    FO->operator()(i1,j1) = CalcActFac(i1,j1,k);
	    FO->operator()(j1,i1) = FO->operator()(i1,j1);
	  }
	}
	//printf("... Blue->FixInp() FO is changed \n"); FO->Print();
      }

      // 2) Do the scaling
      if(IsRhoFacUnc(k) == 1){
	for(Int_t i1 = 0; i1<InpEstOrig; i1++){
	  for(Int_t j1 = i1+1; j1<InpEstOrig; j1++){
	    FO->operator()(i1,j1) = FO->operator()(i1,j1)*CalcActFac(i1,j1,k);
	    FO->operator()(j1,i1) = FO->operator()(i1,j1);
	  }
	}
	//printf("... Blue->FixInp() FO is scaled \n"); FO->Print();
      }

      // 3) Do the reduced correlations
      //UncOrig->Print();
      if(IsRhoRedUnc(k) == 1){
	UncOrig->GetSub(IRowLow, IRowHig, IColLow, IColHig, *FS, "S");
	//printf("... Blue->FixInp(): FS \n"); FS->Print();
	//printf("... Blue->FixInp(): FO \n"); FO->Print();
	Double_t RedCor, Limit = 0.00001, FullCor = 0.99;
	for(Int_t i1 = 0; i1<InpEstOrig; i1++){
	  for(Int_t j1 = i1+1; j1<InpEstOrig; j1++){
	    if(FS->operator()(i1,i1) < FS->operator()(j1,j1)){
	      RedCor = FS->operator()(i1,i1) / FS->operator()(j1,j1);
	    }else{
	      if(FS->operator()(i1,i1) > Limit){
		RedCor = FS->operator()(j1,j1) / FS->operator()(i1,i1);
	      }else{
		RedCor = FS->operator()(j1,j1) / Limit;
	      }
	    }
	    //printf("... Blue->FixInp(): RC %2i, %2i,  %5.3f,  \n",i1,j1,RedCor); 
	    if(FO->operator()(i1,j1) > FullCor){
	      FO->operator()(i1,j1) = RedCor;
	      FO->operator()(j1,i1) = FO->operator()(i1,j1);
	    }
	  }
	}
	//printf("... Blue->FixInp() FO is reduced \n"); FO->Print();
      }

      // Futch that is needed for SolveMaxVar(2)
      if(IsFutchCor() == 1 && k > 0){
	//printf("... Blue->FixInp(): Futch correlations %2i \n",k);
	for(Int_t i1 = 0; i1<InpEstOrig; i1++){
	  for(Int_t j1 = i1+1; j1<InpEstOrig; j1++){
	    FO->operator()(i1,j1) = FO->operator()(i1,j1)*CalcFutchCor(i1,j1);
	    FO->operator()(j1,i1) = FO->operator()(i1,j1);
	  }
	}
	//printf("... Blue->FixInp() FO is futched \n"); FO->Print();
      }

      // Fill result
      for(Int_t iii = 0; iii<InpEstOrig; iii++){
	if(IsActiveEst(iii) == 1){	  
	  for(Int_t jjj = 0; jjj<InpEstOrig; jjj++){
	    if(IsActiveEst(jjj) == 1){
	      //printf("... Blue->FixInp() fill i=%2i, j=%2i, ii=%2i, jj=%2i \n",iii,jjj,ii,jj);
	      FN->operator()(ii,jj) = FO->operator()(iii,jjj);
	      jj = jj + 1;
	    }
	  }
	  ii = ii + 1;
	  jj = 0;
	}
      }
      IRowLow = kk*InpEst; 
      IColLow = IRowLow;
      CC->SetSub(IRowLow, IColLow, *FN);
      ii = 0;
      jj = 0;
      kk = kk + 1;
    }
  }
  Cor->Delete(); Cor = new TMatrixD(InpEst*InpUnc,InpEst*InpUnc);
  Cor->operator=(*CC);
  FN->Delete(); FN = NULL;
  FO->Delete(); FO = NULL;
  FS->Delete(); FS = NULL;
  CC->Delete(); CC = NULL;
  //printf("... Blue->FixInp(): Cor \n"); Cor->Print();

  // Calculate covariance
  FillCov();
  //printf("... Blue->FixInp(): Cov \n"); Cov->Print();

  // Calculate total uncertainties
  FillSig();
  //printf("... Blue->FixInp(): Sig \n"); Sig->Print();

  // Check for --zero-- total uncertainties of estimates, exit if this exists
  for(Int_t i = 0; i<InpEst; i++){     
    if(Sig->operator()(i) > 0){
    }else{
      printf("... Blue->FixInp(): -------------------------------------------------------- \n");
      printf("... Blue->FixInp(): !!! The total uncertainty for estimate %2i is zero   !!! \n", i);
      printf("... Blue->FixInp(): !!! NO COMBINATION POSSIBLE I will not fix the input !!! \n");
      printf("... Blue->FixInp(): -------------------------------------------------------- \n");
      return;
    }
  }

  // Calculate inverse covariance
  FillCovInvert();
  //printf("... Blue->FixInp(): CovI \n"); CovI->Print();

  // Calculate input correlations
  FillRho();
  //printf("... Blue->FixInp(): Rho \n"); Rho->Print();

  // Check whether at least one each is left
  if(InpEst == 0 || InpUnc == 0 || InpObs == 0){
    printf("... Blue->FixInp(): ----------------------------------------------------------- \n");
    printf("... Blue->FixInp(): !!! At least one estimate and uncertaintainty is needed !!! \n");
    printf("... Blue->FixInp(): !!! NO COMBINATION POSSIBLE I will not fix the input    !!! \n");
    printf("... Blue->FixInp(): ----------------------------------------------------------- \n");
    printf("... Blue->FixInp(): InpEst =%3i, InpUnc =%3i, InpObs =%2i\n",
	   InpEst, InpUnc, InpObs);
    return;
  }

  // Report the situation. How many are left and what type they are
  if(IsQuiet() == 0){
    printf("... Blue->FixInp(): Fixed to InpEst =%3i, InpUnc =%3i, InpObs =%2i\n",
	   InpEst, InpUnc, InpObs);
  }

  // Fill the parameter matrices
  CalcParams();

  // Refresh the Likelihood matrix
  LikRes->Delete(); LikRes = new TMatrixD(InpObs,LikDim);

  // Set the flag 
  SetFixedInp(1);
};

//------------------------------------------------------------------------------

void Blue::ReleaseInp(){
  // Reset the flags
  SetFixedInp(0);
  SetIsSolved(0);
  SetIsSolvedRelUnc(0);
  SetIsSolvedAccImp(0);
  SetIsSolvedScaRho(0);
  SetIsSolvedInfWei(0);
  SetIsSolvedPosWei(0);
  SetIsSolvedMaxVar(0);
  SetIsInspectLike(0);
  SetCalcedParams(0);
  SetFutchCor(0,0,0,1);
};

//------------------------------------------------------------------------------

void Blue::ResetInp(){
  ReleaseInp();

  // Reset the estimate and uncertainty selections
  for(Int_t i = 0; i<InpEstOrig; i++){SetActiveEst(i,1);};
  for(Int_t k = 0; k<InpUncOrig; k++){SetActiveUnc(k,1);};
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsRhoValUnc(k) == 1){SetNotRhoValUnc(k);};
    if(IsRhoFacUnc(k) == 1){SetNotRhoFacUnc(k);};
    if(IsRhoRedUnc(k) == 1){SetNotRhoRedUnc(k);};
    if(IsRelValUnc(k) != 0){SetNotRelUnc(k);};
  }
  ResetScaRho(2);
};

//------------------------------------------------------------------------------

void Blue::Solve(){
  if(IsFilledInp() == 0){
    printf("... Blue->Solve(): IGNORED Not all estimates/uncertainties filled IGNORED \n");
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->Solve(): IGNORED Input not yet fixed call FixInp \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->Solve(): IGNORED Two consecutive calls to Solve() are not supported.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  if(IsPrintLevel() >= 1){
    printf("... Blue->Solve(): The present situation is \n");
    PrintStatus();
  }

  // Reset weight matrix and pull vector
  Lam->Delete(); Lam  = new TMatrixD(InpEst,InpObs);
  Pul->Delete(); Pul  = new TVectorD(InpEst);

  // Reset matrices for the results
  XvaRes->Delete(); XvaRes = new TVectorD(InpObs);
  CorRes->Delete(); CorRes = new TMatrixD(InpObs*InpUnc,InpObs*InpUnc);
  CovRes->Delete(); CovRes = new TMatrixD(InpObs,InpObs);
  RhoRes->Delete(); RhoRes = new TMatrixD(InpObs,InpObs);

  // Get weight matrix
  //-1)      H: (Cov-1 * U)
  //-2) CovRes: (UT * Cov-1 * U)-1
  //-3)    Lam: (Cov-1 * U) * (UT * Cov-1 * U)-11 
  TMatrixD *H = new TMatrixD(InpEst,InpObs);
  H->Mult(*CovI, *Uma);
  CovRes->Mult(*Utr,*H);
  CovRes->Invert();
  Lam->Mult(*H, *CovRes);
  H->Delete(); H = NULL;
  //printf("... Blue->Solve(): Weight matrix\n"); Lam->Print();

  // Calculate covariance matrix of observables per uncertainty source
  TMatrixD *C = new TMatrixD(InpEst,InpEst);
  TMatrixD *E = new TMatrixD(InpEst,InpEst);
  TMatrixD *F = new TMatrixD(InpEst,InpEst);
  TMatrixD *I = new TMatrixD(InpEst,InpEst);
  
  TMatrixD *J = new TMatrixD(InpEst,InpObs);
  TMatrixD *K = new TMatrixD(InpObs,InpEst);
  TMatrixD *L = new TMatrixD(InpObs,InpObs);
  
  Int_t IRowLow, IRowHig, IColLow, IColHig;
  for(Int_t k = 0; k<InpUnc; k++){
    // Get covariance matrix per uncertainty source
    IRowLow = k*InpEst; 
    IRowHig = IRowLow + InpEst-1;
    IColLow = IRowLow;
    IColHig = IRowHig;
    Unc->GetSub(IRowLow, IRowHig, IColLow, IColHig, *E, "S");
    Cor->GetSub(IRowLow, IRowHig, IColLow, IColHig, *C, "S");
    F->Mult(*E, *C);
    I->Mult(*F, *E);
    //printf("... Blue->Solve(): Cov_i \n"); I->Print();
    // Get covariance Eq 18. Cov(n,m)_k = Lam(n,i)*Cov(i,j)_k*Lam(j,m)
    J->Mult(*I, *Lam);
    K->Transpose(*Lam);
    L->Mult(*K,*J);
    L->Sqrt();
    CorRes->SetSub(k*InpObs,k*InpObs,*L);
    //CorRes->Print();     
  }
  C->Delete(); C = NULL;
  E->Delete(); E = NULL;
  F->Delete(); F = NULL;
  I->Delete(); I = NULL;
  J->Delete(); J = NULL;
  K->Delete(); K = NULL;
  L->Delete(); L = NULL;

  // Calculate correlation matrix of results
  RhoRes->SetSub(0,0,*CovRes);
  for(Int_t n = 0; n<InpObs; n++){  
    for(Int_t m = 0; m<InpObs; m++){  
	if(n != m){
	  RhoRes->operator()(n,m) = RhoRes->operator()(n,m)/
	  sqrt(RhoRes->operator()(n,n))/
	  sqrt(RhoRes->operator()(m,m));
      }
    }
  }
  for(Int_t n = 0; n<InpObs; n++){  
    RhoRes->operator()(n,n) = 1.;
  }
  
  // Calculate combined values XvaRes(n) = Sum_i Lam(i,n) * Xva(i)
  Double_t val;
  for (Int_t n = 0; n < InpObs; n++){    
    for(Int_t i = 0; i<InpEst; i++){  
      val = val + Lam->operator()(i,n)*Xva->operator()(i);
    }
    XvaRes->operator()(n) = val;
    val = 0.;
  }
  //XvaRes->Print();

  // Calculate Chi-squared, Nndof and Probability Eq(14)
  // Calculate Pull (yi - y)/sqrt(sig_yi^2-sig_y^2)
  //-1) N: (Xva - U * XvaRes)
  //-2) O: N^T
  //-3) M: (CovI * N)
  //-3) P: chiq: (O * M)
  TMatrixD *M = new TMatrixD(InpEst,1);
  TMatrixD *N = new TMatrixD(InpEst,1);
  TMatrixD *O = new TMatrixD(1,InpEst);
  TMatrixD *P = new TMatrixD(1,1);

  // Fill N, calculate O, use val from above
  // Calculate Pull on the fly
  //Sig->Print(); CovRes->Print();
  Double_t unc;
  for(Int_t i = 0; i<InpEst; i++){
    val = 0;
    unc = 0;
    for(Int_t n = 0; n<InpObs; n++){
      val = val + Uma->operator()(i,n)*XvaRes->operator()(n);
      unc = unc + Uma->operator()(i,n)*CovRes->operator()(n,n);
    }
    N->operator()(i,0) = Xva->operator()(i) - val;
    //printf(" Pul %2i %2i %2i %2i \n", i, EstWhichObs(i), IsWhichEst(i),EstWhichObs(IsWhichEst(i)));
    //    Int_t IObs = EstWhichObs(IsWhichEst(i));
    Pul->operator()(i) = N->operator()(i,0) /
      sqrt(Sig->operator()(i) * Sig->operator()(i) - unc);
  }
  O->Transpose(*N);
  M->Mult(*CovI, *N);
  P->Mult(*O, *M);
  //M->Print(); N->Print(); O->Print(); P->Print(); 
  //Pul->Print();
  
  // Calculate ChiQua and NumDof
  ChiQua = P->operator()(0,0);
  NumDof = InpEst - InpObs;
  ChiPro = TMath::Prob(ChiQua,NumDof);

  // Clean up
  M->Delete(); N->Delete(); O->Delete(); P->Delete();
  
  // Success
  SetIsSolved(1);

  // Print some info if not in quiet mode
  if(IsQuiet() == 0){
    // Check what we got
    if(InspectResult() < 0){
      printf("... Blue->Solve(): Bad solution use SetPrintLevel(1) and InspectResult() to digest \n");
    }
    // Print quality of fit
    PrintChiPro();
  }

  return;
};

//------------------------------------------------------------------------------

void Blue::SolveRelUnc(const Double_t Dx){
  if(IsFilledInp() == 0){
    printf("... Blue->SolveRelUnc(%5.3f%%): IGNORED Not all estimates/uncertainties filled IGNORED \n", Dx);
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolveRelUnc(%5.3f%%): IGNORED Input not yet fixed call FixInp \n", Dx);
    return;
  }else if(IsSolvedRelUnc() == 1){
    printf("... Blue->SolveRelUnc(%5.3f%%: IGNORED Nothing changed since last call to this function.", Dx);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolveRelUnc(%5.3f%%): IGNORED Two consecutive calls to Solve() are not supported.", Dx);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
   }

  // Check precision
  if(Dx <= 0){
    printf("... Blue->SolveRelUnc(%5.3f%%): IGNORED I expect a relative precision in %% \n", Dx);
    return;
  }

  // Solve once 
  // This makes sure that XvaRes exists and can be used in the next FixInp
  Solve();

  // Find out if successfull
  if(IsSolved() != 1){
    printf("... Blue->SolveRelUnc(%5.3f%%): IGNORED Could not solve \n", Dx);
    return;
  }
  //PrintResult();

  // Set relative uncertainty flag
  SetRelValUnc(1);

  // Loop until convergence, or at most ILm times
  Int_t    Iok = 0, ILo = 0, ILm = 200;
  Double_t Old = 0, New = 0, Del = 0;
  TVectorD* XvaResLast = new TVectorD(GetActObs());
  while(Iok == 0){
    Iok = 1;
    ILo = ILo + 1;
    // Save result from last iteration
    for(Int_t n = 0; n<GetActObs(); n++){
      XvaResLast->operator()(n) = XvaRes->operator()(n);
    }
    // Solve next time
    ReleaseInp();
    FixInp();
    Solve();
    if(IsPrintLevel() >= 1){
      printf("... Blue->SolveRelUnc(%5.3f%%): Next Iteration = %2i \n", Dx, ILo);
    }
    // Loop over observables and estimates check difference
    for(Int_t n = 0; n<GetActObs(); n++){
      Old = XvaResLast->operator()(n);
      New = XvaRes->operator()(n);
      Del = 100.*TMath::Abs(1.-New/Old);
      if(Del > Dx)Iok = 0;
      if(IsPrintLevel() >= 1){
	printf("... Blue->SolveRelUnc(%5.3f%%): Old, New, Dif = %7.5f %7.5f %7.5f %% \n", Dx, Old, New, Del);
      }
    }
    if(ILo == ILm && Iok == 0){
      printf("... Blue->SolveRelUnc(%5.3f%%): No convergence after %3i iterations.", Dx, ILo);
      printf(" I give up and take the result from the last iteration. \n");
      Iok = 1;
      ILo = 0;
    }
  }
  if(ILo != 0 && IsQuiet() == 0){
    printf("... Blue->SolveRelUnc(%5.3f%%): Success after %2i iterations \n", Dx, ILo);
  }

  // Clean up
  XvaResLast->Delete(); XvaResLast = NULL; 

  // Clear relative uncertainty flag
  SetRelValUnc(0);

  // Success
  SetIsSolvedRelUnc(1);

  return;
};

//------------------------------------------------------------------------------

void Blue::SolveAccImp(const Double_t Dx){
  printf("... Blue->SolveAccImp(%5.3f%%): I use the new default option \n", Dx);
  SolveAccImp(0, Dx);
  return;
};

//------------------------------------------------------------------------------

void Blue::SolveAccImp(const Int_t ImpFla, const Double_t Dx){
  if(IsFilledInp() == 0){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED Not all estimates/uncertainties filled IGNORED \n", 
	   ImpFla, Dx);
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED Input not yet fixed call FixInp \n", 
	   ImpFla, Dx);
    return;
  }else if(IsSolvedAccImp() == 1){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED Nothing changed since last call to this function.", 
	   ImpFla, Dx);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED Two consecutive calls to Solve() are not supported.", 
	   ImpFla, Dx);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Check the precision
  if(Dx <= 0){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED I expect a relative precision in %% \n", ImpFla, Dx);
    return;
  }

  // Print the method
  if(ImpFla == 0){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): I sort according to relative improvement \n", ImpFla, Dx);
  }else if(ImpFla == 1){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): I sort according to absolute BLUE weights \n", ImpFla, Dx);
  }else if(ImpFla == 2){
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): I sort according to inverse variance of the estimates \n", ImpFla, Dx);
  }else{
    printf("... Blue->SolveAccImp(%2i,%5.3f%%): IGNORED invalid option \n", ImpFla, Dx);
    return;
  }

  // Save precision and Flag for PrintAccImp()
  PreAcc = Dx;
  IntAcc = ImpFla;

  // Solve once store estimates by importance
  Solve();
  
  // Go to quiet modus
  SetQuiet(1);

  // Reset matrices
  LasImp->Delete(); LasImp = new TMatrixD(1,InpObs);
  IndImp->Delete(); IndImp = new TMatrixD(InpEst,InpObs);
  ValImp->Delete(); ValImp = new TMatrixD(InpEst,InpObs);
  UncImp->Delete(); UncImp = new TMatrixD(InpEst,InpObs);
  StaImp->Delete(); StaImp = new TMatrixD(InpEst,InpObs);
  SysImp->Delete(); SysImp = new TMatrixD(InpEst,InpObs);

  // Save weights according to ImpFla in an interim matrix
  // Save the number of estimates per observable
  Int_t    IntPre = -1, ActPre = -1;
  Double_t    rho = -2,    zva = -2.;
  TMatrixD* Wei = new TMatrixD(InpEst,InpObs);
  Int_t     NnrEst[InpObs];
  for (Int_t n = 0; n < InpObs; n++){NnrEst[n] = 0;};
  for (Int_t n = 0; n < InpObs; n++){
    // Get most precise estimate i of observable n and its InpEst index
    IntPre = GetPreEst(IsWhichObs(n));
    ActPre = IsIndexEst(IntPre);
    //printf("... Blue->SolveAccImp(): Most precise estimate is %2i %2i \n", IntPre,ActPre);
    for(Int_t i = 0; i<InpEst; i++){
      if(EstWhichObs(IsWhichEst(i)) == IsWhichObs(n)){
	NnrEst[n] = NnrEst[n] + 1;
	if(ImpFla == 0){
	  if(i != ActPre){
	    rho = Rho->operator()(i,ActPre);
	    zva = Sig->operator()(i)/Sig->operator()(ActPre);
	    Wei->operator()(i,n) = 1.- GetPara(2, rho, zva);
	    //printf("... Blue->SolveAccImp(): %2i %2i  %2i %2i z=%5.2f rho=%5.2f sx/s1=%7.4f \n",
	    //		   i,ActPre,IsWhichEst(i),IntPre,zva,rho,Wei->operator()(i,n));
	  }else{
	    Wei->operator()(i,n) = 1.;
	  }
	}else if(ImpFla == 1){
	  Wei->operator()(i,n) = TMath::Abs(Lam->operator()(i,n));
	}else if(ImpFla == 2){
	  if(i != ActPre){
	    Wei->operator()(i,n) = 1./(Sig->operator()(i)*Sig->operator()(i));
	  }else{
	    Wei->operator()(i,n) = 1.;
	  }
	}
      }
    }
  }
  //printf("... Blue->SolveAccImp(): Lam \n"); Lam->Print();
  //printf("... Blue->SolveAccImp(): Wei \n"); Wei->Print();
  
  // Loop over observables
  Double_t compare = 0;
  Int_t ifound = -1;
  Int_t ilastf = -1;
  for (Int_t n = 0; n < InpObs; n++){    
    //printf("... Blue->SolveAccImp(): Next Obs %2i %2i \n", n,IsWhichObs(n));

    // Get most precise estimate i of observable n and its InpEst index
    IntPre = GetPreEst(IsWhichObs(n));
    ActPre = IsIndexEst(IntPre);
    //printf("... Blue->SolveAccImp(): Most precise estimate is %2i %2i \n", IntPre, ActPre);

    // Save most precise estimate at position zero, clear its weight
    LisImp[0] = IntPre;
    IndImp->operator()(0,n) = static_cast<double>(LisImp[0]);
    Wei->operator()(ActPre,n) = 0.0;

    // Find minimum weight
    Double_t MinWei = 100000.;
    for(Int_t i = 0; i<InpEst; i++){if(Wei->operator()(i,n) < MinWei)MinWei = Wei->operator()(i,n);};
    //printf("... Blue->SolveAccImp(): MinWei = %5.3f \n", MinWei);

    // Find the list of importance
    for(Int_t l = 1; l<NnrEst[n]; l++){
      //printf("... Blue->SolveAccImp(): Wei \n"); Wei->Print();
      compare = MinWei;
      ifound  = -1;
      for(Int_t j = 0; j<InpEst; j++){
	//printf("... Blue->SolveAccImp(): Compare: %8.5f %8.5f \n", compare,Wei->operator()(j,n));
	if(Wei->operator()(j,n) > compare){
	  compare = Wei->operator()(j,n);
	  ifound = j;
	}
      }
      if(ifound > -1 && ifound != ilastf){
	//printf("... Blue->SolveAccImp(): %2i which %2i \n", ifound,IsWhichEst(ifound));
	// Fill current
	LisImp[l] = IsWhichEst(ifound);
	IndImp->operator()(l,n) = static_cast<double>(LisImp[l]);
	// reset current
	Wei->operator()(ifound,n) = 0.0;
	ilastf = ifound;
      }
    }
    for(Int_t l = NnrEst[n]; l<InpEst; l++){
      LisImp[l] = -1;
      IndImp->operator()(l,n) = static_cast<double>(LisImp[l]);
    }
    //printf("... Blue->SolveAccImp(): IndImp \n"); IndImp->Print();
    ValImp->operator()(0,n) = Xva->operator()(ActPre);
    UncImp->operator()(0,n) = Sig->operator()(ActPre);
    StaImp->operator()(0,n) = Unc->operator()(0*InpEst+ActPre,0*InpEst+ActPre);
    SysImp->operator()(0,n) = TMath::Sqrt(UncImp->operator()(0,n)*UncImp->operator()(0,n)-
					  StaImp->operator()(0,n)*StaImp->operator()(0,n));

    // Print the findings
    if(IsPrintLevel() > 1){
      printf("... Blue->SolveAccImp(): The order of importance for observable %2i is: %2i ",
	     IsWhichObs(n),LisImp[0]);
      for(Int_t i = 1; i<InpEst; i++){
	if(LisImp[i] != -1){
	  printf(", %2i", LisImp[i]);
	}
      }
      printf("\n");
    }
  
    // Disable all but the most precise estimate 
    ReleaseInp();
    for(Int_t i = 1; i<NnrEst[n]; i++){SetInActiveEst(LisImp[i]);};
    FixInp();
      
    // Add one at a time for decreasing importance, store results
    for(Int_t i = 1; i<NnrEst[n]; i++){
      ReleaseInp();      
      SetActiveEst(LisImp[i]);
      FixInp();
      Solve();
      if(IsPrintLevel() > 0){
	printf("... Blue->SolveAccImp(): Add %2i, %s \n", LisImp[i],GetNamEst(IsWhichEst(i)).Data());
	PrintEst(ActPre);
	PrintResult();
      }
      ValImp->operator()(i,n) = XvaRes->operator()(n);
      UncImp->operator()(i,n) = TMath::Sqrt(CovRes->operator()(n,n));
      StaImp->operator()(i,n) = CorRes->operator()(0*InpObs+n,0*InpObs+n);
      SysImp->operator()(i,n) = TMath::Sqrt(UncImp->operator()(i,n)*UncImp->operator()(i,n)-
					    StaImp->operator()(i,n)*StaImp->operator()(i,n));
    }

    // Solve once if only one estimate for this observable
    if(NnrEst[n] == 1)Solve();

    //printf("... Blue->SolveAccImp(): XvaRes \n"); XvaRes->Print();
    //printf("... Blue->SolveAccImp(): CovRes \n"); CovRes->Print();
    //printf("... Blue->SolveAccImp(): LasImp \n"); LasImp->Print();
    //printf("... Blue->SolveAccImp(): IndImp \n"); IndImp->Print();
    //printf("... Blue->SolveAccImp(): ValImp \n"); ValImp->Print();
    //printf("... Blue->SolveAccImp(): UncImp \n"); UncImp->Print();
    //printf("... Blue->SolveAccImp(): StaImp \n"); StaImp->Print();
    //printf("... Blue->SolveAccImp(): SysImp \n"); SysImp->Print();
    //printf("... Blue->SolveAccImp():    Wei \n"); Wei->Print();
  }
  Wei->Delete(); Wei = NULL;

  // Find the last estimate to be combined
  Double_t DifVal = 0, DifUnc = 0;
  Int_t    NexImp = 0, NumImp = 0;
  for(Int_t n = 0; n<InpObs; n++){
    NumImp = 0;
    for(Int_t i = 1; i<InpEst; i++){
      NexImp = static_cast<int>(IndImp->operator()(i,n));
      if(NexImp != -1){
	DifUnc = 100*(UncImp->operator()(i,n) - UncImp->operator()(i-1,n))/UncImp->operator()(i-1,n);
	DifVal = 100*(ValImp->operator()(i,n) - ValImp->operator()(i-1,n))/ValImp->operator()(i-1,n);
	//printf("... Blue->SolveAccImp(): Add N = %2i ==> Differences: Value = %+5.3f %%, Uncertainty = %+5.3f %% \n",
	//NexImp,DifVal,DifUnc);
	if(DifUnc <= -PreAcc)NumImp = i;
      }
    }
    LasImp->operator()(0,n) = IndImp->operator()(NumImp,n);
    //printf("... Blue->SolveAccImp(): Only combine the first %2i estimates up to and including N = %2i \n",
    //NumImp+1, static_cast<int>(LasImp->operator()(0,n)));
  }

  // Remove the quiet modus
  SetQuiet(0);

  // Set the flag
  SetIsSolvedAccImp(1);

  // Print out if wanted
  if(IsPrintLevel() > 0)PrintAccImp();

  return;
};

//------------------------------------------------------------------------------

void Blue::SolveScaRho(){

  // Set ranges for rho scans to default values
  ResetScaRho(1);

  // Call the default solver
  FlaFac = 0;
  SolveScaRho(FlaFac);
  return;
};

//------------------------------------------------------------------------------

void Blue::SolveScaRho(const Int_t RhoFla){

  // Check the steering flag
  if(RhoFla < 0 || RhoFla > 1){
    printf("... Blue->SolveScaRho(%2i): IGNORED RhoFla = 0/1 is allowed \n", RhoFla);
    return;
  }
  FlaFac = RhoFla;

  // Set ranges for rho scans to default values
  ResetScaRho(1);

  // Call the solver
  SolveScaRho(FlaFac, MinFac, MaxFac);
  return;
};

//------------------------------------------------------------------------------

void Blue::SolveScaRho(const Int_t RhoFla, const Double_t *MinRho, const Double_t *MaxRho){

  // Check the steering flag
  if(RhoFla < 0 || RhoFla > 1){
    printf("... Blue->SolveScaRho(%2i, ...): IGNORED RhoFla = 0/1 is allowed \n", RhoFla);
    return;
  }
  FlaFac = RhoFla;

  // Take the input
  MaxFac->Delete(); MaxFac = new TMatrixD(InpUncOrig,InpFac,MaxRho,"MaxFac");
  MinFac->Delete(); MinFac = new TMatrixD(InpUncOrig,InpFac,MinRho,"MinFac");
     
  // Check for consistency
  Int_t IC = 1;
  for(Int_t k = 0; k<InpUncOrig; k++){
    for(Int_t l = 0; l<InpFac; l++){
      if(MaxFac->operator()(k,l) < MinFac->operator()(k,l)){
	printf("... Blue->SolveScaRho(%2i, ....): IGNORED inconsistent range Max(%2i,%2i) < Min(%2i,%2i) \n",
	       RhoFla , k, l, k, l);
	IC = 0;
      }
      if(MaxFac->operator()(k,l) > 1.){
	printf("... Blue->SolveScaRho(%2i, ....): IGNORED Max(%2i,%2i) > 1 \n", RhoFla,k,l);
	IC = 0;
      }else if(MaxFac->operator()(k,l) < -1.){
	printf("... Blue->SolveScaRho(%2i, ....): IGNORED Max(%2i,%2i) < -1 \n", RhoFla,k,l);
	IC = 0;
      }
      if(MinFac->operator()(k,l) > 1.){
	printf("... Blue->SolveScaRho(%2i, ....): IGNORED Min(%2i,%2i) > 1 \n", RhoFla,k,l);
	IC = 0;
      }else if(MinFac->operator()(k,l) < -1.){
	printf("... Blue->SolveScaRho(%2i, ....): IGNORED Min(%2i,%2i) < -1 \n", RhoFla,k,l);	
	IC = 0;
      }
    }
  }

  // Exit if failures occured
  if(IC == 0){
    printf("... Blue->SolveScaRho(%2i, ....): IGNORED please fix the input \n", RhoFla);	
    return;
  }

  // Call the solver
  SolveScaRho(FlaFac, MinFac, MaxFac);
  return;
};

//------------------------------------------------------------------------------

void Blue::SolveInfWei(){
  if(IsFilledInp() == 0){
    printf("... Blue->SolveInfWei(): IGNORED Not all estimates/uncertainties filled IGNORED \n");
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolveInfWei(): IGNORED Input not yet fixed call FixInp \n");
    return;
  }else if(IsSolvedInfWei() == 1){
    printf("... Blue->SolveInfWei(): IGNORED Nothing changed since last call to this function.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolveInfWei(): IGNORED Two consecutive calls to Solve() are not supported.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Check whether only one active observable
  if(GetActObs() > 1){
    printf("... Blue->SolveInfWei(): IGNORED more than one active observable. Only one is allowed \n");
    return;
  }

  // Reset vectors for information weights
  VarInd->Delete(); VarInd = new TVectorD(InpEstOrig);
  IntWei->Delete(); IntWei = new TVectorD(InpEstOrig+1);
  MarWei->Delete(); MarWei = new TVectorD(InpEstOrig);
  BluWei->Delete(); BluWei = new TVectorD(InpEstOrig);

  // Go to quiet modus
  SetQuiet(1);

  // Solve for all but one, save the variances
  ReleaseInp();
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i)){
      SetInActiveEst(i);
      FixInp();
      Solve();
      ReleaseInp();
      SetActiveEst(i);
      VarInd->operator()(i) = CovRes->operator()(0,0);
    }
  }
  //printf("... Blue->SolveInfWei(): VarInd\n"); VarInd->Print();
  
  // Solve with all estimates, save 
  // 1) variance of the result
  // 2) variances of the active estimates
  // 3) Blue weights of the active estimates
  FixInp();
  Solve();
  //printf("... Blue->SolveInfWei(): Sig \n"); Sig->Print();  
  for(Int_t i = 0; i<InpEst; i++){
    IntWei->operator()(IsWhichEst(i)) = Sig->operator()(i)* Sig->operator()(i);  
    BluWei->operator()(IsWhichEst(i)) = Lam->operator()(i,0);
  }
  //printf("... Blue->SolveInfWei(): IntWei \n"); IntWei->Print();  
  
  // Calculate the information weights for estimates i and the correlation
  // MarWei(i) = 1 - Var_comb(all) / Var_comb(all-i)
  // IntWei(i) = Var_comb(all) / Var_(i)
  // IntWei(InpEstOrig) = 1 - Var_comb(all) * Sum_i(1/(Var_(i))
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i)){
      MarWei->operator()(i) = 1 - CovRes->operator()(0,0)/ VarInd->operator()(i);
      IntWei->operator()(i) = CovRes->operator()(0,0) / IntWei->operator()(i);
      IntWei->operator()(InpEstOrig) = IntWei->operator()(InpEstOrig) + IntWei->operator()(i);
    }
  }
  IntWei->operator()(InpEstOrig) = 1 - IntWei->operator()(InpEstOrig);

  // Remove the quiet modus
  SetQuiet(0);

  // Set the flag
  SetIsSolvedInfWei(1);

  // Print out if wanted
  if(IsPrintLevel() > 0){PrintInfWei();};

  return;
};

//------------------------------------------------------------------------------

void Blue::SolvePosWei(){
  if(IsFilledInp() == 0){
    printf("... Blue->SolvePosWei(): IGNORED Not all estimates/uncertainties filled IGNORED \n");
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolvePosWei(): IGNORED Input not yet fixed call FixInp \n");
    return;
  }else if(IsSolvedPosWei() == 1){
    printf("... Blue->SolvePosWei(): IGNORED Nothing changed since last call to this function.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolvePosWei(): IGNORED Two consecutive calls to Solve() are not supported.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Solve once
  Solve();

  // Due to de-activation of estimates new estimates can get negative weights
  // Iterate until all are positive
  TVectorD* IDeAct = new TVectorD(InpEstOrig);
  Int_t NDeAct = 1, NumCal = 0, NumDea = 0;

  while(NDeAct > 0){
    NDeAct = 0;
    // Store estimates i of observable n with negative weights
    //printf("... Blue->SolvePosWei(): Lam \n"); Lam->Print();
    for(Int_t n = 0; n<InpObs; n++){    
      //printf("... Blue->SolvePosWei(): Next Obs is %2i \n", IsWhichObs(n));
      for(Int_t i = 0; i<InpEst; i++){     
	if(EstWhichObs(IsWhichEst(i)) == IsWhichObs(n)){
	  //printf("... Blue->SolvePosWei(): Next Est=%2i(%2i) determines Obs=%2i \n", 
	  //	 i,IsWhichEst(i),IsWhichObs(n));
	  if(Lam->operator()(i,n) < 0){
	    IDeAct->operator()(NDeAct) = IsWhichEst(i);
	    NDeAct = NDeAct + 1;
	  }
	}
      }
    }
    // Now deactivate if needed and solve
    if(NDeAct > 0){
      ReleaseInp();
      for(Int_t l = 0; l<NDeAct; l++){
	NumDea = NumDea + 1;
	SetInActiveEst(static_cast<int>(IDeAct->operator()(l)));
      }
      NumCal = NumCal + 1;
      FixInp();
      Solve();
    }
  }
  
  // Report the findings
  if(NumCal == 0){
    printf("... Blue->SolvePosWei(): Not needed all weights are positive \n");
  }else{
    printf("... Blue->SolvePosWei(): Succcess after disabling %2i estimates in %2i iterations \n", NumDea, NumCal);
  }

  // Set the flag
  SetIsSolvedPosWei(1);

  // Clean up
  IDeAct->Delete(); IDeAct = NULL;

  return;
};

//------------------------------------------------------------------------------

void Blue::SolveMaxVar(const Int_t IFuRho){
  if(IsFilledInp() == 0){
    printf("... Blue->SolveMaxVar(%2i): IGNORED Not all estimates/uncertainties filled IGNORED \n",IFuRho);
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolveMaxVar(%2i): IGNORED Input not yet fixed call FixInp \n",IFuRho);
    return;
  }else if(IsSolvedMaxVar() == 1){
    printf("... Blue->SolveMaxVar(%2i): IGNORED Nothing changed since last call to this function.",IFuRho);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolveMaxVar(%2i): IGNORED Two consecutive calls to Solve() are not supported.", IFuRho);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Check whether only one active observable
  if(GetActObs() > 1){
    printf("... Blue->SolveMaxVar(%2i): IGNORED more than one active observable only one is allowed \n",IFuRho);
    return;
  }

  // Check the steering flag
  if(IFuRho < -1 ||IFuRho > 2){
    printf("... Blue->SolveMaxVar(%2i): Not Implemented \n", IFuRho);
    return;
  }

  // Set the flag to signal the algorithm to PrintMaxVar()
  IFuFla = TMath::Abs(IFuRho);

  // Check for other correlation steering if Abs(IFuRho) = 2
  Int_t IntFai = 0;
  if(IFuFla ==  2){
    for(Int_t k = 0; k<InpUncOrig; k++){    
      if(IsRhoRedUnc(k) == 1 || IsRhoFacUnc(k) == 1 || IsRhoValUnc(k) == 1)IntFai = 1;
      if(IsRhoRedUnc(k) == 1)printf("... Blue->SolveMaxVar(%2i): Source %2i has -reduced- correlations \n", IFuRho,k);
      if(IsRhoFacUnc(k) == 1)printf("... Blue->SolveMaxVar(%2i): Source %2i has -scaled - correlations \n", IFuRho,k);
      if(IsRhoValUnc(k) == 1)printf("... Blue->SolveMaxVar(%2i): Source %2i has -changed- correlations \n", IFuRho,k);
    }
  }
  if(IntFai == 1){
    printf("... Blue->SolveMaxVar(%2i): IGNORED Not possible \n", IFuRho);
    return;
  }

  // Set the number of failures 
  for(Int_t ind = 0; ind < 3; ind++)IFuFai[ind] = 0;

  // Solve once, save variance and correlations at initial state
  Solve();

  // The smallest variance of any estimate
  Double_t VarMinEst = Sig->operator()(0);
  for(Int_t i = 1; i < InpEst; i++){
    if(VarMinEst > Sig->operator()(i))VarMinEst = Sig->operator()(i);
  }
  VarMinEst = VarMinEst * VarMinEst;
  //printf("... Blue->SolveMaxVar(%2i): VarMinEst = %5.3f \n",IFuRho,VarMinEst);

  // The variances before and after
  VarMax->Delete(); VarMax = new TVectorD(InpObs);
  VarMin->Delete(); VarMin = new TVectorD(InpObs);
  for(Int_t n = 0; n < InpObs; n++){
    VarMin->operator()(n) = CovRes->operator()(n,n);
    VarMax->operator()(n) = CovRes->operator()(n,n);
  }
  //printf("... Blue->SolveMaxVar(%2i): VarMin \n",IFuRho); VarMin->Print();
  //printf("... Blue->SolveMaxVar(%2i): VarMax \n",IFuRho); VarMax->Print();

  // The correlations before and after
  RhoMax->Delete(); RhoMax = new TMatrixD(InpEst,InpEst);
  RhoMin->Delete(); RhoMin = new TMatrixD(InpEst,InpEst);
  RhoMin->SetSub(0,0,*Rho);
  //printf("... Blue->SolveMaxVar(%2i): RhoMin \n",IFuRho); RhoMin->Print();

  // The correlations that equal 1/z
  RhoOoz->Delete(); RhoOoz = new TMatrixD(InpEst,InpEst);
  for(Int_t i = 0; i < InpEst; i++){
    for(Int_t j = i+1; j < InpEst; j++){
      if(Rho->operator()(i,j) > 0){
	RhoOoz->operator()(i,j) = TMath::Min(Sig->operator()(i)/Sig->operator()(j),
					     Sig->operator()(j)/Sig->operator()(i));
	RhoOoz->operator()(j,i) = RhoOoz->operator()(i,j);
      }
    }
    RhoOoz->operator()(i,i) = 1;
  }
  //printf("... Blue->SolveMaxVar(%2i): RhoOoz \n",IFuRho); RhoOoz->Print();

  // The fijk values to fill for IFuFla=0/1/2
  // Need original dimensions
  RhoFpk->Delete(); RhoFpk = new TVectorD(InpUncOrig);
  RhoFij->Delete(); RhoFij = new TMatrixD(InpEstOrig,InpEstOrig);

  // Set quiet modus
  SetQuiet(1);

  // Loop over observables
  Double_t RhoFac = 1; 
  Double_t RhoLim = 0; 
  //Double_t RhoDif = 0.001;
  Double_t RhoDif = 0.01;
  Double_t RhoSav = 1; 
  Double_t RhoTes = 0; 
  Double_t VarRun = 1; 
  for(Int_t n = 0; n < InpObs; n++){    
    //printf("... Blue->SolveMaxVar(%2i): The next observable is %2i %2i \n",IFuRho,n,IsWhichObs(n));
    VarRun = VarMax->operator()(n);
    // Loop over RhoFac values
    if(IFuFla == 0){
      // Find the Fijk = F with MaxVar
      RhoFac = 1 - RhoDif;
      RhoFco = 1;
      while(RhoFac >= RhoLim){
	ReleaseInp();
	//printf("... Blue->SolveMaxVar(%2i): RhoFac %5.3f \n",IFuRho,RhoFac);
	for(Int_t k = 1; k<InpUncOrig; k++){if(IsActiveUnc(k))SetRhoFacUnc(k,RhoFac);};
	FixInp();
	Solve();
	if(InspectResult() < 0){	  	  
	  IFuFai[IFuFla] = IFuFai[IFuFla] + 1;	  
	  if(IsPrintLevel() > 2){
	    printf("... Blue->SolveMaxVar(%2i): Strange Covariance %2i %5.3f %12.9f \n",
		   IFuRho,InspectResult(),RhoFac,CovRes->operator()(n,n));
	  }
	}else if(CovRes->operator()(n,n) > VarRun){
	  //printf("... Blue->SolveMaxVar(%2i): old %12.9f new %12.9f \n",IFuRho,VarRun,CovRes->operator()(n,n));
	  VarRun = CovRes->operator()(n,n);
	  RhoFco = RhoFac;
	}
        RhoFac = RhoFac - RhoDif; 
      }
      // Set the RhoFco. Solve finally
      ReleaseInp();
      for(Int_t k = 1; k<InpUncOrig; k++){if(IsActiveUnc(k))SetRhoFacUnc(k,RhoFco);};
      FixInp();
      Solve();
    }else if(IFuFla == 1){
      // Find the RhoFac[k] with MaxVar
      RhoFpk->operator()(0) = 1;
      for(Int_t k = 1; k<InpUncOrig; k++){
	// Reset for each fk if wanted 
	if(IFuRho > 0)VarRun = VarMax->operator()(n);
	RhoFpk->operator()(k) = 1;
	if(IsActiveUnc(k)){
	  RhoFac = 1 - RhoDif; 
	  //printf("... Blue->SolveMaxVar(%2i): Next source  %2i \n",IFuRho,k);
	  while(RhoFac >= RhoLim){
	    //printf("... Blue->SolveMaxVar(%2i): RhoFac %5.3f \n",IFuRho,RhoFac);
	    if(IsFixedInp() == 1)ReleaseInp();
	    SetRhoFacUnc(k,RhoFac);
	    FixInp();
	    Solve();
	    if(InspectResult() < 0){	  	  
	      IFuFai[IFuFla] = IFuFai[IFuFla] + 1;	  
	      if(IsPrintLevel() > 2){
		printf("... Blue->SolveMaxVar(%2i): Strange Covariance %2i %2i %5.3f %12.9f \n",
		       IFuRho,k,InspectResult(),RhoFac,CovRes->operator()(n,n));
	      }
	    }else if(CovRes->operator()(n,n) > VarRun){
	      //printf("... Blue->SolveMaxVar(%2i): old %12.9f new %12.9f \n",IFuRho,VarRun,CovRes->operator()(n,n));
	      VarRun = CovRes->operator()(n,n);
	      RhoFpk->operator()(k) = RhoFac;
	    }
	    RhoFac = RhoFac - RhoDif;
	  }
	  // Reset Fk for this source
	  ReleaseInp();
	  if(IFuRho > 0)SetRhoFacUnc(k,1);
	}
      }
      // Set all RhoFpk(k=1...NumUnc-1). Solve finally
      ReleaseInp();
      for(Int_t k = 1; k<InpUncOrig; k++){if(IsActiveUnc(k))SetRhoFacUnc(k,RhoFpk->operator()(k));};
      FixInp();
      Solve();
    }else if(IFuFla == 2){
      //printf("... Blue->SolveMaxVar: In Rho \n"); Rho->Print();
      // Preset Futch factors to unity
      ReleaseInp();
      RhoFac = 1; 
      for(Int_t i = 0; i<InpEstOrig; i++){
	for(Int_t j = i; j<InpEstOrig; j++){
	  SetFutchCor(i,j,1,RhoFac);
	}
      }      
      // The indicies of active objects
      Int_t ia = 0, ja = 0;
      //printf("... Blue->SolveMaxVar(%2i): RhoFij \n",IFuRho); RhoFij->Print();
      // Now optimize one by one
      for(Int_t i = 0; i<InpEstOrig; i++){
	if(IsActiveEst(i) == 1){
	  ja = ia + 1;
	  for(Int_t j = i+1; j<InpEstOrig; j++){
	    if(IsActiveEst(j) == 1){
	      RhoFac = 1 - RhoDif; 
	      RhoSav = 1; 
	      RhoTes = Rho->operator()(ia,ja);
	      while(RhoFac >= RhoLim && RhoTes > 0){
		if(IsFixedInp() == 1)ReleaseInp();
		SetFutchCor(i,j,1,RhoFac);	  
		//printf("... Blue->SolveMaxVar(%2i): RhoFij \n",IFuRho); RhoFij->Print();
		FixInp();
		Solve();
		//printf("... Blue->SolveMaxVar(%2i): old %12.9f new %12.9f \n",IFuRho,VarRun,CovRes->operator()(n,n));
		if(InspectResult() < 0){	  	  
		  IFuFai[IFuFla] = IFuFai[IFuFla] + 1;	  
		  if(IsPrintLevel() > 2){
		    printf("... Blue->SolveMaxVar(%2i): Strange Covariance %2i %2i %2i %5.3f %12.9f \n",
			   IFuRho,i,j,InspectResult(),RhoFac,CovRes->operator()(n,n));
		  }
		}else if(CovRes->operator()(n,n) > VarRun){
		  VarRun = CovRes->operator()(n,n);
		  RhoSav = RhoFac;
		}
		RhoFac = RhoFac - RhoDif;
	      }
	      // Set the next achieved Fijk and Solve
	      ReleaseInp();
	      SetFutchCor(i,j,1,RhoSav);
	      FixInp();
	      Solve();
	      VarRun = CovRes->operator()(n,n);
	      if(InspectResult() < 0){
		printf("... Blue->SolveMaxVar(%2i): Strange Covariance at End %2i %2i %5.3f %12.9f \n",
		       IFuRho,i,j,RhoFac,CovRes->operator()(n,n));
		printf("... Blue->SolveMaxVar: (%2i,%2i)=(%2i,%2i) In=%4.2f, Fac=%5.3f, Out=%4.2f, CovRes=%12.9f \n",
		       i,j,ia,ja,RhoTes,RhoSav,RhoTes*RhoSav,VarRun);
	      }
	      ja = ja + 1;
	    }
	  }
	  ia = ia + 1;
	}
      }
    }
    VarMax->operator()(n) = CovRes->operator()(n,n);
    RhoMax->SetSub(0,0,*Rho);
  }
  //printf("... Blue->SolveMaxVar(%2i): RhoMax \n",IFuRho); RhoMax->Print();
  //printf("... Blue->SolveMaxVar(%2i): VarMax \n",IFuRho); VarMax->Print();
  
  // Remove quiet modus and solve with the latest settings
  SetQuiet(0);
  
  // Set the flag
  SetIsSolvedMaxVar(1);

  // Print out if wanted
  if(IsPrintLevel() > 0)PrintMaxVar();

  return;
};

//------------------------------------------------------------------------------
// Filler
//------------------------------------------------------------------------------

void Blue::FillEst(const Int_t i, const Double_t *x){
  if(IsFilledEst() == 1){return;};
  if(IsAllowedEst(i) == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->FillEst(%2i): Fill estimate \n", i);
    }

    // Estimates
    Xva->operator()(i) = x[0];
    //Xva->Print();

    // Uncertainties
    Int_t IRow;
    Double_t sig = 0, val = 0;
    for(Int_t k = 0; k<InpUncOrig; k++){
      // The individual values
      val = x[k+1];
      if(x[k+1] < 0){
	if(IsPrintLevel() >= 2){
	  printf("... Blue->FillEst: I fill this as percentage uncertainty %5.3f%% \n",x[k+1]);
	}
	val = -1. * x[0] * x[k+1] / 100.;
      }

      double* unc = &val;
      IRow=i+k*InpEst;
      //printf("... Blue->FillEst: IRow = %2i, ICol = %2i, val = %5.3f \n",  IRow, IRow, x[k+1]);
      Unc->InsertRow(IRow,IRow,unc,1);
      // The sum
      sig = sig + val*val;
    }

    // Total uncertainty
    Sig->operator()(i) = sqrt(sig);

    //Set flags
    SetActiveEst(i,1);
    SetFilledEst(i);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::FillCor(const Int_t k, const Double_t *x){
  if(IsFilledUnc() == 1)return;
  Int_t Ilong = 1;
  Int_t kk = k;
  if(kk < 0){
    kk = -kk;
    Ilong = 0;
  }
  Int_t IsPosssible = 1;
  if(IsAllowedUnc(kk) == 1){
    if(IsPrintLevel() >= 1)printf("... Blue->FillCor(%2i): Fill correlation matrix \n", k);
    TMatrixD *F;
    if(Ilong == 1){
      F = new TMatrixD(InpEstOrig,InpEstOrig,x,"F");
      // Check if F is a symmetric matrix. If not this is fatal.
      for(Int_t i = 0; i<InpEstOrig; i++){
	for(Int_t j = i+1; j<InpEstOrig; j++){
	  if(F->operator()(j,i) != F->operator()(i,j)){
	    printf("... Blue->FillCor(%2i): The matrix is not symmetric, CHECK INPUT \n", k);
	    printf("... Blue->FillCor(%2i): This mistake for (%2i, %2i) = %5.3f",
		   k,i,j,F->operator()(j,i));
	    printf(" is fatal ! \n");
	    IsPosssible = 0;
	  };
	}
      }
    }else{
      F = new TMatrixD(InpEstOrig,InpEstOrig);      
      F->UnitMatrix();
      Int_t ind = 0;
      for(Int_t i = 0; i<InpEstOrig; i++){
	for(Int_t j = i+1; j<InpEstOrig; j++){
	  F->operator()(i,j) = x[ind];
	  F->operator()(j,i) = F->operator()(i,j);
	  ind = ind + 1;
	}
      }
      //F->Print();
    }

    // Force all diagonal elements to one
    for(Int_t i = 0; i<InpEstOrig; i++)F->operator()(i,i) = 1.0;

    // Check if all elements are within the allowed range of (-1,1). 
    // If not this is fatal.
    for(Int_t i = 0; i<InpEstOrig; i++){
      for(Int_t j = i+1; j<InpEstOrig; j++){	  
	if(F->operator()(i,j)< -1.0 || F->operator()(i,j) > 1.0){
	  printf("... Blue->FillCor(%2i): One element is out of bounds (-1,1), CHECK INPUT \n",k);
	  printf("... Blue->FillCor(%2i): This mistake for (%2i, %2i) = %5.3f",
		 k,i,j,F->operator()(j,i));
	  printf(" is fatal ! \n");
	  IsPosssible = 0;
	};
      }
    }

    // Exit if false input is detected
    if(IsPosssible == 0){
      printf("... Blue->FillCor(%2i): NO COMBINATION POSSIBLE --> return \n", k);
      return;
    }else{
      //F->Print();
    };

    // Fill into local structures
    Int_t IRow = kk*InpEst;
    Int_t ICol = kk*InpEst;
    Cor->SetSub(IRow,ICol,*F);
    //Cor->Print();

    // Cean up
    F->Delete(); F = NULL;

    // Set Flags
    SetActiveUnc(kk,1);
    SetFilledUnc(kk);
  }
};

//------------------------------------------------------------------------------

void Blue::FillCor(const Int_t k, const Double_t rho){
  if(IsFilledUnc() == 1){
    return;
  }
  if(IsAllowedUnc(k) == 1){
    if(rho < -1.0 || rho > 1.0){
      printf("... Blue->FillCor(%2i,%5.3f): IGNORED Correlation out of bound \n",
	     k,rho);
      printf("... Blue->FillCor(%2i): NO COMBINATION POSSIBLE --> return \n", k);
      return;
    }
    if(IsPrintLevel() >= 1){
      printf("... Blue->FillCor(%2i,%5.3f): Fill correlation matrix \n", k,rho);
    }
    TMatrixD *F = new TMatrixD(InpEstOrig,InpEstOrig);
    F->UnitMatrix();
    for(Int_t i = 0; i<InpEstOrig; i++){
      for(Int_t j = i+1; j<InpEstOrig; j++){
	F->operator()(i,j) = rho;
	F->operator()(j,i) = rho;
      }
    }
    Int_t IRow = k*InpEstOrig;
    Cor->SetSub(IRow,IRow,*F);
    //Cor->Print();

    // Cean up
    F->Delete(); F = NULL;

    // Set Flags
    SetActiveUnc(k,1);
    SetFilledUnc(k);
  }
};

//------------------------------------------------------------------------------

void Blue::FillNamEst(const TString* NamEst){
  if(IsFilledNamEst() == 1){
    printf("... Blue->FillNamObs(): These names have been filled already \n");
    return;
  }
  if(IsFilledInp() == 1){
    printf("... Blue->FillNamEst(): Input is closed, please fill names first \n");
    return;
  }
  for(Int_t i = 0; i<InpEstOrig; i++){
    EstNam[i] = NamEst[i];
    if(IsPrintLevel() >= 2){
      printf("... Blue->FillNamEst(%2i): %s \n",i,EstNam[i].Data());
    }
  }
};

//------------------------------------------------------------------------------

void Blue::FillNamUnc(const TString* NamUnc){
  if(IsFilledNamUnc() == 1){
    printf("... Blue->FillNamObs(): These names have been filled already \n");
    return;
  }
  if(IsFilledInp() == 1){
    printf("... Blue->FillNamUnc(): Input is closed, please fill names first \n");
    return;
  }
  for(Int_t k = 0; k<InpUncOrig; k++){
    UncNam[k] = NamUnc[k];
    if(IsPrintLevel() >= 2){
      printf("... Blue->FillNamUnc(%2i): %s \n",k,UncNam[k].Data());
    }
  }
};

//------------------------------------------------------------------------------

void Blue::FillNamObs(const TString* NamObs){
  if(IsFilledNamObs() == 1){
    printf("... Blue->FillNamObs(): These names have been filled already \n");
    return;
  }
  if(IsFilledInp() == 1){
    printf("... Blue->FillNamObs(): Input is closed, please fill names first \n");
    return;
  }
  for(Int_t n = 0; n<InpObsOrig; n++){
    ObsNam[n] = NamObs[n];
    if(IsPrintLevel() >= 2){
      printf("... Blue->FillNamObs(%2i): %s \n",n,ObsNam[n].Data());
    }
  }
};

//------------------------------------------------------------------------------
// Getters
//------------------------------------------------------------------------------

// After FixInp()
Int_t Blue::GetActEst()const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetActEst(): Number of active     estimates = %3i \n", InpEst);
    }
    return InpEst;
  }else{
    printf("... Blue->GetActEst(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetActEst(const Int_t n)const{
  Int_t InpEstObs = 0;
  if(IsFixedInp() == 1){
    if(IsActiveObs(n) == 1){
      for(Int_t i = 0; i<InpEstOrig; i++){
	if(IsActiveEst(i) == 1 && EstWhichObs(i) == n)InpEstObs = InpEstObs + 1;
      }
      if(IsPrintLevel() >= 1){
	printf("... Blue->GetActEst(%2i): Number of active estimates per observable = %3i \n", 
	       n, InpEstObs);
      }
    }else{
      printf("... Blue->GetActEst(%2i): Not an active observable \n", n);
    }
  }else{
    printf("... Blue->GetActEst(%2i): Presently not available, call FixInp() \n", n);
  }
  return InpEstObs;
};

//------------------------------------------------------------------------------

Int_t Blue::GetActUnc()const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetActUnc: Number of active uncertainties = %3i \n", InpUnc);
    }
    return InpUnc;
  }else{
    printf("... Blue->GetActUnc(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetActObs()const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetActObs: Number of active   observables = %3i \n", InpObs);
    }
    return InpObs;
  }else{
    printf("... Blue->GetActObs(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

TString Blue::GetNamEst(const Int_t i)const{
  if(IsFixedInp() == 1){
    if(IsActiveEst(i) == 1){
      return EstNamOrig[i];
    }else{
      printf("... Blue->GetNamEst(%2i): Not an active estimate \n", i);
    }
  }else{
    printf("... Blue->GetNamEst(%2i): IGNORED Input not yet fixed call FixInp \n", i);
  }
  return "NULL";
};

//------------------------------------------------------------------------------

TString Blue::GetNamUnc(const Int_t k)const{
  if(IsFixedInp() == 1){
    if(IsActiveUnc(k) == 1){
      return UncNamOrig[k];
    }else{
      printf("... Blue->GetNamUnc(%2i): Not an active uncertainty \n", k);
    }
  }else{
    printf("... Blue->GetNamUnc(%2i): IGNORED Input not yet fixed call FixInp \n", k);
  }
  return "NULL";
};

//------------------------------------------------------------------------------

TString Blue::GetNamObs(const Int_t n)const{
  if(IsFixedInp() == 1){
    if(IsActiveObs(n) == 1){
      return ObsNamOrig[n];
    }else{
      printf("... Blue->GetNamObs(%2i): Not an active observable \n", n);
    }
  }else{
    printf("... Blue->GetNamObs(%2i): IGNORED Input not yet fixed call FixInp \n", n);
  }
  return "Null";
};

//------------------------------------------------------------------------------

Int_t Blue::GetIndEst(Int_t *IndEst)const{
  if(IsFixedInp() == 1){
    for(Int_t i = 0; i<InpEst; i++){
      IndEst[i] = LisEst[i];
      if(IsPrintLevel() >= 1){
	if(i == 0){
	  printf("... Blue->GetIndEst: List of active     estimates = %2i", IndEst[i]);
	}else{
	  printf(", %2i", IndEst[i]);
	  if(i == InpEst-1)printf("\n");
	}
      }
    }
    return 1;
  }else{
    printf("... Blue->GetIndEst(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetIndUnc(Int_t *IndUnc)const{
  if(IsFixedInp() == 1){
    for(Int_t k = 0; k<InpUnc; k++){
      IndUnc[k] = LisUnc[k];
      if(IsPrintLevel() >= 1){
	if(k == 0){
	  printf("... Blue->GetIndUnc: List of active uncertainties = %2i", IndUnc[k]);
	}else{
	  printf(", %2i", IndUnc[k]);
	  if(k == InpUnc-1)printf("\n");
	}
      }
    }
    return 1;
  }else{
    printf("... Blue->GetIndUnc(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetIndObs(Int_t *IndObs)const{
  if(IsFixedInp() == 1){
    for(Int_t n = 0; n<InpObs; n++){
      IndObs[n] = LisObs[n];
      if(IsPrintLevel() >= 1){
	if(n == 0){
	  printf("... Blue->GetIndObs: List of active   observables = %2i", IndObs[n]);
	}else{
	  printf(", %2i", IndObs[n]);
	  if(n == InpObs-1)printf("\n");
	}
      }
    }
    return 1;
  }else{
    printf("... Blue->GetIndEst(): Presently not available, call FixInp() \n");
  }
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetPreEst(const Int_t n) const{

  // Set variables
  Int_t    IndMos = -1;
  Double_t ValMos =  0;

  // Check if fixed then loop over active estimates 
  if(IsFixedInp() == 1){
    if(IsAllowedObs(n) == 1 && IsActiveObs(n) == 1){
      for(Int_t i = 0; i<InpEst; i++){
	if(EstWhichObs(IsWhichEst(i)) == n){
	  if(IndMos == -1 || Sig->operator()(i) < ValMos){
	    ValMos = Sig->operator()(i);
	    IndMos = IsWhichEst(i);
	  }
	}
      }
    }
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetPreEst(%2i): The most precise estimate is: %2i \n", n, IndMos);
    }
  }else{
    printf("... Blue->GetPreEst(%2i): IGNORED Input not yet fixed call FixInp \n", n);
  }
  return IndMos;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCov(TMatrixD *UseCov)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCov: Return covariance matrix as TMatrixD \n");
    }
    UseCov->SetSub(0,0,*Cov);
    return 1;
  }
  printf("... Blue->GetCov: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCov(Double_t *RetCov)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCov: Return covariance matrix as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpEst,InpEst);
    Int_t IRet = GetCov(Dumm);
    MatrixtoDouble(Dumm, RetCov);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetCov: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCovInvert(TMatrixD *UseCovI)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCovInvert: Return inverse covariance matrix as TMatrixD \n");
    }
    UseCovI->SetSub(0,0,*CovI);
    return 1;
  }
  printf("... Blue->GetCovInvert: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCovInvert(Double_t *RetCovI)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCovInvert: Return inverse covariance matrix as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpEst,InpEst);
    Int_t IRet = GetCovInvert(Dumm);
    MatrixtoDouble(Dumm, RetCovI);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetCovInvert: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetRho(TMatrixD *UseRho)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetRho: Return correlation matrix as TMatrixD \n");
    }
    UseRho->SetSub(0,0,*Rho);
    return 1;
  }
  printf("... Blue->GetRho: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetRho(Double_t *RetRho)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetRho: Return correlation matrix as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpEst,InpEst);
    Int_t IRet = GetRho(Dumm);
    MatrixtoDouble(Dumm, RetRho);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetRho: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEst(TMatrixD *UseEst)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEst: Return the estimates with their uncertainties as TMatrixD \n");
    }

    // Fill the matrix
    // Xva->Print();
    // Unc->Print();
    for(Int_t i = 0; i<InpEst; i++){
      UseEst->operator()(i,0) = Xva->operator()(i);
      for(Int_t k = 0; k<InpUnc; k++){
	UseEst->operator()(i,k+1) = Unc->operator()(i+k*InpEst,i+k*InpEst);
      }
    }
    //UseEst->Print();
    return 1;
  }
  printf("... Blue->GetEst: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEst(Double_t *RetEst)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEst:  Return the estimates with their uncertainties as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpEst,InpUnc+1);
    Int_t IRet = GetEst(Dumm);
    MatrixtoDouble(Dumm, RetEst);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetEst: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEstVal(TMatrixD *UseEstVal)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEstVal: Return the values of the estimates as TMatrixD \n");
    }

    // Fill the Matrix
    for(Int_t i = 0; i<InpEst; i++)UseEstVal->operator()(i,0) = Xva->operator()(i);
    return 1;
  }
  printf("... Blue->GetEstVal: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEstVal(Double_t *RetEstVal)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEstVal: Return the values of the estimates as Double_t array \n");
    }
    // Fill the matrix
    for(Int_t i = 0; i<InpEst; i++)RetEstVal[i] = Xva->operator()(i);
    return 1;
  }
  printf("... Blue->GetEstVal: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEstUnc(TMatrixD *UseEstUnc)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEstUnc: Return the total uncertainties of the estimates as TMatrixD \n");
    }
    // Fill the vector
    for(Int_t i = 0; i<InpEst; i++)UseEstUnc->operator()(i,0) = Sig->operator()(i);
    return 1;
  }
  printf("... Blue->GetEstUnc: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetEstUnc(Double_t *RetEstUnc)const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetEstUnc: Return the total uncertainties of the estimates as Double_t array \n");
    }
    // Fill the array
    for(Int_t i = 0; i<InpEst; i++)RetEstUnc[i] = Sig->operator()(i);
    return 1;
  }
  printf("... Blue->GetEstUnc: Presently not available, call FixInp() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetParams(const Int_t Ifl, TMatrixD *UseParams) const {
  if(IsCalcedParams() == 1){
    // Fill the matrix we want
    TMatrixD* Dumm = new TMatrixD(InpEstOrig,InpEstOrig);
    if(Ifl == 0){Dumm->operator=(*SRat);
    }else if(Ifl == 1){Dumm->operator=(*Beta);
    }else if(Ifl == 2){Dumm->operator=(*Sigx);
    }else if(Ifl == 3){Dumm->operator=(*DBdr);
    }else if(Ifl == 4){Dumm->operator=(*DSdr);
    }else if(Ifl == 5){Dumm->operator=(*DBdz);
    }else if(Ifl == 6){Dumm->operator=(*DSdz);
    }else{printf("... Blue->GetParams(%2i): Not implemented flag \n",Ifl);
    }
    // Now fill the active estimates
    Int_t ii = 0, jj = 0;
    for(Int_t i = 0; i<InpEstOrig; i++){
      jj = 0;
      if(IsActiveEst(i) == 1){	
	for(Int_t j = 0; j<i; j++){
	  if(IsActiveEst(j) == 1){
	    UseParams->operator()(ii,jj) = Dumm->operator()(i,j);
	    jj = jj + 1;
	  }
	}
	ii = ii + 1;
      }
    }
    Dumm->Delete(); Dumm = NULL;
    return 1;
  }else{
    printf("... Blue->GetParams(%2i): Presently not available, call FixInp() \n",Ifl);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::GetParams(const Int_t Ifl, Double_t *RetParams) const{
  // The matrix version is the master
  TMatrixD* Dumm = new TMatrixD(InpEst,InpEst);
  Int_t IRet = GetParams(Ifl, Dumm);
  MatrixtoDouble(Dumm, RetParams);
  Dumm->Delete(); Dumm = NULL;
  return IRet;
};

//------------------------------------------------------------------------------
// After Solve()

Double_t Blue::GetChiq()const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetChiq: Return the Chi-square of the result \n");
    }
    return ChiQua;
  }
  printf("... Blue->GetChiq: Presently not available, call Solve() \n");
  return 0.;
};

//------------------------------------------------------------------------------

Int_t Blue::GetNdof()const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetNdof: Return the number of degrees of");
      printf(" freedom of the result \n");
    }
    return NumDof;
  }
  printf("... Blue->GetNdof: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Double_t Blue::GetProb()const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetProb: Return the chi-square probability");
      printf(" of the result \n");
    }
    return ChiPro;
  }
  printf("... Blue->GetProb: Presently not available, call Solve() \n");
  return 0.;
};

//------------------------------------------------------------------------------

Double_t Blue::GetPull(const Int_t i)const{
  if(IsSolved() == 1){
    if(IsActiveEst(i) == 1){
      Int_t jj = 0;   
      for(Int_t j = 0; j < i; j++){
	jj = jj + IsActiveEst(j);
      }
      if(IsPrintLevel() >= 1){
	printf("... Blue->GetPull(%2i): Return the pull of this estimate = %5.3f \n",
	       i, Pul->operator()(jj));
      }
      return Pul->operator()(jj);
    }else{
      printf("... Blue->GetPull(%2i): Not an active estimate \n", i);
      return 0.;
    }
  }
  printf("... Blue->GetPull(%2i): Presently not available, call Solve() \n", i);
  return 0.;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCovRes(TMatrixD *UseCovRes)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCovRes: Return covariance matrix of the result as TMatrixD \n");
    }
    UseCovRes->SetSub(0,0,*CovRes);
    return 1;
  }
  printf("... Blue->GetCovRes: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetCovRes(Double_t *RetCovRes)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetCovRes: Return covariance matrix of the result as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpObs,InpObs);
    Int_t IRet = GetCovRes(Dumm);
    MatrixtoDouble(Dumm, RetCovRes);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetCovRes: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetRhoRes(TMatrixD *UseRhoRes)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetRhoRes: Return correlation matrix of the result as TMatrixD \n");
    }
    UseRhoRes->SetSub(0,0,*RhoRes);
    return 1;
  }
  printf("... Blue->GetRhoRes: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetRhoRes(Double_t *RetRhoRes)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetRhoRes: Return correlation matrix of the result as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpObs,InpObs);
    Int_t IRet = GetRhoRes(Dumm);
    MatrixtoDouble(Dumm, RetRhoRes);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetRhoRes: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetWeight(TMatrixD *UseWeight)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetWeight: Return weight matrix as TMatrixD \n");
    }
    UseWeight->SetSub(0,0,*Lam);
    return 1;
  }
  printf("... Blue->GetWeight: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetWeight(Double_t *RetWeight)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetWeight: Return weight matrix as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpEst,InpObs);
    Int_t IRet = GetWeight(Dumm);
    MatrixtoDouble(Dumm, RetWeight);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetWeight: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetResult(TMatrixD *UseResult)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetResult: Return matrix of the result as TMatrixD \n");  
    }
    Int_t m = 0;
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(IsActiveObs(n) == 1){
	UseResult->operator()(m,0) = XvaRes->operator()(m);
	for(Int_t k = 0; k<InpUnc; k++){
	  UseResult->operator()(m,k+1) = CorRes->operator()(k*InpObs+m,k*InpObs+m);
	}
	m = m + 1;
      }
    }
    return 1;
  }
  printf("... Blue->GetResult: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetResult(Double_t *RetResult)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetResult:Return matrix of the result as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpObs,InpUnc+1);
    Int_t IRet = GetResult(Dumm);
    MatrixtoDouble(Dumm, RetResult);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetResult: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetUncert(TMatrixD *UseUncert)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetUncert: Return total uncertainties of the results as TMatrixD \n");  
    }
    Int_t m = 0;
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(IsActiveObs(n) == 1){
	UseUncert->operator()(m,0) = sqrt(CovRes->operator()(m,m));
	m = m + 1;
      }
    }
    return 1;
  }
  printf("... Blue->GetUncert: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetUncert(Double_t *RetUncert)const{
  if(IsSolved() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetUncert:Return total uncertainties of the results as Double_t array \n");
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpObs,1);
    Int_t IRet = GetUncert(Dumm);
    MatrixtoDouble(Dumm, RetUncert);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetUncert: Presently not available, call Solve() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetInspectLike(TMatrixD *UseInsLik)const{
  if(IsInspectedLike() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetInspectLike: Return the results from InspectLike() as TMatrixD \n");  
    }
    UseInsLik->SetSub(0,0,*LikRes);
    return 1;
  }
  printf("... Blue->GetInspectLike: Presently not available, call InspectLike() \n");
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetInspectLike(Double_t *RetInsLik)const{
  if(IsInspectedLike() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetInspectLike: Return the results from InspectLike() Double_t array \n");  
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpObsOrig,LikDim);
    Int_t IRet = GetInspectLike(Dumm);
    MatrixtoDouble(Dumm, RetInsLik);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetInspectLike: Presently not available, call InspectLike() \n");
  return 0;
};

//------------------------------------------------------------------------------
// After SolveXXX()

Int_t Blue::GetAccImpLasEst(const Int_t n) const{
  Int_t IndLas = -1;
  if(IsSolvedAccImp() == 1){
    if(IsActiveObs(n) == 1){
      IndLas = static_cast<int>(LasImp->operator()(0,IsIndexObs(n)));
      if(IsPrintLevel() >= 1){
	printf("... Blue->GetAccImpLasEst(%2i): The last estimate to be used is %2i", n, IndLas);
      }
    }else{
      printf("... Blue->GetAccImpLasEst(%2i): Not an active observable  \n", n);
    }
  }else{
    printf("... Blue->GetAccImpLasEst(%2i): Presently not available, call SolveAccImp() \n", n);
  }
  return IndLas;
};
  
//------------------------------------------------------------------------------

Int_t Blue::GetAccImpIndEst(const Int_t n, Int_t *IndEst) const{
  if(IsSolvedAccImp() == 1){
    if(IsActiveObs(n) == 1){
      Int_t na = IsIndexObs(n);
      Int_t ia = 0;
      for(Int_t i = 0; i<InpEst; i++){	
	ia = static_cast<int>(IndImp->operator()(i,na));
	if(ia > -0.5)IndEst[i] = ia;
	if(IsPrintLevel() >= 1){
	  if(i == 0){
	    printf("... Blue->GetAccImpIndEst(%2i): List of importance = %2i", n, IndEst[i]);
	  }else{
	    printf(", %2i", IndEst[i]);
	  }
	  if(ia < 0)printf("\n");
	}
      }
      return 1;
    }else{
      printf("... Blue->GetAccImpIndEst(%2i): Not an active observable  \n", n);
    }
  }else{
    printf("... Blue->GetAccImpIndEst(%2i): Presently not available, call SolveAccImp() \n", n);
  }
  return 0;
};
  
//------------------------------------------------------------------------------

Int_t Blue::GetNumScaFac()const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetNumScaFac(): Return the number of groups for SolveScaRho() \n");
    }
    return InpFac;
  }else{    
    printf("... Blue->GetNumScaFac(): Presently not available, call FixInp() \n");
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::GetNumScaRho()const{
  if(IsFixedInp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetNumScaRho(): Return the number of correlation values per groups for SolveScaRho() \n");
    }
    return InrFac;
  }else{    
    printf("... Blue->GetNumScaRho(): Presently not available, call FixInp() \n");
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::GetScaVal(const Int_t n, TMatrixD *UseScaVal)const{
  if(IsSolvedScaRho() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetScaVal(%2i): Return the values from SolveScaRho() as TMatrixD \n", n);
    }
    if(IsActiveObs(n) == 1){
      UseScaVal->SetSub(0,0,*ValSca[IsIndexObs(n)]);
      return 1;
    }else{
      printf("... Blue->GetScaVal(%2i): %2i is not an active observable \n", n, n);
      return 0;
    }
  }
  printf("... Blue->GetScaVal(%2i): Presently not available, call SolveScaRho() \n", n);
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetScaVal(const Int_t n, Double_t *RetScaVal)const{
  if(IsSolvedScaRho() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetScaVal(%2i): Return the values from SolveScaRho() as Double_t array \n", n);
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpUnc*InpFac,InrFac);
    Int_t IRet = GetScaVal(n,Dumm);
    MatrixtoDouble(Dumm, RetScaVal);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetScaVal(%2i): Presently not available, call SolveScaRho() \n", n);
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetScaUnc(const Int_t n, TMatrixD *UseScaUnc)const{
  if(IsSolvedScaRho() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetScaUnc(%2i): Return the values from SolveScaRho() as TMatrixD \n", n);
    }
    if(IsActiveObs(n) == 1){
      UseScaUnc->SetSub(0,0,*SigSca[IsIndexObs(n)]);
      return 1;
    }else{
      printf("... Blue->GetScaUnc(%2i): %2i is not an active observable \n", n, n);
      return 0;
    }
  }
  printf("... Blue->GetScaUnc(%2i): Presently not available, call SolveScaRho() \n", n);
  return 0;
};

//------------------------------------------------------------------------------

Int_t Blue::GetScaUnc(const Int_t n, Double_t *RetScaUnc)const{
  if(IsSolvedScaRho() == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->GetScaUnc(%2i): Return the values from SolveScaRho() as Double_t array \n", n);
    }
    // The matrix version is the master
    TMatrixD* Dumm = new TMatrixD(InpUnc*InpFac,InrFac);
    Int_t IRet = GetScaUnc(n,Dumm);
    MatrixtoDouble(Dumm, RetScaUnc);
    Dumm->Delete(); Dumm = NULL;
    return IRet;
  }
  printf("... Blue->GetScaUnc(%2i): Presently not available, call SolveScaRho() \n", n);
  return 0;
};

//------------------------------------------------------------------------------
// Setters
//------------------------------------------------------------------------------

void Blue::SetActiveEst(const Int_t i){
  if(IsFilledInp() == 0){
    printf("... Blue->SetActiveEst(%2i): IGNORED input not yet filled \n",i);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetActiveEst(%2i): IGNORED input not yet released.",i);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedEst(i) == 1){
    printf("... Blue->SetActiveEst(%2i):   Activate estimate \n", i);
    SetActiveEst(i, 1);
  }
};

//------------------------------------------------------------------------------

void Blue::SetActiveUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetActiveUnc(%2i): IGNORED input not yet filled \n",k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetActiveUnc(%2i): IGNORED input not yet released.",k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedUnc(k) == 1){
    printf("... Blue->SetActiveUnc(%2i):   Activate uncertainty \n", k);
    SetActiveUnc(k, 1);
  }
};

//------------------------------------------------------------------------------

void Blue::SetInActiveEst(const Int_t i){
  if(IsFilledInp() == 0){
    printf("... Blue->SetInActiveEst(%2i): IGNORED input not yet filled \n", i);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetInActiveEst(%2i): IGNORED input not yet released.",i);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedEst(i) == 1){
    printf("... Blue->SetInActiveEst(%2i): Deactivate estimate \n", i);
    SetActiveEst(i, 0);
  }
};

//------------------------------------------------------------------------------

void Blue::SetInActiveUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetInActiveUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetInActiveUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedUnc(k) == 1){
    printf("... Blue->SetActiveUnc(%2i): Deactivate uncertainty \n", k);
    SetActiveUnc(k, 0);
  }
};

//------------------------------------------------------------------------------

void Blue::SetRhoValUnc(const Double_t RhoVal){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoValUnc(%5.3f): I set the corrrelation for all active uncertainties \n", RhoVal);
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetRhoValUnc(k, RhoVal);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoValUnc(const Int_t k, const Double_t RhoVal){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoValUnc(%2i,%5.3f): I set the corrrelation for this uncertainty \n", k, RhoVal);
  }
  for(Int_t l = 0; l<InpFac; l++){SetRhoValUnc(k, l, RhoVal);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoValUnc(const Int_t k, const Int_t l, const Double_t RhoVal){
  if(IsFilledInp() == 0){
    printf("... Blue->SetRhoValUnc(%2i,%2i,%5.3f): IGNORED input not yet filled \n", k, l, RhoVal);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetRhoValUnc(%2i,%2i,%5.3f): IGNORED input not yet released.", k, l, RhoVal);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoValUnc(%2i,%2i,%5.3f): I set the corrrelation for a subset of this uncertainty \n",
	   k, l, RhoVal);
  }
  SetRhoValUnc(k, l, RhoVal, 1);
  return;
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoValUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRhoValUnc(): I revert to the initially provided corrrelation for all active uncertainties \n");
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetNotRhoValUnc(k);};
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoValUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetNotRhoValUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetNotRhoValUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRhoValUnc(%2i): I revert to the initially provided corrrelation for this uncertainty \n", k);
  }
  for(Int_t l = 0; l<InpFac; l++){SetRhoValUnc(k, l, 1.0, 0);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoFacUnc(const Double_t RhoFac){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoFacUnc(%5.3f): I scale the corrrelation for all active uncertainties \n", RhoFac);
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetRhoFacUnc(k, RhoFac);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoFacUnc(const Int_t k, const Double_t RhoFac){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoFacUnc(%2i,%5.3f): I scale the corrrelation for this uncertainty \n", k, RhoFac);
  }
  for(Int_t l = 0; l<InpFac; l++){SetRhoFacUnc(k, l, RhoFac);};
  return;
};

//------------------------------------------------------------------------------

void Blue::SetRhoFacUnc(const Int_t k, const Int_t l, const Double_t RhoFac){
  if(IsFilledInp() == 0){
    printf("... Blue->SetRhoFacUnc(%2i,%2i,%5.3f): IGNORED input not yet filled \n", k, l, RhoFac);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetRhoFacUnc(%2i,%2i,%5.3f): IGNORED input not yet released.", k, l, RhoFac);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoFacUnc(%2i,%2i,%5.3f): I scale the corrrelation for a subset of this uncertainty \n", 
	   k, l, RhoFac);
  }
  SetRhoFacUnc(k, l, RhoFac, 1);
  return;
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoFacUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRhoFacUnc(): I revert to the initially provided corrrelation for all uncertainties \n");
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetNotRhoFacUnc(k);};
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoFacUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetNotRhoFacUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetNotRhoFacUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRhoFacUnc(%2i): I revert to the initially provided corrrelation for this uncertainty \n", k);
  }
  for(Int_t l = 0; l<InpFac; l++){SetRhoFacUnc(k, l, 1.0, 0);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoRedUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRhoRedUnc(): I use reduced correlations for all active uncertainties \n");
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetRhoRedUnc(k);};
};

//------------------------------------------------------------------------------

void Blue::SetRhoRedUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetRhoRedUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetRhoRedUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedUnc(k) == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->SetRhoRedUnc(%2i): I use reduced correlations for this uncertainty \n", k);
    }
    SetRhoRedUnc(k, 1);
  }
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoRedUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRhoRedUnc(): I revert to the initially provided corrrelation for all uncertainties \n");
  }
  for(Int_t k = 0; k<InpUncOrig; k++){if(IsActiveUnc(k)==1)SetNotRhoRedUnc(k);};
};

//------------------------------------------------------------------------------

void Blue::SetNotRhoRedUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetNotRhoRedUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetNotRhoRedUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedUnc(k) == 1){
    if(IsPrintLevel() >= 1){
      printf("... Blue->SetNotRhoRedUnc(%2i): I revert to the initially provided corrrelation for this uncertainty \n", k);
    }
    SetRhoRedUnc(k, 0);
  }
};

//------------------------------------------------------------------------------

void Blue::SetRelUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetRelUnc(): I treat all uncertainties as relative \n");
  }
  for(Int_t k = 0; k<InpUnc; k++){
    SetRelUnc(k);
  }
};

//------------------------------------------------------------------------------

void Blue::SetRelUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetRelUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetRelUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Do it for uncertainty source k and all estimates
  Double_t ActCof[MaxCof] = {0};
  Double_t Value = 0, Uncer = 0;
  if(IsAllowedUnc(k) == 1){
    if(IsActiveUnc(k) == 1){
      if(IsPrintLevel() >= 1){
	printf("... Blue->SetRelUnc(%2i): I treat this source as a relative uncertainty \n", k);
      }
      for(Int_t i = 0; i<InpEstOrig; i++){	
	if(IsActiveEst(i) == 1){
	  Value = XvaOrig->operator()(i);
	  Uncer = UncOrig->operator()(k*InpEstOrig+i, k*InpEstOrig+i);
	  if(k == 0){
	    ActCof[1] = (Uncer*Uncer) / (TMath::Abs(Value));
	    ActCof[2] = 0;
	  }else{
	    ActCof[1] = 0;
	    ActCof[2] = (Uncer*Uncer) / (Value*Value);
	  }
	  //printf("... Blue->SetRelUnc(%2i): i = %2i Value = %5.3f Uncertainty = %5.3f \n", k, i, Value,Uncer);
	  //printf("... Blue->SetRelUnc(): ActCof = %5.3f %5.3f  %5.3f \n", ActCof[0],ActCof[1],ActCof[2]);
	  SetRelUnc(i, k, &ActCof[0]);	       
	}else{
	  printf("... Blue->SetRelUnc(%2i): IGNORED Not an active estimate i = %2i \n", k, i);
	}
      }
    }else{
      printf("... Blue->SetRelUnc(%2i): IGNORED Not an active uncertainty k = %2i \n", k, k);
    }
  }else{
    printf("... Blue->SetRelUnc(%2i): IGNORED Not an allowed uncertainty k = %2i \n", k, k);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::SetRelUnc(const Int_t i, const Int_t k, const Double_t *ActCof){
  if(IsFilledInp() == 0){
    printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED input not yet filled \n",
	   i,k,ActCof[0],ActCof[1],ActCof[2]);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED input not yet released.", 
	   i,k,ActCof[0],ActCof[1],ActCof[2]);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Check for allowed and active uncertainties and estimates
  if(IsAllowedUnc(k) == 1){
    if(IsActiveUnc(k) == 1){
      if(IsAllowedEst(i) == 1){
	if(IsActiveEst(i) == 1){
	  if(IsPrintLevel() >= 1){
	    printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): I give estimate %2i, source %2i a function \n",
		   i,k,ActCof[0],ActCof[1],ActCof[2],i,k);
	  }
	  for(Int_t l = 0; l<MaxCof; l++){Cof->operator()(i*InpUncOrig+k,l) = ActCof[l];};
	  IndRel->operator()(i,k) = 1.;
	  UncRel[k] = -1;
	  Double_t AllSet = 1.;
	  for(Int_t j = 0;j<InpEstOrig; j++){
	    if(IsActiveEst(j) == 1)AllSet = AllSet * IndRel->operator()(j,k);
	  }
	  if(AllSet > 0.9)UncRel[k] = 1;
	}else{
	  printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED Not an active estimate %2i \n",
		 i,k,ActCof[0],ActCof[1],ActCof[2],i);
	} 
      }else{
	printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED Not an allowed estimate %2i \n",
	       i,k,ActCof[0],ActCof[1],ActCof[2],i);
      }
    }else{
      printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED Not an active uncertainty %2i \n",
	     i,k,ActCof[0],ActCof[1],ActCof[2],k);
    }
  }else{
    printf("... Blue->SetRelUnc(%2i,%2i,%5.3f,%5.3f,%5.3f): IGNORED Not an allowed uncertainty %2i \n",
	   i,k,ActCof[0],ActCof[1],ActCof[2],k);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::SetNotRelUnc(){
  if(IsPrintLevel() >= 1){
    printf("... Blue->SetNotRelUnc(): I treat all uncertainties as absolute \n");
  }
  for(Int_t k = 0; k<InpUnc; k++){
    SetNotRelUnc(k);
  }
};

//------------------------------------------------------------------------------

void Blue::SetNotRelUnc(const Int_t k){
  if(IsFilledInp() == 0){
    printf("... Blue->SetNotRelUnc(%2i): IGNORED input not yet filled \n", k);
    return;
  }else if(IsFixedInp() == 1){
    printf("... Blue->SetNotRelUnc(%2i): IGNORED input not yet released.", k);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }
  if(IsAllowedUnc(k) == 1){
    if(IsActiveUnc(k) == 1){
      if(IsPrintLevel() >= 1){
	printf("... Blue->SetNotRelUnc(%2i): I treat this source as an absolute uncertainty \n", k);
      }
      UncRel[k] = 0;
      for(Int_t i = 0; i<InpEstOrig; i++){
	IndRel->operator()(i,k) = 0.;
      }
    }else{
      printf("... Blue->SetNotRelUnc(%2i): IGNORED Not an active uncertainty \n", k);
      return;
    }
  }else{
    printf("... Blue->SetNotRelUnc(%2i): IGNORED Not an allowed uncertainty \n", k);
    return;
  }
};

//---------------------------------------------------------------------------
// Control Printout
//---------------------------------------------------------------------------

void Blue::SetPrintLevel(const Int_t p){
  if(p >= 0 || p <= 2){
    IPrint = p;
    if(p > 0){
      printf("... Blue->SetPrintLevel(%2i): Print level is set to %2i \n", p, p);
    }
  }else{
    IPrint = 0;
    printf("... Blue->SetPrintLevel(%2i): Out of range set to 0 \n", p);
  }
};

//------------------------------------------------------------------------------

void Blue::SetQuiet(){
  SetQuiet(1);
};

//------------------------------------------------------------------------------

void Blue::SetNotQuiet(){
  SetQuiet(0);
};

//------------------------------------------------------------------------------
// Print Out
//------------------------------------------------------------------------------

void Blue::PrintMatrix(const TMatrixD *TryMat) const {
  TString ForVal = "%5.2f";
  Int_t NumRow = -1;
  Int_t NumCol = -1;
  PrintMatrix(TryMat,NumRow,NumCol,ForVal);
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintMatrix(const TMatrixD *TryMat, const TString ForVal) const {
  Int_t NumRow = -1;
  Int_t NumCol = -1;
  PrintMatrix(TryMat,NumRow,NumCol,ForVal);
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintMatrix(const TMatrixD *TryMat, const Int_t NumRow, const Int_t NumCol) const {
  TString ForVal = "%5.2f";
  PrintMatrix(TryMat,NumRow,NumCol,ForVal);
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintMatrix(const TMatrixD *TryMat, const Int_t NumRow, const Int_t NumCol,
		       const TString ForVal) const {

  // The format string
  TString Format = "To be filled later";

  // Get the dimensions 
  Int_t NRows = 0, NCols = 0;
  if(NumRow == -1){
    NRows = TryMat->GetNrows();
    NCols = TryMat->GetNcols();
  }else{
    NRows = TMath::Min(NumRow, TryMat->GetNrows());
    NCols = TMath::Min(NumCol, TryMat->GetNcols());
  }

  // Do the print out
  for(Int_t i = 0; i<NRows; i++){
    Format = "... Blue->PrintMatrix(): %3i: " + ForVal;
    printf(Format, i, TryMat->operator()(i,0));
    Format = " " + ForVal;
    for(Int_t j = 1; j<NCols; j++){
      printf(Format, TryMat->operator()(i,j));
    }
    printf("\n");
  } 
  printf("\n");
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintDouble(const Double_t *TryDou, const Int_t NumRow, const Int_t NumCol) const {
  TString ForVal = "%5.2f";
  PrintDouble(&TryDou[0], NumRow, NumCol, ForVal);
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintDouble(const Double_t *TryDou, const Int_t NumRow, const Int_t NumCol,
		       const TString ForVal) const {
  // Transform into TMatrix then use PrintMatrix
  TMatrixD *TryMat = new TMatrixD(NumRow,NumCol);
  TryMat->Use(0,NumRow-1,0,NumCol-1,&TryDou[0]);
  if(IsPrintLevel() >= 1){
    printf("... Blue->PrintDouble(): I will let PrintMatrix() do the job \n");
  }
  PrintMatrix(TryMat,NumRow,NumCol,ForVal);
  TryMat->Delete(); TryMat = NULL;
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintListEst() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintListEst(): The present list of     estimates is: ");
    for(Int_t i = 0; i<InpEst-1; i++){
      printf("%2i,", LisEst[i]);
    }
    printf("%2i \n", LisEst[InpEst-1]);
  }else{
    printf("... Blue->PrintListEst(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintListUnc() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintListUnc(): The present list of uncertainties is: ");
    for(Int_t k = 0; k<InpUnc-1; k++){
      printf("%2i,", LisUnc[k]);
    }
    printf("%2i \n", LisUnc[InpUnc-1]);
  }else{
    printf("... Blue->PrintListUnc(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintNamEst() const {
  if(IsFixedInp() == 1){
    Int_t ICount = 0;
    printf("... Blue->PrintNamEst(): The names of the present list of estimates are: \n");
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i)){
	ICount = ICount + 1;
	if(ICount == 1){
	  printf("... Blue->PrintNamEst(): %2i = %s", i, GetNamEst(i).Data());
	}else{
	  printf(", %2i = %s", i, GetNamEst(i).Data());
	}
	if(ICount == 5){
	  printf("\n");
	  ICount = 0;
	}
      }
    }
    if(ICount != 0)printf("\n");
  }else{
    printf("... Blue->PrintNamEst(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintNamUnc() const {
  if(IsFixedInp() == 1){
    Int_t ICount = 0;
    printf("... Blue->PrintNamUnc(): The names of the present list of uncertainties are: \n");
    for(Int_t k = 0; k<InpUncOrig; k++){
      if(IsActiveUnc(k)){
	ICount = ICount + 1;
	if(ICount == 1){
	  printf("... Blue->PrintNamUnc(): %2i = %s", k, GetNamUnc(k).Data());
	}else{
	  printf(", %2i = %s", k, GetNamUnc(k).Data());
	}
	if(ICount == 5){
	  printf("\n");
	  ICount = 0;
	}
      }
    }
    if(ICount != 0)printf("\n");
  }else{
    printf("... Blue->PrintNamUnc(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintEst(const Int_t i) const {
  if(IsFixedInp() == 1){
    if(IsAllowedEst(i) == 1){
      if(IsActiveEst(i) == 1){
	// Find index in active list
	Int_t j = -1;
	for(Int_t jj = 0; jj<InpEst; jj++){
	  if(j == -1 && IsWhichEst(jj) == i)j = jj;
	}
	printf("... Blue->PrintEst(%2i) %s: ",i,GetNamEst(i).Data());
	printf(" %5.3f",Xva->operator()(j));
	printf(" +- %5.3f",Sig->operator()(j));
	if(InpUnc > 1){
	  for(Int_t k = 0; k<InpUnc; k++){
	    if(k == 0){
	      printf(" (%5.3f",Unc->operator()(j+k*InpEst,j+k*InpEst));
	    }else{
	      printf(" +- %5.3f",Unc->operator()(j+k*InpEst,j+k*InpEst));
	    }
	    if(k == InpUnc-1)printf(")");
	  }
	}
	printf("\n");
      }else{
	printf("... Blue->PrintEst(%2i): Not an active estimate \n",i);
      }    
    }else{
      printf("... Blue->PrintEst(%2i):  Not an allowed estimate \n",i);
    }
  }else{
    printf("... Blue->PrintEst(%2i): Presently not available, call FixInp() \n",i);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintEst() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintEst:                Value");
    if(InpUnc != 1)printf("     Full");
    for(Int_t k = 0; k<InpUnc; k++){
      printf("%s  ",GetNamUnc(IsWhichUnc(k)).Data());
    }
    printf("\n");

    for(Int_t i = 0; i<InpEstOrig; i++){
      PrintEst(i);
    }
  }else{
    printf("... Blue->PrintEst: Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCofRelUnc(const Int_t k) const {
  if(IsFixedInp() == 1){
    if(IsAllowedUnc(k) == 1){
      if(IsActiveUnc(k) == 1){
	if(IsRelValUnc(k) == 1){
	  printf("... Blue->PrintCofRelUnc(%2i):  Coefficients to calculate the variance  \n",k);
	  printf("... Blue->PrintCofRelUnc(%2i): --------- For uncertainty = %2i -------- \n",k,k);
	  printf("... Blue->PrintCofRelUnc(%2i): Est          a0          a1          a2 \n",k);
	  for(Int_t i = 0; i<InpEstOrig; i++){
	    if(IsActiveEst(i) == 1 && IsRelValUnc(k) == 1){
	      Int_t IRow = i*InpUncOrig + k;
	      printf("... Blue->PrintCofRelUnc(%2i): %3i    %4.2e    %4.2e    %4.2e \n",k,i,
		     Cof->operator()(IRow,0),Cof->operator()(IRow,1),
		     Cof->operator()(IRow,2));
	    }
	  }
	  printf("\n");
	}else{
	  printf("... Blue->PrintCofRelUnc(%2i): Not a relative uncertainty \n",k);
	}
      }else{
	printf("... Blue->PrintCofRelUnc(%2i): Not an active uncertainty \n",k);
      }
    }else{
      printf("... Blue->PrintCofRelUnc(%2i): Not an allowed uncertainty \n",k);
    }
  }else{
    printf("... Blue->PrintCofRelUnc(%2i): Presently not available, call FixInp() \n",k);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCofRelUnc() const {
  if(IsFixedInp() == 1){
    //printf("... Blue->PrintCofRelUnc(): Cof \n");Cof->Print();
    for(Int_t k = 0; k<InpUncOrig; k++){
      if(IsActiveUnc(k) == 1){
	PrintCofRelUnc(k);
      }
    }
  }else{
    printf("... Blue->PrintCofRelUnc(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCor(const Int_t k) const {
  if(IsFixedInp() == 1){
    if(IsAllowedUnc(k) == 1){
      if(IsActiveUnc(k) == 1){
	// Find index in active list
	Int_t l = -1;
	for(Int_t ll = 0; ll<InpUncOrig; ll++){
	  if(l == -1 && IsWhichUnc(ll) == k)l = ll;
	}
	printf("... Blue->PrintCor(%2i): Correlation matrix for: %s \n",k,GetNamUnc(k).Data());
	TMatrixD *F = new TMatrixD(InpEst,InpEst);
	Int_t IRowLow = l*InpEst; 
	Int_t IRowHig = IRowLow + InpEst-1;
	Int_t IColLow = IRowLow; 
	Int_t IColHig = IRowHig;
	//  printf("... %3i, %3i, %3i, %3i \n",IRowLow, IColLow, IRowHig, IColHig);
	Cor->GetSub(IRowLow, IRowHig, IColLow, IColHig, *F, "S");
	F->Print();
	F->Delete(); F = NULL;
      }else{
	printf("... Blue->PrintCor(%2i): Not an active uncertainty \n", k);
      }    
    }else{
      printf("... Blue->PrintCor(%2i): Presently not available, call FixInp() \n",k);
    }
  }else{
    printf("... Blue->PrintCor(%2i):  Not an allowed uncertainty \n", k);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCor() const {
  if(IsFixedInp() == 1){
    for(Int_t k = 0; k<InpUncOrig; k++){
      PrintCor(k);
    }
  }else{
    printf("... Blue->PrintCor: Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCov(const Int_t k) const {
  if(IsFixedInp() == 1){
    if(IsAllowedUnc(k) == 1){
      if(IsActiveUnc(k) == 1){
	// Find index in active list
	Int_t l = -1;
	for(Int_t ll = 0; ll<InpUncOrig; ll++){
	  if(l == -1 && IsWhichUnc(ll) == k)l = ll;
	}
	printf("... Blue->PrintCov(%2i): Covariance matrix for: %s \n",k,GetNamUnc(k).Data());
	TMatrixD *F = new TMatrixD(InpEst,InpEst);
	TMatrixD *G = new TMatrixD(InpEst,InpEst);
	TMatrixD *H = new TMatrixD(InpEst,InpEst);
	Int_t IRowLow = l*InpEst; 
	Int_t IRowHig = IRowLow + InpEst-1;
	Int_t IColLow = IRowLow; 
	Int_t IColHig = IRowHig;
	//  printf("... %3i, %3i, %3i, %3i \n",IRowLow, IColLow, IRowHig, IColHig);
	// The correlation matrix
	Cor->GetSub(IRowLow, IRowHig, IColLow, IColHig, *F, "S");
	// The uncertainty matrix
	Unc->GetSub(IRowLow, IRowHig, IColLow, IColHig, *G, "S");
	//Calculate the covariance for this source
	Unc->GetSub(IRowLow, IRowHig, IColLow, IColHig, *H, "S");
	H->operator*=(*F);
	H->operator*=(*G);
	H->Print();
	// Clean up
	F->Delete(); F = NULL;
	G->Delete(); G = NULL;
	H->Delete(); H = NULL;
      }else{
	printf("... Blue->PrintCov(%2i): Not an active uncertainty \n", k);
      }    
    }else{
      printf("... Blue->PrintCov(%2i): Presently not available, call FixInp() \n",k);
    }
  }else{
    printf("... Blue->PrintCov(%2i):  Not an allowed uncertainty \n", k);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCov() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintCov: The covariance matrix of the input \n");
    Cov->Print();
  }else{
    printf("... Blue->PrintCov: Presently not available, call FixInp() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintCovInvert() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintCovInvert: The inverse covariance matrix of the input \n");
    CovI->Print();
  }else{
    printf("... Blue->PrintCovInvert: Presently not available, call FixInp() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintRho() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintRho: The correlation matrix of the input \n");
    Rho->Print();
  }else{
    printf("... Blue->PrintRho: Presently not available, call FixInp() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintCompatEst(const TString FilNam) const {

  // Return if not fixed
  if(IsFixedInp() != 1){
    printf("... Blue->PrintCompatEst(): Presently not available, call FixInp() \n");
    return;
  }

  // Check number of estimates, return if not sensible
  if(GetActEst() == 1){
    printf("... Blue->PrintCompatEst(): There is only one active estimate \n");
    return;
  }else if(GetActEst() == GetActObs()){
    printf("... Blue->PrintCompatEst(): There are as many active estimates as observables \n");
    return;
  }
  
  // List the estimates
  if(IsPrintLevel() >= 1){
    printf("... Blue->PrintCompatEst(): The follwing estimates were used \n");
    for(Int_t i = 0; i<InpEstOrig; i++)if(IsActiveEst(i) == 1)PrintEst(i);
  }
  
  // Define maximum chisquare for histogram
  Double_t ChiMax = 5;
  
  // Find maximum number of estimates of any observable
  Int_t MaxNum = 0;
  Int_t ActNum = 0;
  for(Int_t n = 0; n<InpObsOrig; n++){
    if(IsActiveObs(n)){
      ActNum = 0;
      for(Int_t i = 0; i<InpEst; i++){
	if(EstWhichObs(IsWhichEst(i)) == n)ActNum = ActNum + 1;
      }
      if(ActNum > MaxNum)MaxNum = ActNum;
    }
  }
  MaxNum = 0.5*(MaxNum-1)*MaxNum;
  //printf("... Blue->PrintCompatEst(): MaxNum = %6i \n",MaxNum);

  // Book Matricees to store the chisquare and probability values
  TMatrixD *C = new TMatrixD(InpObs,MaxNum);
  TMatrixD *P = new TMatrixD(InpObs,MaxNum);
  
  // Check compatibility of the estimates per observable
  //Xva->Print();  Rho->Print(); Cov->Print();
  printf("... Blue->PrintCompatEst(): The compatibility of the estimates for Ndof = 1 \n");
  Int_t iok[4] = {0};
  Int_t npa = 0;
  Double_t rho = 0., vao = 0., vat = 0., vac = 0., chi = 0.;
  Int_t na = 0, ia = 0;
  for(Int_t n = 0; n<InpObsOrig; n++){
    if(IsActiveObs(n)){
      printf("... Blue->PrintCompatEst(): The next observable is: %s \n", 
	     GetNamObs(n).Data());
      npa = 0;
      ia  = 0;
      for(Int_t i = 0; i<InpEst-1; i++){
	for(Int_t j = i+1; j<InpEst; j++){
	  if(EstWhichObs(IsWhichEst(i)) == EstWhichObs(IsWhichEst(j)) && 
	     EstWhichObs(IsWhichEst(i)) == n){
	    npa = 1;
	    rho = Rho->operator()(i,j);
	    vao = Cov->operator()(i,i);
	    vat = Cov->operator()(j,j);
	    vac = vao + vat - 2 * rho * TMath::Sqrt(vao) * TMath::Sqrt(vat);
	    chi = (Xva->operator()(i) - Xva->operator()(j)) *
	      (Xva->operator()(i) - Xva->operator()(j)) / vac;
	    if(chi >= 1)iok[0] = iok[0] + 1;
	    if(chi >= 2)iok[1] = iok[1] + 1;
	    if(chi >= 3)iok[2] = iok[2] + 1;
	    if(chi >= 5)iok[3] = iok[3] + 1;
	    //printf("... Blue->PrintCompatEst(): i = %2i , j = %2i, rho= %5.3f, sig1 = %5.3f,  sig2 = %5.3f \n", 
	    //i, j, rho, TMath::Sqrt(vao), TMath::Sqrt(vat));
	    printf("... Blue->PrintCompatEst(): (%s,%s) Sigma(Delta) = %5.3f,", 
		   GetNamEst(IsWhichEst(i)).Data(), 
		   GetNamEst(IsWhichEst(j)).Data(), TMath::Sqrt(vac));
	    printf(" Chiq = %5.3f, Prob = %5.3f \n", chi, TMath::Prob(chi,1));
	    C->operator()(na,ia) = TMath::Min(chi,ChiMax);
	    P->operator()(na,ia) = TMath::Prob(chi,1);
	    ia = ia + 1;
	  }
	}
      }
      if(npa == 0)printf("... Blue->PrintCompatEst(): There is only one estimate for this observable \n");
      na = na + 1;
      printf("\n");
    }
  }
  
  // Report large chisquared if spotted
  if(iok[3] > 0)printf("... Blue->PrintCompatEst(): Pairs of estimates with Chiq > 5 = %2i \n", iok[3]);
  if(iok[2] > 0)printf("... Blue->PrintCompatEst(): Pairs of estimates with Chiq > 3 = %2i \n", iok[2]);
  if(iok[1] > 0)printf("... Blue->PrintCompatEst(): Pairs of estimates with Chiq > 2 = %2i \n", iok[1]);
  if(iok[0] > 0){
    printf("... Blue->PrintCompatEst(): Pairs of estimates with Chiq > 1 = %2i \n", iok[0]);
  }else{
    printf("... Blue->PrintCompatEst(): The compatibility of the estimates is good");
    printf(" (all differences agree within 1 sigma) \n");
  }
  
  // Clean up and return if no plots wanted
  if(FilNam == "NoGraphWanted"){
    C->Delete(); C = NULL;
    P->Delete(); P = NULL;
    return;
  }

  // Report matrices
  //printf("... Blue->PrintCompatEst(): C: \n"); C->Print();
  //printf("... Blue->PrintCompatEst(): P: \n"); P->Print();
  
  // Set pointers and book histos
  Int_t NumBin = 5;
  TH1F* HisDel;
  TH1F* HisChi[GetActObs()];
  TH1F* HisPro[GetActObs()];
  char DumName[50];
  na = 0;
  for(Int_t n = 0; n<InpObsOrig; n++){
    if(IsActiveObs(n)){
      // The chi squared
      sprintf(DumName,"ChiQua Obs%i",n);
      HisDel = (TH1F*) gROOT->FindObject(DumName);
      if(HisDel)delete HisDel;
      HisChi[na] = new TH1F(DumName,GetNamObs(n).Data(),NumBin,0.,ChiMax);
      HisChi[na]->GetXaxis()->SetTitleOffset(1.0);
      HisChi[na]->GetXaxis()->SetTitleSize(0.05);
      HisChi[na]->GetXaxis()->SetLabelSize(0.04);
      HisChi[na]->GetXaxis()->SetNdivisions(NumBin);
      HisChi[na]->GetYaxis()->SetTitleOffset(1.4);
      HisChi[na]->GetYaxis()->SetTitleSize(0.05);
      HisChi[na]->GetYaxis()->SetLabelSize(0.04);
      HisChi[na]->SetFillColor(kRed);
      HisChi[na]->GetXaxis()->SetTitle("#chi^{2}");
      HisChi[na]->GetYaxis()->SetTitle("Entries");
      
      // The chi squared probabilities
      sprintf(DumName,"ChiPro Obs%i",n);
      HisDel = (TH1F*) gROOT->FindObject(DumName);
      if(HisDel)delete HisDel;
      HisPro[na] = new TH1F(DumName,GetNamObs(n).Data(),NumBin,0.,1.);
      HisPro[na]->GetXaxis()->SetTitleOffset(1.0);
      HisPro[na]->GetXaxis()->SetTitleSize(0.05);
      HisPro[na]->GetXaxis()->SetLabelSize(0.04);
      HisPro[na]->GetXaxis()->SetNdivisions(NumBin);
      HisPro[na]->GetYaxis()->SetTitleOffset(1.4);
      HisPro[na]->GetYaxis()->SetTitleSize(0.05);
      HisPro[na]->GetYaxis()->SetLabelSize(0.04);
      HisPro[na]->SetFillColor(kRed);
      HisPro[na]->GetXaxis()->SetTitle("P(#chi^{2},1)");
      HisPro[na]->GetYaxis()->SetTitle("Entries");
    
      // Fill histos from matrices
      for(Int_t i = 0; i<MaxNum; i++){
	if(C->operator()(na,i) > 0){
	  HisChi[na]->Fill(C->operator()(na,i));
	  HisPro[na]->Fill(P->operator()(na,i));
	}
      }
      // Set histo scales
      HisChi[na]->SetMinimum(0.);
      HisChi[na]->SetMaximum(1.2*HisChi[na]->GetMaximum()); 
      HisPro[na]->SetMinimum(0.);
      HisPro[na]->SetMaximum(1.2*HisPro[na]->GetMaximum()); 
      na = na + 1;
    }
  }

  // Book Canvases with at most five columns
  Int_t NumRow = 1 + GetActObs()/5;
  Int_t NumCol = TMath::Min(GetActObs(),5);
  //printf("... Blue->PrintCompatEst(): Row = %2i Col = %2i \n", NumRow, NumCol);
  TCanvas* CanDel;
  TCanvas* CanPoi[2];
  TString  CanNam[2] = {"Canv ChiQua", "Canv ChiPro"};
  Int_t CanDix = NumCol*300;
  Int_t CanDiy = TMath::Min(1200,NumRow*400);
  for(Int_t l = 0; l<2; l++){
    CanDel = (TCanvas*) gROOT->FindObject(CanNam[l]);
    if(CanDel) delete CanDel;
    CanPoi[l] = new TCanvas(CanNam[l],CanNam[l],0,0,CanDix,CanDiy);
    CanPoi[l]->Divide(NumCol,NumRow);
    for(Int_t i = 0; i<NumCol*NumRow; i++){
      CanPoi[l]->GetPad(i+1)->SetTopMargin(0.08);
      CanPoi[l]->GetPad(i+1)->SetBottomMargin(0.10);
      CanPoi[l]->GetPad(i+1)->SetLeftMargin(0.15);
      CanPoi[l]->GetPad(i+1)->SetRightMargin(0.02);
    }
  }

  // Draw histos
  for(Int_t n = 0; n<InpObs; n++){
    CanPoi[0]->cd(n+1); HisChi[n]->Draw();
    CanPoi[1]->cd(n+1); HisPro[n]->Draw();
  }

  // Write Files
  TString OutFil = " ";
  OutFil = FilNam + "_ComEst_ChiQua.pdf"; CanPoi[0]->Print(OutFil); 
  OutFil = FilNam + "_ComEst_ChiPro.pdf"; CanPoi[1]->Print(OutFil); 

  // Delete matrices
  C->Delete(); C = NULL;
  P->Delete(); P = NULL;

  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCompatEst() const {
  const TString FilNam = "NoGraphWanted";
  PrintCompatEst(FilNam);
};

//------------------------------------------------------------------------------

void Blue::PrintParams(const Int_t Ifl) const {
  if(IsCalcedParams() == 1){

    // Define dummy matrix
    TMatrixD* Dumm = new TMatrixD(InpEstOrig,InpEstOrig);

    // Write the header get the matrix 
    printf("\n");
    if(Ifl == 0){
      printf("... Blue->PrintParams(%2i): Sigma_i/Sigma_j i<j \n", Ifl);
      Dumm->operator=(*SRat);
    }else if(Ifl == 1){
      printf("... Blue->PrintParams(%2i): Beta \n", Ifl);
      Dumm->operator=(*Beta);
    }else if(Ifl == 2){
      printf("... Blue->PrintParams(%2i): Sigma_x/Sigma_i \n", Ifl);
      Dumm->operator=(*Sigx);
    }else if(Ifl == 3){
      printf("... Blue->PrintParams(%2i): dBeta/drho \n", Ifl);
      Dumm->operator=(*DBdr);
    }else if(Ifl == 4){
      printf("... Blue->PrintParams(%2i): 1/Sigma_i dSigma_x/dRho \n", Ifl);
      Dumm->operator=(*DSdr);
    }else if(Ifl == 5){
      printf("... Blue->PrintParams(%2i): dBeta/dz \n", Ifl);
      Dumm->operator=(*DBdz);
    }else if(Ifl == 6){
      printf("... Blue->PrintParams(%2i): 1/Sigma_i dSigma_x/dz \n", Ifl);
      Dumm->operator=(*DSdz);
    }else{
      printf("... Blue->PrintParams(%2i): Not implemented flag \n",Ifl);
    }
    
    // Print out the matrix
    Int_t IFound = 0;
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i) == 1){
	if(IFound == 0){
	  printf("... Blue->PrintParams(%2i):         %s",Ifl,GetNamEst(i).Data());
	}else{
	  printf(" %s",GetNamEst(i).Data());
	}
	IFound = IFound  + 1;
      }
    }
    printf("\n");
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i) == 1){
	printf("... Blue->PrintParams(%2i): %s",Ifl,GetNamEst(i).Data());
	for(Int_t j = 0; j<i; j++){
	  if(IsActiveEst(j) == 1)printf(" %+7.2f", Dumm->operator()(i,j));
	}
	printf("\n");
      }
    }
    // Delete dummy matrix
    //printf("... Blue->PrintParams(%2i): Dumm \n",Ifl); Dumm->Print();
    Dumm->Delete(); Dumm = NULL;
  }else{
    printf("... Blue->PrintParams(%2i): Presently not available, call FixInp() \n",Ifl);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintParams() const {
  for(Int_t k = 0; k<7; k++)PrintParams(k);
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintPull(const Int_t i) const {
  if(IsSolved() == 1){
    if(IsActiveEst(i) == 1){
      Int_t jj = 0;   
      for(Int_t j = 0; j < i; j++){
	jj = jj + IsActiveEst(j);
      }
      printf("... Blue->PrintPull(%2i): The pull for estimate %s is: %5.3f \n",
	     i, GetNamEst(i).Data(), Pul->operator()(jj));
    }else{
      if(IsPrintLevel() >= 1){
	printf("... Blue->PrintPull(%2i): Not an active estimate \n",i);
      }
    }
  }else{
    printf("... Blue->PrintPull: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintPull() const {
  if(IsSolved() == 1){
    for(Int_t i = 0; i<InpEstOrig; i++){
      PrintPull(i);
    }
  }else{
    printf("... Blue->PrintPull: Presently not available, call Solve() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintListObs() const {
  if(IsFixedInp() == 1){
    printf("... Blue->PrintListObs(): The present list of   observables is: ");
    for(Int_t n = 0; n<InpObs-1; n++){
      printf("%2i,", LisObs[n]);
    }
    printf("%2i \n", LisObs[InpObs-1]);
  }else{
    printf("... Blue->PrintListObs(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintNamObs() const {
  if(IsFixedInp() == 1){
    Int_t ICount = 0;
    printf("... Blue->PrintNamObs(): The names of the present list of observables are: \n");
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(IsActiveObs(n)){	
	ICount = ICount + 1;
	if(ICount == 1){
	  printf("... Blue->PrintNamObs(): %2i = %s", n, GetNamObs(n).Data());
	}else{
	  printf(", %2i = %s", n, GetNamObs(n).Data());
	}
	if(ICount == 5){
	  printf("\n");
	  ICount = 0;
	}
      }
    }
    if(ICount != 0)printf("\n");
  }else{
    printf("... Blue->PrintNamObs(): Presently not available, call FixInp() \n");
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::PrintCovRes() const {
  if(IsSolved() == 1){
    if(InpObs == 1){
      printf("... Blue->PrintCovRes: Single observable sig^2 = %5.3f \n",CovRes->operator()(0,0));
    }else{     
      printf("... Blue->PrintCovRes: The covariance matrix of the results \n");
      CovRes->Print();
    }
  }else{
    printf("... Blue->PrintCovRes: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintRhoRes() const {
  if(IsSolved() == 1){
    if(InpObs == 1){
      printf("... Blue->PrintCovRes: Single observable rho = %5.3f \n",RhoRes->operator()(0,0));
    }else{     
      printf("... Blue->PrintRhoRes: The correlation matrix of the results \n");
      RhoRes->Print();      
    }
  }else{
    printf("... Blue->PrintRhoRes: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintWeight() const {
  if(IsSolved() == 1){
    printf("... Blue->PrintWeight: The weight matrix of the combination \n");
    Lam->Print();
  }else{
    printf("... Blue->PrintWeight: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintResult() const {
  if(IsSolved() == 1){
    // Write out the linear combination
    printf("... Blue->PrintResult: Linear combination of estimates for observables \n");
    printf("               xi =  ");
    for(Int_t i = 0; i<InpEst; i++){
      printf("   %s    ",GetNamEst(IsWhichEst(i)).Data());
    }
    printf("\n");

    Int_t m = 0;
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(IsActiveObs(n) == 1){
	printf("%s = %5.4f = ",GetNamObs(n).Data(),XvaRes->operator()(m));
	for(Int_t i = 0; i<InpEst; i++){
	  if(i == 0){	  
	    printf("%+5.4f * x%1i ",Lam->operator()(i,m),IsWhichEst(i));
	  }else if(i<10){
	    printf("+ %+5.4f * x%1i ",Lam->operator()(i,m),IsWhichEst(i));
	  }else{
	    printf("+ %+5.4f * x%2i ",Lam->operator()(i,m),IsWhichEst(i));
	  }
	}
	m = m + 1;
	printf("\n");
      }
    }
    printf("\n");
    
    // Write out the uncertainty breakdown
    printf("... Blue->PrintResult: Breakdown of uncertainties \n");      

    printf(" Uncert =                 ");
    for(Int_t k = 0; k<InpUnc; k++){
      printf("%s  ",GetNamUnc(IsWhichUnc(k)).Data());
    }
    printf("\n");

    m = 0;
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(IsActiveObs(n) == 1){
	printf("%s = %5.4f +- %5.4f ",GetNamObs(n).Data(),XvaRes->operator()(m),
	       sqrt(CovRes->operator()(m,m)));
	if(InpUnc >= 2){     
	  for(Int_t k = 0; k<InpUnc; k++){
	    if(k == 0){
	      printf("(%5.4f",CorRes->operator()(k*InpObs+m,k*InpObs+m));
	    }else{
	      printf(" +- %5.4f",CorRes->operator()(k*InpObs+m,k*InpObs+m));
	    }
	    if(k == InpUnc-1){printf(")\n");}	  	   
	  }
	}else{
	  printf("\n");
	}
	m = m + 1;
      }
    } 
    printf("\n");

    // Write out stat + sys uncertainty
    printf("... Blue->PrintResult: Breakdown in stat (k==0)");
    if(InpUnc > 1)printf("+ sys (all k>0)");
    printf("\n");
    printf("\n\n**************************************************************************************** \n");
    printf(" Result =   Value (+-  stat");
    if(InpUnc > 1){printf(" +-  syst)");
    }else{printf(")");};
    printf(" = +-  full \n");
    m = 0;
    Double_t SysUnc = 0;
    for(Int_t n = 0; n<InpObsOrig; n++){
      SysUnc = 0;
      if(IsActiveObs(n) == 1){
	printf("%s = %5.4f",GetNamObs(n).Data(),XvaRes->operator()(m));
	for(Int_t k = 0; k<InpUnc; k++){
	  if(k == 0){
	    printf(" (+- %5.4f", CorRes->operator()(k*InpObs+m,k*InpObs+m));
	    if(InpUnc == 1)printf(")");
	  }else{
	    SysUnc = SysUnc + CorRes->operator()(k*InpObs+m,k*InpObs+m) * CorRes->operator()(k*InpObs+m,k*InpObs+m);
	    if(k == InpUnc-1){
	      printf(" +- %5.4f)", TMath::Sqrt(SysUnc));
	    }
	  }
	}
	printf(" = +- %5.4f",TMath::Sqrt(CovRes->operator()(m,m)));
	printf("\n");
        printf("**************************************************************************************** \n\n\n");
	m = m + 1;
      }
    }     
    printf("\n");
  }else{
    printf("... Blue->PrintResult: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintCompatObs() const {
  if(IsSolved() == 1){
    // Check number of observables
    if(GetActObs() == 1){
      printf("... Blue->PrintCompatObs: There is only one active observable \n");
      return;
    }
    // List the observables
    if(IsPrintLevel() >= 1){
      printf("... Blue->PrintCompatObs: The following observables were determined \n");
      PrintResult();
    }

    // Check compatibility of the observables
    // XvaRes->Print();  RhoRes->Print(); CovRes->Print();
    printf("\n");
    printf("... Blue->PrintCompatObs: The compatibility of the observables for Ndof = 1 \n");
    Int_t iok1 = 0, iok2 = 0;
    Double_t rho = 0., vao = 0., vat = 0., vac = 0., chi = 0.;
    for(Int_t n = 0; n<InpObs-1; n++){
      for(Int_t m = n+1; m<InpObs; m++){
	rho = RhoRes->operator()(n,m);
	vao = CovRes->operator()(n,n);
	vat = CovRes->operator()(m,m);
	vac = vao + vat - 2. * rho * TMath::Sqrt(vao) * TMath::Sqrt(vat);
	chi = (XvaRes->operator()(n) - XvaRes->operator()(m)) *
	  (XvaRes->operator()(n) - XvaRes->operator()(m)) / vac;
	if(chi >= 1)iok1 = iok1 + 1;
	if(chi >= 2)iok2 = iok2 + 1;
	//printf("... Blue->PrintCompatObs: i = %2i, j = %2i, rho= %5.3f, sig1 = %5.3f,  sig2 = %5.3f \n", 
	//n, m, rho, TMath::Sqrt(vao), TMath::Sqrt(vat));
	printf("... Blue->PrintCompatObs: (%s,%s): Sigma(Delta) = %5.3f,", 
	       GetNamObs(IsWhichObs(n)).Data(), GetNamObs(IsWhichObs(m)).Data(),
	       TMath::Sqrt(vac));
	printf(" Chiq = %5.3f, Prob = %5.3f \n", chi, TMath::Prob(chi,1));
      }
    }
    // Report problems if spotted
    if(iok2 > 0){
      printf("... Blue->PrintCompatObs: Pairs of observables with Chiq > 2 = %2i \n", iok2);
    }
    if(iok1 > 0){
      printf("... Blue->PrintCompatObs: Pairs of observables with Chiq > 1 = %2i \n", iok1);
    }else{
      printf("... Blue->PrintCompatObs: The compatibility of the observables is good");
      printf(" (all differences agree within 1 sigma) \n");
    }
    printf("\n");
  }else{
    printf("... Blue->PrintCompatObs: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintChiPro() const {
  if(IsSolved() == 1){
    printf("... Blue->PrintChiPro(): ChiQua = %5.3f for NDof = %2i", ChiQua, NumDof);
    printf(" Probability = %5.2f %% \n", 100*ChiPro);
  }else{
    printf("... Blue->PrintChiPro: Presently not available, call Solve() \n");
  }
};

//------------------------------------------------------------------------------

void Blue::PrintInspectLike() const {
  if(IsInspectedLike() == 1){
    Int_t LocFla = 0;
    for(Int_t n = 0; n<InpObs; n++){
      if(LikRes->operator()(n, 6) > 0.5){
	printf("... Blue->PrintInspectLike(): The next observable is: %s \n", 
	       GetNamObs(IsWhichObs(n)).Data());

	// The Flag of what was caculated
	LocFla = static_cast<Int_t>(LikRes->operator()(n, 6));
	printf("... Blue->PrintInspectLike(): The likelihood is");
	if(LocFla > 9){
	  LocFla = LocFla - 10;
	  printf(" Gaussian and");
	}
	if(LocFla == 1){printf(" exact.");
	}else if(LocFla == 2){printf(" only an approximation.");
	}
	printf("\n");
       
	printf("... Blue->PrintInspectLike(): Likelihood: x = %5.2f %+5.2f %+5.2f ==> %5.2f(+-%5.2f) \n",
	       LikRes->operator()(n,0), LikRes->operator()(n,1), LikRes->operator()(n,2),
	       LikRes->operator()(n,0), 0.5*(LikRes->operator()(n,2)-LikRes->operator()(n,1)));
	printf("... Blue->PrintInspectLike():       BLUE: x = %5.2f %+5.2f %+5.2f ==> %5.2f(+-%5.2f) \n",
	       LikRes->operator()(n,3), LikRes->operator()(n,4), LikRes->operator()(n,5),
	       LikRes->operator()(n,3), 0.5*(LikRes->operator()(n,5)-LikRes->operator()(n,4)));
      }
    }
  }else{
    printf("... Blue->PrintInspectLike(): Presently not available, call InspectLike() \n");
  }
};

//------------------------------------------------------------------------------
// After SolveXXX()

void Blue::PrintAccImp() const {
  if(IsSolvedAccImp() == 1){
    //printf("... Blue->PrintAccImpp(): IndImp\n"); IndImp->Print();
    //printf("... Blue->PrintAccImp(): ValImp\n"); ValImp->Print();
    //printf("... Blue->PrintAccImp(): UncImp\n"); UncImp->Print();

    // Print which algorithm and what cut was chosen
    printf("... Blue->PrintAccImp():\n");
    if(IntAcc == 0){
      printf("... Blue->PrintAccImp(): Sorted according to relative improvement, the default \n", IntAcc);
    }else if(IntAcc == 1){
      printf("... Blue->PrintAccImp(): Sorted according to absolute BLUE weights \n", IntAcc);
    }else if(IntAcc == 2){
      printf("... Blue->PrintAccImp(): Sorted according to inverse variance of the estimates \n", IntAcc);
    }
    printf("... Blue->PrintAccImp(): The decision was to stop when improvement is always less than %5.1f %% \n", PreAcc);
    printf("... Blue->PrintAccImp():\n");

    //printf("... Blue->PrintAccImp(): ValImp \n"); ValImp->Print();
    //printf("... Blue->PrintAccImp(): UncImp \n"); UncImp->Print();
    //printf("... Blue->PrintAccImp(): IndImp \n"); IndImp->Print();

    // Do the loop
    Double_t Values =  0, Uncert =  0;
    Double_t DifVal =  0, DifUnc =  0;
    Int_t    NumImp =  0;
    Int_t    IstImp =  1;
    Int_t    NexImp = -1;
    Int_t    IndNex = -1; 
    Int_t    MosImp = -1;
    Int_t    IndMos = -1;
    for(Int_t n = 0; n<InpObs; n++){
      NumImp = 0;
      IstImp = 1;
      printf("... Blue->PrintAccImp(): The list of importance for observable %2i = %s is: %2i",
	     IsWhichObs(n), GetNamObs(IsWhichObs(n)).Data(), static_cast<int>(IndImp->operator()(0,n))); 
      for(Int_t i = 1; i<InpEst; i++){
	NexImp = static_cast<int>(IndImp->operator()(i,n));
	if(NexImp != -1){
	  printf(", %2i",NexImp);
	  IstImp = IstImp + 1;
	}
      }
      printf(" \n");
      Values = ValImp->operator()(0,n);
      Uncert = UncImp->operator()(0,n);
      MosImp = static_cast<int>(IndImp->operator()(0,n));
      IndMos = IsIndexEst(MosImp);
      printf("... Blue->PrintAccImp(): The most precise estimate '%s' yields: %5.3f +- %5.3f \n", 
	     GetNamEst(MosImp).Data(),Values, Uncert);

      // Print the pair wise info
      Double_t dumm = 0;
      //      PrintListEst();
      if(IstImp > 1){
	TMatrixD* Dumm = new TMatrixD(InpEstOrig,InpEstOrig);
	printf("... Blue->PrintAccImp(): %2i = %s:", MosImp, GetNamEst(MosImp).Data());
	printf("   Value");
	printf(" +- Stat");
	printf("   +- Syst");
	printf("  rho");
	printf("          z");
	printf("       Beta");
	printf("    S_x/S_i");
	printf(" dBeta/drho");
	printf(" 1/S_i*dS_x/drho");
	printf("   dBeta/dz");
	printf("   1/S_i*dS_x/dz\n");
	for(Int_t i = 1; i<InpEst; i++){
	  NexImp = static_cast<int>(IndImp->operator()(i,n));
	  IndNex = IsIndexEst(NexImp);
	  if(NexImp != -1){
	    printf("... Blue->PrintAccImp(): %2i = %s:", NexImp, GetNamEst(NexImp).Data());
	    // Print the individual estimates
	    printf(" %5.3f", Xva->operator()(IndNex));
	    printf(" +- %5.3f", Unc->operator()(IndNex,IndNex));
	    dumm = TMath::Sqrt(TMath::Power(Sig->operator()(IndNex),2.0)-
			       TMath::Power(Unc->operator()(IndNex,IndNex),2.0));
	    printf(" +- %5.3f", dumm);

	    // Print the correlation
	    printf(" %4.2f", Rho->operator()(IndMos,IndNex));

	    // Print the parameters of the pair wise combinations
	    for(Int_t Ifl = 0; Ifl<7; Ifl++){
	      GetParams(Ifl, Dumm);
	      if(Ifl == 4 || Ifl == 6)printf("     ");
	      dumm = Dumm->operator()(IndMos,IndNex);
	      if(IndMos < IndNex)dumm = Dumm->operator()(IndNex,IndMos);
	      if(Ifl == 0 && dumm < 1.0)dumm = 1./dumm;
	      printf("     %+5.3f", dumm);
	    }
	    printf("\n");
	  }
	}
	Dumm->Delete(); Dumm = NULL;
      }
      printf("... Blue->PrintAccImp():\n");

      // Print the differences
      for(Int_t i = 1; i<InpEst; i++){
	NexImp = static_cast<int>(IndImp->operator()(i,n));
	if(NexImp != -1){
	  DifUnc = 100*(UncImp->operator()(i,n) - UncImp->operator()(i-1,n))/UncImp->operator()(i-1,n);
	  DifVal = 100*(ValImp->operator()(i,n) - ValImp->operator()(i-1,n))/ValImp->operator()(i-1,n);
	  printf("... Blue->PrintAccImp(): Add N = %2i ==> Differences: Value = %+5.3f %%, Uncertainty = %+5.3f %% \n",
		 NexImp,DifVal,DifUnc);
	  if(DifUnc <= -PreAcc)NumImp = i;
	}
      }

      // List the findings
      if(NumImp == 0){
	printf("... Blue->PrintAccImp(): There is only one estimate for this observable \n");
      }else{
	if(NumImp == IstImp){
	  printf("... Blue->PrintAccImp(): Combine all %2i estimates \n",IstImp);
	}else{
	  printf("... Blue->PrintAccImp(): Only combine the first %2i estimates up to and including N = %2i \n",
		 NumImp+1, static_cast<int>(IndImp->operator()(NumImp,n)));
	}
      }
      printf("... Blue->PrintAccImp():\n");
    }
    PrintChiPro();
  }else{
    printf("... Blue->PrintAccImp(): Presently not available, call SolveAccImp() \n");
  }
  return;  
};

//------------------------------------------------------------------------------

void Blue::PrintScaRho(const TString FilNam) const {

  // Check whether solved
  if(IsSolvedScaRho() == 0)return;

  // Print the steering
  TString ModCor = "to be filled later";
  if(FlaFac == 0){
    printf("... Blue->PrintScaRho(): Independent (recommended) variation of the uncertainties sources \n");
    ModCor = "_Ind";
  }else{
    printf("... Blue->PrintScaRho(): Combined variation of the uncertainties sources \n");
    ModCor = "_Sim";
  }
  if(InpFac == 1){ModCor = ModCor + "One";
  }else{ModCor = ModCor + "Mor";};

  // Print the remaining matrix
  printf("... Blue->PrintScaRho():\n");
  printf("... Blue->PrintScaRho(): The remaining matrix of groups of correlations is: \n");
  printf("... Blue->PrintScaRho():     ");
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i) == 1){
      printf(" %2i", i);
    }
  }
  printf("\n");
  printf("... Blue->PrintScaRho():    ---------------- \n");
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(IsActiveEst(i) == 1){
      printf("... Blue->PrintScaRho(): %2i |", i);
      for(Int_t j = 0; j<InpEstOrig; j++){
	if(IsActiveEst(j) == 1){
	  printf(" %2i", static_cast<int>(MatFac->operator()(i,j)));
	}
      }
      printf("\n");
    }
  }

  // Print the rho ranges
  printf("... Blue->PrintScaRho():\n");
  printf("... Blue->PrintScaRho(): The scanned ranges for the groups are:\n");
  printf("... Blue->PrintScaRho():  Source:");
  for(Int_t l = 0; l<InpFac; l++){printf("     l = %2i    |", l);};
  printf("\n");
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      printf("... Blue->PrintScaRho(): %s:", GetNamUnc(k).Data());
      for(Int_t l = 0; l<InpFac; l++){
	printf(" %5.3f - %5.3f |", MinFac->operator()(k,l), MaxFac->operator()(k,l));
      }
      printf("\n");
    }
  }
  printf("... Blue->PrintScaRho(): \n");
  printf("... Blue->PrintScaRho(): The values are given for");
  printf(" %2i decreasing correlation assumptions, using the above ranges \n", InrFac);

  // Print the failures if needed
  if(FaiFac > 0){
    printf("... Blue->PrintScaRho(): The matrix inversion failed %3i times \n", FaiFac);
  }

  // Loop over observables
  Int_t ka = 0;
  for(Int_t n = 0; n<InpObs; n++){
    printf("... Blue->PrintScaRho():\n");
    printf("... Blue->PrintScaRho(): Next observable: %s\n", GetNamObs(IsWhichObs(n)).Data());
    for(Int_t l = 0; l<InpFac; l++){
      printf("... Blue->PrintScaRho():\n");
      printf("... Blue->PrintScaRho(): Next group l =  %2i\n", l);
      // The values
      printf("... Blue->PrintScaRho(): -------------------- The differences in the        values");
      printf(" ----------------------------- \n");
      ka = 0;
      for(Int_t k = 0; k<InpUncOrig; k++){
	if(IsActiveUnc(k) == 1){
	  printf("... Blue->PrintScaRho(): %s:", GetNamUnc(k).Data());
	  for(Int_t ll = 0; ll<InrFac; ll++){
	    printf(" %+5.3f ", ValSca[n]->operator()(ka+l*InpUnc,ll));
	  }
	  printf("\n");
	  ka = ka + 1;
	}
      }
      printf("... Blue->PrintScaRho(): ----------------------------------------------");
      printf("------------------------------------------ \n");
      printf("... Blue->PrintScaRho():   Total:");
      for(Int_t ll = 0; ll<InrFac; ll++){
	printf(" %+5.3f ", VtoSca[n]->operator()(l,ll));
      }
      printf("\n");

      // The uncertainties
      printf("... Blue->PrintScaRho():\n");
      printf("... Blue->PrintScaRho(): -------------------- The differences in the uncertainties");
      printf(" ----------------------------- \n");
      ka = 0;
      for(Int_t k = 0; k<InpUncOrig; k++){
	if(IsActiveUnc(k) == 1){
	  printf("... Blue->PrintScaRho(): %s:", GetNamUnc(k).Data());
	  for(Int_t ll = 0; ll<InrFac; ll++){
	    printf(" %+5.3f ", SigSca[n]->operator()(ka+l*InpUnc,ll));
	  }
	  printf("\n");
	  ka = ka + 1;
	}
      }
      printf("... Blue->PrintScaRho(): ----------------------------------------------");
      printf("------------------------------------------ \n");
      printf("... Blue->PrintScaRho():   Total:");
      for(Int_t ll = 0; ll<InrFac; ll++){
	printf(" %+5.3f ", StoSca[n]->operator()(l,ll));
      }
      printf("\n");
    }
  }

  // Return if no plots wanted
  if(FilNam == "NoGraphWanted")return;

  //
  // ------------- Now do the figures
  //
  // Book Canvases with at most five columns
  char  DumName[50];
  Int_t NumRow = 1 + GetNumScaFac()/5;
  Int_t NumCol = TMath::Min(GetNumScaFac(),5);
  //printf("... Blue->PrintScaRho(): Row = %2i Col = %2i \n", NumRow, NumCol);
  
  TCanvas* CanDel;
  TCanvas* CanVal[InpObs];
  TCanvas* CanUnc[InpObs];
  Int_t CanDix = NumCol*400;
  Int_t CanDiy = TMath::Min(1200,NumRow*400);
  for(Int_t n = 0; n<InpObs; n++){   
    sprintf(DumName,"Canv:%s_Result",GetNamObs(IsWhichObs(n)).Data());
    CanDel = (TCanvas*) gROOT->FindObject(DumName);
    if(CanDel) delete CanDel;
    CanVal[n] = new TCanvas(DumName,DumName,0,0,CanDix,CanDiy);
    CanVal[n]->Divide(NumCol,NumRow);
    for(Int_t i = 1; i < NumCol*NumRow+1; i++){
      CanVal[n]->GetPad(i)->SetTopMargin(0.10);
      CanVal[n]->GetPad(i)->SetBottomMargin(0.17);
      CanVal[n]->GetPad(i)->SetLeftMargin(0.12);
      CanVal[n]->GetPad(i)->SetRightMargin(0.02);
    }

    sprintf(DumName,"Canv:%s_Uncert",GetNamObs(IsWhichObs(n)).Data());
    CanDel = (TCanvas*) gROOT->FindObject(DumName);
    if(CanDel) delete CanDel;
    CanUnc[n] = new TCanvas(DumName,DumName,0,0,CanDix,CanDiy);
    CanUnc[n]->Divide(NumCol,NumRow);
    for(Int_t i = 1; i < NumCol*NumRow+1; i++){
      CanUnc[n]->GetPad(i)->SetTopMargin(0.10);
      CanUnc[n]->GetPad(i)->SetBottomMargin(0.17);
      CanUnc[n]->GetPad(i)->SetLeftMargin(0.12);
      CanUnc[n]->GetPad(i)->SetRightMargin(0.02);
    }
  }
  
  // The Dummy Histos and the TGraphs that hold the data
  Int_t NumHis = InpObs * InpFac;
  TH1F* DumVal[NumHis];
  TH1F* DumUnc[NumHis];
  const Int_t NumTgr = 3;
  TGraph* GraDel;
  TGraph* ScaVal[NumHis][NumTgr];
  TGraph* ScaUnc[NumHis][NumTgr];
  Int_t   NumBin = GetActUnc();
  Int_t   IdxCol[NumTgr] = {kBlue, kGreen, kRed};
  //printf("... Blue->PrintScaRho(): NumHis = %2i \n", NumHis);

  Int_t   IdxFac[NumTgr] = {3, 6, 9};
  sprintf(DumName,"Steps_%i/%i/%i=B/G/R",IdxFac[0]+1,IdxFac[1]+1,IdxFac[2]+1);
  TString IndSte = &DumName[0];
  TString DumNam = "To be filled later";

  // A horizontal line at zero
  Double_t MinHis = -0.5, MaxHis = InpUnc - 0.5;
  TLine* l0 = new TLine(MinHis,0,MaxHis,0);
  l0->SetLineColor(kBlack);
  l0->SetLineWidth(2);

  //printf("... Blue->PrintScaRho(): MaxFac \n"); MaxFac->Print();
  //printf("... Blue->PrintScaRho(): MinFac \n"); MinFac->Print();

  // Loop over observables
  for(Int_t n = 0; n<InpObs; n++){
    
    // Loop over rho factors
    for(Int_t l = 0; l<InpFac; l++){

      // Book Dummy Histos
      sprintf(DumName,"Res:%s_Group_%i_%s",GetNamObs(IsWhichObs(n)).Data(),l,IndSte.Data());
      //gDirectory->ls();
      gDirectory->Delete(DumName);
      DumVal[n*InpFac+l] = new TH1F(DumName,DumName,InpUnc,MinHis,MaxHis);
      DumVal[n*InpFac+l]->GetXaxis()->SetNdivisions(InpUnc);
      
      sprintf(DumName,"Unc:%s_Group_%i_%s",GetNamObs(IsWhichObs(n)).Data(),l,IndSte.Data());
      gDirectory->Delete(DumName);
      DumUnc[n*InpFac+l] = new TH1F(DumName,DumName,InpUnc,MinHis,MaxHis);
      DumUnc[n*InpFac+l]->GetXaxis()->SetNdivisions(InpUnc);

      // Book TGraphs
      for(Int_t ij = 0; ij<NumTgr; ij++){
	ScaVal[n*InpFac+l][ij] = new TGraph(NumBin);
	ScaVal[n*InpFac+l][ij]->SetFillColor(IdxCol[ij]);

     	ScaUnc[n*InpFac+l][ij] = new TGraph(NumBin);
	ScaUnc[n*InpFac+l][ij]->SetFillColor(IdxCol[ij]);
      }

      // Loop over uncertainties
      ka = 0;
      for(Int_t k = 0; k<InpUncOrig; k++){
	if(IsActiveUnc(k) == 1){
	  //printf("... Blue->PrintScaRho(): Next uncertainty %2i \n",ka);
	  DumVal[n*InpFac+l]->GetXaxis()->SetBinLabel(ka+1, GetNamUnc(k).Data());
	  DumUnc[n*InpFac+l]->GetXaxis()->SetBinLabel(ka+1, GetNamUnc(k).Data());
	  //for(Int_t ll = 0; ll<InrFac; ll++){
	  //printf("%2i, v=%+9.5f, s=%+9.5f \n",ka, ValSca[n]->operator()(ka+l*InpUnc,ll), 
	  //SigSca[n]->operator()(ka+l*InpUnc,ll));
	  //}
	  // Fill the TGraphs
	  for(Int_t ij = 0; ij<NumTgr; ij++){	    
	    ScaVal[n*InpFac+l][ij]->SetPoint(ka,ka,ValSca[n]->operator()(ka+l*InpUnc,IdxFac[ij]));
            ScaUnc[n*InpFac+l][ij]->SetPoint(ka,ka,SigSca[n]->operator()(ka+l*InpUnc,IdxFac[ij]));
	  }
	  ka = ka + 1;
	}
      }

      // Get maximum per histogram      
      Double_t MaxVal = 0, MinVal = 0;
      Double_t MaxUnc = 0, MinUnc = 0;
      Double_t MaxAct = 0, MinAct = 0;
      for(Int_t ij = 0; ij<NumTgr; ij++){
	MaxAct = TMath::MaxElement(NumBin,ScaVal[n*InpFac+l][ij]->GetY());
	MinAct = TMath::MinElement(NumBin,ScaVal[n*InpFac+l][ij]->GetY());
	//printf("... Blue->PrintScaRho(): Val max = %5.3f, min = %5.3f \n", MaxAct, MinAct);
	if(ij == 0){
	  MaxVal = MaxAct;
	  MinVal = MinAct;
	}else{
	  MaxVal = TMath::Max(MaxVal, MaxAct);
	  MinVal = TMath::Min(MinVal, MinAct);
	}
	if(ij == NumTgr-1){
	  MaxAct = TMath::Abs(MinVal);
	  MaxVal = 1.2*TMath::Max(MaxVal, MaxAct);
	  MinVal = -1.*MaxVal;
	}
  
	MaxAct = TMath::MaxElement(NumBin,ScaUnc[n*InpFac+l][ij]->GetY());
	MinAct = TMath::MinElement(NumBin,ScaUnc[n*InpFac+l][ij]->GetY());
	//printf("... Blue->PrintScaRho(): Unc max = %5.3f, min = %5.3f \n", MaxAct, MinAct);
	if(ij == 0){
	  MaxUnc = MaxAct;
	  MinUnc = MinAct;
	}else{
	  MaxUnc = TMath::Max(MaxUnc, MaxAct);
	  MinUnc = TMath::Min(MinUnc, MinAct);
	}	
	if(ij == NumTgr-1){
	  MaxAct = TMath::Abs(MinUnc);
	  MaxUnc = 1.2*TMath::Max(MaxUnc, MaxAct);
	  MinUnc = -1.*MaxUnc;
	}
      }

      // Now set max and min
      DumVal[n*InpFac+l]->SetMaximum(MaxVal);
      DumVal[n*InpFac+l]->SetMinimum(MinVal);
      DumVal[n*InpFac+l]->LabelsOption("v");
      //printf("... Blue->PrintScaRho(): Val: max = %5.3f, min = %5.3f \n", MaxVal, MinVal);

      DumUnc[n*InpFac+l]->SetMaximum(MaxUnc);
      DumUnc[n*InpFac+l]->SetMinimum(MinUnc);
      DumUnc[n*InpFac+l]->LabelsOption("v");
      //printf("... Blue->PrintScaRho(): Unc: max = %5.3f, min = %5.3f \n", MaxUnc, MinUnc);

      // Set the Yaxis and offsets
      DumNam = "#Delta" + GetNamObs(IsWhichObs(n));
      DumVal[n*InpFac+l]->GetXaxis()->SetTitleOffset(1.8);
      DumVal[n*InpFac+l]->GetXaxis()->SetLabelSize(0.05);
      DumVal[n*InpFac+l]->GetYaxis()->SetTitle(DumNam);
      DumVal[n*InpFac+l]->GetYaxis()->SetTitleOffset(1.8);
      DumVal[n*InpFac+l]->GetYaxis()->SetTitleSize(0.03);
      DumVal[n*InpFac+l]->GetYaxis()->SetLabelSize(0.03);
     
      DumNam = "#Delta #sigma(" + GetNamObs(IsWhichObs(n)) + ")";
      DumUnc[n*InpFac+l]->GetXaxis()->SetTitleOffset(1.8);
      DumUnc[n*InpFac+l]->GetXaxis()->SetLabelSize(0.05);
      DumUnc[n*InpFac+l]->GetYaxis()->SetTitle(DumNam);
      DumUnc[n*InpFac+l]->GetYaxis()->SetTitleOffset(1.8);
      DumUnc[n*InpFac+l]->GetYaxis()->SetTitleSize(0.03);
      DumUnc[n*InpFac+l]->GetYaxis()->SetLabelSize(0.03);   

      CanVal[n]->cd(l+1); 
      DumVal[n*InpFac+l]->Draw();
      for(Int_t ij = NumTgr-1; ij>-1; ij--)ScaVal[n*InpFac+l][ij]->Draw("Bsame");
      l0->Draw("same");

      CanUnc[n]->cd(l+1);
      DumUnc[n*InpFac+l]->Draw();
      for(Int_t ij = NumTgr-1; ij>-1; ij--)ScaUnc[n*InpFac+l][ij]->Draw("Bsame");
      l0->Draw("same");
    }
    
    // Write Files
    char ObsNum[8];
    sprintf(ObsNum,"_%i",n);
    TString OutFil = " ";
    OutFil = FilNam + "_ScaRho" + ModCor + "_Val_Obs" + ObsNum + ".pdf";
    CanVal[n]->Print(OutFil);
    OutFil = FilNam + "_ScaRho" + ModCor + "_Unc_Obs" + ObsNum + ".pdf";
    CanUnc[n]->Print(OutFil);
  }

  return;
};

//------------------------------------------------------------------------------

void Blue::PrintScaRho() const {
  const TString FilNam = "NoGraphWanted";
  PrintScaRho(FilNam);
};

//------------------------------------------------------------------------------

void Blue::PrintInfWei() const {
  if(IsSolvedInfWei() == 1){
    // Print the information weights
    //    Marginal weights == MarWei(i)
    //   Intrinsic weights == IntWei(i)
    // Correlation weights == IntWei(InpEstOrig)
    //         Blue weights == BluWei(i)
    // Sum up Blue weights
    Double_t RelWei = 0;
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i)){RelWei = RelWei + TMath::Abs(BluWei->operator()(i));};
    }
    if(RelWei == 0){
      printf("... Blue->PrintInfWei: IGNORED: Sum of weights is zero. \n");
      return;  
    }
    
    printf("... Blue->PrintInfWei():     Estimate   Blue Intrinsic Marginal Relative \n");
    printf("... Blue->PrintInfWei(): ----------------------------------------------- \n");
    for(Int_t i = 0; i<InpEstOrig; i++){
      if(IsActiveEst(i)){
	printf("... Blue->PrintInfWei(): %12i %+5.1f%%    %+5.1f%%   %+5.1f%%   %+5.1f%% \n", i, 
	       100*BluWei->operator()(i),100*IntWei->operator()(i),
	       100*MarWei->operator()(i),100*TMath::Abs(BluWei->operator()(i))/RelWei);
      }
    }
    printf("... Blue->PrintInfWei(): ----------------------------------------------- \n");
    printf("... Blue->PrintInfWei(): Correlations          %+5.1f%% \n", 
	   100*IntWei->operator()(InpEstOrig));
    printf("... Blue->PrintInfWei: \n");    
    PrintChiPro();
  }else{
    printf("... Blue->PrintInfWei: Presently not available, call SolveInfWei() \n");
  }
  return;  
};

//------------------------------------------------------------------------------

void Blue::PrintMaxVar() const {
  if(IsSolvedMaxVar() == 1){
    printf("... Blue->PrintMaxVar: The variance at the start = %10.8f \n",
	   VarMin->operator()(0));
    printf("... Blue->PrintMaxVar: The correlation at the start \n"); RhoMin->Print();

    if(IFuFai[IFuFla]>0){
      printf("... Blue->PrintMaxVar: The number of Matrix-Inversion problems encountered are %5i \n",
	     IFuFai[IFuFla]);
    }
    if(IFuFla == 0){
      printf("... Blue->PrintMaxVar: The scale factors f_ijk = f = %5.3f \n", RhoFco);
    }else if(IFuFla == 1){
      printf("... Blue->PrintMaxVar: The scale factors f_ijk = f_k (k>0) \n"); RhoFpk->Print();
    }else if(IFuFla == 2){
      TMatrixD* RhoAct = new TMatrixD(InpEst,InpEst);
      Int_t ia = 0, ja = 0;
      for(Int_t i = 0; i<InpEstOrig; i++){
	if(IsActiveEst(i) == 1){
	  RhoAct->operator()(ia,ia) = 1;
	  ja = ia + 1;	  
	  for(Int_t j = i+1; j<InpEstOrig; j++){
	    if(IsActiveEst(j) == 1){
	      RhoAct->operator()(ia,ja) = RhoFij->operator()( i, j);
	      RhoAct->operator()(ja,ia) = RhoAct->operator()(ia,ja);
	      ja = ja + 1;
	    }
	  }
	  ia = ia + 1;
	}
      }
      printf("... Blue->PrintMaxVar: The scale factors f_ijk = f_ij for  (k>0) \n"); 
      RhoAct->Print();
      RhoAct->Delete(); RhoAct = NULL;
    }   
    printf("... Blue->PrintMaxVar: The variance after convergence = %10.8f \n",
	   VarMax->operator()(0));
    printf("... Blue->PrintMaxVar: The correlation after convergence \n"); RhoMax->Print();
    //--not yet
    //TMatrixD* A = new TMatrixD(InpEst,InpEst);
    //A->SetSub(0,0,*RhoOoz);
    //A->operator-=(*RhoMax);
    //printf("... Blue->PrintMaxVar(): The correlations for rho_ij = 1/zij \n"); RhoOoz->Print();
    //printf("... Blue->PrintMaxVar(): The difference of the maximisation from rho_ij = 1/zij \n"); A->Print();
    PrintChiPro();
  }else{
    printf("... Blue->PrintMaxVar: Presently not available, call SolveMaxVar() \n");
    return;
  }
  return;
};


//------------------------------------------------------------------------------

void Blue::PrintStatus() const {
  printf("\n");

  //------------------------- Available after constructor
  // Print version and date
  printf("... Blue->PrintStatus(): This is the BLUE Software Version %s run on %s \n",
	 Versio.Data(),Today.Data());
  // Print version and date
  printf("... Blue->PrintStatus(): Constructed InpEst =%3i, InpUnc =%3i, InpObs =%2i\n",
	 InpEstOrig, InpUncOrig, InpObsOrig);

  if(InpFac >1){
    printf("... Blue->PrintStatus(): The matrix of the %2i groups of correlations: \n", InpFac);
    for(Int_t i = 0; i<InpEstOrig; i++){
      printf("... Blue->PrintStatus(): ");
      for(Int_t j = 0; j<InpEstOrig; j++){
	printf("%3i",static_cast<int>(MatFac->operator()(i,j)));
	if(j < InpEst-1)printf(",");
      }
      printf("\n");
    }
  }

  //------------------------- Available after fixed
  if(IsFixedInp() == 1){
    // Changed correlations
    Int_t iok = 0;
    for(Int_t k = 0; k<InpUncOrig; k++){ 
      if(IsActiveUnc(k) == 1 && IsRhoValUnc(k) == 1){
	if(iok == 0){
	  printf("... Blue->Status(): Changed correlations for uncertainties: %2i", k);
	  iok = 1;
	}else{
	  printf(", %2i",k);
	}
      }
    }
    if(iok == 1){printf("\n");};
    
    // Scaled correlations
    iok = 0;
    for(Int_t k = 0; k<InpUncOrig; k++){ 
      if(IsActiveUnc(k) == 1 && IsRhoFacUnc(k) == 1){
	if(iok == 0){
	  printf("... Blue->Status(): Scaled  correlations for uncertainties: %2i", k);
	  iok = 1;
	}else{
	  printf(", %2i",k);
	}
      }
    }
    if(iok == 1){printf("\n");};
    
    // Reduced correlations
    iok = 0;
    for(Int_t k = 0; k<InpUncOrig; k++){ 
      if(IsActiveUnc(k) == 1 && IsRhoRedUnc(k) == 1){
	if(iok == 0){
	  printf("... Blue->Status(): Reduced correlations for uncertainties: %2i", k);
	  iok = 1;
	}else{
	  printf(", %2i",k);
	}
      }
    }
    if(iok == 1){printf("\n");};
    
    // Relative uncertainties
    iok = 0;
    for(Int_t k = 0; k<InpUncOrig; k++){ 
      if(IsActiveUnc(k) == 1 && IsRelValUnc(k) != 0){
	if(iok == 0){
	  printf("\n");
	  printf("... Blue->Status(): Relative uncertainties for source: %2i", k);
	  iok = 1;
	}else{
	  printf(", %2i",k);
	}
      }
    }
    if(iok == 1){ printf("\n");};
  }

  //------------------------- Available after fixed

  // The active lists
  if(IsFixedInp() == 1){
    printf("\n");
    PrintListEst();
    PrintListUnc();
    PrintNamUnc();
    PrintListObs();
    PrintNamObs();
    printf("\n");
    PrintEst();

    // Details of the input  
    PrintParams();
    if(IsPrintLevel() >= 1){
      PrintRho();
      PrintCov();
      if(IsPrintLevel() >= 2){
	PrintCovInvert();
	PrintCompatEst();
      }
    }
  }
  
  //------------------------- Available after solved

  if(IsSolved() == 1){
    printf("\n");
    PrintResult();
    PrintChiPro();
    PrintRhoRes();
    PrintWeight();
    if(IsPrintLevel() >= 1){
      PrintCovRes();
      if(IsPrintLevel() >= 2){
	PrintCompatObs();
      }
    }
  }
  return;
};

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

Double_t Blue::GetPara(const Int_t ifl, const Double_t rho, const Double_t zva) const {

  // Set variables
  char what[20];
  Double_t zael;
  Double_t nenn = 1 - 2*rho*zva + pow(zva,2);
  Double_t resu = -99;

  Double_t tolzva = TMath::Abs(1 - zva);
  //Double_t absrho = TMath::Abs(rho);

  Int_t lfl = ifl;
  if(lfl > 6){lfl = lfl - 6;
  }else if(lfl < -6){lfl = lfl + 6;
  }
  if((rho >= 0.9999 && tolzva <= 0.0001) || 
     zva <=  0.9999 || 
     rho <= -1.0001 || 
     rho >=  1.0001 || 
     TMath::Abs(lfl) >= 7){
    printf("... GetPara(%2i,%5.3f,%5.3f): IGNORED Outside of validity \n", ifl,rho,zva);
    resu = -99.;
  }else{
    if(lfl == 1 || lfl == -1){
      // beta
      sprintf(what,"%1i:beta",lfl);
      zael = (1 - rho*zva);
      resu = zael/nenn;
    }else if(lfl == 2 || lfl == -2){
      // sigx/sig1
      sprintf(what,"%1i:sigx",lfl);
      zael = zva*zva * (1 - rho*rho);
      resu = sqrt(zael/nenn);
    }else if(lfl == 3|| lfl == -3){
      // dbeta/drho 
      sprintf(what,"%1i:dbdr",lfl);
      zael = zva*(1-zva*zva);
      resu = zael/(nenn*nenn);
    }else if(lfl == 4|| lfl == -4){
      // dsigx/drho
      sprintf(what,"%1i:dsdr",lfl);
      // RN:10.01.2014, should read z instead of z^2
      // resu = zva * zva * (zva - rho) * (1-rho*zva) / sqrt((1-rho*rho)*pow(nenn,3));
      resu = zva * (zva - rho) * (1-rho*zva) / sqrt((1-rho*rho)*pow(nenn,3));
    }else if(lfl == 5|| lfl == -5){
      // dbeta/dzva
      sprintf(what,"%1i:dbdz",lfl);
      zael = rho * (1 + zva*zva) - 2*zva;
      resu = zael/(nenn*nenn);
    }else if(lfl == 6|| lfl == -6){
      // dsigx/dzva
      sprintf(what,"%1i:dsdz",lfl);
      //RN20.11.13 if(absrho >= 0.999999)rho = 0.999999*rho;
      resu = (1-rho*zva) * sqrt((1-rho*rho)/pow(nenn,3));
    }
  }
  if (lfl <= 0){
    printf("... Blue->GetPara: %s(rho=%5.3f,z=%5.3f) = %5.3f \n",
	   what,rho,zva,resu);
  }
  return resu;
};

//------------------------------------------------------------------------------

Double_t Blue::FunPara(const Double_t *x, const Double_t *par) const {
  Int_t ifl = static_cast<Int_t>(par[1]);
  Double_t zva = 0;
  Double_t rho = 0;
  Double_t val = 0;
  if(TMath::Abs(ifl)< 5){
    zva = static_cast<Double_t>(par[0]);
    rho = static_cast<Double_t>(x[0]);
  }else{
    zva = static_cast<Double_t>(x[0]);
    rho = static_cast<Double_t>(par[0]);
  }
  val = GetPara(ifl, rho, zva);
  //RN 10.01.14
  //if(ifl == 4)val = TMath::Abs(val);
  return val;
};

//------------------------------------------------------------------------------

void Blue::InspectPair(const Int_t i,  const Int_t j) const {
  TString FilNam = "NoGraphWanted";
  Int_t IndFig = 0;
  InspectPair(i, j, FilNam, IndFig);
};

//------------------------------------------------------------------------------

void Blue::InspectPair(const Int_t i, const Int_t j, const TString FilNam) const {
  Int_t IndFig = 0;
  InspectPair(i, j, FilNam, IndFig);
};

//------------------------------------------------------------------------------

void Blue::InspectPair(const Int_t i, const Int_t j, const TString FilNam, const Int_t IndFig) const {
  if(IsFixedInp() == 0){
    printf("... Blue->InspectPair(%2i,%2i,FilNam,IndFig): IGNORED Please first fix the input \n",
	   i,j,IndFig);
    return;
  }
  if(IsActiveEst(i) == 0 || IsActiveEst(j) == 0){
    printf("... Blue->InspectPair(%2i,%2i,FilNam,IndFig): IGNORED Works only for active estimates \n",
	   i,j,IndFig);
    return;
  }
  if(i == j){
    printf("... Blue->InspectPair(%2i,%2i,FilNam,)IndFig: IGNORED Works only for two different estimates \n",
	   i,j,IndFig);
    return;
  }
  if(EstWhichObs(i) != EstWhichObs(j)){
    printf("... Blue->InspectPair(%2i,%2i,FilNam,IndFig): IGNORED Estimate of a different observables \n",
	   i,j,IndFig);
    return;
  }
  if(IndFig < 0 || IndFig > 1){
    printf("... Blue->InspectPair(%2i,%2i,FilNam): IGNORED Invalid IndFig only (0,1) is allowed \n",i,j,IndFig);
    return;
  }

  // Save indices in -all- estimates
  Int_t ii = i, ij = j;

  // Get indices in -active- estimates
  Int_t i1 = -1, i2 = -1;
  for(Int_t jj = 0; jj<InpEst; jj++){
    if(IsWhichEst(jj) == i)i1 = jj;
    if(IsWhichEst(jj) == j)i2 = jj;
  }
  
  // Get values and correlation
  Double_t xv1, xv2, sv1, sv2, rho, zva;
  // 1) Check which has larger variance exchange if needed
  sv1 = Sig->operator()(i1);
  sv2 = Sig->operator()(i2);
  if(sv1 > sv2){
    Int_t i3 = i2;
    i2 = i1;
    i1 = i3;
    sv1 = Sig->operator()(i1);
    sv2 = Sig->operator()(i2);
    i3 = ii;
    ii = ij;
    ij = i3;
  }
  printf("... Blue->InspectPair: Estimate '%s' is the more precise result \n", GetNamEst(ii).Data());
  zva = sv2 / sv1;
  xv1 = Xva->operator()(i1);
  xv2 = Xva->operator()(i2);
  rho = Cov->operator()(i1,i2)/(sv1*sv2);
  printf("... Blue->InspectPair: x1(%s) = %5.3f +- %5.3f \n", GetNamEst(ii).Data(), xv1, sv1);
  printf("... Blue->InspectPair: x2(%s) = %5.3f +- %5.3f \n", GetNamEst(ij).Data(), xv2, sv2);
  printf("... Blue->InspectPair: rho = %5.3f, z = %5.3f \n", rho, zva);

  // Look for consistency of the input
  Double_t dx, sd, adx, nsd;
  dx  = xv1 - xv2;
  adx = TMath::Abs(dx);
  sd  = sqrt(pow(sv1,2) + pow(sv2,2) - 2*rho*sv1*sv2);
  if(sd == 0){
    printf("... Blue->InspectPair: Vanishing variance(%5.3f) of the difference \n", sd);
  }else{
    printf("... Blue->InspectPair: Compatibility: dx +- sd = %5.3f +- %5.3f (= %5.3f sigma)\n",
	   dx, sd, adx/sd);
  }
  nsd = adx/sd;
  if(nsd <= 1){
    printf("... Blue->InspectPair: Compatibility is good. \n");
  }else if(nsd <= 2){
    printf("... Blue->InspectPair: Compatibility is about ok. \n");
  }else if(nsd <= 3){
    printf("... Blue->InspectPair: Compatibility is bad! \n");
  }else if(nsd <= 4){
    printf("... Blue->InspectPair: Compatibility is very bad!! \n");
  }

  // Get parameters for combination
  Double_t beta = GetPara(1, rho, zva);
  Double_t dbdr = GetPara(3, rho, zva);
  Double_t dbdz = GetPara(5, rho, zva);
  Double_t sigx = GetPara(2, rho, zva);
  Double_t dsdr = GetPara(4, rho, zva);
  Double_t dsdz = GetPara(6, rho, zva);

  // Do the combination
  Double_t xv, sv;
  xv = (1-beta) * xv1 + beta*xv2;
  sv = sigx*sv1;
  printf("... Blue->InspectPair: The combined result: x +- sx = %5.3f +- %5.3f \n", xv, sv);

  // Print out parameters
  printf("... Blue->InspectPair: The parameters of Eq.1 - Eq.6 \n");
  printf("... Blue->InspectPair: Par: beta = %+5.3f, dbeta/drho = %+5.3f, dbeta/dz = %+5.3f \n",
	 beta, dbdr, dbdz);
  printf("... Blue->InspectPair: Par: sigx = %+5.3f, dsigx/drho = %+5.3f, dsigx/dz = %+5.3f \n",
	 sigx, dsdr, dsdz);

  // Look at sensitivity 
  // A) change rho by at most +-10% absolute
  // B) change sigma2 by at most +-10% relative
  // Keep same in DrawSens
  Double_t r[3] = {TMath::Max(rho-0.1, -0.99), rho, TMath::Min(rho+0.1, 0.99)};
  Double_t z[3] = {TMath::Min(zva,TMath::Max(zva*0.9, 1.01)), zva, zva*1.1};
  printf("... Blue->InspectPair: The  variation is performed for: %5.3f < rho < %5.3f and %5.3f < z < %5.3f \n",
	 r[0], r[2], z[0], z[2]);

  // Loop over all combinations, find range in x and sigmax
  Double_t x[9], s[9];
  Double_t xmi = 0, xma = 0, smi = 0, sma = 0;
  for(Int_t iii = 0;iii<3; iii++){
    for(Int_t jjj = 0; jjj<3; jjj++){
      beta = GetPara(1, r[iii], z[jjj]);
      x[iii+jjj] = (1-beta) * xv1 + beta*xv2;
      s[iii+jjj] = GetPara(2, r[iii], z[jjj]) * sv1;
      if(iii == 0 && jjj == 0){
	xmi = x[iii+jjj];
	xma = xmi;
	smi = s[iii+jjj];
	sma = smi;
      }else{
	xmi = TMath::Min(xmi, x[iii+jjj]);
	xma = TMath::Max(xma, x[iii+jjj]);
	smi = TMath::Min(smi, s[iii+jjj]);
	sma = TMath::Max(sma, s[iii+jjj]);
      }
      printf("... Blue->InspectPair: rho = %5.3f, z = %5.3f, x +- sx = %5.3f +- %5.3f \n", 
	   r[iii], z[jjj], x[iii+jjj], s[iii+jjj]);
    }
  }
  // Get percentage deviations
  Double_t pxmi, pxma, psmi, psma;
  pxmi = 100*(xv-xmi)/xv;
  pxma = 100*(xma-xv)/xv;
  psmi = 100*(sv-smi)/sv;
  psma = 100*(sma-sv)/sv;

  // Print the result
  printf("... Blue->InspectPair: Range for central Value: %5.3f < %5.3f < %5.3f (-%4.1f + %4.1f)%% \n",
	 xmi, xv, xma, pxmi, pxma);
  printf("... Blue->InspectPair: Range for   Uncertainty: %5.3f < %5.3f < %5.3f (-%4.1f + %4.1f)%% \n",
	 smi, sv, sma, psmi, psma);

  // Graphics if wanted
  if(FilNam == "NoGraphWanted"){return;};

  // Concatenate filename
  char EstOne[8]; sprintf(EstOne,"%i",i);
  char EstTwo[8]; sprintf(EstTwo,"%i",j);
  TString FilInd = FilNam + "_x" + EstOne +"_x" + EstTwo;
  DrawSens(xv1, xv2, sv1, sv2, rho, FilInd, IndFig);
};

//------------------------------------------------------------------------------

void Blue::DrawSens(const Double_t xv1, const Double_t xv2, 
		    const Double_t sv1, const Double_t sv2, 
		    const Double_t rho, const TString FilNam) const {
  Int_t IndFig = 0;
  DrawSens(xv1, xv2, sv1, sv2, rho, FilNam, IndFig);
};

//------------------------------------------------------------------------------

void Blue::DrawSens(const Double_t xv1, const Double_t xv2, 
		    const Double_t sv1, const Double_t sv2, 
		    const Double_t rho, const TString FilNam, const Int_t IndFig) const {
  Double_t zva = sv2/sv1;
  if(zva < 1){
    printf("... Blue->DrawSens: IGNORED x2 is the more precise result, ");
    printf("please exchange the values! return \n");
    return;
  }
  printf("... Blue->DrawSens: x1 = %5.2f +- %5.2f, x2 = %5.2f +- %5.2f with (rho = %5.2f, z = %5.2f) \n", 
	 xv1,sv1,xv2,sv2, rho,zva);
  
  // Keep same in InspectPair
  Double_t r[3] = {TMath::Max(rho-0.1, -0.99), rho, TMath::Min(rho+0.1, 0.99)};
  Double_t z[3] = {TMath::Min(zva,TMath::Max(zva*0.9, 1.01)), zva, zva*1.1};

  // Get parameters for combination
  Double_t beta = GetPara(1, rho, zva);
  Double_t sigx = GetPara(2, rho, zva);
  Double_t    x = (1-beta) * xv1 + beta * xv2;
  Double_t   sx = sigx * sv1;

  // Parameters for the functions, 8 parameters and three values
  const Int_t NPoint = 500;
  const Int_t NumPar = 8;
  const Int_t NumVal = 3;
  const Int_t NumFunc = NumPar * NumVal;

  // Get description names
  // Dynmamic ones per function (function and derivative)
  char funName[120];
  char poiValu[120];
  Double_t der = 0;
  Double_t dhi = 0;
  Double_t dlo = 0;
  Double_t ddu = 0;
  Double_t dzz = 1 + 0.20*(2*zva-1);
  char relName[120];
  Double_t rel = 0;
  Double_t relx[NumPar] = {-0.9, -0.80, -0.9, -0.9,  dzz,  dzz,  dzz,  dzz};
  Double_t rely[NumPar] = {-0.5,  0.08, -2.3, -0.4, -2.2, -0.9, -0.5, 0.08};

  // Colors
  const Int_t IcoFla[4] = {kBlue, kBlack, kRed};

  // Static ones only once
  char xcoName[50]; sprintf(xcoName,"x = %5.2f +- %4.2f",x,sx);
  TLatex *Txco = new TLatex(-0.85,-0.9,xcoName);
  char xv1Name[50]; sprintf(xv1Name,"x_{1} = %5.2f +- %4.2f",xv1,sv1);
  TLatex *Txv1 = new TLatex(-0.90,-1.1,xv1Name);
  char xv2Name[50]; sprintf(xv2Name,"x_{2} = %5.2f +- %4.2f",xv2,sv2);
  TLatex *Txv2 = new TLatex(-0.90,-1.3,xv2Name);

  // The ranges in z
  char dzrsName[50]; sprintf(dzrsName,"f(#rho): z =");
  char dzr0Name[50]; sprintf(dzr0Name," %4.2f,",z[0]);
  char dzr1Name[50]; sprintf(dzr1Name," %4.2f,",z[1]);
  char dzr2Name[50]; sprintf(dzr2Name," %4.2f",z[2]);
  TLatex *Tdzrs = new TLatex(-0.90-0.03,-2.7,dzrsName);
  TLatex *Tdzr0 = new TLatex(-0.90+0.42,-2.7,dzr0Name);
  TLatex *Tdzr1 = new TLatex(-0.90+0.76,-2.7,dzr1Name);
  TLatex *Tdzr2 = new TLatex(-0.90+1.10,-2.7,dzr2Name);
  Tdzr0->SetTextColor(IcoFla[0]);
  Tdzr1->SetTextColor(IcoFla[1]);
  Tdzr2->SetTextColor(IcoFla[2]);

  // The ranges in rho
  Double_t length = 2*zva-1;
  char drrsName[50]; sprintf(drrsName,"f(z): #rho =");
  char drr0Name[50]; sprintf(drr0Name," %4.2f,",r[0]);
  char drr1Name[50]; sprintf(drr1Name," %4.2f,",r[1]);
  char drr2Name[50]; sprintf(drr2Name," %4.2f",r[2]);
  TLatex *Tdrrs = new TLatex(dzz            ,-2.6,drrsName);
  TLatex *Tdrr0 = new TLatex(dzz+0.23*length,-2.6,drr0Name);
  TLatex *Tdrr1 = new TLatex(dzz+0.40*length,-2.6,drr1Name);
  TLatex *Tdrr2 = new TLatex(dzz+0.57*length,-2.6,drr2Name);
  Tdrr0->SetTextColor(IcoFla[0]);
  Tdrr1->SetTextColor(IcoFla[1]);
  Tdrr2->SetTextColor(IcoFla[2]);

  // Set up the combined canvas
  TCanvas *CanDel;
  const Int_t CanDix = 4*300;
  const Int_t CanDiy = 2*400;
  CanDel = (TCanvas*) gROOT->FindObject("myCanvas");
  if(CanDel) delete CanDel;
  TCanvas *myCanvas = new TCanvas("myCanvas","Combine",0,0,CanDix,CanDiy);
  myCanvas->Divide(4,2);
  for(Int_t i = 0; i < NumPar; i++){
    myCanvas->GetPad(i+1)->SetTopMargin(0.02);
    myCanvas->GetPad(i+1)->SetBottomMargin(0.12);
    myCanvas->GetPad(i+1)->SetLeftMargin(0.15);
    myCanvas->GetPad(i+1)->SetRightMargin(0.02);
  }
  
  // The individaul ones
  TCanvas* Can[NumPar];
  Int_t lmax = 1;
  if(IndFig == 1){
    lmax = 2;
    char CanName[50];
    for(Int_t i = 0; i < NumPar; i++){
      sprintf(CanName,"Canvas%i",i);
      CanDel = (TCanvas*) gROOT->FindObject(CanName);
      if(CanDel) delete CanDel;
      Can[i] = new TCanvas(CanName,CanName,0,0,300,400);
      Can[i]->SetTopMargin(0.02);
      Can[i]->SetBottomMargin(0.12);
      Can[i]->SetLeftMargin(0.15);
      Can[i]->SetRightMargin(0.02);
    }	
  }

  // Pointer for the marker of the actual pair, plus x value
  TMarker *TP[NumFunc];
  Double_t xval[NumPar] = {rho, rho ,rho, rho, zva, zva, zva, zva};

  // Pointer to functions and dummy histograms plus names
  TF1* Func[NumFunc];
  char FuncName[50];
  TH2F* Dum[NumFunc];
  char DumName[50];

  // Dimensions for dummy histos and list of subfigures vs i
  // subfig = the canvas number 
  // numfig = the number of the figure
  Int_t  numfig[NumPar] = {  1,   2,  3,    4,     7,     8,     5,     6};
  Int_t  subfig[NumPar] = {  1,   5,  2,    6,     4,     8,     3,     7};
  Double_t xhig[NumPar] = {   1,  1,  1,    1, 2*zva, 2*zva, 2*zva, 2*zva};
  Double_t xlow[NumPar] = {  -1, -1, -1,   -1,     1,     1,     1,     1};
  Double_t yhig[NumPar] = { 0.6,  1,  0,  1.5,   0.5,   0.5,   0.6,     1};
  Double_t ylow[NumPar] = {-1.5,  0, -3, -0.5,    -3,    -1,  -1.5,     0};

  // Loop over functions(i) and values(j)  
  Int_t jj = 0;
  Int_t ifl;
  for(Int_t i = 0; i < NumPar; i++){
    ifl = i + 1;

    // Dummy plot, to be filled in second loop
    sprintf(DumName,"Dummy%i",i);
    gDirectory->Delete(DumName);
    Dum[i] = new TH2F(DumName, "", 20, xlow[i], xhig[i], 20, ylow[i], yhig[i]);
    Dum[i]->GetXaxis()->SetTitleOffset(0.75);
    Dum[i]->GetXaxis()->SetTitleSize(0.07);
    Dum[i]->GetXaxis()->SetLabelSize(0.043);
    Dum[i]->GetYaxis()->SetTitleOffset(1.00);
    Dum[i]->GetYaxis()->SetTitleSize(0.07);
    Dum[i]->GetYaxis()->SetLabelSize(0.043);
    //    Dum[i]->Draw();
    
    // Prepare text for each subfigure
    // 1) X-axis
    if(ifl < 5){Dum[i]->GetXaxis()->SetTitle("#rho");
    }else{Dum[i]->GetXaxis()->SetTitle("z");
    }

    // 2) Y-axis, function names and derivative calculation
    der = GetPara(ifl, rho, zva);
    if(ifl < 3 || ifl > 6){
      dlo = der;
      dhi = der;
      for(Int_t j = 0; j < 3; j++){
	for(Int_t k = 0; k < 3; k++){
	  ddu = GetPara(ifl, r[j], z[k]);
	  //printf("... Blue->DrawSens: ifl = %2i, ddu(%2i,%2i) = %7.5f \n", ifl, j, k, ddu);
	  if(dlo > ddu)dlo = ddu;
	  if(dhi < ddu)dhi = ddu;
	}
      }
    }else if(ifl < 5){
      dlo = GetPara(ifl, rho, z[0]);
      dhi = GetPara(ifl, rho, z[2]);
    }else{
      dlo = GetPara(ifl, r[0], zva);
      dhi = GetPara(ifl, r[2], zva);
    }
    if(dhi < dlo){
      ddu = dhi;
      dhi = dlo;
      dlo = ddu;
    }
    //printf("... Blue->DrawSens: ifl = %2i, min = %5.2f, max = %5.2f \n", ifl, dlo, dhi);

    // Now fill the names per sub-figure
    if(ifl == 1 || ifl == 7){
      sprintf(poiValu,"%4.2f",der);
      Dum[i]->GetYaxis()->SetTitle("#beta");
      if(ifl == 1)sprintf(relName,"%4.2f < #beta < %4.2f",dlo,dhi);
    }else if(ifl == 2 || ifl == 8){
      sprintf(poiValu,"%4.2f",der);      
      Dum[i]->GetYaxis()->SetTitle("#sigma_{x}/#sigma_{1}");
      if(ifl == 2)sprintf(relName,"%4.2f < #sigma_{x}/#sigma_{1} < %4.2f",dlo,dhi);
    }else if(ifl == 3){
      sprintf(poiValu,"%4.2f",der);
      Dum[i]->GetYaxis()->SetTitle("#partial#beta/#partial#rho");
      sprintf(relName,"%4.2f < #partial#beta/#partial#rho < %4.2f ",dlo,dhi);
    }else if(ifl == 4){
      sprintf(poiValu,"%4.2f",der);
      Dum[i]->GetYaxis()->SetTitle("1/#sigma_{1} #partial#sigma_{x}/#partial#rho");
      sprintf(relName,"%4.2f < 1/#sigma_{1}#partial#sigma_{x}/#partial#rho < %4.2f",dlo,dhi);
    }else if(ifl == 5){
      sprintf(poiValu,"%4.2f",der);
      Dum[i]->GetYaxis()->SetTitle("#partial#beta/#partialz");
      sprintf(relName,"%4.2f < #partial#beta/#partialz < %4.2f ",dlo,dhi);
    }else if(ifl == 6){
      sprintf(poiValu,"%4.2f",der);
      Dum[i]->GetYaxis()->SetTitle("1/#sigma_{1} #partial#sigma_{x}/#partialz");
      sprintf(relName,"%4.2f < 1/#sigma_{1}#partial#sigma_{x}/#partialz < %4.2f",dlo,dhi);
    }
    
    // Now draw on the individual and the combined Canvas
    for(Int_t l = 0; l < lmax; l++){
      if(l == 0){myCanvas->cd(subfig[i]);
      }else{Can[i]->cd();
      }
      Dum[i]->Draw();
      
      if(ifl == 1){
	Txco->Draw("same");
	Txv1->Draw("same");
	Txv2->Draw("same");
      }else if(ifl == 3){
	Tdzrs->Draw("same");
	Tdzr0->Draw("same");
	Tdzr1->Draw("same");
	Tdzr2->Draw("same");
      }else if(ifl == 5){
	Tdrrs->Draw("same");
	Tdrr0->Draw("same");
	Tdrr1->Draw("same");
	Tdrr2->Draw("same");
      }

      // The three functions per plot
      for(Int_t j = 0; j < NumVal; j++){
	sprintf(FuncName,"FunPara%i",i*NumPar+j);
	Func[jj] = new TF1(FuncName,this,&Blue::FunPara,xlow[i],xhig[i],2,"Blue","FunPara");
	if(ifl < 5){Func[jj]->SetParameter(0,z[j]);
	}else{Func[jj]->SetParameter(0,r[j]);
	}
	Func[jj]->SetParameter(1,ifl);
	Func[jj]->SetNpx(NPoint);
	Func[jj]->SetLineColor(IcoFla[j]);
	if(j == 1){Func[jj]->SetLineStyle(kDashed);};
	Func[jj]->Draw("same");
	
	// Increment counter
	jj = jj +1;
      }
      
      // Draw the function, the actual point and the ranges
      TLatex *Trel = new TLatex(relx[i],rely[i],relName);
      TP[i] = new TMarker(xval[i],GetPara(ifl, rho, zva),20);
      TP[i]->SetMarkerColor(IcoFla[1]);
      TP[i]->SetMarkerSize(1.5);
      TP[i]->Draw("same");
      if(ifl < 7)Trel->Draw("same");
    }
  }

  // Save canvas
  TString pdf = FilNam + "_InsPai.pdf";
  myCanvas->Print(pdf); 

  if(IndFig == 1){
    char FileName[80];
    TString IflNam[8] = {"a", "e", "b", "f", "d", "h", "c", "g"};
    for(Int_t i = 0; i < NumPar; i++){
      sprintf(FileName,"%s_InsPai_%s.pdf",FilNam.Data(),IflNam[i].Data());
      Can[i]->Print(FileName);
    }
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::InspectLike(const Int_t n) {
  TString FilNam = "NoGraphWanted";
  InspectLike(n, FilNam);
};

//------------------------------------------------------------------------------

void Blue::InspectLike(const Int_t n, const TString FilNam) {

  // Check for consistency 
  if(IsSolved() == 0){
    // Solved ?
    printf("... Blue->InspectLike(%2i): IGNORED Input not yet Solved call Solve \n", n);
    return;
  }else if(IsActiveObs(n) == 0){
    // Active ?
    printf("... Blue->InspectLike(%2i): IGNORED Not an active observable \n", n);
    return;
  }

  // Get the index within the present list
  Int_t na = IsIndexObs(n);

  // Check if filled already
  if(LikRes->operator()(na,6) != 0){
    // Already done ?
    printf("... Blue->InspectLike(%2i): IGNORED Already called for this observable \n", n);
    return;
  }

  // Check if at least one relative uncertainty exist
  LikFla = 0;
  Int_t IsRela = 0;
  for(Int_t k = 0; k<InpUncOrig; k++)if(IsRelValUnc(k) == 1)IsRela = 1;
  if(IsRela == 0){
    printf("... Blue->InspectLike(%2i): The likelihood is Gaussian",n);
    printf(" (There is no relative uncertainty.) \n");
    LikFla = 10;
  }else if(IsSolvedRelUnc() == 0){
    printf("... Blue->InspectLike(%2i): WARNING: Relative uncertainties exist,",n);
    printf(" but not solved by SolveRelUnc(). Hope this was intended. \n");
  }

  // Check the number of observables
  if(InpObs == 1){
    printf("... Blue->InspectLike(%2i): The likelihood is   exact",n);
    printf(" (There is only one observable) \n");
    LikFla = LikFla + 1;
  }else{
    printf("... Blue->InspectLike(%2i): The likelihood is only an approximation",n);
    printf(" (There are more than one observable) \n");
    LikFla = LikFla + 2;
  }

  // Fill the type of likelihood
  LikRes->operator()(na,6) = static_cast<Double_t>(LikFla);

  // Variable for the names
  char AnyName[80];

  // Variables for the functions
  Double_t sigx = 0, xlow =     0, xhig = 0, nsig = 5;
  Double_t xste = 0, nste = 10000, xact = 0, fact = 0;
  Double_t par0 = static_cast<Double_t>(n);
  Double_t par1 = static_cast<Double_t>(InpEst);

  // Variables for the likelihood scan
  Double_t LikMax = 0, LogLik = 0, LxxMax = 0, BluMax = 0, BxxMax = 0;
  Double_t LikLow = 0, LogBlu = 0, LikHig = 0, BluLow = 0, BluHig = 0;

  // Set the x-axis limits and the step size
  sigx = TMath::Sqrt(CovRes->operator()(na,na));
  xlow = XvaRes->operator()(na) - nsig*sigx;
  xhig = XvaRes->operator()(na) + nsig*sigx;
  xste = (xhig-xlow)/nste;

  // The garbage collector
  TF1* TF1Del;

  // Define the likelihood function
  TF1Del = (TF1*) gROOT->GetListOfFunctions()->FindObject("FuncLike");
  if(TF1Del)delete TF1Del;
  TF1* fLike = new TF1("FuncLike",this,&Blue::Likelihood,xlow,xhig,3,"Blue","Likelihood");
  fLike->SetParameter(0, par0);
  fLike->SetParameter(1, par1);
  fLike->SetParameter(2,   1.);
  fLike->SetNpx(100);
  fLike->SetLineColor(kRed);
  sprintf(AnyName,"Obs_{%i}",n);
  fLike->GetXaxis()->SetTitle(AnyName);
  fLike->GetYaxis()->SetTitleOffset(1.4);
  sprintf(AnyName,"pdf(Est_{i}, Obs_{%i})",n); 
  fLike->GetYaxis()->SetTitle(AnyName);

  // Get the maximum of the likelihood and its x value
  LxxMax = fLike->GetMaximumX(xlow,xhig);
  LikMax = fLike->GetMaximum(xlow,xhig);
  LogLik = -2.*TMath::Log(LikMax);

  // Get the aymmetric uncertainties of the likelihood
  LikLow = xhig;
  xact   = LxxMax;
  while(LikLow == xhig){
    xact = xact - xste;
    fact = -2.*TMath::Log(fLike->Eval(xact));
    if(fact > LogLik + 1)LikLow = xact;
    if(xact < xlow){
      LikLow = xlow;
      printf("... Blue->InspectLike(%2i): Failed: No lower uncertainty found", n);
      printf(" for the likelihood within plot range. Set to %5.3f \n", LikLow);
    }
  }
  LikHig = xlow;
  xact   = LxxMax;
  while(LikHig == xlow){
    xact = xact + xste;
    fact = -2.*TMath::Log(fLike->Eval(xact));
    if(fact > LogLik + 1)LikHig = xact;
    if(xact > xhig){
      LikHig = xhig;
      printf("... Blue->InspectLike(%2i): Failed: No upper uncertainty found", n);
      printf(" for the likelihood within plot range. Set to %5.3f \n", LikHig);
    }
  }

  // Fill the matrix
  LikRes->operator()(na,0) = LxxMax;
  LikRes->operator()(na,1) = LikLow-LxxMax;
  LikRes->operator()(na,2) = LikHig-LxxMax;
  printf("... Blue->InspectLike(%2i): L_Like(%5.3f) = %6.3e, %+5.3f < x < %+5.3f ==> %5.3f(+-%5.3f) \n",
	 n,LxxMax,LikMax,LikLow-LxxMax,LikHig-LxxMax,LxxMax,0.5*(LikHig-LikLow));

  // Define the BLUE function
  TF1Del = (TF1*) gROOT->GetListOfFunctions()->FindObject("FuncBlue");
  if(TF1Del)delete TF1Del;
  TF1* fBlue = new TF1("FuncBlue",this,&Blue::Likelihood,xlow,xhig,3,"Blue","Likelihood");
  fBlue->SetParameter(0, par0);
  fBlue->SetParameter(1, par1);
  fBlue->SetParameter(2,   0.);
  fBlue->SetNpx(100);
  fBlue->SetLineColor(kBlue);

  // Get the maximum of the BLUE function and its x value
  BxxMax = fBlue->GetMaximumX(xlow,xhig);
  BluMax = fBlue->GetMaximum(xlow,xhig);
  LogBlu = -2.*TMath::Log(BluMax);

  // Get the uncertainties of the BLUE function
  BluLow = xhig; 
  xact   = BxxMax;
  while(BluLow == xhig){
    xact = xact - xste;
    fact = -2.*TMath::Log(fBlue->Eval(xact));
    if(fact > LogBlu + 1)BluLow = xact;
    if(xact < xlow){
      BluLow = xlow;
      printf("... Blue->InspectLike(%2i): Failed: No lower uncertainty found", n);
      printf(" for BLUE within plot range. Set to %5.3f \n", BluLow);
    }
  }
  BluHig = xlow;
  xact   = BxxMax;
  while(BluHig == xlow){
    xact = xact + xste;
    fact = -2.*TMath::Log(fBlue->Eval(xact));
    if(fact > LogBlu + 1)BluHig = xact;
    if(xact > xhig){
      BluHig = xhig;
      printf("... Blue->InspectLike(%2i): Failed: No upper uncertainty found", n);
      printf(" for BLUE within plot range. Set to %5.3f \n", BluHig);
    }
  }

  // Fill the matrix
  LikRes->operator()(na,3) = BxxMax;
  LikRes->operator()(na,4) = BluLow-BxxMax;
  LikRes->operator()(na,5) = BluHig-BxxMax;
  printf("... Blue->InspectLike(%2i): L_Blue(%5.3f) = %6.3e, %+5.3f < x < %+5.3f ==> %5.3f(+-%5.3f) \n",
	 n,BxxMax,BluMax,BluLow-BxxMax,BluHig-BxxMax,BxxMax,0.5*(BluHig-BluLow));
  
  // Do the plots if wanted
  if(FilNam != "NoGraphWanted"){

    // Define the canvas     
    sprintf(AnyName,"%s_CanvLike_%i",FilNam.Data(),n);
    TCanvas* CanDel;
    CanDel = (TCanvas*) gROOT->FindObject(AnyName);
    if(CanDel) delete CanDel;
    TCanvas* CanLik;
    CanLik = new TCanvas(AnyName,AnyName,0,0,600,600);
    CanLik->SetBottomMargin(0.10);
    CanLik->SetTopMargin(0.05);
    CanLik->SetLeftMargin(0.15);
    CanLik->SetRightMargin(0.02);
    CanLik->SetBorderSize(4);

    // Set the maximum
    Double_t ymax = 1.1*TMath::Max(BluMax,LikMax);
    fLike->SetMaximum(ymax);

    // Get lines for the likelihood +- 1 sigma    
    TLine* lLikel = new TLine(LikLow,0,LikLow,fLike->Eval(LikLow));
    lLikel->SetLineColor(kRed); lLikel->SetLineWidth(2); lLikel->SetLineStyle(kDashed);
    TLine* lLikec = new TLine(LxxMax,0,LxxMax,fLike->Eval(LxxMax));
    lLikec->SetLineColor(kRed); lLikec->SetLineWidth(2);
    TLine* lLikeh = new TLine(LikHig,0,LikHig,fLike->Eval(LikHig));
    lLikeh->SetLineColor(kRed); lLikeh->SetLineWidth(2); lLikeh->SetLineStyle(kDashed);

    // Get lines for BLUE +- 1 sigma    
    TLine* lBluel = new TLine(BluLow,0,BluLow,fBlue->Eval(BluLow));
    lBluel->SetLineColor(kBlue); lBluel->SetLineWidth(2); lBluel->SetLineStyle(kDashed);
    TLine* lBluec = new TLine(BxxMax,0,BxxMax,fBlue->Eval(BxxMax));
    lBluec->SetLineColor(kBlue); lBluec->SetLineWidth(2);
    TLine* lBlueh = new TLine(BluHig,0,BluHig,fBlue->Eval(BluHig));
    lBlueh->SetLineColor(kBlue); lBlueh->SetLineWidth(2); lBlueh->SetLineStyle(kDashed);

    // Get values lines for likelihood +- 1 sigma    
    Double_t xtxt[2] = {xlow + 0.063*(xhig-xlow), xlow + 0.040*(xhig-xlow)};
    Double_t ytxt[2] = {0.96*ymax,0.92*ymax};

    sprintf(AnyName,"Like = %5.3f (%-5.3f %+5.3f) (+-%5.3f)", LxxMax, LikLow-LxxMax, LikHig-LxxMax, 0.5*(LikHig-LikLow)); 
    TLatex *TLike = new TLatex(xtxt[0],ytxt[0],AnyName);
    TLike->SetTextColor(kRed);
    TLike->SetTextSize(0.032);

    // Get values lines for BLUE +- 1 sigma    
    sprintf(AnyName,"BLUE = %5.3f (%-5.3f %+5.3f)", BxxMax, BluLow-BxxMax, BluHig-BxxMax); 
    TLatex *TBlue = new TLatex(xtxt[1],ytxt[1],AnyName);
    TBlue->SetTextColor(kBlue);
    TBlue->SetTextSize(0.032);

    // Set Style and draw
    CanLik->cd(1); 
    fLike->Draw(); 
    lLikel->Draw("same"); lLikec->Draw("same"); lLikeh->Draw("same"); TLike->Draw("same");
    fBlue->Draw("same"); 
    lBluel->Draw("same"); lBluec->Draw("same"); lBlueh->Draw("same"); TBlue->Draw("same");
    
    // Write the figure
    char ObsNum[8];
    sprintf(ObsNum,"%i",n);
    TString OutFil = FilNam + "_InsLik_Obs_" + ObsNum + ".pdf";
    CanLik->Print(OutFil);
  }

  // Clean up    
  // gROOT->GetListOfFunctions()->Print();
  delete fLike; fLike = NULL;
  delete fBlue; fBlue = NULL;

  // Set the Flag
  SetIsInspectLike(1);

  return;
};

//------------------------------------------------------------------------------

Int_t Blue::InspectResult() const {

  // Check whether solved
  if(IsSolved() == 0){
    printf("... Blue->InspectResult(): Presently not available, call Solve() \n");
    return 0;
  }

  // Check the results for negative or -nan total variance
  Int_t CovNeg = 0;
  for(Int_t n = 0; n<InpObs; n++){
    if(CovRes->operator()(n,n) < 0 || TMath::IsNaN(CovRes->operator()(n,n)) == 1)CovNeg = 1;
  }

  // Check the results for larger variance than the best estimate
  Int_t CovLar = 0, CouEst = 0;
  for(Int_t n = 0; n<InpObs; n++){
    CouEst = 0;
    for(Int_t i = 0; i<InpEst; i++){
      if(EstWhichObs(IsWhichEst(i)) == IsWhichObs(n))CouEst = CouEst + 1;
    }
    if(CouEst > 1){
      for(Int_t i = 0; i<InpEst; i++){
	if(IsWhichEst(i) == GetPreEst(IsWhichObs(n))){
	  if(TMath::Sqrt(CovRes->operator()(n,n)) > Sig->operator()(i))CovLar = 1;
	}
      }
    }
  }

  // Check the results for negative or -nan individual uncertainties
  Int_t CorNeg = 0;
  Int_t m = 0;
  for(Int_t n = 0; n<InpObsOrig; n++){
    if(IsActiveObs(n) == 1){
      //printf("... Blue->InspectResult():");
      for(Int_t k = 0; k<InpUnc; k++){
	//printf(" %7.5f", CorRes->operator()(k*InpObs+m,k*InpObs+m));
	if(CorRes->operator()(k*InpObs+m,k*InpObs+m) < 0 ||
	   TMath::IsNaN(CorRes->operator()(k*InpObs+m,k*InpObs+m)) == 1)CorNeg = 1;
      }
      m = m + 1;
      //      printf("\n");
    }
  }

  // Show the findings if wanted
  if(CovNeg == 1 || CorNeg == 1 || CovLar == 1){
    if(IsPrintLevel() >= 1){
      Int_t NegEig = GetMatEigen(Cov);
      if(CorNeg == 1)printf("... Blue->InspectResult(): At least one uncertainty is negative or -nan \n");
      if(CovNeg == 1)printf("... Blue->InspectResult(): The total uncertainty is negative or -nan \n");
      if(CovLar == 1){
	printf("... Blue->InspectResult(): The total uncertainty is larger than");
	printf(" the one of the most precise estimate\n");
      }
      if(NegEig == 0)printf("... Blue->InspectResult(): Covariance with negative Eigenvalue(s) \n");
      PrintStatus();
  }
    return -1*(100*CovLar + 10*CovNeg + CorNeg);
  }
  return 1;
};

//------------------------------------------------------------------------------

void Blue::LatexResult(const TString FilNam) const {
  TString ForVal = "%5.2f";
  TString ForUnc = ForVal;
  TString ForWei = "%4.2f";
  TString ForRho = ForWei;
  TString ForPul = ForVal;
  LatexResult(FilNam, ForVal, ForUnc, ForWei, ForRho, ForPul);
  return;
};

//------------------------------------------------------------------------------

void Blue::LatexResult(const TString FilNam, const TString ForVal, const TString ForUnc,
		       const TString ForWei, const TString ForRho, const TString ForPul) const {
  if(IsSolved() == 0){
    printf("... Blue->LatexResult: Presently not available, call Solve() \n");
    return;
  }

  // Decide on sidewaystable based on NumMax
  const Int_t NumMax = 10;

  // Open file
  TString outfile = FilNam + ".tex";
  std::ofstream ofs (outfile, std::ofstream::out);

  // Declare variables
  char c[100];
  TVectorD* SysEst = new TVectorD(InpEst);
  TVectorD* SysObs = new TVectorD(InpObs);
  Double_t Sys = 0;
  Int_t NumCol = 0, Isidew = 0;
  TString Format = "To be filled later";

  // Get an idea of the dimension
  NumCol = InpEst + InpObs + 1;
  if(NumCol > NumMax)Isidew = 1;

  // Dummy return variable for sprintf
  Int_t IRet = 0;
  if(IRet == 1)return;

  // Define header
  IRet = sprintf(c,"\\documentclass[11pt,a4paper]{article}\n"); ofs<<c;
  IRet = sprintf(c,"\\usepackage{underscore}\n"); ofs<<c;
  IRet = sprintf(c,"\\usepackage{rotating}\n"); ofs<<c;
  IRet = sprintf(c,"\\setlength{\\textheight} {23.5 true cm} \n"); ofs<<c;
  IRet = sprintf(c,"\\setlength{\\textwidth} {17 true cm} \n"); ofs<<c;
  IRet = sprintf(c,"\\setlength{\\oddsidemargin} {0 mm} \n"); ofs<<c;
  IRet = sprintf(c,"\\setlength{\\evensidemargin} {0 mm} \n"); ofs<<c;

  // Start document
  IRet = sprintf(c,"\\begin{document}\n"); ofs<<c;

  // Write some header text
  IRet = sprintf(c," \\section{Results} \n"); ofs<<c;
  IRet = sprintf(c," \\label{sec:result} \n"); ofs<<c;
  IRet = sprintf(c," The name of the file was %s \n",FilNam.Data()); ofs<<c;
  IRet = sprintf(c," The results have been obtained with the BLUE Software Version %s \n",
		 Versio.Data()); ofs<<c;
  IRet = sprintf(c," on %s \n",Today.Data()); ofs<<c;
  IRet = sprintf(c," \\par \n"); ofs<<c;

  // Print the Table of estimates and observables
  IRet = sprintf(c,"%%\n"); ofs<<c;
  IRet = sprintf(c,"%% -------------------------------------------------------\n"); ofs<<c;
  IRet = sprintf(c,"%%\n"); ofs<<c;
  if(Isidew == 1){IRet = sprintf(c,"\\begin{sidewaystable}[tbp!]\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\begin{table}[tbp!]\n"); ofs<<c;};
  IRet = sprintf(c,"\\begin{center}\n"); ofs<<c;
  if(NumCol > 15){IRet = sprintf(c,"\\small\n"); ofs<<c;};

  // The tabular statement
  IRet = sprintf(c,"\\begin{tabular}{|l|"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,"r|"); ofs<<c;
  }
  IRet = sprintf(c,"|"); ofs<<c;
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c,"r|"); ofs<<c;
  }
  IRet = sprintf(c,"}\n"); ofs<<c;

  // The first line
  IRet = sprintf(c,"\\cline{2-%2i}\n",NumCol); ofs<<c;
  IRet = sprintf(c,"  \\multicolumn{1}{c}{} \n"); ofs<<c;
  IRet = sprintf(c,"& \\multicolumn{%2i}{|c||}{Estimates} \n",InpEst); ofs<<c;
  if(InpObs > 1){IRet = sprintf(c,"& \\multicolumn{%2i}{|c|}{Observables} \\\\ \n",InpObs);
  }else{IRet = sprintf(c,"& \\multicolumn{%2i}{|c|}{Observable} \\\\ \n",InpObs);
  }
  ofs<<c;
  IRet = sprintf(c,"\\hline\n"); ofs<<c;

  // The Estimate line
  //  IRet = sprintf(c,"Estimate  "); ofs<<c;
  IRet = sprintf(c,"Uncertainty"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c," & %s",GetNamEst(IsWhichEst(i)).Data()); ofs<<c;
  }
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c," & %s",GetNamObs(IsWhichObs(n)).Data()); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \n"); ofs<<c;

  if(InpObs > 1){
    // The Observable line
    IRet = sprintf(c,"Observable"); ofs<<c;
    for(Int_t i = 0; i<InpEst; i++){
      IRet = sprintf(c," & %s",GetNamObs(EstWhichObs(IsWhichEst(i))).Data()); ofs<<c;
    }
    for(Int_t n = 0; n<InpObs; n++){
      IRet = sprintf(c," &  "); ofs<<c;
    }
    IRet = sprintf(c," \\\\ \n"); ofs<<c;
  }
  IRet = sprintf(c,"\\hline\n"); ofs<<c;

  // The Values
  IRet = sprintf(c,"Values "); ofs<<c;
  Format = " & " + ForVal;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,Format,Xva->operator()(i)); ofs<<c;
  }
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c,Format,XvaRes->operator()(n)); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \n"); ofs<<c;

  // The Stat Unc k=0 
  // The Syst Unc k>0
  Format = " & " + ForUnc;
  for(Int_t k = 0; k<InpUnc; k++){
    IRet = sprintf(c,"%s",GetNamUnc(IsWhichUnc(k)).Data()); ofs<<c;
    for(Int_t i = 0; i<InpEst; i++){
      Sys = Unc->operator()(i+k*InpEst,i+k*InpEst);
      if(k > 0)SysEst->operator()(i) = SysEst->operator()(i) + Sys*Sys;
      IRet = sprintf(c,Format,Sys); ofs<<c;
    }
    for(Int_t n = 0; n<InpObs; n++){
      Sys = CorRes->operator()(k*InpObs+n,k*InpObs+n);
      if(k > 0)SysObs->operator()(n) = SysObs->operator()(n) + Sys*Sys;
      IRet = sprintf(c,Format,Sys); ofs<<c;
    }
    IRet = sprintf(c," \\\\ \n"); ofs<<c;
  }
  IRet = sprintf(c,"\\hline\n"); ofs<<c;

  // The Total Syst 
  IRet = sprintf(c,"Total Syst"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,Format,TMath::Sqrt(SysEst->operator()(i))); ofs<<c;
  }
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c,Format,TMath::Sqrt(SysObs->operator()(n))); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \n"); ofs<<c;

  // The Full Unc
  IRet = sprintf(c,"Full      "); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,Format,Sig->operator()(i)); ofs<<c;
  }
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c,Format,sqrt(CovRes->operator()(n,n))); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \n"); ofs<<c;
  IRet = sprintf(c,"\\hline\n"); ofs<<c;

  // Finish up the table
  IRet = sprintf(c,"\\end{tabular}\n"); ofs<<c;
  IRet = sprintf(c,"\\end{center}\n"); ofs<<c;
  if(InpObs>1){IRet = sprintf(c,"\\caption{Combination of %2i observables", InpObs);
  }else{IRet = sprintf(c,"\\caption{Combination of %2i observable", InpObs);
  }
  ofs<<c;
  IRet = sprintf(c," from %2i correlated estimates",InpEst); ofs<<c;
  IRet = sprintf(c," using the BLUE software.}\n"); ofs<<c;
  IRet = sprintf(c,"\\label{tab:BlueRes}.\n"); ofs<<c;
  if(Isidew == 1){IRet = sprintf(c,"\\end{sidewaystable}\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\end{table}\n"); ofs<<c;};

  // --------------------------------------------------------------------------------------------

  // Print the correlation matrix of the estimates
  NumCol = InpEst + 1;
  if(NumCol > NumMax)Isidew = 1;
  IRet = sprintf(c,"%%\n"); ofs<<c;
  IRet = sprintf(c,"%% -------------------------------------------------------\n"); ofs<<c;
  IRet = sprintf(c,"%%\n"); ofs<<c;

  if(Isidew == 1){IRet = sprintf(c,"\\begin{sidewaystable}[tbp!]\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\begin{table}[tbp!]\n"); ofs<<c;};
  IRet = sprintf(c,"\\begin{center}\n"); ofs<<c;

  // The tabular statement
  IRet = sprintf(c,"\\begin{tabular}{|l|"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,"r|"); ofs<<c;
  }
  IRet = sprintf(c,"}\\hline \n"); ofs<<c;

  // The first line
  IRet = sprintf(c,"Estimates"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c," & \\rotatebox{90}{%s \\,}",GetNamEst(IsWhichEst(i)).Data()); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \\hline\n"); ofs<<c;

  // The Matrix
  Format = " & " + ForRho;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,"%s",GetNamEst(IsWhichEst(i)).Data()); ofs<<c;
    for(Int_t j = 0; j<InpEst; j++){
      if(j <= i){IRet = sprintf(c,Format,Rho->operator()(i,j)); ofs<<c;
      }else if(j < InpEst-1){sprintf(c," &      "); ofs<<c;
      }else{IRet = sprintf(c," &       \\\\\n"); ofs<<c;
      }
    }
  }
  IRet = sprintf(c," \\\\ \\hline \n"); ofs<<c;

  // Finish up the table
  IRet = sprintf(c,"\\end{tabular}\n"); ofs<<c;
  IRet = sprintf(c,"\\end{center}\n"); ofs<<c;
  IRet = sprintf(c,"\\caption{Correlation matrix for the %2i correlated estimates",InpEst); ofs<<c;
  IRet = sprintf(c," using the BLUE software.}\n"); ofs<<c;
  IRet = sprintf(c,"\\label{tab:BlueCor}.\n"); ofs<<c;
  if(Isidew == 1){IRet = sprintf(c,"\\end{sidewaystable}\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\end{table}\n"); ofs<<c;};

  // --------------------------------------------------------------------------------------------

  // Print the blue weights
  NumCol = InpEst + 1;
  if(NumCol > NumMax)Isidew = 1;
  IRet = sprintf(c,"%%\n"); ofs<<c;
  IRet = sprintf(c,"%% -------------------------------------------------------\n"); ofs<<c;
  IRet = sprintf(c,"%%\n"); ofs<<c;

  if(Isidew == 1){IRet = sprintf(c,"\\begin{sidewaystable}[tbp!]\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\begin{table}[tbp!]\n"); ofs<<c;};
  IRet = sprintf(c,"\\begin{center}\n"); ofs<<c;

  // The tabular statement
  IRet = sprintf(c,"\\begin{tabular}{|l|"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c,"r|"); ofs<<c;
  }
  IRet = sprintf(c,"}\\hline \n"); ofs<<c;

  // The first line
  IRet = sprintf(c,"Estimates  "); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c," & %s",GetNamEst(IsWhichEst(i)).Data()); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \n"); ofs<<c;

  // The second line
  IRet = sprintf(c,"Observables"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    IRet = sprintf(c," & %s",GetNamObs(EstWhichObs(IsWhichEst(i))).Data()); ofs<<c;    
  }
  IRet = sprintf(c," \\\\ \\hline\n"); ofs<<c;

  // The Weight Matrix
  Format = " & " + ForWei;
  for(Int_t n = 0; n<InpObs; n++){
    IRet = sprintf(c,"Weights (%s)",GetNamObs(IsWhichObs(n)).Data()); ofs<<c;
    for(Int_t i = 0; i<InpEst; i++){
      IRet = sprintf(c,Format,Lam->operator()(i,n)); ofs<<c;
    }
    IRet = sprintf(c," \\\\\n"); ofs<<c;
  }
  IRet = sprintf(c," \\hline\n"); ofs<<c;

  // The pull
  Format = " & " + ForPul;
  IRet = sprintf(c,"Pull"); ofs<<c;
  for(Int_t i = 0; i<InpEst; i++){
    //201213 IRet = sprintf(c,Format,Pul->operator()(IsWhichEst(i))); ofs<<c;
    IRet = sprintf(c,Format,Pul->operator()(i)); ofs<<c;
  }
  IRet = sprintf(c," \\\\ \\hline \n"); ofs<<c;

  // Finish up the table
  IRet = sprintf(c,"\\end{tabular}\n"); ofs<<c;
  IRet = sprintf(c,"\\end{center}\n"); ofs<<c;
  IRet = sprintf(c,"\\caption{Weights and pulls of %2i correlated estimates \n",InpEst); ofs<<c;
  if(InpObs >1){IRet = sprintf(c,"for %2i observables",InpObs);
  }else{IRet = sprintf(c,"for %2i observable",InpObs);
  }
  ofs<<c;
  IRet = sprintf(c," using the BLUE software.}\n"); ofs<<c;
  IRet = sprintf(c,"\\label{tab:BlueWeight}.\n"); ofs<<c;
  if(Isidew == 1){IRet = sprintf(c,"\\end{sidewaystable}\n"); ofs<<c;
  }else{IRet = sprintf(c,"\\end{table}\n"); ofs<<c;};

  // --------------------------------------------------------------------------------------------

  // Print the correlation matrix of the observables
  if(InpObs > 1){
    NumCol = InpObs + 1;
    if(NumCol > NumMax)Isidew = 1;
    IRet = sprintf(c,"%%\n"); ofs<<c;
    IRet = sprintf(c,"%% -------------------------------------------------------\n"); ofs<<c;
    IRet = sprintf(c,"%%\n"); ofs<<c;
    
    if(Isidew == 1){IRet = sprintf(c,"\\begin{sidewaystable}[tbp!]\n"); ofs<<c;
    }else{IRet = sprintf(c,"\\begin{table}[tbp!]\n"); ofs<<c;};
    IRet = sprintf(c,"\\begin{center}\n"); ofs<<c;
    
    // The tabular statement
    IRet = sprintf(c,"\\begin{tabular}{|l|"); ofs<<c;
    for(Int_t n = 0; n<InpObs; n++){
      IRet = sprintf(c,"r|"); ofs<<c;
    }
    IRet = sprintf(c,"}\\hline \n"); ofs<<c;
    
    // The first line
    IRet = sprintf(c,"Observables"); ofs<<c;
    for(Int_t n = 0; n<InpObs; n++){
      IRet = sprintf(c," & \\rotatebox{90}{%s \\,}",GetNamObs(IsWhichObs(n)).Data()); ofs<<c;
    }
    IRet = sprintf(c," \\\\ \\hline\n"); ofs<<c;
    
    // The Matrix
    // RhoRes->Print();
    Format = " & " + ForRho;
    for(Int_t n = 0; n<InpObs; n++){
      IRet = sprintf(c,"%s",GetNamObs(IsWhichObs(n)).Data()); ofs<<c;
      for(Int_t m = 0; m<InpObs; m++){
	if(m <= n){IRet = sprintf(c,Format,RhoRes->operator()(n,m)); ofs<<c;
	}else if(m < InpObs-1){IRet = sprintf(c," &      "); ofs<<c;
	}else{IRet = sprintf(c," &       \\\\\n"); ofs<<c;
	}
      }
    }
    IRet = sprintf(c," \\\\ \\hline \n"); ofs<<c;

    // Finish up the table
    IRet = sprintf(c,"\\end{tabular}\n"); ofs<<c;
    IRet = sprintf(c,"\\end{center}\n"); ofs<<c;
    IRet = sprintf(c,"\\caption{Correlation matrix for the %2i",InpObs); ofs<<c;
    IRet = sprintf(c," correlated observables using the BLUE software.}\n"); ofs<<c;
    IRet = sprintf(c,"\\label{tab:BlueCorRes}.\n"); ofs<<c;
    if(Isidew == 1){IRet = sprintf(c,"\\end{sidewaystable}\n"); ofs<<c;
    }else{IRet = sprintf(c,"\\end{table}\n"); ofs<<c;};
  }

  // End of tex file
  IRet = sprintf(c,"\\end{document}\n"); ofs<<c;

  // Close file
  ofs.close();

  printf("... Blue->LatexResult(): The Latex File %s has been produced. \n", 
	 FilNam.Data());

  // Clean up
  SysEst->Delete(); SysEst = NULL;
  SysObs->Delete(); SysObs = NULL;

  return;
};

//------------------------------------------------------------------------------

void Blue::DisplayResult(const Int_t n, const TString FilNam) const {
  TString ForVal = "%5.2f";
  TString ForUnc = ForVal;
  DisplayResult(n, FilNam, ForVal, ForUnc);
  return;
};

//------------------------------------------------------------------------------

void Blue::DisplayResult(const Int_t n, const TString FilNam, 
			 const TString ForVal, const TString ForUnc) const {
  if(IsSolved() == 0){
    printf("... Blue->DisplayResult(%2i): Presently not available, call Solve() \n", n);
    return;
  }

  // Check whether n is an active observable
  Int_t na = 0;
  if(IsActiveObs(n) == 1){
    printf("... Blue->DisplayResult(%2i): The next observable is: %s \n", 
	   n, GetNamObs(n).Data());
    na = IsIndexObs(n);
  }else{
    printf("... Blue->DisplayResult(%2i): Not an active observable \n", n);
    return;
  }

  // Reset the arrays for PlotRes
  for(Int_t i = 0; i<InpEst+1; i++){
    Indx[i] = 0;
    Colo[i] = 0;
    Name[i] = "";
    Valu[i] = 0.;
    Stat[i] = 0.;
    Syst[i] = 0.;
    Full[i] = 0.;
  }

  // Declare variables
  TVectorD* SysEst = new TVectorD(InpEst);
  TVectorD* SysObs = new TVectorD(InpObs);
  Double_t Sys = 0;
  Int_t ind = 0;

  // The dimension
  Int_t N = 1;
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      N = N + 1;
    }
  }

  // The Indices
  for(Int_t i = 0; i<N; i++)Indx[i] = i;

  // The Colors
  Colo[0] = kRed;
  for(Int_t i = 1; i<N; i++)Colo[i] = kBlack;

  // The Names
  ind = 0;
  Name[ind] = GetNamObs(n);
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      ind = ind + 1;
      Name[ind] = GetNamEst(IsWhichEst(i));
    }
  }

  // The Values
  ind = 0;
  Valu[ind] = XvaRes->operator()(na);
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      ind = ind + 1;
      Valu[ind] = Xva->operator()(i);
    }
  }
  
  // The Stat uncertainties
  ind = 0;
  Stat[ind] = CorRes->operator()(0*InpObs+na,0*InpObs+na);
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      ind = ind + 1;
      Stat[ind] = Unc->operator()(i+0*InpEst,i+0*InpEst);
    }
  }

  // Calculate the syst unc k>0
  for(Int_t k = 1; k<InpUnc; k++){
    for(Int_t i = 0; i<InpEst; i++){
      if(EstWhichObs(IsWhichEst(i)) == n){
	Sys = Unc->operator()(i+k*InpEst,i+k*InpEst);
	SysEst->operator()(i) = SysEst->operator()(i) + Sys*Sys;
      }
    }
    Sys = CorRes->operator()(k*InpObs+na,k*InpObs+na);
    SysObs->operator()(na) = SysObs->operator()(na) + Sys*Sys;
  }

  // Fill the syst 
  ind = 0;
  Syst[ind] = TMath::Sqrt(SysObs->operator()(na));
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      ind = ind + 1;
      Syst[ind] = TMath::Sqrt(SysEst->operator()(i));
    }
  }

  // The Full Unc
  ind = 0;
  Full[ind] = TMath::Sqrt(CovRes->operator()(na,na));
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      ind = ind + 1;
      Full[ind] = Sig->operator()(i);
    }
  }

  // Prepare filename
  char ObsNum[8];
  sprintf(ObsNum,"%i",n);
  TString FulNam = FilNam + "_DisRes_Obs_" + ObsNum;

  // Produce the file
  PlotRes(N, FulNam, ForVal, ForUnc);

  // Clean up
  SysEst->Delete(); SysEst = NULL;
  SysObs->Delete(); SysObs = NULL;
  return;
};

//------------------------------------------------------------------------------
// After SolveXXX()

void Blue::DisplayAccImp(const Int_t n, const TString FilNam) const {
  TString ForVal = "%5.2f";
  TString ForUnc = ForVal;
  DisplayAccImp(n, FilNam, ForVal, ForUnc);
  return;
};

//------------------------------------------------------------------------------

void Blue::DisplayAccImp(const Int_t n, const TString FilNam, 
			 const TString ForVal, const TString ForUnc) const {

  if(IsSolvedAccImp() == 1){
    if(IsPrintLevel() >= 1){
      printf("... DisplayAccImp(%s): \n", FilNam.Data());
      printf("... DisplayAccImp(): LasImp \n"); LasImp->Print();
      printf("... DisplayAccImp(): IndImp \n"); IndImp->Print();
      printf("... DisplayAccImp(): ValImp \n"); ValImp->Print();
      printf("... DisplayAccImp(): UncImp \n"); UncImp->Print();
      printf("... DisplayAccImp(): StaImp \n"); StaImp->Print();
      printf("... DisplayAccImp(): SysImp \n"); SysImp->Print();
    }

    // Check whether n is an active observable
    Int_t na = 0;
    if(IsActiveObs(n) == 1){
      printf("... Blue->DisplayAccImp(%2i): The next observable is: %s \n", 
	     n, GetNamObs(n).Data());
      na = IsIndexObs(n);
    }else{
      printf("... Blue->DisplayAccImp(%2i): Not an active observable \n", n);
      return;
    }
    
    // Reset the arrays for PlotRes
    for(Int_t i = 0; i<InpEst+1; i++){
      Indx[i] = 0;
      Colo[i] = 0;
      Name[i] = "";
      Valu[i] = 0.;
      Stat[i] = 0.;
      Syst[i] = 0.;
      Full[i] = 0.;
    }
    
    // Get the dimension and the index of the last one to be combined for the
    // wanted precision
    Int_t NumCom = 0, LasCom = 0;
    for(Int_t i = 0; i<InpEst; i++){
      if(IndImp->operator()(i,na) > -0.5)NumCom = NumCom+1;
      if(IndImp->operator()(i,na) == LasImp->operator()(0,na))LasCom = i;
    }
    //printf("... Blue->DisplayAccImp(%2i): NumCom = %2i, LasCom = %2i\n", 
    //	   n, NumCom,LasCom);

    // Fill the values
    TString ComNam = "";
    Int_t ind = NumCom - 1;
    for(Int_t i = 0; i<NumCom; i++){
      Indx[i] = i;
      Colo[ind] = kBlack;
      if(i == LasCom)Colo[ind] = kRed;

      if(i > 0)ComNam = "+";
      ComNam = ComNam + GetNamEst(static_cast<int>(IndImp->operator()(i,na)));
      Name[ind] = ComNam.Data();

      Valu[ind] = ValImp->operator()(i,na);
      Stat[ind] = StaImp->operator()(i,na);
      Syst[ind] = SysImp->operator()(i,na);
      Full[ind] = UncImp->operator()(i,na);
      ind = ind - 1;
    }

    // Prepare the filename
    char ObsNum[8];
    sprintf(ObsNum,"%i",n);
    TString FulNam = FilNam + "_AccImp_Obs_" + ObsNum;
    
    // Produce the file
    PlotRes(NumCom, FulNam, ForVal, ForUnc);    
  }else{
    printf("... Blue->DisplayAccImp(%2i): Presently not available, call SolveAccImp() \n",n);
  }
  return;
};

// ---- Now the private member functions

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

void Blue::Construct(const Int_t NumEst, const Int_t NumUnc, const Int_t NumObs, 
		     const Int_t *IntObs, const Int_t *IntFac){

  // Set the version number
  Versio   = "2.1.0";

  // Set print level to low
  SetPrintLevel(0);

  // Set the date
  Date   = new TDatime();
  Months = new TString[12];
  Months[0] = "Jan"; Months[1] = "Feb"; Months[ 2] = "Mar"; Months[ 3] = "Apr";
  Months[4] = "May"; Months[5] = "Jun"; Months[ 6] = "Jul"; Months[ 7] = "Aug";
  Months[8] = "Sep"; Months[9] = "Oct"; Months[10] = "Nov"; Months[11] = "Dec";

  char Buffer[100];
  sprintf(Buffer, "%s %2i, %4i", Months[Date->GetMonth()-1].Data(),
	   Date->GetDay(), Date->GetYear());
  Today = &Buffer[0]; 

  // Set actual dimensions
  InpEst = NumEst;
  InpUnc = NumUnc;
  InpObs = NumObs;

  // Save original dimensions
  InpEstOrig = InpEst;
  InpUncOrig = InpUnc;
  InpObsOrig = InpObs;

  // Set names of estimates, uncertainties and observables and filling flags
  TString Format = "To be filled later";
  EstNam = new TString[InpEstOrig];
  UncNam = new TString[InpUncOrig];
  ObsNam = new TString[InpObsOrig];
  Format = "Est %3i";
  for(Int_t i = 0; i<InpEstOrig; i++){
    //if(i>=10)Format = "Est %2i";
    sprintf(Buffer,Format,i);
    EstNam[i] = &Buffer[0];
  }
  UncNam[0] = "   Stat";
  Format = "Syst%3i";
  for(Int_t k = 1; k<InpUncOrig; k++){
    //if(k>=10)Format = "Syst %2i";
    sprintf(Buffer,Format,k);
    UncNam[k] = &Buffer[0];
  }
  Format = "Obs %3i";
  for(Int_t n = 0; n<InpObsOrig; n++){
    //if(n>=10)Format = "Obs %2i";
    sprintf(Buffer,Format,n);
    ObsNam[n] = &Buffer[0];
  }
  IsFillEstNam = 0;
  IsFillUncNam = 0;
  IsFillObsNam = 0;

  // Save harbour for names
  EstNamOrig = new TString[InpEstOrig];
  UncNamOrig = new TString[InpUncOrig];
  ObsNamOrig = new TString[InpObsOrig];

  // Set which estimates determines which observable
  EstObs = new Int_t[InpEstOrig];
  for(Int_t i = 0; i<InpEstOrig; i++){
    if(InpObsOrig == 1){
      EstObs[i] = 0;
    }else{
      // Check consistency
      EstObs[i] = *(IntObs+i);
      if(EstObs[i] < 0 ){
	printf("... Blue->Blue: Inconsistent input negative observable %2i --> stop \n", EstObs[i]);
	printf("... Blue->Blue: Constructor FAILED! \n");
	return;
      }
      if(EstObs[i] >= InpObsOrig){
	printf("... Blue->Blue: Inconsistent input only %2i observables", InpObsOrig);
	printf(" but estimate %2i determines observable %2i \n", i,EstObs[i]);
	printf("... Blue->Blue: Constructor FAILED! \n");
	return;
      }
    }
    // printf("... Blue->Blue: Next %2i,\n", EstObs[i]);
  }

  // Set the matrix of scale factors
  InpFac = 0;
  MatFac = new TMatrixD(InpEstOrig,InpEstOrig);
  Int_t IntTes = 0;
  for(Int_t i = 0; i<InpEstOrig; i++){
    for(Int_t j = 0; j<InpEstOrig; j++){
      IntTes = *(IntFac+i*InpEstOrig+j);
      if(i != j){
	// Check off-diagonal < 0, do not care about diagonal
	if(IntTes < 0){
	  printf("... Blue->Blue: Inconsistent input,");
	  printf(" negative value in IWhichFac = %2i \n", IntTes);
	  printf("... Blue->Blue: Constructor FAILED! \n");
	  return;
	}
	// Get maximum value
	if(IntTes > InpFac)InpFac = IntTes;
      }
      MatFac->operator()(i,j) = static_cast<double>(IntTes);
    }
  }
  InpFac = InpFac + 1;
  //printf("... Blue->Blue: MatFac \n");  MatFac->Print();

  // Check for symmetry
  if(IsWhichMatrix(MatFac) == 3){
    printf("... Blue->Blue: Inconsistent input, IWhichFac is not symmetric \n");
    printf("... Blue->Blue: Constructor FAILED! \n");
    return;
  }

  // Check for consecutive values up to InpFac-1
  Int_t DumFac[InpFac];
  for(Int_t l = 0; l<InpFac; l++){DumFac[l] = 0;};
  for(Int_t i = 0; i<InpEstOrig; i++){
    for(Int_t j = i+1; j<InpEstOrig; j++){
      //printf("... Blue->Blue: %5.2f \n", MatFac->operator()(i,j));
      DumFac[static_cast<int>(MatFac->operator()(i,j))] = 1;
    }
  }
  for(Int_t l = 0; l<InpFac; l++){
    //printf("... Blue->Blue: DumFac %2i \n", DumFac[l]);
    if(DumFac[l] == 0){
      printf("... Blue->Blue: Inconsistent input, IWhichFac does not have");
      printf(" consecutive values, %2i is missing. \n", l);
      printf("... Blue->Blue: Constructor FAILED! \n");
      return;
    } 
  }

  // Set the actual factors and ranges for SolveScaRho
  ActFac = new TMatrixD(InpUncOrig,InpFac);
  MinFac = new TMatrixD(InpUncOrig,InpFac);
  MaxFac = new TMatrixD(InpUncOrig,InpFac);
  ResetScaRho(2);
  FlaFac = 0;
  FaiFac = 0;

  // Set the number of scan values and the result matrices
  InrFac = 10;
  ValSca = new TMatrixD* [InpObsOrig];
  SigSca = new TMatrixD* [InpObsOrig];
  VtoSca = new TMatrixD* [InpObsOrig];
  StoSca = new TMatrixD* [InpObsOrig];
  for(Int_t n = 0; n<InpObsOrig; n++){
    ValSca[n] = new TMatrixD(InpUncOrig*InpFac,InrFac);
    SigSca[n] = new TMatrixD(InpUncOrig*InpFac,InrFac);
    VtoSca[n] = new TMatrixD(InpFac,InrFac);
    StoSca[n] = new TMatrixD(InpFac,InrFac);
  }

  // Fill lists of active estimates and uncertainties
  LisEst = new Int_t[InpEstOrig];
  LisUnc = new Int_t[InpUncOrig];
  for(Int_t i = 0; i<InpEstOrig; i++){LisEst[i]=i;};
  for(Int_t k = 0; k<InpUncOrig; k++){LisUnc[k]=k;};

  // Set up U-Matrices for the combination
  Uma = new TMatrixD(InpEstOrig,InpObsOrig);
  Utr = new TMatrixD(InpObsOrig,InpEstOrig);
  for(Int_t i = 0; i<InpEstOrig; i++){
    for(Int_t n = 0; n<InpObsOrig; n++){
      if(n == EstObs[i])Uma->operator()(i,n) = 1;
    }
  }
  Utr->Transpose(*Uma);
  //Uma->Print();
  //Utr->Print();

  // Reset arrays for controlling the estimates
  EstAct = new Int_t[InpEstOrig];
  EstFil = new Int_t[InpEstOrig];
  for (Int_t i = 0; i < InpEstOrig; i++){
    EstAct[i] = 0;
    EstFil[i] = 0;
  }
  UncAct = new Int_t[InpUncOrig];
  UncFil = new Int_t[InpUncOrig];
  UncCha = new Int_t[InpUncOrig];
  UncFac = new Int_t[InpUncOrig];
  UncRed = new Int_t[InpUncOrig];
  UncRel = new Int_t[InpUncOrig];
  for (Int_t k = 0; k < InpUncOrig; k++){
    UncAct[k] = 0;
    UncFil[k] = 0;
    UncCha[k] = 0;
    UncFac[k] = 0;
    UncRed[k] = 0;
    UncRel[k] = 0;
  }

  // Fill active observables and list of observables
  LisObs = new Int_t[InpObsOrig];
  ObsAct = new Int_t[InpObsOrig];
  for (Int_t n = 0; n < InpObsOrig; n++){
    ObsAct[n] = 1;
    LisObs[n] = n;
  }

  // Reset global steering variables
  InpFil = 0;
  SetFixedInp(0);
  SetIsSolved(0);
  SetIsSolvedRelUnc(0);
  SetIsSolvedAccImp(0);
  SetIsSolvedScaRho(0);
  SetIsSolvedInfWei(0);
  SetIsSolvedPosWei(0);
  SetIsSolvedMaxVar(0);

  IsRelUncMode  = 0;

  // Array of coefficients forBLUE with relative uncertainties
  Cof  = new TMatrixD(InpEstOrig*InpUncOrig,MaxCof);

  // Reset list of importance and the steering flag for SolveAccImp
  PreAcc = 1.;
  IntAcc = 0;
  LisImp = new Int_t[InpEstOrig];
  for (Int_t i = 0; i < InpEstOrig; i++){LisImp[i] = -1;};
  LasImp = new TMatrixD(1,InpObsOrig);
  IndImp = new TMatrixD(InpEstOrig,InpObsOrig);
  ValImp = new TMatrixD(InpEstOrig,InpObsOrig);
  UncImp = new TMatrixD(InpEstOrig,InpObsOrig);
  StaImp = new TMatrixD(InpEstOrig,InpObsOrig);
  SysImp = new TMatrixD(InpEstOrig,InpObsOrig);

  // The arrays for PlotRes
  Indx = new Int_t[InpEstOrig+1];
  Colo = new Int_t[InpEstOrig+1];
  Name = new TString[InpEstOrig+1];
  Valu = new Double_t[InpEstOrig+1];
  Stat = new Double_t[InpEstOrig+1];
  Syst = new Double_t[InpEstOrig+1];
  Full = new Double_t[InpEstOrig+1];
  for(Int_t i = 0; i<InpEstOrig+1; i++){
    Indx[i] = 0;
    Colo[i] = 0;
    Name[i] = "";
    Valu[i] = 0.;
    Stat[i] = 0.;
    Syst[i] = 0.;
    Full[i] = 0.;
  }

  // The list of initial/final variances/correlations for SolveMaxVar
  VarMax = new TVectorD(InpObsOrig);
  VarMin = new TVectorD(InpObsOrig);
  RhoMax = new TMatrixD(InpEstOrig,InpEstOrig);
  RhoMin = new TMatrixD(InpEstOrig,InpEstOrig);
  RhoOoz = new TMatrixD(InpEstOrig,InpEstOrig);

  RhoFco = 1.;
  RhoFpk = new TVectorD(InpUncOrig);
  RhoFij = new TMatrixD(InpEstOrig,InpEstOrig);
  IFuCor = 0;
  IFuFla = 0;
  IFuFai = new Int_t[3];

  // Save harbor for all inputs
  XvaOrig = new TVectorD(InpEstOrig);
  SigOrig = new TVectorD(InpEstOrig);
  UncOrig = new TMatrixD(InpEstOrig*InpUncOrig,InpEstOrig*InpUncOrig);
  CorOrig = new TMatrixD(InpEstOrig*InpUncOrig,InpEstOrig*InpUncOrig);
  UmaOrig = new TMatrixD(InpEstOrig,InpObsOrig);
  UtrOrig = new TMatrixD(InpObsOrig,InpEstOrig);

  // Vectors to do the job
  Xva = new TVectorD(InpEstOrig);
  Sig = new TVectorD(InpEstOrig);

  // Matrices to do the job
  Unc  = new TMatrixD(InpEstOrig*InpUncOrig,InpEstOrig*InpUncOrig);
  Cor  = new TMatrixD(InpEstOrig*InpUncOrig,InpEstOrig*InpUncOrig);
  Cov  = new TMatrixD(InpEstOrig,InpEstOrig);
  CovI = new TMatrixD(InpEstOrig,InpEstOrig);
  Rho  = new TMatrixD(InpEstOrig,InpEstOrig);
  Lam  = new TMatrixD(InpEstOrig,InpObsOrig);
  Pul  = new TVectorD(InpEstOrig);

  // Matrices for the results
  XvaRes = new TVectorD(InpObsOrig);
  CorRes = new TMatrixD(InpObsOrig*InpUncOrig,InpObsOrig*InpUncOrig);
  CovRes = new TMatrixD(InpObsOrig,InpObsOrig);
  RhoRes = new TMatrixD(InpObsOrig,InpObsOrig);

  // Chiq information of the results
  ChiQua = 0;
  NumDof = 0;
  ChiPro = 0;

  // Vectors for the information weights
  VarInd = new TVectorD(InpEstOrig);
  IntWei = new TVectorD(InpEstOrig+1);
  MarWei = new TVectorD(InpEstOrig);
  BluWei = new TVectorD(InpEstOrig);

  // Matrices for the parameters
  SRat = new TMatrixD(InpEstOrig,InpEstOrig);
  Beta = new TMatrixD(InpEstOrig,InpEstOrig);
  Sigx = new TMatrixD(InpEstOrig,InpEstOrig);
  DBdr = new TMatrixD(InpEstOrig,InpEstOrig);
  DSdr = new TMatrixD(InpEstOrig,InpEstOrig);
  DBdz = new TMatrixD(InpEstOrig,InpEstOrig);
  DSdz = new TMatrixD(InpEstOrig,InpEstOrig);

  // Reset the flag
  SetCalcedParams(0);

  // Matrix for control of the relative uncertainties 
  IndRel = new TMatrixD(InpEstOrig,InpUncOrig);

  // Matrices for the InspectLike results
  // ChiRes holds Chiq and det(V) per point
  ChiRes = new TMatrixD(2,1);
  // LikRes holds result x, dxlow, dxhig
  LikRes = new TMatrixD(InpObsOrig,LikDim);
  LikFla = 0;
  SetIsInspectLike(0);

  // Switch off quiet mode
  SetQuiet(0);

  //Setup ROOT style
  SetupRoot();
};

//------------------------------------------------------------------------------

void Blue::SetupRoot(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
};

//------------------------------------------------------------------------------
// Solver
//------------------------------------------------------------------------------

void Blue::SolveScaRho(const Int_t FlaFac, const TMatrixD* MinFac, const TMatrixD* MaxFac){

  if(IsFilledInp() == 0){
    printf("... Blue->SolveScaRho(): IGNORED Not all estimates/uncertainties filled IGNORED \n");
    return;
  }else if(IsFixedInp() == 0){
    printf("... Blue->SolveScaRho(): IGNORED Input not yet fixed call FixInp \n");
    return;
  }else if(IsSolvedScaRho() == 1){
    printf("... Blue->SolveScaRho(): IGNORED Nothing changed since last call to this function.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }else if(IsSolved() == 1){
    printf("... Blue->SolveScaRho(): IGNORED Two consecutive calls to Solve() are not supported.");
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Check for changed or reduced correlations
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      if(IsRhoValUnc(k) == 1){
	printf("... Blue->SolveScaRho(%2i): IGNORED I found a changed correlation. Please fix. \n", FlaFac,k);    
	return;
      }
      if(IsRhoRedUnc(k) == 1){
	printf("... Blue->SolveScaRho(%2i): IGNORED I found a reduced correlation. Please fix. \n", FlaFac,k);    
	return;
      }
    }
  }

  // Check for values != 1 and save the actual factors
  Int_t IFound = 0;
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      for(Int_t l = 0; l<InpFac; l++){	
	if(ActFac->operator()(k,l) != 1)IFound = 1; 
      }
    }
  }
  if(IFound == 1){
    printf("... Blue->SolveScaRho(%2i): WARNING: At least one factor initially is not 1. \n", FlaFac);
    printf("... Blue->SolveScaRho(%2i): WARNING: If this is not what you want reset to 1. \n", FlaFac);
  }
  TMatrixD* SavFac = new TMatrixD(InpUncOrig,InpFac);
  SavFac->SetSub(0,0,*ActFac);

  // Reset the failure counter
  FaiFac = 0;

  // Set the quiet modus
  SetQuiet(1);

  // Solve once store the result
  Solve();
  TVectorD* ValDef = new TVectorD(InpObs);
  TVectorD* SigDef = new TVectorD(InpObs);
  for(Int_t n = 0; n<InpObs; n++){
    ValDef->operator()(n) = XvaRes->operator()(n);
    SigDef->operator()(n) = TMath::Sqrt(CovRes->operator()(n,n));
    //printf("... Blue->SolveScaRho(): Next Obs n = %2i = %5.3f +- %5.3f \n",
    //IsWhichObs(n), ValDef->operator()(n), SigDef->operator()(n));
  }
  //printf("... Blue->SolveScaRho(): ActFac\n"); PrintMatrix(ActFac);
  //printf("... Blue->SolveScaRho(): MinFac\n"); PrintMatrix(MinFac);
  //printf("... Blue->SolveScaRho(): MaxFac\n"); PrintMatrix(MaxFac);

  // Set the matrices for the results of the scan
  for(Int_t n = 0; n<InpObs; n++){  
    ValSca[n]->Delete(); ValSca[n] = new TMatrixD(InpUnc*InpFac,InrFac);
    SigSca[n]->Delete(); SigSca[n] = new TMatrixD(InpUnc*InpFac,InrFac);
    VtoSca[n]->Delete(); VtoSca[n] = new TMatrixD(InpFac,InrFac);
    StoSca[n]->Delete(); StoSca[n] = new TMatrixD(InpFac,InrFac);
  }

  // Set the first value and the step size
  TMatrixD* RhoLoo = new TMatrixD(InpUncOrig,InpFac);
  TMatrixD* RhoSte = new TMatrixD(InpUncOrig,InpFac);
  for(Int_t k = 0; k<InpUncOrig; k++){
    if(IsActiveUnc(k) == 1){
      for(Int_t l = 0; l<InpFac; l++){	
	RhoSte->operator()(k,l) = (MaxFac->operator()(k,l) - MinFac->operator()(k,l))/InrFac;
	RhoLoo->operator()(k,l) =  MaxFac->operator()(k,l) - RhoSte->operator()(k,l);
      }
    }
  }
    
  // Perform the loop
  Int_t ka = 0;
  Int_t FaiAct = 0;
  ReleaseInp();  
  for(Int_t ll = 0; ll<InrFac; ll++){
    //printf("... Blue->SolveScaRho(): \n");
    for(Int_t l = 0; l<InpFac; l++){
      ka = 0;    
      for(Int_t k = 0; k<InpUncOrig; k++){
	FaiAct = 0;
	if(IsActiveUnc(k) == 1){
	  SetRhoFacUnc(k, l, RhoLoo->operator()(k,l));
	  //printf("... Blue->SolveScaRho(): k = %2i, l = %2i, rho = %5.2f(%5.2f,%5.2f) \n", 
	  //	 k, l, ActFac->operator()(k,l), MinFac->operator()(k,l), MaxFac->operator()(k,l));
	  //PrintMatrix(ActFac);
	  FixInp();
	  //PrintCor(k);
	  Solve();
	  FaiAct = InspectResult();
	  for(Int_t n = 0; n<InpObs; n++){
	    // Check for failures
	    if(FaiAct < 0){
	      FaiFac = FaiFac + 1;
	      //printf("... Blue->SolveScaRho(): Failure = %2i \n", FaiFac);
	      ValSca[n]->operator()(ka+l*InpUnc,ll) = -1.00;
	      SigSca[n]->operator()(ka+l*InpUnc,ll) = -1.00;
	    }else{
	      ValSca[n]->operator()(ka+l*InpUnc,ll) = XvaRes->operator()(n) - ValDef->operator()(n);
	      SigSca[n]->operator()(ka+l*InpUnc,ll) = TMath::Sqrt(CovRes->operator()(n,n))-SigDef->operator()(n);
	    }
	  }
	  //printf("... Blue->SolveScaRho(): ValSca \n"); ValSca[0]->Print();
	  //printf("... Blue->SolveScaRho(): SigSca \n"); SigSca[0]->Print();
	  ReleaseInp();
	  // Keep scaling independent per group l
	  if(FlaFac == 0)SetRhoFacUnc(k, l, 1.0);
	  //printf("... Blue->SolveScaRho(): k = %2i, l = %2i, rho = %5.2f(%5.2f,%5.2f) \n", 
	  //k, l, ActFac->operator()(k,l), MinFac->operator()(k,l), MaxFac->operator()(k,l));
	  //PrintMatrix(ActFac);
	  ka = ka + 1;
	}
      }
      for(Int_t k = 0; k<InpUncOrig; k++){
	if(IsActiveUnc(k) == 1){SetRhoFacUnc(k, l, 1.0);};
      }
    }
    for(Int_t k = 0; k<InpUncOrig; k++){
      if(IsActiveUnc(k) == 1){
	for(Int_t l = 0; l<InpFac; l++){
	  RhoLoo->operator()(k,l) = RhoLoo->operator()(k,l) - RhoSte->operator()(k,l);
	}
      }
    }
  }

  // Fill the total values depending on FlaFac
  Double_t ValSum = 0, SigSum = 0; 
  for(Int_t n = 0; n<InpObs; n++){  
    for(Int_t ll = 0; ll<InrFac; ll++){
      for(Int_t l = 0; l<InpFac; l++){
	ValSum = 0;
	SigSum = 0; 
	if(FlaFac == 0){
	  for(Int_t k = 0; k<InpUnc; k++){
	    if(ValSca[n]->operator()(k+l*InpUnc,ll) != -1.0){
	      ValSum = ValSum + ValSca[n]->operator()(k+l*InpUnc,ll)*ValSca[n]->operator()(k+l*InpUnc,ll);
	    }
	    if(SigSca[n]->operator()(k+l*InpUnc,ll) != -1.0){
	      SigSum = SigSum + SigSca[n]->operator()(k+l*InpUnc,ll)*SigSca[n]->operator()(k+l*InpUnc,ll);
	    }
	  }
	  VtoSca[n]->operator()(l,ll) = TMath::Sqrt(ValSum);
	  StoSca[n]->operator()(l,ll) = TMath::Sqrt(SigSum);
	}else{
	  VtoSca[n]->operator()(l,ll) = ValSca[n]->operator()(InpUnc*(l+1)-1,ll);
	  StoSca[n]->operator()(l,ll) = SigSca[n]->operator()(InpUnc*(l+1)-1,ll);
	}
      }
    }
  }

  // Print out of result
  //for(Int_t n = 0; n<InpObs; n++){  
  //printf("... Blue->SolveScaRho(): ValSca \n"); ValSca[n]->Print();
  //printf("... Blue->SolveScaRho(): VtoSca \n"); VtoSca[n]->Print();
  //printf("... Blue->SolveScaRho(): SigSca \n"); SigSca[n]->Print();
  //printf("... Blue->SolveScaRho(): StoSca \n"); StoSca[n]->Print();
  //}

  // Restore the initial factors
  ActFac->SetSub(0,0,*SavFac);
  FixInp();
  Solve();

  // Remove the quiet modus
  SetQuiet(0);

  // Clean up
  SavFac->Delete(); SavFac = NULL;
  ValDef->Delete(); ValDef = NULL;
  SigDef->Delete(); SigDef = NULL;
  RhoLoo->Delete(); RhoLoo = NULL;
  RhoSte->Delete(); RhoSte = NULL;

  // Set the flag
  SetIsSolvedScaRho(1);

  // Print out if wanted
  if(IsPrintLevel() > 0){PrintScaRho();};

  return;
};

//------------------------------------------------------------------------------

void Blue::SetIsSolved(const Int_t l){
  IsSolve = l;
  if(l == 1 && IsQuiet() == 0){
      printf("... Blue->SetIsSolved(l): Input was solved! \n");
      printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedRelUnc(const Int_t l){
  IsSolveRelUnc = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedRelUnc(l): Input was solved with Rel-ative Unc-ertainties! \n");
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedAccImp(const Int_t l){
  IsSolveAccImp = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedAccImp(l): Input was solved Acc-ording to the estimate Imp-ortance! \n");
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedScaRho(const Int_t l){
  IsSolveScaRho = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedScaRho(l): Input was solved Sca-nning the correlations Rho! \n");
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedInfWei(const Int_t l){
  IsSolveInfWei = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedInfWei(l): Input was solved with Inf-ormation Wei-ghts! \n");
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedPosWei(const Int_t l){
  IsSolvePosWei = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedPosWei(l): Input was solved with only Pos-itive Wei-ghts! \n");
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsSolvedMaxVar(const Int_t l){
  IsSolveMaxVar = l;
  if(l == 1 && IsQuiet() == 0){
    printf("... Blue->SetIsSolvedMaxVar(%2i): Input was solved to Max-imise its Var-iance! \n",l);
    printf("\n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetIsInspectLike(const Int_t l){
  IsInspectLike = l;
  if(l == 1){
    printf("... Blue->SetIsInspectLike(%2i): The likelihood was inspected! \n",l);
    printf("\n");
  }
};

//------------------------------------------------------------------------------
// Filler
//------------------------------------------------------------------------------

void Blue::FillCov(){
  // Reset the Covariance Matrix
  Cov->Delete(); Cov  = new TMatrixD(InpEst,InpEst);

  TMatrixD *C = new TMatrixD(InpEst,InpEst);
  TMatrixD *E = new TMatrixD(InpEst,InpEst);
  TMatrixD *H = new TMatrixD(InpEst,InpEst);
  TMatrixD *I = new TMatrixD(InpEst,InpEst);

  Int_t IRowLow, IRowHig, IColLow, IColHig;
  for(Int_t k = 0; k<InpUnc; k++){
    IRowLow = k*InpEst; 
    IRowHig = IRowLow + InpEst-1;
    IColLow = IRowLow;
    IColHig = IRowHig;
    //printf("... Blue->FillCov(): %3i, %3i, %3i, %3i \n",IRowLow, IColLow, IRowHig, IColHig);

    // Sigmas per source
    Unc->GetSub(IRowLow, IRowHig, IColLow, IColHig, *E, "S");
    //E->Print();
    
    // Correlation per source
    Cor->GetSub(IRowLow, IRowHig, IColLow, IColHig, *C, "S");
    //C->Print();

    // Sum it up = sum E*C*E
    H->Mult(*E, *C);
    I->Mult(*H, *E);
    //I->Print();

    Cov->operator+=(*I);
    //printf("... Blue->FillCov(): Covariance \n"); Cov->Print();
  }

  // Enable to check for positive Eigenvalues
  // ICheck = GetMatEigen(Cov);

  // Clean up
  C->Delete(); C = NULL;
  E->Delete(); E = NULL;
  H->Delete(); H = NULL;
  I->Delete(); I = NULL;
  return;
};

//------------------------------------------------------------------------------

void Blue::FillCovInvert(){
  // Reset inverse covariance
  CovI->Delete(); CovI = new TMatrixD(InpEst,InpEst);
  
  // Get covarianca and invert
  CovI->SetSub(0,0,*Cov);
  CovI->Invert();
  //printf("... Blue->FillCovInvert(): Inverse covariance \n"); CovI->Print();

  return;
};
 
//------------------------------------------------------------------------------

void Blue::FillSig(){
  // Reset the uncertainty vector
  Sig->Delete(); Sig = new TVectorD(InpEst);
  for(Int_t i = 0; i<InpEst; i++){
    Sig->operator()(i) = sqrt(Cov->operator()(i,i));
  }
  if(IsPrintLevel() >= 2){
    printf("... Blue->FillSig() \n"); 
    Sig->Print();
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::FillRho(){
  // Reset the correlation matrix
  Rho->Delete(); Rho = new TMatrixD(InpEst,InpEst);
  Rho->operator=(*Cov);
  //Rho->Print();
  for(Int_t i = 0; i<InpEst; i++){
    for(Int_t j = i+1; j<InpEst; j++){
      Rho->operator()(i,j) = Rho->operator()(i,j) /
	sqrt(Rho->operator()(i,i) * Rho->operator()(j,j));
      Rho->operator()(j,i) = Rho->operator()(i,j);
    }
  }
  for(Int_t i = 0; i<InpEst; i++){
    Rho->operator()(i,i) = 1.0;
  }
  //Rho->Print();
  return;
};

//------------------------------------------------------------------------------
// Getters
//------------------------------------------------------------------------------

Int_t Blue::IsAllowedEst(const Int_t i) const {
  if(i >= 0 && i < InpEstOrig){
    return 1;
  }else{
    printf("... Blue->IsAllowedEst(%2i): Not an allowed estimate \n", i);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsActiveEst(const Int_t i) const {
  if(IsAllowedEst(i) == 1){
    return EstAct[i];
  }else{
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledEst(const Int_t i) const {
  if(IsAllowedEst(i) == 1){
    return EstFil[i];
  }else{
    printf("... Blue->IsFilledEst(%2i): Not a filled estimate \n", i);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledEst() const {
  Int_t iok = 1;
  for(Int_t i = 0; i<InpEstOrig; i++){
    iok = iok * IsFilledEst(i);
  }
  if(iok == 1){
    return 1;
  }else{
    //printf("... Blue->IsFilledEst(): Not all estimates are filled \n");
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsWhichEst(const Int_t ia) const {
  // Returns index of original list
  if(ia < InpEst){
    return LisEst[ia];
  }else{
    printf("... Blue->IsWhichEst(%2i): Not a defined estimate \n", ia);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsIndexEst(const Int_t i) const {
  // Returns index of actual list
  if(i < InpEstOrig){
    for(Int_t j = 0; j<InpEst; j++){
      if(LisEst[j] == i)return j;
    }
  }else{
    printf("... Blue->IsIndexEst(%2i): Not a defined estimate \n", i);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::EstWhichObs(const Int_t i) const {
  if(IsAllowedEst(i) == 1){
    return EstObs[i];
  }else{
    printf("... Blue->EstWhichObs(%2i): Not an allowed estimate \n", i);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsAllowedUnc(const Int_t k) const {
  if(k >= 0 && k < InpUncOrig){
    return 1 ;
  }else{
    printf("... Blue->IsAllowedUnc(%2i): Not an allowed uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsActiveUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncAct[k];
  }else{
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncFil[k];
  }else{
    printf("... Blue->IsFilledUnc(%2i): Not a filled uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledUnc() const {
  Int_t iok = 1;
  for(Int_t k = 0; k<InpUncOrig; k++){
    iok = iok * IsFilledUnc(k);
  }
  if(iok == 1){
    return 1;
  }else{
    //printf("... Blue->IsFilledUnc(): Not all uncertainties filled \n");
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsAllowedFac(const Int_t l) const {
  if(l >= 0 && l < InpFac){
    return 1;
  }else{
    printf("... Blue->IsAllowedFac(%2i): Not an allowed group \n", l);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsRhoValUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncCha[k];
  }else{
    printf("... Blue->IsRhoValUnc(%2i): Not an allowed uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsRhoFacUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncFac[k];
  }else{
    printf("... Blue->IsRhoFacUnc(%2i): Not an allowed uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsRhoRedUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncRed[k];
  }else{
    printf("... Blue->IsRhoRedUnc(%2i): Not an allowed uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsRelValUnc(const Int_t k) const {
  if(IsAllowedUnc(k) == 1){
    return UncRel[k];
  }else{
    printf("... Blue->IsRelValUnc(%2i): Not an allowed uncertainty \n", k);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsWhichUnc(const Int_t ka) const {
  // Returns index of original list
  if(ka < InpUnc){
    return LisUnc[ka];
  }else{
    printf("... Blue->IsWhichUnc(%2i): Not a defined uncertainty \n", ka);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsIndexUnc(const Int_t k) const {
  // Returns index of actual list
  if(k < InpUncOrig){
    for(Int_t kk = 0; kk<InpUnc; kk++){
      if(LisUnc[kk] == k)return kk;
    }
  }else{
    printf("... Blue->IsIndexUnc(%2i): Not a defined uncertainty \n", k);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsAllowedObs(const Int_t n) const {
  if(n > -1 && n < InpObsOrig){
    return 1 ;
  }else{
    printf("... Blue->IsAllowedObs(%2i): Not an allowed observable \n", n);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsActiveObs(const Int_t n) const {
  if(IsAllowedObs(n) == 1){
    return ObsAct[n];
  }else{
    printf("... Blue->IsActiveObs(%2i): Not an allowed observable \n", n);
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsWhichObs(const Int_t na) const {
  // Returns index of original list
  if(na < InpObs){
    return LisObs[na];
  }else{
    printf("... Blue->IsWhichObs(%2i): Not a defined observable \n", na);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsIndexObs(const Int_t n) const {
  // Returns index of actual list
  if(n < InpObsOrig){
    for(Int_t m = 0; m<InpObs; m++){
      if(LisObs[m] == n)return m;
    }
  }else{
    printf("... Blue->IsIndexObs(%2i): Not a defined observable \n", n);
    return -1;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledNamEst() const {
  if(IsFillEstNam == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsFilledNamEst(): Estimator names not yet filled \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledNamUnc() const {
  if(IsFillUncNam == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsFilledNamUnc(): Uncertainty names not yet filled \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledNamObs() const {
  if(IsFillObsNam == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsFilledNamObs(): Observable names not yet filled \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFilledInp() const {
  if(IsFilledEst()*IsFilledUnc() == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsFilledInp(): Not all input filled \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsFixedInp() const {
  if(InpFix == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsFixedInp(): Not fixed input \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsCalcedParams() const {
  return IsCalcParams;
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolved() const {
  if(IsSolve == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolved(): Not yet solved Call Solve() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedRelUnc() const {
  if(IsSolveRelUnc == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedRelUnc(): Not yet solved Call SolveRelUnc() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedAccImp() const {
  if(IsSolveAccImp == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedAccImp(): Not yet solved Call SolveAccImp() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedScaRho() const {
  if(IsSolveScaRho == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedScaRho(): Not yet solved Call SolveScaRho() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedInfWei() const {
  if(IsSolveInfWei == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedInfWei(): Not yet solved Call SolveInfWei() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedPosWei() const {
  if(IsSolvePosWei == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedPosWei(): Not yet solved Call SolvePosWei() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedMaxVar() const {
  if(IsSolveMaxVar == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedMaxUnc(): Not yet solved Call SolveMaxVar() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsSolvedAnyWay() const {
  if(IsSolve == 1 || IsSolveRelUnc == 1 || IsSolveAccImp == 1 || 
     IsSolveInfWei == 1 || IsSolvePosWei == 1 || IsSolveMaxVar == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsSolvedAnyWay(): Not yet solved Call any Solver() \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsInspectedLike() const {
  if(IsInspectLike == 0 && IsPrintLevel() >= 2){
    printf("... Blueinspectl->IsInspectLike(): Not yet done Call InspectLike() \n");
  }
  return IsInspectLike;
};

//------------------------------------------------------------------------------

Int_t Blue::IsRelValUnc() const {
  if(IsRelUncMode == 1){
    return 1;
  }else{
    if(IsPrintLevel() >= 2){
      printf("... Blue->IsRelValUnc(): Not in the mode with relative uncertainties \n");
    }
    return 0;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsPrintLevel() const {
  return IPrint;
};

//------------------------------------------------------------------------------

Int_t Blue::IsQuiet() const {
  return IQuiet;
};

//------------------------------------------------------------------------------

Int_t Blue::IsFutchCor() const {
  return IFuCor;
};

//------------------------------------------------------------------------------

Int_t Blue::GetMatEigen(const TMatrixD *TryMat) const {

  // The return value  1 / 0 = only positive / some negative Eigenvalues
  Int_t ISign = 1;

  // Get the dimensions define the structures
  Int_t          NRows = TryMat->GetNrows();
  TMatrixDEigen* EMatrix = new TMatrixDEigen(*Cov);
  TMatrixD*      EValues = new TMatrixD(NRows,NRows);
  TMatrixD*      EVector = new TMatrixD(NRows,NRows);
  
  // Fill the Eigenvectors and Eigenvalues
  *EVector = EMatrix->GetEigenVectors();
  *EValues = EMatrix->GetEigenValues();
  //printf("... Blue->GetMatEigen(): EVector \n"); EVector->Print();
  //printf("... Blue->GetMatEigen(): EValues \n"); EValues->Print();

  // Watch out for negative Eigenvalues and report
  Int_t INega = 0;
  for(Int_t i = 0; i<InpEst; i++){
    if(INega == 0 && EValues->operator()(i,i) < 0){
      INega = 1;
      if(IsPrintLevel() >= 2){
	printf("... Blue->GetMatEigen(): Covariance with negative Eigenvalue(s) \n");
	EValues->Print();
	ISign = 0;
      }
    }
  }

  // Clean up
  delete EMatrix; EMatrix = NULL;
  EValues->Delete(); EValues = NULL;
  EVector->Delete(); EVector = NULL;
  
  return ISign;
};

//------------------------------------------------------------------------------
// Setters
//------------------------------------------------------------------------------

void Blue::SetActiveEst(const Int_t i, const Int_t l){
  if(IsAllowedEst(i) == 1){
    if(l == 0){
      if(IsPrintLevel() >= 2){
	printf("... Blue->SetActiveEst(%2i,0): Deactivate estimate \n", i);
      }
      EstAct[i] = 0;
    }else if(l == 1){
      if(IsPrintLevel() >= 2){
	printf("... Blue->SetActiveEst(%2i,l):   Activate estimate \n", i);
      }
      EstAct[i] = 1;
    }else{
      printf("... Blue->SetActiveEst(%2i,%2i): IGNORED Only 1/0 = active/inactive are allowed \n", i,l);
    }
  }
};

//------------------------------------------------------------------------------

void Blue::SetFilledEst(const Int_t i) {
  EstFil[i] = 1;
  if(IsFilledEst() == 1){SetFilledEst();};
};

//------------------------------------------------------------------------------

void Blue::SetFilledEst() {
  if(IsQuiet() == 0)printf("... Blue->SetFilledEst: All estimates (%2i) filled \n", InpEstOrig);
  if(IsFilledUnc() == 1){SetFilledInp();};
};

//------------------------------------------------------------------------------

void Blue::SetActiveUnc(const Int_t k, const Int_t l){
  if(IsAllowedUnc(k) == 1){
    if(l == 0){
      if(IsPrintLevel() >= 2){
	printf("... Blue->SetActiveUnc(%2i,0): Deactivate uncertainty \n", k);
      }
      UncAct[k] = 0;
    }else if(l == 1){
      if(IsPrintLevel() >= 2){
	printf("... Blue->SetActiveUnc(%2i,l):   Activate uncertainty \n", k);
      }
      UncAct[k] = 1;
    }else{
      printf("... Blue->SetActiveUnc(%2i,%2i): IGNORED Only 1/0 = active/inactive allowed \n",k,l);
    }
  }
};

//------------------------------------------------------------------------------

void Blue::SetFilledUnc(const Int_t k) {
  UncFil[k] = 1;
  if(IsFilledUnc() == 1){SetFilledUnc();};
};

//------------------------------------------------------------------------------

void Blue::SetFilledUnc() {
  if(IsQuiet() == 0)printf("... Blue->SetFilledUnc: All uncertainties (%2i) filled \n", InpUncOrig);
  if(IsFilledEst() == 1){SetFilledInp();};
};

//------------------------------------------------------------------------------

void Blue::SetRhoValUnc(const Int_t k, const Int_t l, const Double_t RhoVal, const Int_t m){
  if(IsAllowedUnc(k) == 1 && IsAllowedFac(l) == 1){
    if(IsActiveUnc(k) == 1){
      if(RhoVal < -1 ||  RhoVal > 1){
	printf("... Blue->SetRhoValUnc(%2i,%5.3f,%2i): IGNORED correlation outside (-1, 1) \n", k, RhoVal, m);
	return;
      }else if(m < 0 ||  m>1){
	printf("... Blue->SetRhoValUnc(%2i,%5.3f,%2i): IGNORED m outside validity. 0,1 \n", 
	       k, RhoVal, m);
	return;
      }
      if(m == 0){
	if(IsPrintLevel() >= 1){
	  printf("... Blue->SetRhoValUnc(%2i,%5.3f,%2i): I revert to the original correlation for this uncertainty \n",  k, RhoVal, m);
	}
	UncCha[k] = m;
        ActFac->operator()(k,l) = 1.0;
      }else if(m == 1){
	if(IsRhoFacUnc(k) == 1 || IsRhoRedUnc(k) == 1){
	  printf("... Blue->SetRhoValUnc(%2i,%5.3f,%2i): IGNORED Only one of -Changed- -Scaled =%2i- or -Reduced =%2i- possible \n", 
		 k, RhoVal, m, IsRhoFacUnc(k), IsRhoRedUnc(k));
	  return;
	}
	if(IsPrintLevel() >= 1){
	  printf("... Blue->SetRhoValUnc(%2i,%5.3f,%2i): I change the correlation for the uncertainty \n",  k, RhoVal, m);
	}
	UncCha[k] = m;
        ActFac->operator()(k,l) = RhoVal;
      }
    }else{
      printf("... Blue->SetRhoValUnc(%2i): IGNORED not an active uncertainty \n", k);
    }
  }else{
    printf("... Blue->SetRhoValUnc(%2i): IGNORED not an allowed uncertainty or group %2i \n", k, l);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::SetRhoFacUnc(const Int_t k, const Int_t l, const Double_t RhoFac, const Int_t m){
  if(IsAllowedUnc(k) == 1 && IsAllowedFac(l) == 1){
    if(IsActiveUnc(k) == 1){
      if(RhoFac < -1 ||  RhoFac > 1){
	printf("... Blue->SetRhoFacUnc(%2i,%5.3f,%2i): IGNORED scale factor outside (-1, 1) \n", k, RhoFac, m);
	return;
      }else if(m < 0 ||  m>1){
	printf("... Blue->SetRhoFacUnc(%2i,%5.3f,%2i): IGNORED m outside validity. 0,1 \n", k, RhoFac, m);
	return;
      }
      if(m == 0){
	if(IsPrintLevel() >= 1){
	  printf("... Blue->SetRhoFacUnc(%2i,%5.3f,%2i): I revert to the original correlation for this uncertainty \n",  k, RhoFac, m);
	}
	UncFac[k] = m;
        ActFac->operator()(k,l) = 1.0;
      }else if(m == 1){
	if(IsRhoValUnc(k) == 1 || IsRhoRedUnc(k) == 1){
	  printf("... Blue->SetRhoFacUnc(%2i,%5.3f,%2i): IGNORED Only one of -Changed =%2i- -Scaled- or -Reduced =%2i- possible \n", 
		 k, RhoFac, m, IsRhoValUnc(k), IsRhoRedUnc(k));
	  return;
	}
	if(IsPrintLevel() >= 1){
	  printf("... Blue->SetRhoFacUnc(%2i,%5.3f,%2i): I scale the correlation for this uncertainty \n",  k, RhoFac, m);
	}
	UncFac[k] = m;
        ActFac->operator()(k,l) = RhoFac;
      }
    }else{
      printf("... Blue->SetRhoFacUnc(%2i): IGNORED not an active uncertainty \n", k);
    }
  }else{
    printf("... Blue->SetRhoFacUnc(%2i): IGNORED not an allowed uncertainty our group %2i \n", k,l);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::SetRhoRedUnc(const Int_t k, const Int_t l){
  if(IsAllowedUnc(k) == 1 && IsActiveUnc(k) == 1){
    if(l == 0){
      if(IsPrintLevel() >= 1){
	printf("... Blue->SetRhoRedUnc(%2i,0): I remove the reduced correlations for this uncertainty \n", k);
      }
      UncRed[k] = l;
    }else if(l == 1){
      if(IsRhoValUnc(k) == 1 || IsRhoFacUnc(k) == 1){
	printf("... Blue->SetRhoRedUnc(%2i,%2i): IGNORED Only one of -Changed =%2i- -Scaled =%2i- or -Reduced- possible \n", 
	       k, l, IsRhoValUnc(k), IsRhoFacUnc(k));
	return;
      }
      if(IsPrintLevel() >= 1){
	printf("... Blue->SetRhoRedUnc(%2i,l): Set reduced correlations for this uncertainty \n", k);
      }
      UncRed[k] = l;
    }else{
      printf("... Blue->SetRhoRedUnc(%2i,%2i): IGNORED Only 1/0 = active/inactive allowed \n",k,l);
    }
  }
};

//------------------------------------------------------------------------------

void Blue::SetRelValUnc(const Int_t l){
  if(l == 0){
    if(IsPrintLevel() >= 2){
      printf("... Blue->SetRelValUnc(%2i): I switch off the Blue with relative uncertainties \n", l);
    }
  }else if(l == 1){
    if(IsPrintLevel() >= 2){
      printf("... Blue->SetRelValUnc(%2i): I use the Blue with relative uncertainties \n", l);
    }
  }else{
    printf("... Blue->SetRelValUnc(%2i): IGNORED Only 1/0 = active/inactive allowed \n", l);
    return;
  }
  IsRelUncMode = l;
  return;
};

//------------------------------------------------------------------------------

void Blue::SetFilledInp(){
  // Fix input for ever
  if(IsFilledEst()*IsFilledUnc() == 1){
    // Not liked by .L Blue.cxx
    XvaOrig->operator=(*Xva);
    SigOrig->operator=(*Sig);
    UncOrig->operator=(*Unc);
    CorOrig->operator=(*Cor);
    UmaOrig->operator=(*Uma);
    UtrOrig->operator=(*Utr);

    // Save names
    for(Int_t i = 0; i<InpEstOrig; i++)EstNamOrig[i] = EstNam[i];
    for(Int_t k = 0; k<InpUncOrig; k++)UncNamOrig[k] = UncNam[k];
    for(Int_t n = 0; n<InpObsOrig; n++)ObsNamOrig[n] = ObsNam[n];

    // Much simpler see above
    //Xva->GetSub(0, InpEstOrig-1, *XvaOrig, "S");
    //Sig->GetSub(0, InpEstOrig-1, *SigOrig, "S");
    //Unc->GetSub(0, InpEstOrig*InpUncOrig-1, 0, InpEstOrig*InpUncOrig-1, *UncOrig, "S");
    //Cor->GetSub(0, InpEstOrig*InpUncOrig-1, 0, InpEstOrig*InpUncOrig-1, *CorOrig, "S");
    //Uma->GetSub(0, InpEstOrig-1, 0, InpObsOrig-1, *UmaOrig, "S");
    //Utr->GetSub(0, InpObsOrig-1, 0, InpEstOrig-1, *UtrOrig, "S");

    InpFil = 1;
    if(IsQuiet() == 0)printf("... Blue->SetFilledInp: All input filled! \n");
  }
};

//------------------------------------------------------------------------------

void Blue::SetFixedInp(const Int_t l){
  if(l == 0){
    InpFix = 0;
  }else if(l == 1){
    InpFix = 1;
    if(IsQuiet() == 0){
      printf("... Blue->SetFixedInp(l): Input fixed for solving! \n");
    }
  }else{
    printf("... Blue->SetFixedInp(%2i): IGNORED Only 1/0 = active/inactive allowed \n",l);
  }
  return;
};

//------------------------------------------------------------------------------

void Blue::SetQuiet(const Int_t l){
  IQuiet = l;
};

//------------------------------------------------------------------------------

void Blue::SetCalcedParams(const Int_t l){
  IsCalcParams = l;
};

//------------------------------------------------------------------------------

void Blue::SetFutchCor(const Int_t i, const Int_t j, const Int_t l, const Double_t RhoFut){
  if(IsFixedInp() == 1){
    printf("... Blue->SetFutchCor(%2i): IGNORED input not yet released.", l);
    printf(" Call ReleaseInp() or ResetInp() \n");
    return;
  }

  // Fill the matrix and set the flag
  if(l == 1){
    if(i == j){
      RhoFij->operator()(i,j) = 1;
    }else{
      RhoFij->operator()(i,j) = RhoFut;
      RhoFij->operator()(j,i) = RhoFij->operator()(i,j);
    }
  }else{
    RhoFij->operator()(i,j) = 1;
  }
  IFuCor = l;

  return;
};

//----------------------------------------------------------------------------
// PrintOut
//----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Calculate
//------------------------------------------------------------------------------

Double_t Blue::CalcRelUnc(const Int_t i, const Int_t k, const Double_t ValCom) const{
  // The return value
  Double_t UncCom = 0;
    
  // The default scaling
  if(IsRelValUnc(k) == 1){
    Int_t IRow = i*InpUncOrig + k;
    UncCom = Cof->operator()(IRow,0) + Cof->operator()(IRow,1)*TMath::Abs(ValCom) + Cof->operator()(IRow,2)*ValCom*ValCom;
    if(UncCom < 0){
      printf("... Blue->CalcRelUnc(%2i, %2i, %5.3f): Unphysical Uncertainty = %5.3f is smaller than zero", i, k, ValCom,UncCom);
      printf("==> No combination will be performed \n");
      UncCom = -1.0;
    }else{
      UncCom = TMath::Sqrt(UncCom);
    }
    if(IsPrintLevel() >= 2){
      printf("... Blue->CalcRelUnc(%2i, %2i, %5.3f): Uncertainty = %5.3f", i, k, ValCom,UncCom);
      printf(" 0 = %8.6f, 1 =  %8.6f, 2 =  %8.6f \n",
	     Cof->operator()(IRow,0),
	     Cof->operator()(IRow,1)*TMath::Abs(ValCom),
	     Cof->operator()(IRow,2)*ValCom*ValCom);
    }
  }
  return UncCom;
};

//------------------------------------------------------------------------------

Double_t Blue::Likelihood(const Double_t *x, const Double_t *par) const {

  //------------------------------------------------------------------------------
  //               1               1                    1
  // Like = ---------------- -----------------  exp [- --- (xi-xx)^T V^-1 (xi-xx)] 
  //        sqrt(2pi)^InpEst sqrt(Abs(Det(V)))          2
  //
  //        with: xx = x[0] <==>     true x for xi determining Obs == n
  //              xx = xm   <==> combined x for xi determining Obs == m != n
  //
  // par[0] = Index of the observable n in question 
  // par[1] = Total number of estimates
  // par[2] = 0/1 for the BLUE Gaussian / Likelihood function
  //------------------------------------------------------------------------------

  // The observable in question
  Int_t n = static_cast<Int_t>(par[0]);
  
  // The dimension of the Gauss == The number of estimates
  Double_t d = par[1];
  
  // Get the Matrix that holds:
  // 0,0 = the chi squared for the likelihood 
  // 1,0 = the determinant of the covariance matrix
  Int_t IsLike = static_cast<Int_t>(par[2]);
  CalcChiRes(n, IsLike, x[0]);
  //PrintMatrix(ChiRes,"%8.4f");

  // Calculate the factors in front of the Gauss
  Double_t PiiFac = 1. / TMath::Power(2*TMath::ACos(-1.),d/2);
  Double_t DetFac = 1. / TMath::Sqrt(ChiRes->operator()(1,0));
  
  // Return the Gauss
  return PiiFac*DetFac*TMath::Exp(-0.5*ChiRes->operator()(0,0));
};

//------------------------------------------------------------------------------

void Blue::CalcChiRes(const Int_t n, const Int_t IsLike, const Double_t xt) const {

  // If IsLike == 0/1 <==> BLUE Gaussian / Likelihood
  // If IsRela == 0/1 <==> Absolute / relative uncertainties
  // Calculate covariance for IsLike == IsRela == 1 and IsSolvedRelUnc
  // else use the existing one

  // Check if at least one relative uncertainty exist
  Int_t IsRela = 0;
  for(Int_t k = 0; k<InpUncOrig; k++)if(IsRelValUnc(k) == 1)IsRela = 1;

  // Define the matrices to work with
  TMatrixD *EE = new TMatrixD(InpEst,1);
  TMatrixD *ET = new TMatrixD(1,InpEst);
  TMatrixD *CO = new TMatrixD(InpEst, InpEst);
  TMatrixD *CI = new TMatrixD(InpEst, InpEst);
  TMatrixD *DU = new TMatrixD(InpEst, 1);
  TMatrixD *VA = new TMatrixD(1, 1);

  // Fill the vectors with difference of:
  // Estimate and         xt for obs == n
  // Estimate and observable for obs != n
  Int_t na;
  for(Int_t i = 0; i<InpEst; i++){
    if(EstWhichObs(IsWhichEst(i)) == n){
      EE->operator()(i,0) = Xva->operator()(i) - xt;
    }else{
      na = IsIndexObs(EstWhichObs(IsWhichEst(i)));
      EE->operator()(i,0) = Xva->operator()(i) - XvaRes->operator()(na);      
    }
    ET->operator()(0,i) = EE->operator()(i,0);
  }
  // printf("... Blue->CalcChiRes(): EE, ET \n"); EE->Print(); ET->Print();

  // Calculate covariance as a function of xt
  Double_t Sigik = 0, Sigjk = 0, Rhova = 0, xi = 0, xj = 0;
  Int_t    io = 0, jo = 0, ko = 0, no = 0;
  if(IsLike == 1 && IsRela == 1){
    for(Int_t i = 0; i<InpEst; i++){
      io = IsWhichEst(i);
      no = EstWhichObs(io);
      xi = xt;
      if(no != n)xi = XvaRes->operator()(IsIndexObs(EstWhichObs(IsWhichEst(i))));
      for(Int_t j = i; j<InpEst; j++){
	jo = IsWhichEst(j);
	no = EstWhichObs(jo);	
	xj = xt;
	if(no != n)xj = XvaRes->operator()(IsIndexObs(EstWhichObs(IsWhichEst(j))));
	for(Int_t k = 0; k<InpUnc; k++){
	  ko = IsWhichUnc(k);
	  //printf("... Blue->CalcChiRes(): i=%2i(%2i), j=%2i(%2i), k=%2i(%2i)\n",i,io,j,jo,k,ko);
	  if(i == j){
	    if(IsRelValUnc(ko) == 1){
	      Sigik = CalcRelUnc(io, ko, xi);
	    }else{
	      Sigik = Unc->operator()(i+k*InpEst,i+k*InpEst);
	    }
	    CO->operator()(i,i) = CO->operator()(i,i) + Sigik * Sigik;
	  }else{
	    if(IsRelValUnc(ko) == 1){
	      Sigik = CalcRelUnc(io, ko, xi);
	      Sigjk = CalcRelUnc(jo, ko, xj);
	    }else{
	      Sigik = Unc->operator()(i+k*InpEst,i+k*InpEst);
	      Sigjk = Unc->operator()(j+k*InpEst,j+k*InpEst);
	    }
	    Rhova = Cor->operator()(i+k*InpEst,j+k*InpEst);
	    CO->operator()(i,j) = CO->operator()(i,j) + Rhova * Sigik * Sigjk;
	  }
	}
	CO->operator()(j,i) = CO->operator()(i,j);
      }
    }
  }else{
    // Take the existing covariance matrix
    CO->SetSub(0,0,*Cov);
  }
  //printf("... Blue->CalcChiRes(): CovI, CI \n"); CovI->Print(); CI->Print();

  // Invert the covariance
  CI->SetSub(0,0,*CO);
  CI->Invert();
  //printf("... Blue->CalcChiRes(): Cov, CO \n");  Cov->Print(); CO->Print();

  // Calculate ET * CI * EE 
  DU->Mult(*CI, *EE);
  VA->Mult(*ET, *DU);

  // The results
  ChiRes->operator()(0,0) = VA->operator()(0, 0);
  ChiRes->operator()(1,0) = TMath::Abs(CO->Determinant());

  // Clean up matrices
  EE->Delete(); EE = NULL;
  ET->Delete(); ET = NULL;
  CO->Delete(); CO = NULL;
  CI->Delete(); CI = NULL;
  DU->Delete(); DU = NULL;
  VA->Delete(); VA = NULL;

  return;
};

//------------------------------------------------------------------------------

void Blue::CalcParams(){
  //printf("... Blue->CalcParams(): Rho \n"); Rho->Print()
  //printf("... Blue->CalcParams(): Sig \n"); Sig->Print();
  
  // Reset matrices
  SRat->Delete(); SRat = new TMatrixD(InpEstOrig,InpEstOrig);
  Beta->Delete(); Beta = new TMatrixD(InpEstOrig,InpEstOrig);
  Sigx->Delete(); Sigx = new TMatrixD(InpEstOrig,InpEstOrig);
  DBdr->Delete(); DBdr = new TMatrixD(InpEstOrig,InpEstOrig);
  DSdr->Delete(); DSdr = new TMatrixD(InpEstOrig,InpEstOrig);
  DBdz->Delete(); DBdz = new TMatrixD(InpEstOrig,InpEstOrig);
  DSdz->Delete(); DSdz = new TMatrixD(InpEstOrig,InpEstOrig);
  
  // Loop and fill lower half
  Double_t rho = 0, zva = 0;
  Int_t ii = 0, jj = 0;
  for(Int_t i = 0; i<InpEst; i++){
    for(Int_t j = i+1; j<InpEst; j++){
      if(EstWhichObs(IsWhichEst(i)) == EstWhichObs(IsWhichEst(j))){
	rho = Rho->operator()(i,j);
	zva = Sig->operator()(i)/Sig->operator()(j);
	ii = IsWhichEst(i);
	jj = IsWhichEst(j);
	SRat->operator()(jj,ii) = zva;
	if(zva < 1)zva = 1/zva;
	//printf("... Blue->CalcParams(): rho=%5.3f zva=%5.3f \n",rho,zva);
	Beta->operator()(jj,ii) = GetPara(1, rho, zva);
	Sigx->operator()(jj,ii) = GetPara(2, rho, zva);
	DBdr->operator()(jj,ii) = GetPara(3, rho, zva);
	DSdr->operator()(jj,ii) = GetPara(4, rho, zva);
	DBdz->operator()(jj,ii) = GetPara(5, rho, zva);
	DSdz->operator()(jj,ii) = GetPara(6, rho, zva);	 
      }
    }
  }
  SetCalcedParams(1);
  return;
};

//------------------------------------------------------------------------------

void Blue::ResetScaRho(const Int_t IFlag){

  // Set the actual factors to unity
  if(IFlag == 0 || IFlag == 2){ 
    for(Int_t k = 0; k<InpUncOrig; k++){
      for(Int_t l = 0; l<InpFac; l++){
	ActFac->operator()(k,l) = 1;
      }
    }
  }

  // Set the scan ranges to 0 -> 1
  if(IFlag == 0 || IFlag == 2){ 
    for(Int_t k = 0; k<InpUncOrig; k++){
      for(Int_t l = 0; l<InpFac; l++){
	MinFac->operator()(k,l) = 0;
	MaxFac->operator()(k,l) = 1;
      }
    }
  }
  return;
};

//------------------------------------------------------------------------------

Double_t Blue::CalcActFac(const Int_t i, const Int_t j, const Int_t k){

  // The subgroup index
  Int_t l = static_cast<int>(MatFac->operator()(i,j));

  // Return the factor
  return ActFac->operator()(k,l);
};

//------------------------------------------------------------------------------

Double_t Blue::CalcFutchCor(const Int_t i, const Int_t j){
  if(IsFutchCor() == 1){
    return RhoFij->operator()(i,j);
  }else{
    return 1.;
  }
};

//------------------------------------------------------------------------------

Int_t Blue::IsWhichMatrix(const TMatrixD *TryMat) const {
  // The flag 0 / 1 / 2 / 3 = uncorrelated / fully correlated / symmetric / not symmetric
  Int_t IsWhichM = 3;

  // Get the dimensions, exit if not symmetric matrix
  Int_t NRows = TryMat->GetNrows();
  Int_t NCols = TryMat->GetNcols();
  if(NRows != NCols){
    printf("... Blue->IsWhichMatrix(): IGNORED Only square matrices are allowed \n");
    return 3;
  }

  //
  // Try to see what matrix it is
  //
  TMatrixD* NCor = new TMatrixD(NRows,NRows);
  TMatrixD* FCor = new TMatrixD(NRows,NRows);
  TMatrixD* Null = new TMatrixD(NRows,NRows);
  TMatrixD* Comp = new TMatrixD(NRows,NRows);
  NCor->UnitMatrix();
  for(Int_t i = 0; i < NRows; i++){
    for(Int_t j = 0; j < NRows; j++){
      FCor->operator()(i,j) = 1;
      Null->operator()(i,j) = 0;
    }
  }
  //printf("... Blue->IsWhichMatrix(): NCor \n"); NCor->Print();
  //printf("... Blue->IsWhichMatrix(): FCor \n"); FCor->Print();
  //printf("... Blue->IsWhichMatrix(): Null \n"); Null->Print();

  // The testing flags
  Int_t ITest = 0;
  Double_t Small = 0.000000001;

  // Test for Un-correlated or Fully-correlated
  for(Int_t io = 0; io < 2; io++){    
    ITest = 0;
    Comp->operator=(*Null);
    if(io == 0){Comp->operator+=(*NCor);
    }else{Comp->operator+=(*FCor);};
    Comp->operator-=(*TryMat);
    //printf("... Blue->IsWhichMatrix(): Comp \n"); Comp->Print();
    for(Int_t i = 0; i < NRows; i++){
      for(Int_t j = i; j < NRows; j++){
	if(Comp->operator()(j,i) != 0)ITest = 1;
      }     
    }
    if(ITest == 0)IsWhichM = io;    
  }

  // Test for symmetric if the above failed
  ITest = 0;
  if(IsWhichM == 3){
    for(Int_t i = 0; i < NRows; i++){
      for(Int_t j = i; j < NRows; j++){
	if(TMath::Abs(TryMat->operator()(j,i)-TryMat->operator()(i,j)) > Small){
	  if(IsPrintLevel() >= 2){
	    printf("... Blue->IsWhichMatrix(): %2i %2i %10.8f  %10.8f \n",
		   i,j,TryMat->operator()(i,j),TryMat->operator()(j,i));
	  }
	  ITest = 1;
	}
      }     
    }
    if(ITest == 0)IsWhichM = 2;
  }

  // Clean up
  NCor->Delete(); NCor = NULL;
  FCor->Delete(); FCor = NULL;
  Null->Delete(); Null = NULL;
  Comp->Delete(); Comp = NULL;

  // We are done
  return IsWhichM;
};

//------------------------------------------------------------------------------

void Blue::MatrixtoDouble(const TMatrixD *InpMat, Double_t *OutDou) const {

  // Get Dimensions
  Int_t NRows = InpMat->GetNrows();
  Int_t NCols = InpMat->GetNcols();

  // Fill array from TMatrix
  Int_t ind = 0;
  for(Int_t i = 0; i<NRows; i++){
    for(Int_t j = 0; j<NCols; j++){
      OutDou[ind] = InpMat->operator()(i,j);
      ind = ind + 1;
    }
  }  
  return;
};

//------------------------------------------------------------------------------

void Blue::WriteInput(const Int_t IWhat, const TString FilNam, 
		      const TString ForVal, const TString ForRho) const {

  // Check whether input is fixed
  if(IsFixedInp() == 0){
    printf("... Blue->WriteInput(%2i): IGNORED Input not yet fixed call FixInp \n",IWhat);
    return;
  }
  
  // Check whether Flag is allowed
  if(IWhat != 1 && IWhat != 2){
    printf("... Blue->WriteInput(%2i): Not implemented Flag \n",IWhat);
    return;
  }
  
  // Dummy return variable for sprintf
  Int_t IRet = 0;
  if(IRet == 1)return;

  // Open file
  std::ofstream ofs (FilNam, std::ofstream::out);

  // Declare variables
  char c[100];
  TString Format = "To be filled later";
  
  // Write the header
  IRet = sprintf(c,"./combine<<!\n"); ofs<<c;
  IRet = sprintf(c,"'Input generated by Blue Version %s on %s' \n",Versio.Data(),Today.Data()); ofs<<c;
  IRet = sprintf(c,"-1 -1 0 internal & minuit debug level, dependency flag\n"); ofs<<c;
  Int_t LocUnc = InpUnc;
  if(IWhat == 1)LocUnc = 1;
  IRet = sprintf(c,"%2i %2i %2i # of observables, measurements, error classes\n",InpObs,InpEst,LocUnc); ofs<<c;
  for(Int_t n = 0; n < InpObsOrig; n++){
    if(IsActiveObs(n) == 1)IRet = sprintf(c," '%s'",ObsNam[n].Data()); ofs<<c;
  }  
  IRet = sprintf(c," name of observable \n"); ofs<<c;
  if(IWhat == 1){
    IRet = sprintf(c,"                       '%s'","Full"); ofs<<c;
  }else{
    for(Int_t k = 0; k < InpUncOrig; k++){
      if(IsActiveUnc(k) == 1)IRet = sprintf(c," '%s'",UncNam[k].Data()); ofs<<c;
    }    
  }
  IRet = sprintf(c,"\n"); ofs<<c;

  // The estimates and uncertainties
  Int_t ii = 0;
  for(Int_t i = 0; i < InpEst; i++){
    ii = IsWhichEst(i);
    IRet = sprintf(c," '%s' '%s' ", EstNam[ii].Data(), ObsNam[EstWhichObs(ii)].Data()); ofs<<c;
    Format = " " + ForVal;
    IRet = sprintf(c,Format,Xva->operator()(i)); ofs<<c;
    if(IWhat == 1){	
      IRet = sprintf(c,Format,Sig->operator()(i)); ofs<<c;
    }else{
      for(Int_t k = 0; k < InpUnc; k++){
	IRet = sprintf(c,Format,Unc->operator()(i+k*InpEst,i+k*InpEst)); ofs<<c;
      }
    }
    IRet = sprintf(c,"\n"); ofs<<c;
  }
  IRet = sprintf(c,"\n"); ofs<<c;

  // The correlations
  Format = " " + ForRho;
  if(IWhat == 1){	
    for(Int_t i = 0; i < InpEst; i++){
      for(Int_t j = 0; j < InpEst; j++){
	IRet = sprintf(c,Format,Rho->operator()(i,j)); ofs<<c;
      }
      if(i == 0){IRet = sprintf(c," %s ","Full"); ofs<<c;};
      IRet = sprintf(c,"\n"); ofs<<c;
    }
  }else{
    for(Int_t k = 0; k < InpUnc; k++){
      // Get the matrix per source
      TMatrixD *F = new TMatrixD(InpEst,InpEst);
      Int_t IRowLow = k*InpEst; 
      Int_t IRowHig = IRowLow + InpEst-1;
      Int_t IColLow = IRowLow; 
      Int_t IColHig = IRowHig;
      //  printf("... %3i, %3i, %3i, %3i \n",IRowLow, IColLow, IRowHig, IColHig);
      Cor->GetSub(IRowLow, IRowHig, IColLow, IColHig, *F, "S");

      // Write the findings
      for(Int_t i = 0; i < InpEst; i++){
	for(Int_t j = 0; j < InpEst; j++){
	  IRet = sprintf(c,Format,F->operator()(i,j)); ofs<<c;
	}
	if(i == 0){IRet = sprintf(c," '%s'",UncNam[IsWhichUnc(k)].Data()); ofs<<c;};
	IRet = sprintf(c,"\n"); ofs<<c;
      }
      IRet = sprintf(c,"\n"); ofs<<c;	
      // clean up
      F->Delete(); F = NULL;
    }
  }
  
  //Fill tail
  IRet = sprintf(c,"!\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  // Close file
  ofs.close();

  return;
};

//----------------------------------------------------------------------------
// Display
//----------------------------------------------------------------------------

void Blue::PlotRes(const Int_t N, const TString FilNam, 
		   const TString ForVal, const TString ForUnc) const {

  // Declare variables
  char c[100];
  TString Format = "To be filled later";

  // Dummy return variable for sprintf
  Int_t IRet = 0;
  if(IRet == 1)return;

  // Open file
  TString pdffile = FilNam + ".pdf";
  TString namfile = "void " + FilNam + "(TString pdffile=\"" + pdffile +"\"){";
  TString outfile = FilNam + ".cxx";
  std::ofstream ofs (outfile, std::ofstream::out);

  //Fill head
  IRet = sprintf(c,"#include \"TROOT.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TStyle.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TCanvas.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TMarker.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TH2F.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TLatex.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TPave.h\" \n"); ofs<<c;
  IRet = sprintf(c,"#include \"TGraphAsymmErrors.h\" \n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  // Start function
  ofs<<namfile;
  IRet = sprintf(c,"\n"); ofs<<c;

  // The dimensions
  IRet = sprintf(c,"\t const Int_t N = %2i; \n",N); ofs<<c;

  // The array of indices
  IRet = sprintf(c,"\t Float_t Indx[N] = {"); ofs<<c;
  IRet = sprintf(c,"%2i",Indx[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,",%2i",Indx[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // The Colors
  IRet = sprintf(c,"\t Int_t   Colo[N] = {"); ofs<<c;
  IRet = sprintf(c,"%2i",Colo[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,",%2i",Colo[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // The Names
  IRet = sprintf(c,"\t TString Name[N] = {"); ofs<<c;
  IRet = sprintf(c,"\" %s \"",Name[0].Data()); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,", \"%s\"",Name[i].Data()); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // The Values
  Format = "," + ForVal;
  IRet = sprintf(c,"\t Float_t Valu[N] = {"); ofs<<c;
  IRet = sprintf(c,ForVal, Valu[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
      IRet = sprintf(c,Format,Valu[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // The Stat uncertainties
  Format = "," + ForUnc;
  IRet = sprintf(c,"\t Float_t Stat[N] = {"); ofs<<c;
  IRet = sprintf(c,ForUnc,Stat[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,Format,Stat[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // Fill the syst 
  IRet = sprintf(c,"\t Float_t Syst[N] = {"); ofs<<c;
  IRet = sprintf(c,ForUnc, Syst[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,Format,Syst[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  // The Full Unc
  IRet = sprintf(c,"\t Float_t Full[N] = {"); ofs<<c;
  IRet = sprintf(c,ForUnc,Full[0]); ofs<<c;
  for(Int_t i = 1; i<N; i++){
    IRet = sprintf(c,Format,Full[i]); ofs<<c;
  }
  IRet = sprintf(c,"};\n"); ofs<<c;

  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t // Text Strings \n"); ofs<<c;
  IRet = sprintf(c,"\t TString CanTit = \"The combination You always wanted to do\";\n"); ofs<<c;
  IRet = sprintf(c,"\t TString AXIS = \"Obs_{yours} [units]\";\n"); ofs<<c;
  IRet = sprintf(c,"\t TString HEAD = \"Combined by BLUE %s on %s \";\n",
		 Versio.Data(),Today.Data()); ofs<<c;
  IRet = sprintf(c,"\t TString UNCB = \"(stat)        (syst)\";\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Canvas parameters\n"); ofs<<c;
  IRet = sprintf(c,"\t Int_t Canx = 800, Cany = 650;\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Histogram x-axis is booked depending on Nxmi and Nxma with\n"); ofs<<c;
  IRet = sprintf(c,"\t // x=(Valu[0]-Fxmi*largest extension of unc bar to tge left,\n"); ofs<<c;
  IRet = sprintf(c,"\t //    Valu[0]+Fxma*largest extension of unc bar to the right) and y(-2, N+2)\n"); ofs<<c;
  IRet = sprintf(c,"\t // The offset of values (unc text) from combined valu can be steered with Fxva(Fxun)\n"); ofs<<c;
  IRet = sprintf(c,"\t //---------------------------------------------------------------------- \n"); ofs<<c;
  IRet = sprintf(c,"\t // MOST CASES SHOULD BE ADJUSTABLE BY CHANGING THE FOLLOWING FOUR VALUES \n"); ofs<<c;
  IRet = sprintf(c,"\t //---------------------------------------------------------------------- \n"); ofs<<c;
  IRet = sprintf(c,"\t Float_t Fxmi = 2.0, Fxma = 3.0;\n"); ofs<<c;
  IRet = sprintf(c,"\t Float_t Fxva = 1.1, Fxun = 2.0;\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\n"); ofs<<c;
  IRet = sprintf(c,"\t //----------------------------------------------------------------------------\n"); ofs<<c;
  IRet = sprintf(c,"\t // No number or text needed below this line\n"); ofs<<c;
  IRet = sprintf(c,"\t //----------------------------------------------------------------------------\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Get min and max in x\n"); ofs<<c;
  IRet = sprintf(c,"\t Float_t xmi = 0, xma = 0, xfu = 0;\n"); ofs<<c;
  IRet = sprintf(c,"\t Float_t xText = 0, xValu = 0, xUnco = 0;\n"); ofs<<c;
  IRet = sprintf(c,"\t Float_t xmin = Valu[0], xmax = xmin;\n"); ofs<<c;
  IRet = sprintf(c,"\t for(Int_t j = 0; j < N; j++) {\n"); ofs<<c;
  IRet = sprintf(c,"\t\t if(Valu[j]-Full[j] < xmin)xmin = Valu[j]-Full[j];\n"); ofs<<c;
  IRet = sprintf(c,"\t\t if(Valu[j]+Full[j] > xmax)xmax = Valu[j]+Full[j];\n"); ofs<<c;
  IRet = sprintf(c,"\t }\n"); ofs<<c;
  IRet = sprintf(c,"\t xfu = (xmax - xmin);\n"); ofs<<c;
  IRet = sprintf(c,"\t xmi = Valu[0] - Fxmi*xfu;\n"); ofs<<c;
  IRet = sprintf(c,"\t xma = Valu[0] + Fxma*xfu;\n"); ofs<<c;
  IRet = sprintf(c,"\t xText = xmi + xfu/10.;\n"); ofs<<c;
  IRet = sprintf(c,"\t xValu = Valu[0] + Fxva*xfu;\n"); ofs<<c;
  IRet = sprintf(c,"\t xUnco = Valu[0] + Fxun*xfu;\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Set up graphs with stat and full uncertainties\n"); ofs<<c;
  IRet = sprintf(c,"\t TGraphAsymmErrors *grStat = new TGraphAsymmErrors(N, Valu, Indx, Stat, Stat, 0, 0);\n"); ofs<<c;
  IRet = sprintf(c,"\t TGraphAsymmErrors *grFull = new TGraphAsymmErrors(N, Valu, Indx, Full, Full, 0, 0);\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Now draw everything\n"); ofs<<c;
  IRet = sprintf(c,"\t gROOT->Reset(); \n"); ofs<<c;
  IRet = sprintf(c,"\t gROOT->SetStyle(\"Plain\");\n"); ofs<<c;
  IRet = sprintf(c,"\t gStyle->SetOptStat(0);\n"); ofs<<c;
  IRet = sprintf(c,"\t gStyle->SetEndErrorSize(6);\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;
  
  IRet = sprintf(c,"\t // Make canvas\n"); ofs<<c;
  IRet = sprintf(c,"\t TCanvas *CanDel;\n"); ofs<<c;
  IRet = sprintf(c,"\t CanDel = (TCanvas*) gROOT->FindObject(\"myCanvas\");\n"); ofs<<c;
  IRet = sprintf(c,"\t if(CanDel) delete CanDel;\n"); ofs<<c;
  IRet = sprintf(c,"\t TCanvas *myCanvas = new TCanvas(\"myCanvas\",CanTit, Canx, Cany);\n"); ofs<<c;
  IRet = sprintf(c,"\t myCanvas->SetTopMargin(0.025);\n"); ofs<<c;
  IRet = sprintf(c,"\t myCanvas->SetBottomMargin(0.14);\n"); ofs<<c;
  IRet = sprintf(c,"\t myCanvas->SetLeftMargin(0.04);\n"); ofs<<c;
  IRet = sprintf(c,"\t myCanvas->SetRightMargin(0.04);\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;
  
  IRet = sprintf(c,"\t // Book the TH2F plot\n"); ofs<<c;
  IRet = sprintf(c,"\t TH2F *H2FDel = (TH2F*) gROOT->FindObject(\"myPlot\");\n"); ofs<<c;
  IRet = sprintf(c,"\t if(H2FDel)delete H2FDel;\n"); ofs<<c;
  IRet = sprintf(c,"\t TH2F *myPlot = new TH2F(\"myPlot\", \"\", 20, xmi, xma, 20, -2.0, N + 1.0);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->SetTitle(AXIS);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->CenterTitle();\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->SetTitleFont(42);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->SetTitleSize(0.06);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->SetNdivisions(6);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetXaxis()->SetLabelFont(42);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetYaxis()->SetLabelColor(0);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetYaxis()->SetNdivisions(1);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->GetYaxis()->SetLabelFont(42);\n"); ofs<<c;
  IRet = sprintf(c,"\t myPlot->Draw();  \n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;
  
  IRet = sprintf(c,"\t // Put the vertical line around the combined value\n"); ofs<<c;
  IRet = sprintf(c,"\t TPave *Comb = new TPave(Valu[0]-Full[0], Indx[0]-0.5, Valu[0]+Full[0], Indx[N-1]+0.5, 0, \"br\");\n"); ofs<<c;
  IRet = sprintf(c,"\t Comb->SetFillColor(17);\n"); ofs<<c;
  IRet = sprintf(c,"\t Comb->Draw();\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;
  
  IRet = sprintf(c,"\t // Draw results as points\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->SetLineColor(4);\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->SetLineWidth(3);\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->SetMarkerColor(2);\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->SetMarkerStyle(20);\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->SetMarkerSize(1.3);\n"); ofs<<c;
  IRet = sprintf(c,"\t grFull->Draw(\"P\");\n"); ofs<<c;  
  IRet = sprintf(c,"\t grStat->SetLineColor(2);\n"); ofs<<c;
  IRet = sprintf(c,"\t grStat->SetLineWidth(3);\n"); ofs<<c;
  IRet = sprintf(c,"\t grStat->SetMarkerColor(2);\n"); ofs<<c;
  IRet = sprintf(c,"\t grStat->SetMarkerStyle(20);\n"); ofs<<c;
  IRet = sprintf(c,"\t grStat->SetMarkerSize(1.3);\n"); ofs<<c;
  IRet = sprintf(c,"\t grStat->Draw(\"Psame\");\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Store numerical values in a TString\n"); ofs<<c;
  IRet = sprintf(c,"\t // Print points and write values\n"); ofs<<c;
  IRet = sprintf(c,"\t TLatex *Text = new TLatex();\n"); ofs<<c;
  IRet = sprintf(c,"\t char cx[200];\n"); ofs<<c;
  IRet = sprintf(c,"\t for(Int_t j = 0; j < N; j++) {\n"); ofs<<c;
  //  IRet = sprintf(c,"\t\t sprintf(cx,\" %%5.2F #pm %%5.2F #pm %%5.2F\",Valu[j],Stat[j],Syst[j]);\n"); ofs<<c;
    
  IRet = sprintf(c,"\t\t sprintf(cx,\""); ofs<<c; 
  IRet = sprintf(c,"%s",ForVal.Data()); ofs<<c;
  IRet = sprintf(c," #pm "); ofs<<c;
  IRet = sprintf(c,"%s",ForUnc.Data()); ofs<<c;
  IRet = sprintf(c," #pm "); ofs<<c;
  IRet = sprintf(c,"%s",ForUnc.Data()); ofs<<c;
  IRet = sprintf(c,"\",Valu[j],Stat[j],Syst[j]);\n"); ofs<<c;

  IRet = sprintf(c,"\t\t TString Resu = cx;\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->SetTextSize(0.015);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->SetTextColor(Colo[j]);\n"); ofs<<c;
  //  IRet = sprintf(c,"\t\t if(j==0)Text->SetTextColor(kRed);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->SetTextFont(42);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->SetTextSize(0.035);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->DrawLatex(xText, Indx[j]-0.12, Name[j]);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->SetTextSize(0.045);\n"); ofs<<c;
  IRet = sprintf(c,"\t\t Text->DrawLatex(xValu, Indx[j]-0.12, Resu);\n"); ofs<<c;
  IRet = sprintf(c,"\t}\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Print the header and the uncertainty breakdown\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->SetTextSize(0.04);\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->SetTextFont(42);\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->SetTextColor(kBlack);\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->DrawLatex(xText,Indx[N-1]+1.0, HEAD);\n"); ofs<<c;  
  IRet = sprintf(c,"\t Text->SetTextSize(0.03);\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->SetTextColor(kBlack);\n"); ofs<<c;
  IRet = sprintf(c,"\t Text->DrawLatex(xUnco,Indx[0]-1.0, UNCB);\n"); ofs<<c;
  IRet = sprintf(c,"\n"); ofs<<c;

  IRet = sprintf(c,"\t // Write the file\n"); ofs<<c;  
  IRet = sprintf(c,"\t myCanvas->Print(pdffile);\n"); ofs<<c;
  
  //Fill tail
  IRet = sprintf(c,"}\n"); ofs<<c;

  // Close file
  ofs.close();
};
