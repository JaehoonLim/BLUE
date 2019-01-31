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
#include "Blue.h"

void B_AN2015_240(Int_t Flag = 0){

  //-----------------------------------------------------------------------------
  // Original File : B_NIMA_500_391.cxx
  // Modified by JaeHoon Lim, Seungkyu Ha for CMS-AN2015-240
  //
  // Flag steers which of the results should be calculated
  //
  // 1-7
  // 1: un-correlated stat. uncertainties
  // 2: un-correlated stat. and dilution uncertainties
  // 3: un-correlated stat. uncertainties and 100% correlated dilution uncertainties
  // 4: un-correlated stat. and sys. uncertainties
  // 5: un-correlated stat. uncertainties and 100% correlated sys. uncertainties
  // 6: un-correlated all uncertainties
  // 7: un-correlated stat. uncertainties and 100% correlated dilution and sys. uncertainties
  //
  // 10: Asymmetry
  // 20: Asymmetry (total uncertainties)
  // 30: dtG
  // 
  // 100: O_1
  // 300: O_3
  //
  //-----------------------------------------------------------------------------

  // Print log option
  bool PrintLog = false;

  // Print Flag
  const TString FlagName1[7] = {
    " with un-correlated stat. uncertainties",
    " with un-correlated stat. and dilution uncertainties",
    " with un-correlated stat. uncertainties and 100% correlated dilution uncertainties",
    " with un-correlated stat. and sys. uncertainties",
    " with un-correlated stat. uncertainties and 100% correlated sys. uncertainties",
    " with un-correlated all uncertainties",
    " with un-correlated stat. uncertainties and 100% correlated dilution and sys. uncertainties"
  };

  const TString FlagName2[3] = {
    " Asymmetry",
    " Asymmetry (root sum squared)",
    " dtG"
  };

  const TString FlagName3[3] = {
    "O_1",
    "",
    "O_3"
  };

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst   = 3;
  static const Int_t StatUnc  = 1;
  static const Int_t ADiluUnc = 1;
  static const Int_t EDiluUnc = 4;
  static const Int_t ASystUnc = 1;
  static const Int_t ESystUnc = 24;
  static const Int_t NumObs   = 1;
  Int_t DiluUnc = 0;
  Int_t SystUnc = 0;
  Int_t NumUnc = 0;

  // Index for which estimate determines which observable
  Int_t IWhichObs[NumEst] = {
    1, 1, 1
  };

  TString NamEst[NumEst] = {
    "MuMu", "ElEl", "ElMu"
  };

  Int_t Flag1 = Flag%10; 
  Int_t Flag2 = (Flag%100)/10; 
  Int_t Flag3 = Flag/100; 

  // Steer the input depending on Flag
  if(Flag3 != 1 && Flag3 != 3){
    printf("B_AN2015_240: Not implemented Flag IGNORED %2i \n", Flag);
    return;
  };

  if(Flag2 == 1){
    DiluUnc = EDiluUnc;
    SystUnc = ESystUnc;
  }else if(Flag2 == 2 || Flag2 == 3){
    DiluUnc = ADiluUnc;
    SystUnc = ASystUnc;
  }else{
    printf("B_AN2015_240: Not implemented Flag IGNORED %2i \n", Flag);
    return;
  };

  if(Flag1 == 1){
    NumUnc = StatUnc;
  }else if(Flag1 ==  2 || Flag1 ==  3){
    NumUnc = DiluUnc;
  }else if(Flag1 ==  4 || Flag1 ==  5){
    NumUnc = SystUnc;
  }else if(Flag1 ==  6 || Flag1 ==  7){
    NumUnc = StatUnc + DiluUnc + SystUnc;
  }else{
    printf("B_AN2015_240: Not implemented Flag IGNORED %2i \n", Flag);
    return;
  };


  // Fill estimates for Asymmetry stat. uncertainties
  static const Int_t LenXEstStat = 2 * NumEst * (StatUnc + 1);
  Double_t XEstStat[LenXEstStat] = {
    -0.0060, 0.0058,
    -0.0012, 0.0083,
     0.0041, 0.0037,
    -0.0047, 0.0058,
    -0.0085, 0.0083,
    -0.0051, 0.0037
  };

  // Fill estimates for dtG stat. uncertainty
  Double_t XEstStatdtG[LenXEstStat] = {
    -0.4132, 0.4594,
    -0.1756, 0.5861,
     0.2086, 0.2955,
    -0.3778, 0.5482,
    -0.8321, 0.7496,
    -0.4261, 0.3187
  };

  // Fill estimates for Asymmetry dilution uncertainties - Table 36
  static const Int_t LenXEstEDilu = 2 * NumEst * (StatUnc + EDiluUnc + 1);
  TString NamUncEDilu[EDiluUnc+1] = {
    "stat.",
    "jet pT",
    "jet ang.",
    "charge",
    "b/bbar"
  };
  Double_t XEstEDilu[LenXEstEDilu]= {
    -0.0060, 0.0058, 0.0024, 0.0015, 0.0001, 0.0001, 
    -0.0012, 0.0083, 0.0065, 0.0029, 0.0010, 0.0000,
     0.0041, 0.0037, 0.0002, 0.0001, 0.0004, 0.0000,
    -0.0047, 0.0058, 0.0010, 0.0017, 0.0001, 0.0000,
    -0.0085, 0.0083, 0.0025, 0.0019, 0.0011, 0.0000,
    -0.0051, 0.0037, 0.0006, 0.0005, 0.0005, 0.0000
  // Asym.   stat.   jet pT  jet ang.charge  b/bbar 
  };

  // Fill estimates for Asymmetry dilution uncertainties - Table 45, 46, 47
  static const Int_t LenXEstADilu = 2 * NumEst * (StatUnc + ADiluUnc + 1);
  TString NamUncADilu[ADiluUnc + 1] = {
    "stat.",
    "dilution"
  };
  Double_t XEstADilu[LenXEstADilu]= {
    -0.0060, 0.0058, 0.0028, 
    -0.0012, 0.0083, 0.0066,
     0.0041, 0.0037, 0.0005,
    -0.0047, 0.0058, 0.0020,
    -0.0085, 0.0083, 0.0029,
    -0.0051, 0.0037, 0.0007
  // Asym.   stat.   dilution
  };

  // Fill estimates for dtG dilution uncertainties
  Double_t XEstDiludtG[LenXEstADilu]= {
    -0.4132, 0.4594, 0.2218,
    -0.1756, 0.5861, 0.4546,
     0.2086, 0.2955, 0.0398,
    -0.3778, 0.5482, 0.1836,
    -0.8321, 0.7496, 0.2556,
    -0.4261, 0.3187, 0.0625
  // dtG     stat.   dilution
  };

  // Fill estimates for Asymmetry sys. uncertainties - Table 42, 42, 44
  static const Int_t LenXEstESyst = 2 * NumEst * (StatUnc + ESystUnc + 1);
  TString NamUncESyst[ESystUnc + 1] = {
    "stat.",
    "Lumi.",
    "PileUp",
    "PDF",
    "Fac/Re",
    "ME-PS",
    "ISR",
    "FSR",
    "Unerly",
    "Color",
    "pT rew",
    "Trigger",
    "b-frag",
    "B had",
    "Had. model",
    "Muon",
    "Electron",
    "LepID",
    "LepIso",
    "JES",
    "JER",
    "MET",
    "BTagEff.",
    "BTagSF",
    "BG models"
  };
  Double_t XEstESyst[LenXEstESyst]= {
    -0.0060, 0.0058, 0.0000, 0.0002, 0.0019, 0.0002, 0.0035, 0.0033, 0.0019, 0.0055, 0.0045, 0.0000, 0.0000 ,0.0003, 0.0000, 0.0022, 0.0001, 0.0000, 0.0000, 0.0000, 0.0021, 0.0035, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0012, 0.0083, 0.0000, 0.0017, 0.0003, 0.0002, 0.0029, 0.0022, 0.0024, 0.0047, 0.0021, 0.0000, 0.0005, 0.0000, 0.0011, 0.0001, 0.0000, 0.0013, 0.0003, 0.0000, 0.0039, 0.0033, 0.0000, 0.0000, 0.0005, 0.0063,
     0.0041, 0.0037, 0.0000, 0.0001, 0.0003, 0.0002, 0.0023, 0.0011, 0.0005, 0.0011, 0.0013, 0.0000, 0.0003, 0.0000, 0.0006, 0.0000, 0.0002, 0.0007, 0.0000, 0.0000, 0.0005, 0.0000, 0.0000, 0.0000, 0.0000, 0.0009,
    -0.0047, 0.0058, 0.0000, 0.0013, 0.0003, 0.0006, 0.0008, 0.0012, 0.0007, 0.0015, 0.0039, 0.0000, 0.0003, 0.0001, 0.0005, 0.0000, 0.0002, 0.0000, 0.0000, 0.0000, 0.0002, 0.0023, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0085, 0.0083, 0.0000, 0.0002, 0.0003, 0.0002, 0.0028, 0.0015, 0.0033, 0.0030, 0.0054, 0.0001, 0.0004, 0.0001, 0.0028, 0.0001, 0.0000, 0.0029, 0.0002, 0.0000, 0.0036, 0.0019, 0.0000, 0.0000, 0.0001, 0.0063,
    -0.0051, 0.0037, 0.0000, 0.0001, 0.0004, 0.0002, 0.0015, 0.0003, 0.0012, 0.0023, 0.0009, 0.0000, 0.0001, 0.0001, 0.0002, 0.0000, 0.0003, 0.0004, 0.0000, 0.0000, 0.0006, 0.0001, 0.0000, 0.0000, 0.0000, 0.0009
  // Asym.   stat.   Lumi    PileUp  PDF     Fac/Re  ME-PS   ISR     FSR     Under   Color   pT rew  Trigger b-frag  B-had   Had.    Muon    Elec    LepID   LepIso  JES     JER     MET     bTagEff bTagSF  BG
  };

  // Fill estimates for Asymmetry sys. uncertainties - Table 45, 46, 47
  static const Int_t LenXEstASyst = 2 * NumEst * (StatUnc + ASystUnc + 1);
  TString NamUncASyst[ASystUnc + 1] = {
    "stat.",
    "sys."
  };
  Double_t XEstASyst[LenXEstASyst]= {
    -0.0060, 0.0058, 0.0111, 
    -0.0012, 0.0083, 0.0095,
     0.0041, 0.0037, 0.0028,
    -0.0047, 0.0058, 0.0064,
    -0.0085, 0.0083, 0.0108,
    -0.0051, 0.0037, 0.0029
  // Asym.   stat.   sys.
  };

  // Fill estimates for dtG sys. uncertainties
  Double_t XEstSystdtG[LenXEstASyst]= {
    -0.4132, 0.4594, 0.8646,
    -0.1756, 0.5861, 0.6543,
     0.2086, 0.2955, 0.2204,
    -0.3778, 0.5482, 0.5895,
    -0.8321, 0.7496, 0.9455,
    -0.4261, 0.3187, 0.2515
  // dtg     stat.   sys.
  };

  // Fill estimates for Asymmetry all uncertainties
  static const Int_t LenXEstEAll = 2 * NumEst * (StatUnc + EDiluUnc + ESystUnc + 1);
  TString NamUncEAll[StatUnc + EDiluUnc + ESystUnc] = {
    "stat.",
    "jet pT",
    "jet ang.",
    "charge",
    "b/bbar",
    "Lumi.",
    "PileUp",
    "PDF",
    "Fac/Re",
    "ME-PS",
    "ISR",
    "FSR",
    "Unerly",
    "Color",
    "pT rew",
    "Trigger",
    "b-frag",
    "B had",
    "Had. model",
    "Muon",
    "Electron",
    "LepID",
    "LepIso",
    "JES",
    "JER",
    "MET",
    "BTagEff.",
    "BTagSF",
    "BG models"
  };
  Double_t XEstEAll[LenXEstEAll]= {
    -0.0060, 0.0058, 0.0024, 0.0015, 0.0001, 0.0001, 0.0000, 0.0002, 0.0019, 0.0002, 0.0035, 0.0033, 0.0019, 0.0055, 0.0045, 0.0000, 0.0000 ,0.0003, 0.0000, 0.0022, 0.0001, 0.0000, 0.0000, 0.0000, 0.0021, 0.0035, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0012, 0.0083, 0.0065, 0.0029, 0.0010, 0.0000, 0.0000, 0.0017, 0.0003, 0.0002, 0.0029, 0.0022, 0.0024, 0.0047, 0.0021, 0.0000, 0.0005, 0.0000, 0.0011, 0.0001, 0.0000, 0.0013, 0.0003, 0.0000, 0.0039, 0.0033, 0.0000, 0.0000, 0.0005, 0.0063,
     0.0041, 0.0037, 0.0002, 0.0001, 0.0004, 0.0000, 0.0000, 0.0001, 0.0003, 0.0002, 0.0023, 0.0011, 0.0005, 0.0011, 0.0013, 0.0000, 0.0003, 0.0000, 0.0006, 0.0000, 0.0002, 0.0007, 0.0000, 0.0000, 0.0005, 0.0000, 0.0000, 0.0000, 0.0000, 0.0009,
    -0.0047, 0.0058, 0.0010, 0.0017, 0.0001, 0.0000, 0.0000, 0.0013, 0.0003, 0.0006, 0.0008, 0.0012, 0.0007, 0.0015, 0.0039, 0.0000, 0.0003, 0.0001, 0.0005, 0.0000, 0.0002, 0.0000, 0.0000, 0.0000, 0.0002, 0.0023, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0085, 0.0083, 0.0025, 0.0019, 0.0011, 0.0000, 0.0000, 0.0002, 0.0003, 0.0002, 0.0028, 0.0015, 0.0033, 0.0030, 0.0054, 0.0001, 0.0004, 0.0001, 0.0028, 0.0001, 0.0000, 0.0029, 0.0002, 0.0000, 0.0036, 0.0019, 0.0000, 0.0000, 0.0001, 0.0063,
    -0.0051, 0.0037, 0.0006, 0.0005, 0.0005, 0.0000, 0.0000, 0.0001, 0.0004, 0.0002, 0.0015, 0.0003, 0.0012, 0.0023, 0.0009, 0.0000, 0.0001, 0.0001, 0.0002, 0.0000, 0.0003, 0.0004, 0.0000, 0.0000, 0.0006, 0.0001, 0.0000, 0.0000, 0.0000, 0.0009
  };

  // Fill estimates for Asymmetry all uncertainties - Table 45, 46, 47
  static const Int_t LenXEstAAll = 2 * NumEst * (StatUnc + ADiluUnc + ASystUnc + 1);
  TString NamUncAAll[StatUnc + ADiluUnc + ASystUnc] = {
    "stat.",
    "dilution",
    "sys."
  };
  Double_t XEstAAll[LenXEstAAll]= {
    -0.0060, 0.0058, 0.0028, 0.0111, 
    -0.0012, 0.0083, 0.0066, 0.0095,
     0.0041, 0.0037, 0.0005, 0.0028,
    -0.0047, 0.0058, 0.0020, 0.0064,
    -0.0085, 0.0083, 0.0029, 0.0108,
    -0.0051, 0.0037, 0.0007, 0.0029
  // Asym.   stat.   dilu    sys.
  };

  // Fill estimates for dtG all uncertainties
  Double_t XEstAlldtG[LenXEstAAll]= {
    -0.4132, 0.4594, 0.2218, 0.8646,
    -0.1756, 0.5861, 0.4546, 0.6543,
     0.2086, 0.2955, 0.0398, 0.2204,
    -0.3778, 0.5482, 0.1836, 0.5895,
    -0.8321, 0.7496, 0.2556, 0.9455,
    -0.4261, 0.3187, 0.0625, 0.2515
  // dtG     stat.   dilu    sys.
  };

  // Correlation Matrix for un-correlated uncertainty
  static const Int_t LenCor = NumEst * NumEst;
  Double_t UnCor[LenCor] = {
    1.00, 0.00, 0.00,
    0.00, 1.00, 0.00,
    0.00, 0.00, 1.00
  };

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  if(PrintLog) myBlue->PrintStatus();

  // Fill all estimates
  Int_t ind = 0;
  if(Flag3 == 3) ind = 3;
  for(Int_t i = 0; i<NumEst; i++){
    if(Flag1 == 1){
      if(Flag2 == 3){
        myBlue->FillEst(i,&XEstStatdtG[ind]);
      }else{
        myBlue->FillEst(i,&XEstStat[ind]);
      }
    }else if(Flag1 == 2 || Flag1 == 3){
      if(Flag2 == 1){
        myBlue->FillEst(i,&XEstEDilu[ind]);
        myBlue->FillNamUnc(NamUncEDilu);
      }else if(Flag2 == 2){
        myBlue->FillEst(i,&XEstADilu[ind]);
        myBlue->FillNamUnc(NamUncADilu);
      }else{
        myBlue->FillEst(i,&XEstDiludtG[ind]);
        myBlue->FillNamUnc(NamUncADilu);
      }
    }else if(Flag1 == 4 || Flag1 == 5){
      if(Flag2 == 1){
        myBlue->FillEst(i,&XEstESyst[ind]);
        myBlue->FillNamUnc(NamUncESyst);
      }else if(Flag2 == 2){
        myBlue->FillEst(i,&XEstASyst[ind]);
        myBlue->FillNamUnc(NamUncASyst);
      }else{
        myBlue->FillEst(i,&XEstSystdtG[ind]);
        myBlue->FillNamUnc(NamUncASyst);
      }
    }else{
      if(Flag2 == 1){
        myBlue->FillEst(i,&XEstEAll[ind]);
        myBlue->FillNamUnc(NamUncEAll);
      }else if(Flag2 == 2){
        myBlue->FillEst(i,&XEstAAll[ind]);
        myBlue->FillNamUnc(NamUncAAll);
      }else{
        myBlue->FillEst(i,&XEstAlldtG[ind]);
        myBlue->FillNamUnc(NamUncAAll);
      }
    }
    ind = ind + NumUnc + 1;
  };

  // Set names
  myBlue->FillNamEst(NamEst);
  myBlue->FillNamObs(&FlagName2[Flag2-1]);

  // Fill all Correlations
  myBlue->FillCor(0, &UnCor[0]);
  for(Int_t k = 1; k<NumUnc; k++){
    if(Flag1 == 2 || Flag1 == 4 || Flag1 == 6){
      myBlue->FillCor(k, &UnCor[0]);
    }else if(Flag1 == 3 || Flag1 == 5 || Flag1 == 7){
      myBlue->FillCor(k, 1.0);
    }
  };

  // Fix input and inspect input
  myBlue->FixInp();
  if(PrintLog){
    myBlue->PrintEst();
    myBlue->PrintCor();
    myBlue->PrintCov();
    myBlue->PrintRho();
  };

  // Solve and inspect result
  myBlue->Solve();
  printf("\n\n**************************************************************************************** \n");
  printf("... B_AN2015_240: Calculate %s%s%s \n", FlagName3[Flag3-1].Data(),FlagName2[Flag2-1].Data(),FlagName1[Flag1-1].Data());
  printf("**************************************************************************************** \n\n\n");
  myBlue->PrintResult();
  if(PrintLog){
    myBlue->PrintRhoRes();
  }

  // Delete Object
  delete myBlue;
  myBlue = NULL;
  printf("\n\n\n\n");
  return;
}
