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
  // 0 == CEDM
  // 0: all uncertainties
  //
  // 1-7 == Asymmetry
  // 1: un-correlated stat. uncertainties
  // 2: un-correlated dilution uncertainties
  // 3: 100% correlated dilution uncertainties
  // 4: un-correlated sys. uncertainties
  // 5: 100% correlated sys. uncertainties
  // 6: un-correlated all uncertainties
  // 7: un-correlated stat. uncertainties and 100% correlated dilution and sys. uncertainties
  //-----------------------------------------------------------------------------

  // Print log option
  bool PrintLog = false;

  // Print Flag
  const string FlagName[8] = {
    "CEDM with all uncertianties",
    "Asymmetry with un-correlated stat. uncertainties",
    "Asymmetry with un-correlated dilution uncertainties",
    "Asymmetry with 100% correlated dilution uncertainties",
    "Asymmetry with un-correlated sys. uncertainties",
    "Asymmetry with 100% correlated sys. uncertainties",
    "Asymmetry with un-correlated all uncertainties",
    "Asymmetry with un-correlated stat. uncertainties and 100% correlated dilution and sys. uncertainties"
  };

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst  = 6;
  static const Int_t StatUnc = 1;
  static const Int_t DiluUnc = 4;
  static const Int_t SystUnc = 24;
  static const Int_t NumObs  = 2;
  Int_t NumUnc;

  // Index for which estimate determines which observable
  Int_t IWhichObs[NumEst] = {
    0, 0, 0,  // O_1 mm, ee, em 
    1, 1, 1,  // O_3 mm, ee, em
  };

  TString NamEst[NumEst] = {
    "O_1 MuMu", "O_1 ElEl", "O_1 ElMu",
    "O_3 MuMu", "O_3 ElEl", "O_3 ElMu"
  };

  TString NamObsCEDM[NumObs] = {"CEDM O_1",      "CEDM O_3"};
  TString NamObsAsym[NumObs] = {"Asymmetry O_1", "Asymmetry O_3"};

  // Steer the input depending on Flag
  if(Flag == 0 || Flag == 1){
    NumUnc = StatUnc;
  }else if(Flag ==  2 || Flag ==  3){
    NumUnc = DiluUnc;
  }else if(Flag ==  4 || Flag ==  5){
    NumUnc = SystUnc;
  }else if(Flag ==  6 || Flag ==  7){
    NumUnc = StatUnc + DiluUnc + SystUnc;
  }else{
    printf("B_AN2015_240: Not implemented Flag IGNORED %2i \n", Flag);
    return;
  }

  // Fill estimates for CEDM stat. uncertainty
  static const Int_t LenXEstStat = NumEst * (StatUnc+1);
  Double_t XEstCEDM[LenXEstStat] = {
    -1.20, 2.91,
    -0.51, 2.75,
     0.60, 1.07,
    -1.10, 2.35,
    -2.41, 3.58,
    -1.23, 1.19
  };

  // Fill estimates for Asymmetry stat. uncertainties
  Double_t XEstStat[LenXEstStat] = {
    -0.0060, 0.0058,
    -0.0012, 0.0083,
     0.0041, 0.0037,
    -0.0047, 0.0058,
    -0.0085, 0.0083,
    -0.0051, 0.0037
  };

  // Fill estimates for Asymmetry dilution uncertainties - Table 36
  static const Int_t LenXEstDilu = NumEst * (DiluUnc+1);
  TString NamUncDilu[DiluUnc] = {
    "jet pT resolution",
    "jet angular resolution",
    "charge mis-id",
    "b/bbar response"
  };
  Double_t XEstDilu[LenXEstDilu]= {
    -0.0060, 0.0024, 0.0015, 0.0001, 0.0001, 
    -0.0012, 0.0065, 0.0029, 0.0010, 0.0000,
     0.0041, 0.0002, 0.0001, 0.0004, 0.0000,
    -0.0047, 0.0010, 0.0017, 0.0001, 0.0000,
    -0.0085, 0.0025, 0.0019, 0.0011, 0.0000,
    -0.0051, 0.0006, 0.0005, 0.0005, 0.0000
  // Asym.   jet pT  angular charge  b/bbar 
  };

  // Fill estimates for Asymmetry syst.  uncertainties
  static const Int_t LenXEstSyst = NumEst * (SystUnc+1);
  TString NamUncSyst[SystUnc] = {
    "Luminosity",
    "Pile Up",
    "PDF",
    "Factorization and Renormalization",
    "ME-PS matching",
    "ISR",
    "FSR",
    "Unerlying event",
    "Color reconnection",
    "top pT reweighting",
    "Trigger SF",
    "b-fragmentation",
    "Semileptonic BR of B hadron",
    "Hadronization model",
    "Muon correction",
    "Electron energy scale and smearing",
    "Lepton ID",
    "Lepton Iso",
    "JES",
    "JER",
    "Unclustered energy of MET",
    "BTagging Eff.",
    "BTagging SF",
    "Limited statistics of BG models"
  };
  Double_t XEstSyst[LenXEstSyst]= {
    -0.0060, 0.0000, 0.0002, 0.0019, 0.0002, 0.0035, 0.0033, 0.0019, 0.0055, 0.0045, 0.0000, 0.0000 ,0.0003, 0.0000, 0.0022, 0.0001, 0.0000, 0.0000, 0.0000, 0.0021, 0.0035, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0012, 0.0000, 0.0017, 0.0003, 0.0002, 0.0029, 0.0022, 0.0024, 0.0047, 0.0021, 0.0000, 0.0005, 0.0000, 0.0011, 0.0001, 0.0000, 0.0013, 0.0003, 0.0000, 0.0039, 0.0033, 0.0000, 0.0000, 0.0005, 0.0063,
     0.0041, 0.0000, 0.0001, 0.0003, 0.0002, 0.0023, 0.0011, 0.0005, 0.0011, 0.0013, 0.0000, 0.0003, 0.0000, 0.0006, 0.0000, 0.0002, 0.0007, 0.0000, 0.0000, 0.0005, 0.0000, 0.0000, 0.0000, 0.0000, 0.0009,
    -0.0047, 0.0000, 0.0013, 0.0003, 0.0006, 0.0008, 0.0012, 0.0007, 0.0015, 0.0039, 0.0000, 0.0003, 0.0001, 0.0005, 0.0000, 0.0002, 0.0000, 0.0000, 0.0000, 0.0002, 0.0023, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0085, 0.0000, 0.0002, 0.0003, 0.0002, 0.0028, 0.0015, 0.0033, 0.0030, 0.0054, 0.0001, 0.0004, 0.0001, 0.0028, 0.0001, 0.0000, 0.0029, 0.0002, 0.0000, 0.0036, 0.0019, 0.0000, 0.0000, 0.0001, 0.0063,
    -0.0051, 0.0000, 0.0001, 0.0004, 0.0002, 0.0015, 0.0003, 0.0012, 0.0023, 0.0009, 0.0000, 0.0001, 0.0001, 0.0002, 0.0000, 0.0003, 0.0004, 0.0000, 0.0000, 0.0006, 0.0001, 0.0000, 0.0000, 0.0000, 0.0009
  // Asym.   Lumi    PileUp  PDF     Fac/Re  ME-PS   ISR     FSR     Under   Color   pT rew  Trigger b-frag  B-had   Had.    Muon    Elec    LepID   LepIso  JES     JER     MET     bTagEff bTagSF  BG
  };

  // Fill estimates for Asymmetry all uncertainties
  static const Int_t LenXEstAll = NumEst * (StatUnc + DiluUnc + SystUnc + 1);
  TString NamUncAll[29] = {
    "statistic",
    "jet pT resolution",
    "jet angular resolution",
    "charge mis-id",
    "b/bbar response"
    "Luminosity",
    "Pile Up",
    "PDF",
    "Factorization and Renormalization",
    "ME-PS matching",
    "ISR",
    "FSR",
    "Unerlying event",
    "Color reconnection",
    "top pT reweighting",
    "Trigger SF",
    "b-fragmentation",
    "Semileptonic BR of B hadron",
    "Hadronization model",
    "Muon correction",
    "Electron energy scale and smearing",
    "Lepton ID",
    "Lepton Iso",
    "JES",
    "JER",
    "Unclustered energy of MET",
    "BTagging Eff.",
    "BTagging SF",
    "Limited statistics of BG models"
  };
  Double_t XEstAll[LenXEstAll]= {
    -0.0060, 0.0058, 0.0024, 0.0015, 0.0001, 0.0001, 0.0000, 0.0002, 0.0019, 0.0002, 0.0035, 0.0033, 0.0019, 0.0055, 0.0045, 0.0000, 0.0000 ,0.0003, 0.0000, 0.0022, 0.0001, 0.0000, 0.0000, 0.0000, 0.0021, 0.0035, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0012, 0.0083, 0.0065, 0.0029, 0.0010, 0.0000, 0.0000, 0.0017, 0.0003, 0.0002, 0.0029, 0.0022, 0.0024, 0.0047, 0.0021, 0.0000, 0.0005, 0.0000, 0.0011, 0.0001, 0.0000, 0.0013, 0.0003, 0.0000, 0.0039, 0.0033, 0.0000, 0.0000, 0.0005, 0.0063,
     0.0041, 0.0037, 0.0002, 0.0001, 0.0004, 0.0000, 0.0000, 0.0001, 0.0003, 0.0002, 0.0023, 0.0011, 0.0005, 0.0011, 0.0013, 0.0000, 0.0003, 0.0000, 0.0006, 0.0000, 0.0002, 0.0007, 0.0000, 0.0000, 0.0005, 0.0000, 0.0000, 0.0000, 0.0000, 0.0009,
    -0.0047, 0.0058, 0.0010, 0.0017, 0.0001, 0.0000, 0.0000, 0.0013, 0.0003, 0.0006, 0.0008, 0.0012, 0.0007, 0.0015, 0.0039, 0.0000, 0.0003, 0.0001, 0.0005, 0.0000, 0.0002, 0.0000, 0.0000, 0.0000, 0.0002, 0.0023, 0.0000, 0.0000, 0.0003, 0.0045,
    -0.0085, 0.0083, 0.0025, 0.0019, 0.0011, 0.0000, 0.0000, 0.0002, 0.0003, 0.0002, 0.0028, 0.0015, 0.0033, 0.0030, 0.0054, 0.0001, 0.0004, 0.0001, 0.0028, 0.0001, 0.0000, 0.0029, 0.0002, 0.0000, 0.0036, 0.0019, 0.0000, 0.0000, 0.0001, 0.0063,
    -0.0051, 0.0037, 0.0006, 0.0005, 0.0005, 0.0000, 0.0000, 0.0001, 0.0004, 0.0002, 0.0015, 0.0003, 0.0012, 0.0023, 0.0009, 0.0000, 0.0001, 0.0001, 0.0002, 0.0000, 0.0003, 0.0004, 0.0000, 0.0000, 0.0006, 0.0001, 0.0000, 0.0000, 0.0000, 0.0009
  };

  // Correlation Matrix for un-correlated uncertainty
  static const Int_t LenCor = NumEst * NumEst;
  Double_t UnCor[LenCor] = {
    1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 1.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 1.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 1.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 1.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 1.00
  };

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  if(PrintLog) myBlue->PrintStatus();

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    if(Flag == 0){
      myBlue->FillEst(i,&XEstCEDM[ind]);
    }else if(Flag == 1){
      myBlue->FillEst(i,&XEstStat[ind]);
    }else if(Flag == 2 || Flag == 3){
      myBlue->FillEst(i,&XEstDilu[ind]);
      myBlue->FillNamUnc(NamUncDilu);
    }else if(Flag == 4 || Flag == 5){
      myBlue->FillEst(i,&XEstSyst[ind]);
      myBlue->FillNamUnc(NamUncSyst);
    }else{
      myBlue->FillEst(i,&XEstAll[ind]);
      myBlue->FillNamUnc(NamUncAll);
    }
    ind = ind + NumUnc + 1;
  }

  // Set names
  myBlue->FillNamEst(NamEst);
  if(Flag == 0){
    myBlue->FillNamObs(NamObsCEDM);
  }else{
    myBlue->FillNamObs(NamObsAsym);
  }

  // Fill all Correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(Flag == 0 || Flag == 1 || Flag == 2 || Flag == 4 || Flag == 6){
      myBlue->FillCor(k, &UnCor[0]);
    }else if(Flag == 3 || Flag == 5){
      myBlue->FillCor(k, 1.0);
    }else{
      if(k == 1){
        myBlue->FillCor(k, &UnCor[0]);
      }else{
        myBlue->FillCor(k, 1.0);
      }
    }
  }

  // Fix input and inspect input
  myBlue->FixInp();
  if(PrintLog){
    myBlue->PrintEst();
    if(Flag == 0){
      myBlue->PrintCor();
      myBlue->PrintCov();
      myBlue->PrintRho();
    };
  };

  // Solve and inspect result
  myBlue->Solve();
  printf("**************************************************************************************** \n");
  printf("... B_AN2015_240: Calculate %s \n", FlagName[Flag].c_str());
  printf("**************************************************************************************** \n");
  myBlue->PrintResult();
  printf("**************************************************************************************** \n");
  if(Flag == 0 && PrintLog){
    myBlue->PrintRhoRes();
  }

  // Delete Object
  delete myBlue;
  myBlue = NULL;
  printf("\n\n\n\n");
  return;
}
