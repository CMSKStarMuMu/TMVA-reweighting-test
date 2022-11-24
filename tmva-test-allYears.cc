//g++ -O3 -o tmva-test-allYears tmva-test-allYears.cc  `root-config --cflags --libs` -lRooFit -lGenVector -lMathCore -lTMVA -l TMVAGui 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include "Riostream.h"
#include <map>
#include <string>
#include "RooGlobalFunc.h"
#include <vector>
#include <math.h>
#include <TGenericClassInfo.h> 
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TAttLine.h>
#include "TRatioPlot.h"
#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/BDT.h"
#include "TMVA/BDT_Reg.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/PDF.h"
#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
void tmva_test_bdt();
void tmva_evaluate_bdt();
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;

double mc_sigma = 0.040;
double mc_mass  = 5.27783;
double JPsiMass = 3.096916;
double nSigma_psiRej = 3.;
double tagged_mass=0;
double nsig_sw=0;
double mumuMass=-99;
double mumuMassE=-99;
double kstTrk1Pt;
double kstTrk2Pt;
double kstTrk1Eta;
double kstTrk2Eta;
Float_t weight=-99;

double kstTrk1PtD=-99;
double kstTrk2PtD=-99;
double kstTrk1EtaD=-99;
double kstTrk2EtaD=-99;
Double_t truthMatchMum=-99 ;
Double_t truthMatchMup=-99 ; 
Double_t truthMatchTrkm =-99;
Double_t truthMatchTrkp=-99 ;
Double_t trig=-99 ;
Double_t eventN=-99;
int pass_preselection=-99;
int num_threads=40;
std::string year=""; 
std::string year_default="2016";
std::string NameFileMCp0 = \
"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/MC_JPSI_2016_preBDT_Nov21.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/MC_JPSI_2016_preBDT_scale_add_vars_median_no_mass_window.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_MC_JPSI_scale_and_preselection_forreweighting_p0.root";
std::string NameFileMCp1 = \
"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/MC_JPSI_2016_preBDT_Nov21.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/MC_JPSI_2016_preBDT_scale_add_vars_median_no_mass_window.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_MC_JPSI_scale_and_preselection_forreweighting_p1.root";
std::string NameFileDatap0 = \
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_data_beforsel.root";
"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016Data_passPreselection_passSPlotCuts_mergeSweights.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_data_aftersel_p0.root";
std::string NameFileDatap1 = \
"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016Data_passPreselection_passSPlotCuts_mergeSweights.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_data_beforsel.root";
//"/gwpool/users/dini/FittoneNew3/DataMCreweighting-XGBV5_new/B0KstMuMu/reweight/Tree/final/XGBV5/2016/2016_data_aftersel_p1.root";
std::string datasetYear;
std::string dirfilexml;
std::string datasetname = "dataset";
std::string filexml= "/weights/TMVAClassification_BDT.weights.xml";
std::string OutputFileName="TMVA_ClassificationOutput.root";
std::vector<std::string> features = {"bVtxCL", "bLBS", "bLBSE" , "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt",\
                                     "kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE",\
				     "kstTrk2DCABS", "kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
// std::vector<std::string> features = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt",\
//                                      "kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE",\
// 				     "kstTrk2DCABS", "kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
//std::vector<std::string> features = {"kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta"};
std::vector<TH1D*> HistData;
std::vector<TH1D*> HistMC;
std::vector<TH1D*> HistMCW;
std::vector<TCanvas*> cstudies;			     
std::vector<TRatioPlot*> RatiosDataMC;			     
std::vector<TRatioPlot*> RatiosDataMCW;			     
int main (int argc, char** argv) {


  if( argc<=1 ){
    std::cout<<Form("Usage: %s [year=2016,2017,2018] {bdt}",argv[0])<<std::endl;
    std::cout<<Form("Please, set the year (at least)")<<std::endl;
    std::cout<<Form("example: %s 2016 bdt; if you want to train a new bdt & plots\n",argv[0])<<std::endl;
    std::cout<<Form("         %s 2016	 ; if you want just produce the  plots\n",argv[0])<<std::endl;
    exit(0);
  }else{
   if ((strcmp(argv[1],"2016") == 0 || \
        strcmp(argv[1],"2017") == 0 || \
	strcmp(argv[1],"2018") == 0)){
    year=argv[1];
    replaceAll( NameFileMCp0,  year_default, year);
    replaceAll( NameFileMCp1,  year_default, year);
    replaceAll( NameFileDatap0,  year_default, year);
    replaceAll( NameFileDatap1,  year_default, year);
    datasetYear=datasetname+year;
    dirfilexml=datasetYear+filexml;
    ROOT::EnableImplicitMT(num_threads);
    TMVA::Tools::Instance();
    std::cout<<Form("Setting year=%s",year.c_str())<<std::endl;
    std::cout<<Form("Setting dataset year=%s",datasetYear.c_str())<<std::endl;
    }else{
     std::cout<<Form("not recognize year=",year.c_str())<<std::endl;
     exit(0);
    }
	
  }  
   
  if( argc>2 ){
   if ((strcmp(argv[2],"bdt") == 0)){
    std::cout<<Form("Start training the BDT \n")<<std::endl;
    tmva_test_bdt();
    tmva_evaluate_bdt();
   }
  }else{
    tmva_evaluate_bdt();
  }
}
 void tmva_test_bdt(){
  gROOT ->Reset();
  gROOT->SetStyle("Plain");
//   ROOT::EnableImplicitMT(num_threads); 
//   TMVA::Tools::Instance();



//  TH1D* HxMassWn     = new TH1D( "HxMassWn"          , "B^{0} Mass sPlot <S0",100, 5., 5.6);
  
  TFile *fData = new TFile(NameFileDatap0.c_str(),"READ");
  TTree *TreeData     = (TTree*)fData->Get("ntuple");
  std::cout<<Form("Open Data File :%s \n",NameFileDatap0.c_str())<<std::endl;
  int nentriesData = (int)TreeData->GetEntries();
  std::cout<<Form("Found Data entries= = %d",nentriesData)<<std::endl;
  TFile *fMC   = new TFile(NameFileMCp0.c_str(),"READ");
  std::cout<<Form("Open MC File :%s \n",NameFileMCp0.c_str())<<std::endl;
  TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
  int nentriesMC = (int)TreeMC->GetEntries();
  std::cout<<Form("Found MC entries= = %d",nentriesMC)<<std::endl;
  DIR* dir = opendir(datasetYear.c_str());
  
  if(!dir) {
   gSystem->Exec(Form("mkdir -p %s",datasetYear.c_str()));
   std::cout<<Form("Create Dir  %s",datasetYear.c_str())<<std::endl;
  } 
  auto outputFile = TFile::Open(Form("%s/%s",datasetYear.c_str(),OutputFileName.c_str()), "RECREATE");

//  TMVA::Factory factory("TMVAClassification", outputFile,"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 
  TMVA::Factory factory =  TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
//  TMVA::Factory factory =  TMVA::Factory( "TMVAMulticlass", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
//  TMVA::Factory factory("TMVAClassification", outputFile,"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" ); 
  TMVA::DataLoader * loader = new TMVA::DataLoader(datasetYear.c_str());
//  TMVA::DataLoader * loader = new TMVA::DataLoader(Form("dataset");
//global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
   
// You can add an arbitrary number of signal or background trees
  TCut mycuts  = "pass_preselection==1&& abs(mumuMass - 3.096916) < 3*mumuMassE && eventN%2==0&&!TMath::IsNaN(bDCABSE)&&!TMath::IsNaN(bLBSE)&&abs(bCosAlphaBS)<=1&&sum_isopt_04<10"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
//  TCut mycuts  = "pass_preselection==1&& abs(mumuMass - 3.096916) < 3*mumuMassE &&!TMath::IsNaN(bDCABSE)&&!TMath::IsNaN(bLBSE)"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
//   TCut mycutb1 = "pass_preselection==1&& abs(mumuMass - 3.096916) < 3*mumuMassE && (trig==1)";
//   TCut mycutb2 = "(truthMatchMum==1)&&(truthMatchMup==1) && truthMatchTrkm==1 && truthMatchTrkp==1"; // for example: TCut mycutb = "abs(var1)<0.5";
//   TCut mycutb=mycutb1&&mycutb2;
  TCut mycutb="pass_preselection==1&& abs(mumuMass - 3.096916) < 3*mumuMassE && eventN%2==0 && (trig==1)&&(truthMatchMum==1)&&(truthMatchMup==1) && truthMatchTrkm==1 && truthMatchTrkp==1&&!TMath::IsNaN(bDCABSE)&&!TMath::IsNaN(bLBSE)&&abs(bCosAlphaBS)<=1&&sum_isopt_04<10";
//  TCut mycutb="pass_preselection==1&& abs(mumuMass - 3.096916) < 3*mumuMassE && (trig==1)&&(truthMatchMum==1)&&(truthMatchMup==1) && truthMatchTrkm==1 && truthMatchTrkp==1&&!TMath::IsNaN(bDCABSE)&&!TMath::IsNaN(bLBSE)";
// loader->AddSignalTree    ( TreeData,     signalWeight );
// loader->AddBackgroundTree( TreeMC  , backgroundWeight );
  loader->AddTree( TreeData,"Signal"	, signalWeight , mycuts   );
  loader->AddTree( TreeMC, "Background",backgroundWeight ,mycutb );
  for (Int_t i=0;i<=features.size()-1;i++) { 
   loader->AddVariable( features[i].c_str(), 'F' );
  }
//   
//   loader->AddVariable( "kstTrk1Pt", 'D' );
//   loader->AddVariable( "kstTrk2Pt", 'D' );
//   loader->AddVariable( "kstTrk1Eta", 'D' );
//   loader->AddVariable( "kstTrk2Eta", 'D' );
  
//   loader->AddVariable( "kstTrk1Pt", 'F' );
//   loader->AddVariable( "kstTrk2Pt", 'F' );
//   loader->AddVariable( "kstTrk1Eta", 'F' );
//   loader->AddVariable( "kstTrk2Eta", 'F' );
//   loader->AddSpectator( "truthMatchMum",  'F' );
//   loader->AddSpectator( "truthMatchMup",  'F' );
  loader->SetSignalWeightExpression("nsig_sw");
  loader->SetBackgroundWeightExpression("weight");
//  loader->AddVariable( "nsig_sw", 'D' );
// Apply additional cuts on the signal and background samples (can be different)
// Tell the factory how to use the training and testing events
//
// If no numbers of events are given, half of the events in the tree are used 
// for training, and the other half for testing:
//    loader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
// To also specify the number of testing events, use:
//    loader->PrepareTrainingAndTestTree( mycut,
//                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   loader->PrepareTrainingAndTestTree( "", "",
//   loader->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                  "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );
//                                    "nTrain_Signal=500000:nTrain_Background=500000:SplitMode=Random:NormMode=NumEvents:!V" );
                                      "SplitMode=Random:NormMode=NumEvents:!V" );
//Boosted Decision Trees
   factory.BookMethod( loader,  TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
//   factory.BookMethod( loader,  TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");

//   factory.BookMethod(loader,TMVA::Types::kBDT,"BDT","!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
//   factory.BookMethod(loader,TMVA::Types::kBDT,"BDT","!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
//   factory.BookMethod(loader,TMVA::Types::kBDT,"BDT","!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs:NegWeightTreatment=Pray" );
//   factory.BookMethod(loader,TMVA::Types::kBDT, "BDT","!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs:IgnoreNegWeightsInTraining" );
   factory.TrainAllMethods();

   factory.TestAllMethods();  

   factory.EvaluateAllMethods();  
   TCanvas *TRoc = factory.GetROCCurve(loader);
//    c1->Draw();
//    c1->Print("c1-tmva-test-allYears-2016.pdf");
   outputFile->Close();   
//   TCanvas *TRoc = factory.GetROCCurve (loader);
   TRoc->Print(Form("%s/roc-tmva-test-allYears-%s.pdf",datasetYear.c_str(),year.c_str()));
 } 
//
//=================================================================================================================================================================
//   
 void tmva_evaluate_bdt(){
 
   gROOT ->Reset();
   gROOT->SetStyle("Plain");

   TH1D* HxBDTData   	 = new TH1D( "HxBDTData"	, Form("Bdt Output Data %s",year.c_str()) ,100, 0., 1.0);
   TH1D* HxBDTMC   	 = new TH1D( "HxBDTMC"		, Form("Bdt Output MC %s",year.c_str()) ,100, 0., 1.0);
   TH1D* HxBDTMCW   	 = new TH1D( "HxBDTMCW"		, Form("Bdt Output MC %s reweighted",year.c_str()) ,100, 0., 1.0);
   std::vector<double> VarD;			     
   std::vector<float> VarF;

//    ROOT::EnableImplicitMT(num_threads);
//    TMVA::Tools::Instance();
//
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   for (Int_t i=0;i<=features.size()-1;i++) { 
    VarD.push_back(0);
    VarF.push_back(0);
   }
   for (Int_t i=0;i<=features.size()-1;i++) { 
    reader->AddVariable( features[i].c_str(), &VarF[i]);
    cstudies.push_back(new TCanvas(Form("c_%s",features[i].c_str()),Form("MC %s reweighting studies %s",features[i].c_str(),year.c_str()),200,10,900,780));
   }

//   TString weightfile = dir + prefix  + TString(".weights.xml");
   FILE *file = fopen(dirfilexml.c_str(), "r");
   if (file){
    std::cout<<"Reading BDT file "<<dirfilexml.c_str()<<std::endl;
    reader->BookMVA( "BDT",dirfilexml.c_str() );
   }else{
    std::cout<<dirfilexml.c_str()<<" not found; exit"<<std::endl;
    exit(1);
   } 

//=============== MC TTree ================== 
   TFile *fMC   = new TFile(NameFileMCp1.c_str(),"READ");
   std::cout<<Form("Opening MC File :%s \n",NameFileMCp1.c_str())<<std::endl;
   TTree *TreeMC     = (TTree*)fMC->Get("ntuple");
   for (Int_t i=0;i<=features.size()-1;i++) { 
     TreeMC->SetBranchAddress(features[i].c_str()     ,&VarD[i] );
     float hMin = TreeMC->GetMinimum(features[i].c_str());
     float hMax = TreeMC->GetMaximum(features[i].c_str());
     
     if(features[i]=="kstTrk1Pt"|| features[i]=="kstTrk2Pt") hMin=0;
     if(features[i]=="kstTrk1Pt"|| features[i]=="kstTrk2Pt") hMax=40;
     if(features[i]=="kstTrk1Eta"|| features[i]=="kstTrk2Eta") hMin=-3.14159;
     if(features[i]=="kstTrk1Eta"|| features[i]=="kstTrk2Eta") hMax= 3.14159;
     if(features[i]=="sum_isopt_04") hMin=0;
     if(features[i]=="sum_isopt_04") hMax=30;
     if(features[i]=="bDCABS") hMin=-0.015;
     if(features[i]=="bDCABS") hMax=0.015;
     if(features[i]=="bDCABSE") hMin=0.;
     if(features[i]=="bDCABSE") hMax=0.04;
     if(features[i]=="bCosAlphaBS") hMin=0.998;
     if(features[i]=="bCosAlphaBS") hMax=1.;
     if(features[i]=="bVtxCL") hMin=0.;
     if(features[i]=="bVtxCL") hMax=1.;
     if(features[i]=="kstTrk1DCABS") hMin=-1.;
     if(features[i]=="kstTrk1DCABS") hMax= 1.;
     if(features[i]=="kstTrk2DCABS") hMin=-1.;
     if(features[i]=="kstTrk2DCABS") hMax= 1.;
     if(features[i]=="kstTrk1DCABSE") hMin= 0.;
     if(features[i]=="kstTrk1DCABSE") hMax= 0.025;
     if(features[i]=="kstTrk2DCABSE") hMin= 0.;
     if(features[i]=="kstTrk2DCABSE") hMax= 0.025;
//   float hMin = 0;
//   float hMax = 35;
     HistData.push_back(new TH1D(  Form("Hx%s_Data" ,features[i].c_str()), Form("%s Data %s"	        ,features[i].c_str(),year.c_str()),100,hMin,hMax));
     HistMC.push_back(  new TH1D(  Form("Hx%s_MC"   ,features[i].c_str()), Form("%s MC %s" 	        ,features[i].c_str(),year.c_str()),100,hMin,hMax));
     HistMCW.push_back( new TH1D(  Form("Hx%s_MCW"  ,features[i].c_str()), Form("%s MC %s %s Reweighted",features[i].c_str(),year.c_str()),100,hMin,hMax));
   }
//   fMC->cd();
   TreeMC->SetBranchAddress("pass_preselection"     ,&pass_preselection );
   TreeMC->SetBranchAddress("tagged_mass"	    ,&tagged_mass );
   TreeMC->SetBranchAddress("mumuMass"  	    ,&mumuMass );
   TreeMC->SetBranchAddress("mumuMassE" 	    ,&mumuMassE );
   TreeMC->SetBranchAddress("truthMatchMum"	    ,&truthMatchMum );
   TreeMC->SetBranchAddress("truthMatchMup"	    ,&truthMatchMup );
   TreeMC->SetBranchAddress("truthMatchTrkm"	    ,&truthMatchTrkm );
   TreeMC->SetBranchAddress("truthMatchTrkp"	    ,&truthMatchTrkp );
   TreeMC->SetBranchAddress("trig"	            ,&trig );
   TreeMC->SetBranchAddress("weight"	            ,&weight );
   float val_mva=0;
   float val_mva_S=0;
   float val_mva_B=0;
   float rewe=0;
   int nentriesMC = (int)TreeMC->GetEntries(); 
   std::cout<<Form("Found MC entries= = %d",nentriesMC)<<std::endl;
//   for (Int_t i=0;i<100000;i++) {
   for (Int_t i=0;i<nentriesMC;i++) {
       TreeMC->GetEntry(i);
       if ( i%100000==0 ) {
          std::cout<<Form("Event %d",i)<<std::endl;
       }
       if( int(eventN)%2==0 ) continue; // parity!!!!!!!!!!!!!!!!!!!!!!111
       if(pass_preselection!=1) continue;
       if(trig!=1) continue;
       if(truthMatchMum !=1) continue;
       if(truthMatchMup !=1) continue;
       if(truthMatchTrkm!=1) continue; 
       if(truthMatchTrkp!=1) continue;
       if(fabs(mumuMass - JPsiMass) < nSigma_psiRej*mumuMassE){
        for (Int_t i=0;i<=features.size()-1;i++) { 
 	 if(TMath::IsNaN(VarD[i])) continue;
         VarF[i]=float(VarD[i]);
        }
 
 	val_mva_S = reader->EvaluateMulticlass( "BDT"	 )[0];
 	val_mva_B = reader->EvaluateMulticlass( "BDT"	 )[1];
 	rewe=val_mva_S/val_mva_B;
 	HxBDTMC ->Fill(val_mva_S);
 	HxBDTMCW->Fill(val_mva_S,rewe);
 	for (Int_t i=0;i<=features.size()-1;i++) {
 	 HistMC[i] ->Fill(VarD[i],weight);
 	 HistMCW[i]->Fill(VarD[i],rewe*weight);
 	}
      } 
    }  
//=============== Data TTree ================== 
    TFile *fData   = new TFile(NameFileDatap1.c_str(),"READ");
    std::cout<<Form("Opening Data File :%s \n",NameFileDatap1.c_str())<<std::endl;
//    fData->cd();
    TTree *TreeData     = (TTree*)fData->Get("ntuple");
    for (Int_t i=0;i<=features.size()-1;i++) { 
     TreeData->SetBranchAddress(features[i].c_str()     ,&VarD[i] );
    }
    TreeData->SetBranchAddress("pass_preselection" ,&pass_preselection );
    TreeData->SetBranchAddress("tagged_mass"	   ,&tagged_mass );
    TreeData->SetBranchAddress("nsig_sw"	   ,&nsig_sw );
    TreeData->SetBranchAddress("mumuMass"	   ,&mumuMass );
    TreeData->SetBranchAddress("mumuMassE"	   ,&mumuMassE );
    TreeData->SetBranchAddress("eventN"	           ,&eventN );
    int nentriesData = (int)TreeData->GetEntries();
    std::cout<<Form("Found Data entries= = %d",nentriesData)<<std::endl;
    for (Int_t i=0;i<nentriesData;i++) {
//    for (Int_t i=0;i<100000;i++) {
    	   TreeData->GetEntry(i);
 	   if ( i%100000==0 ) {
 	      std::cout<<Form("Event %d",i)<<std::endl;
 	    }
	   if( int(eventN)%2==0 ) continue; // parity!!!!!!!!!!!!!!!!!!!!!!111
    	   if(pass_preselection!=1) continue;
    	   if(fabs(mumuMass - JPsiMass) < nSigma_psiRej*mumuMassE){
 	    for (Int_t i=0;i<=features.size()-1;i++) {
 	     if(TMath::IsNaN(VarD[i])) continue;
 	     VarF[i]=float(VarD[i]);
 	    }
 	    val_mva_S = reader->EvaluateMulticlass( "BDT"    )[0];
    	    HxBDTData->Fill(val_mva_S,nsig_sw);
 	    for (Int_t i=0;i<=features.size()-1;i++) {
 	     HistData[i]->Fill(VarD[i],nsig_sw);
 	    }
     	   }
     }
 TCanvas* c_bdt = new TCanvas("c_bdt",Form("BDT score %s",year.c_str()),200,10,900,780);
 c_bdt->cd();
 HxBDTData->Sumw2();
 HxBDTMC->Sumw2();
 HxBDTMCW->Sumw2();
 HxBDTData->SetMaximum(HxBDTData->GetMaximum()*1.2);
 HxBDTData->Draw();
 HxBDTMC->Scale(HxBDTData->Integral()/HxBDTMC->Integral());
 HxBDTMC->SetLineColor(kRed);
 HxBDTMC->SetMarkerColor(kRed);
 HxBDTMC->Draw("same,Hist");
 HxBDTMCW->Scale(HxBDTData->Integral()/HxBDTMCW->Integral());
 HxBDTMCW->SetLineColor(kBlue);
 HxBDTMCW->SetMarkerColor(kBlue);
 HxBDTMCW->Draw("same,Hist");
 c_bdt->Print(Form("%s/tmva-bdt-score-%s.pdf",datasetYear.c_str(),year.c_str()));
 for (Int_t i=0;i<=features.size()-1;i++) {
   cstudies[i]->cd();
   HistData[i]->Sumw2();
   HistMC[i]->Sumw2();
   HistMCW[i]->Sumw2();
   HistData[i]->Draw();
   HistData[i]->SetMaximum(HistData[i]->GetMaximum()*1.2);
   HistMC[i]->Scale(HistData[i]->Integral()/HistMC[i]->Integral());
   HistMC[i]->SetLineColor(kRed);
   HistMC[i]->SetMarkerColor(kRed);
   HistMC[i]->Draw("same,Hist");
   HistMCW[i]->Scale(HistData[i]->Integral()/HistMCW[i]->Integral());
   HistMCW[i]->SetLineColor(kBlue);
   HistMCW[i]->SetMarkerColor(kBlue);
   HistMCW[i]->Draw("same,Hist");
   RatiosDataMC.push_back(new TRatioPlot(HistData[i],HistMC[i],"diffsig"));
   RatiosDataMC[i]->SetGraphDrawOpt("L");
   RatiosDataMC[i]->SetSeparationMargin(0.0);
   RatiosDataMC[i]->SetH1DrawOpt("E1");
   RatiosDataMC[i]->SetH2DrawOpt("HIST");
   RatiosDataMC[i]->Draw();
//    RatiosDataMC[i]->GetLowerRefGraph()->SetMinimum(RatiosDataMC[i]->GetLowerRefGraph()->GetMinimum()*2.);
//    RatiosDataMC[i]->GetLowerRefGraph()->SetMaximum(RatiosDataMC[i]->GetLowerRefGraph()->GetMaximum()*2.);
//   TGraph *g = RatiosDataMC[i]->GetLowerRefGraph();
   RatiosDataMCW.push_back(new TRatioPlot(HistData[i],HistMCW[i],"divsym"));
   RatiosDataMCW[i]->SetGraphDrawOpt("L");
   RatiosDataMCW[i]->SetSeparationMargin(0.0);
   RatiosDataMCW[i]->SetH1DrawOpt("E1");
   RatiosDataMCW[i]->SetH2DrawOpt("HIST");
//   RatiosDataMCW[i]->Draw("same");
//   RatiosDataMCW[i]->GetLowerRefGraph()->
   RatiosDataMCW[i]->Draw();
   RatiosDataMCW[i]->GetLowerRefGraph()->SetLineColor(kBlue);
   RatiosDataMCW[i]->GetLowerRefGraph()->SetMarkerColor(kBlue);
   RatiosDataMCW[i]->GetUpperPad()->cd();
   HistMC[i]->SetLineColor(kRed);
   HistMC[i]->SetMarkerColor(kRed);
//    HistMCW[i]->SetLineColor(kBlue);
//    HistMCW[i]->SetMarkerColor(kBlue);
//    HistData[i]->Draw("same");
   HistMC[i]->Draw("same,Hist");
   RatiosDataMCW[i]->GetLowerPad()->cd();
   RatiosDataMC[i]->GetLowerRefGraph()->SetLineColor(kRed);
   RatiosDataMC[i]->GetLowerRefGraph()->SetMarkerColor(kRed);
   RatiosDataMCW[i]->GetLowerRefGraph()->SetLineColor(kBlue);
   RatiosDataMCW[i]->GetLowerRefGraph()->SetMarkerColor(kBlue);
//    RatiosDataMC[i]->GetLowerRefGraph()->SetMinimum(RatiosDataMC[i]->GetLowerRefGraph()->GetMinimum()*1.4);
//    RatiosDataMC[i]->GetLowerRefGraph()->SetMaximum(RatiosDataMC[i]->GetLowerRefGraph()->GetMaximum()*1.4);
//    RatiosDataMCW[i]->GetLowerRefGraph()->SetMinimum(RatiosDataMC[i]->GetLowerRefGraph()->GetMinimum()*1.4);
//    RatiosDataMCW[i]->GetLowerRefGraph()->SetMaximum(RatiosDataMC[i]->GetLowerRefGraph()->GetMaximum()*1.4);
   RatiosDataMC[i]->GetLowerRefGraph()->Draw("same");
//   RatiosDataMCW[i]->GetLowerPad()->cd();
//   g->Draw();
//   HistMCW[i]->Draw("same,Hist");
//   cstudies[i]->Update();
   cstudies[i]->Print( Form("%s/tmva-studies-%s-%s.pdf",datasetYear.c_str()  ,features[i].c_str(),year.c_str()));

 }

}
//===============================================================================================================
void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}
