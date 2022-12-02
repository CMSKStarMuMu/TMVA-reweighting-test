void produce_DRweight(int year, bool usePDG = true)
{

  double Fl  = 0.5534;
  double P1  = -0.0159;
  double P2  = -0.0012;
  double P3  = 0.237;
  double P4p = -0.9557;
  double P5p = -0.0064;
  double P6p = 0.0018;
  double P8p = -0.218;

  double FlMC =  5.9999e-01;
  double P1MC = -1.9816e-01;
  double P2MC = -2.7908e-04;
  double P3MC = -4.4279e-04;
  double P4pMC = -8.7036e-01;
  double P5pMC =  1.1506e-03;
  double P6pMC = -3.2764e-04;
  double P8pMC =  4.2577e-04;

  double FlPDG = 0.571;

  double FlVal = Fl;
  double P2Val = P2;
  if (usePDG) {
    FlVal = FlPDG;
    P2Val = 0;
  }

  auto fin = TFile::Open(Form("/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/%iMC_JPSI_noIP2D_addxcutvariable.root",year));
  auto tin = (TTree*)fin->Get("ntuple");

  double ctK, ctL, phi, DRweight;
  Long64_t eventN;
  tin->SetBranchAddress("eventN",&eventN);
  tin->SetBranchAddress("gen_cos_theta_k",&ctK);
  tin->SetBranchAddress("gen_cos_theta_l",&ctL);
  tin->SetBranchAddress("gen_phi_kst_mumu",&phi);

  auto fout = new TFile(Form("%iMC_JPSI%s.root",year,usePDG?"_PDGinputs":""),"RECREATE");
  auto tout = new TTree("ntuple_DRw", "MC decay-rate weights");
  tout->Branch("eventN", &eventN,"eventN/L");
  tout->Branch("DRweight", &DRweight,"DRweight/D");

  for (int i=0; i<tin->GetEntries(); ++i) {

    tin->GetEntry(i);

    double decMC = ( 0.75 * (1-FlMC) * (1-ctK*ctK) +
		     FlMC * ctK*ctK +
		     ( 0.25 * (1-FlMC) * (1-ctK*ctK) - FlMC * ctK*ctK ) * ( 2 * ctL*ctL -1 ) +
		     0.5 * P1MC * (1-FlMC) * (1-ctK*ctK) * (1-ctL*ctL) * cos(2*phi) +
		     2 * cos(phi) * ctK * sqrt(FlMC * (1-FlMC) * (1-ctK*ctK)) * ( P4pMC * ctL * sqrt(1-ctL*ctL) + P5pMC * sqrt(1-ctL*ctL) ) +
		     2 * sin(phi) * ctK * sqrt(FlMC * (1-FlMC) * (1-ctK*ctK)) * ( P8pMC * ctL * sqrt(1-ctL*ctL) - P6pMC * sqrt(1-ctL*ctL) ) +
		     2 * P2MC * (1-FlMC) * (1-ctK*ctK) * ctL -
		     P3MC * (1-FlMC) * (1-ctK*ctK) * (1-ctL*ctL) * sin(2*phi) );

    double decTarget = ( 0.75 * (1-FlVal) * (1-ctK*ctK) +
			 FlVal * ctK*ctK +
			 ( 0.25 * (1-FlVal) * (1-ctK*ctK) - FlVal * ctK*ctK ) * ( 2 * ctL*ctL -1 ) +
			 0.5 * P1 * (1-FlVal) * (1-ctK*ctK) * (1-ctL*ctL) * cos(2*phi) +
			 2 * cos(phi) * ctK * sqrt(FlVal * (1-FlVal) * (1-ctK*ctK)) * ( P4p * ctL * sqrt(1-ctL*ctL) + P5p * sqrt(1-ctL*ctL) ) +
			 2 * sin(phi) * ctK * sqrt(FlVal * (1-FlVal) * (1-ctK*ctK)) * ( P8p * ctL * sqrt(1-ctL*ctL) - P6p * sqrt(1-ctL*ctL) ) +
			 2 * P2Val * (1-FlVal) * (1-ctK*ctK) * ctL -
			 P3 * (1-FlVal) * (1-ctK*ctK) * (1-ctL*ctL) * sin(2*phi) );

    DRweight=decTarget/decMC;

    tout->Fill();

  }

  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();

  return;

}
