{
gStyle->SetNumberContours(99);
gStyle->SetOptStat(0);
TFile fileOutput1("Kaon_Lambda_Output_test.root", "recreate");
TFile file1("KL_flux_p_Pythia_24m.root"); //Macro for the momentum distribution of Kaon Long Beam

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TH1F* h24p = (TH1F*)file1.Get("h24p");

TH1F* h_K_Long = new TH1F("h_K_Long", "#Lambda + #pi Event Rate Over 100 Days; Momentum [GeV/c];Produced #Lambda + #pi+ Events", 12000, 0, 12);
//TH1F* h_K_Long = new TH1F("h_K_Long", "Convolution of Beam Flux and Cross Section; Momentum [GeV/c]; Normalized Yield [Arbitrary Units]", 12000, 0, 12); //normalised to cross section we care about distrbution not the y axis (where are the higher statistics)

TH1F* h_K_Long1 = new TH1F("h_K_Long1", "Momentum of Kaon L lambda pi +; Beam Momentum [GeV/c]; #sigma (#Lambda #pi+) [mb]", 12000, 0, 12);
TH1F* h_K_Long2 = new TH1F("h_K_Long2", "Momentum of Kaon L vs Lambda proton cross section; p [GeV/c]; #sigma (#Lambda p) [mb]", 12000, 0, 12);
TH1F* h_K_Long3 = new TH1F("h_K_Long3", "#Lambda + Proton Event Rate Over 100 Days; Momentum [GeV/c];Produced #Lambda + Proton Events", 12000, 0, 12);

TH1F* h_K_Long5 = new TH1F("h_K_Long5", "Event Rate; p [GeV/c]; Events [number of events]", 12000, 0, 12);
TH1F* h_K_Long6 = new TH1F("h_K_Long6", "Event Rate; p [GeV/c]; Events [number of events]", 12000, 0, 12);
TH1F* h_K_Long7 = new TH1F("h_K_Long7", "Event Rate; Beam Momentum [GeV/c]; Acceptance", 12000, 0, 12);
TH1F* h_K_Long8 = new TH1F("h_K_Long8", "Accepted #Lambda + #pi Event Rate Over 100 days; Momentum [GeV/c]; Accepted Events", 12000, 0, 12);
TH1F* h_K_Long9 = new TH1F("h_K_Long9", "lambda rho; #Lambda Momentum [GeV/c]; Accepted #Lambda", 12000, 0, 12);
TH1F* h_K_Long10 = new TH1F("h_K_Long10", "lambda rho ; #Lambda Momentum [GeV/c]; Produced #Lambda", 12000, 0, 12);
TH1F* h_K_Long11 = new TH1F("h_K_Long11", "Accepted; #Lambda Momentum [GeV/c]; Acceptance", 12000, 0, 12);
TH1F* h_K_Long12 = new TH1F("h_K_Long12", "Accepted Rescatered #Lambda + Proton Event Rate Over 100 Days; Momentum [GeV/c]; Measured #Lambda + Proton Events", 12000, 0, 12);

TH1F* h_Nick = new TH1F("h_Nick", "Event Rate; p [GeV/c]; Events [number of events]", 100000, 0, 12);


for (int k = 1; k < 12001; k++) {
    if (h_K_Long1->GetBinCenter(k) > 1) {
        h_K_Long1->SetBinContent(k, 3 / (h_K_Long1->GetBinCenter(k) * h_K_Long1->GetBinCenter(k)));
    }
    else {
        h_K_Long1->SetBinContent(k, 3);
    }
    h_K_Long2->SetBinContent(k, 10);
}

h_K_Long->Multiply(h_K_Long1, h24p);

h_K_Long->Scale(40 * (0.602 / 1.007947) * 0.071 * 100 * 24 * 60 * 60 / 1000);

h_K_Long3->Multiply(h_K_Long, h_K_Long2);

//Path length histograms
TH2D* h_t_path = new TH2D("h_t_path", "togethertransverse path length vs Phi; #Phi [rad]; Path Length [m]", 1000, 0, TMath::Pi(), 1000, 0, 0.04);
TH2D* h_t_path1 = new TH2D("h_t_path1", "exit from end cylinder path length vs Phi; #Theta [rad]; Path Length [m]", 1000, 0, 2, 1000, 0, 0.4);
TH2D* h_t_path2 = new TH2D("h_t_path2", "exit from sides path length vs Phi; #Theta [rad]; Path Length [m]", 1000, 0, 2, 1000, 0, 0.4);

TH2D* h_l_p = new TH2D("h_l_p", "longitudinal path length vs theta; #Theta [rad]; Longitudinal Path Length [m]", 1000, 0, 2, 1000, 0, 0.4);
TH2D* h_L = new TH2D("h_L", "Grand path length vs theta; #Theta [rad]; Total Path Length [m]", 1000, 0, 2/*TMath::Pi()*/, 1000, 0, 0.4);
TH2D* h_L_D = new TH2D("h_L_D", "Decay Grand path length vs theta; #Theta [rad]; Total Path Length (After Decay) [m]", 1000, 0, 2/*TMath::Pi()*/, 1000, 0, 0.4);
TH2D* h_KL_rho_L = new TH2D("h_KL_rho_L", "Grand path length vs momentum of kaon; momentum [GeV/c]; path length", 12000, 0, 12, 2000, 0, 0.4);

//TH1F* h_L = new TH1F("h_L", "p length y projection; p length; Counts", 200, 0, 1);
TF1* decay_path = new TF1("decay_path", "TMath::Exp(-x/([0]*[1]*30*0.26)); [cm]; ", 0, 41);

/*
TF1* decay_path_1 = new TF1("decay_path_1", "TMath::Exp(-x/(0.8*[1]*30*0.26)); [cm]; ", 0, 41);
TF1* decay_path_2 = new TF1("decay_path_2", "TMath::Exp(-x/(0.8*[1]*30*0.26)); [cm]; ", 0, 41);
TF1* decay_path_3 = new TF1("decay_path_3", "TMath::Exp(-x/(0.8*[1]*30*0.26)); [cm]; ", 0, 41);
TF1* decay_path_4 = new TF1("decay_path_4", "TMath::Exp(-x/(0.8*[1]*30*0.26)); [cm]; ", 0, 41);
TF1* decay_path_5 = new TF1("decay_path_5", "TMath::Exp(-x/(0.8*[1]*30*0.26)); [cm]; ", 0, 41);
*/

TH2D* h_xsection = new TH2D("h_xsection", "Cross Section Distribution; p [GeV/c]; Cross Section [mb];", 2000, 0, 15, 2000, 0, 1E9);
TH2D* h_N_total = new TH2D("h_N_total", "Total Number of #Lambda over 100 days; p [GeV/c]; Number of lambda events;", 2000, 0, 15, 2000, 0, 1E9);
TH2D* h_N_resc_total = new TH2D("h_N_resc_total", "Total Number of Rescattered #Lambda over 100 days; p [GeV/c]; Number of lambda events;", 2000, 0, 15, 2000, 0, 6E6);

//Luminosity
//TH2D* h_Luminosity = new TH2D("h_Luminosity", "Luminosity;  Path Length; Luminosity", 2000, 0, 0.4, 1000, 0, 1E30);
TH1F* h_Luminosity = new TH1F("h_Luminosity", "Luminosity; Kaon Momentum; Luminosity excluding flux", 1000, 0, 1);
//TH2D* myvertex = new TH2D("myvertex", "Generation of x-y Coordinates;x [cm];y [cm]", 100, -2, 2, 100, -2, 2);

TH1F* h_Lambda_N = new TH1F("h_Lambda_N", "#Lambda Events;  not sure; Number of #Lambda", 10000, 0, 1E8);

//Acceptance
TH2D* h_pion_plus_Ac = new TH2D("h_pion_plus_Ac", "Acceptance of #pi+;#pi+ Momentum [GeV/c];#pi+ #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());
TH2D* h_resc_Proton_Ac = new TH2D("h_resc_Proton_Ac", "Acceptance of rescattered proton;Re-scattered Proton Momentum [GeV/c];Re-scattered Proton #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());
TH2D* h_decay_Pi_Minus_Ac = new TH2D("h_decay_Pi_Minus_Ac", "Acceptance of decayed #pi-;Decay #pi- Momentum [GeV/c];Decay #pi- #Theta [rad];", 1000, 0, 2.5, 1000, 0, TMath::Pi());
TH2D* h_decay_Proton_Ac = new TH2D("h_decay_Proton_Ac", "Acceptance of the Decayed proton;Decay #pi- Momentum [GeV/c];Decay #pi- #Theta [rad];", 1000, 0, 10, 1000, 0, TMath::Pi());

//(x,y) histograms
TH2D* myvertex = new TH2D("myvertex", "Generation of x-y Coordinates;x [m];y [m]", 100, -2, 2, 100, -2, 2);
TH1F* myvertex_z = new TH1F("myvertex_z", "decay probabilty; z; Counts", 200, 0, 10);

//Lambda Histrograms
TH1F* h_Lambda = new TH1F("h_Lambda", "Momentum of #Lambda; #Lambda Momentum [GeV/c]; Counts", 200, 0, 10);
TH2D* h_Lambda_Angular = new TH2D("h_Lambda_Angular", "Angular dependence of #Lambda;#Lambda Momentum [GeV/c]; #Lambda #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());
TH2D* h_Lambda_Phi = new TH2D("h_Lambda_Phi", "Angular dependence of #Lambda; Momentum [GeV/c]; #Phi [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());

//Pion Plus Histrograms
TH1F* h_pion_plus = new TH1F("h_pion_plus", "Momentum of #pi+; p [GeV/c]; Counts", 1000, 0, 10);
TH1F* h_pion_plus_C = new TH1F("h_pion_plus_C", "Momentum of #pi+; p [GeV/c]; Counts", 1000, 0, 10);

TH2D* h_pion_plus_Angular = new TH2D("h_pion_plus_Angular", "Angular dependence of #pi+;#pi+ Momentum [GeV/c]; #pi+ #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi()); //the right 1000 was at 200 for this one and below
TH2D* h_pion_plus_Angular_C = new TH2D("h_pion_plus_Angular_C", "Angular dependence of #pi+ (Constrained); Momentum of Particle [GeV/c]; #Theta Depedence of Particle [radians];", 1000, 0, 12, 1000, 0, TMath::Pi());//^

TH2D* h_pion_plus_Phi = new TH2D("h_pion_plus_Phi", "Angular dependence of #pi+;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());
TH2D* h_pion_plus_Phi_C = new TH2D("h_pion_plus_Phi_C", "Angular dependence of #pi+;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());

//Rescattered lambda histograms
TH1F* h_resc_lambda_p = new TH1F("h_resc_lambda_p", "Momentum of Rescattered #Lambda; p [GeV/c]; Counts", 200, 0, 10);
TH2D* h_resc_lambda_angular = new TH2D("h_resc_lambda_angular", "Angular dependence of Rescattered #Lambda;Re-scattered #Lambda Momentum [GeV/c];Re-scattered #Lambda #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());
TH2D* h_resc_lambda_Phi = new TH2D("h_resc_lambda_Phi", "Angular dependence of  Rescattered #Lambda; Momentum [GeV/c]; #Phi Dependence [radians];", 200, 0, 12, 200, 0, TMath::Pi());

//rescattered proton histograms
TH1F* h_resc_Proton = new TH1F("h_resc_Proton", "Momentum of rescattered proton; p [GeV/c]; Counts", 1000, 0, 10);
TH1F* h_resc_Proton_C = new TH1F("h_resc_Proton_C", "Momentum of rescattered proton; p [GeV/c]; Counts", 1000, 0, 10);

TH2D* h_resc_Proton_Angular = new TH2D("h_resc_Proton_Angular", "Angular dependence of rescattered proton;Re-scattered Proton Momentum [GeV/c];Re-scattered Proton #Theta [rad];", 1000, 0, 12, 1000, 0, TMath::Pi());
TH2D* h_resc_Proton_Angular_C = new TH2D("h_resc_Proton_Angular_C", "Angular dependence of rescattered proton (Constrained); Momentum of Particle [GeV/c]; #Theta Depedence of Particle [radians];", 1000, 0, 12, 1000, 0, TMath::Pi());

TH2D* h_resc_Proton_Phi = new TH2D("h_resc_Proton_Phi", "Angular dependence of Rescattered Proton;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());
TH2D* h_resc_Proton_Phi_C = new TH2D("h_resc_Proton_Phi_C", "Angular dependence of Rescattered Proton;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());

//Decay Proton Histograms
TH1F* h_decay_Proton = new TH1F("h_decay_Proton", "Momentum of the decayed proton; p [GeV/c]; Counts", 1000, 0, 10);
TH1F* h_decay_Proton_C = new TH1F("h_decay_Proton_C", "Momentum of the decayed proton; p [GeV/c]; Counts", 1000, 0, 10);

TH2D* h_decay_Proton_Angular = new TH2D("h_decay_Proton_Angular", "Angular dependence of the Decayed proton;Decay Proton Momentum [GeV/c];Decay Proton #Theta [rad];", 1000, 0, 10, 1000, 0, TMath::Pi());
TH2D* h_decay_Proton_Angular_C = new TH2D("h_decay_Proton_Angular_C", "Angular dependence of the Decayed proton (Constrained); Momentum [GeV/c]; #Theta [rad];", 1000, 0, 10, 1000, 0, TMath::Pi());

TH2D* h_decay_Proton_Phi = new TH2D("h_decay_Proton_Phi", "Angular dependence of Decayed Proton;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());
TH2D* h_decay_Proton_Phi_C = new TH2D("h_decay_Proton_Phi_C", "Angular dependence of Decayed Proton;p [GeV/c]; #Phi [radians];", 200, 0, 10, 200, 0, TMath::Pi());

//Decay Pi minus Histograms
TH1F* h_decay_Pi_Minus = new TH1F("h_decay_Pi_Minus", "Momentum of decayed #pi-; p [GeV/c]; Counts", 200, 0, 4);
TH1F* h_decay_Pi_Minus_C = new TH1F("h_decay_Pi_Minus_C", "Momentum of decayed #pi-; p [GeV/c]; Counts", 200, 0, 4);

TH2D* h_decay_Pi_Minus_Angular = new TH2D("h_decay_Pi_Minus_Angular", "Angular dependence of decayed #pi-;Decay #pi- Momentum [GeV/c];Decay #pi- #Theta [rad];", 1000, 0, 2.5, 1000, 0, TMath::Pi());
TH2D* h_decay_Pi_Minus_Angular_C = new TH2D("h_decay_Pi_Minus_Angular_C", "Angular dependence of decayed #pi- (Constrained); Momentum [GeV/c]; #Theta [rad];", 1000, 0, 2.5, 1000, 0, TMath::Pi());

TH2D* h_decay_Pi_Minus_Phi = new TH2D("h_decay_Pi_Minus_Phi", "Angular dependence of Decayed #pi-;p [GeV/c]; #Phi [radians];", 200, 0, 4, 200, 0, TMath::Pi());
TH2D* h_decay_Pi_Minus_Phi_C = new TH2D("h_decay_Pi_Minus_Phi_C", "Angular dependence of Decayed #pi-;p [GeV/c]; #Phi [radians];", 200, 0, 4, 200, 0, TMath::Pi());


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creating particles and terms

// How many events to simulate and percentage completed

Int_t nevents = 1000000;
Int_t Percentage = nevents / 100;

// Creating TLorentzVectors of particles
TLorentzVector Beam, Target; // Beam and target

TLorentzVector* pion_plus, * Lambda;// First vertex particles
TLorentzVector* resc_Proton, * resc_Lambda; // Second vertex particles

TLorentzVector q;
TRandom3* mygen = new TRandom3();
TRandom3* mygenz = new TRandom3(0);
// Making Weights
Double_t Phasespace_Weight_1;
Double_t qWeight;

// Setting TLorentzVectors for beam and target in GeV (Px,Py,Pz,M)

Target.SetXYZM(0, 0, 0, 0.93827);

// Defining initial vertex and masses of particles
TLorentzVector V1 = Beam + Target;
Double_t Masses_1[2] = { 1.11568, 0.139570 }; // Lambda(0) + pi+(1) (primary vertex)
Double_t Masses_2[2] = { 1.11568, 0.93827 }; //  lambda + proton -> lambda(0) + proton(1)                            
Double_t Masses_3[2] = { 0.93827,0.139570 } // Lambda -> proton(0) and pi-(1)
TLorentzVector Inv_Lambda; // 4-vector of Lambda from proton + pion
//double W = V1.M();


// Creating decay vertices
TGenPhaseSpace Vertex_1, Vertex_2, Vertex_3;
Double_t BeamRandom;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t zverte;

double radius = 0.03;
double T_L = 0.4;
double T_R = 0.03;
double xverte;
double yverte;
double Lambda_Phi;
double Lambda_Theta;
double Lambda_Phi_V;
double alpha;
double t_path;
double l_p;
double z_exit;
double angle_counter;
double L;
double L_Choice;
double vertex_d;
double Lambda_Beta;
double Lambda_Gamma;
double GDecay;
double Luminosity;
double MeanPL;
double lambda_events;
double resc_events;
double resc_events_lumi;
double xsection;
double MeanDerivedPL;

// Loop over events
// Looping over simulated events
for (Long64_t i = 0; i < nevents; i++) {

    // Counter, shows percentage of simulated events completed
    if (i % Percentage == 0) {
        fprintf(stderr, "%lld\r", i / Percentage);
        fflush(stderr);
    }

    zverte = mygenz->Rndm()*0.4;
    //std::cout << zverte << endl;
    
    double Radius_F = 0.04;
    while (Radius_F > 0.03) {
    xverte = mygen->Gaus(0, 0.002);
    yverte = mygen->Gaus(0, 0.002);
    Radius_F = TMath::Sqrt(xverte * xverte + yverte * yverte);
    }

   myvertex->Fill(xverte, yverte);

   BeamRandom = h24p->GetRandom();

   Beam.SetXYZM(0, 0, BeamRandom, 0.497611);
   V1 = Beam + Target;
    // SetDecay(total energy, no. particles, mass array)
    // Setting the first decay
    h_K_Long->Fill(BeamRandom);
    if (!Vertex_1.SetDecay(V1, 2, Masses_1)) continue;
    h_Nick->Fill(BeamRandom);
    // Generating event and defining the phasespace weight for this decay
    Phasespace_Weight_1 = Vertex_1.Generate();

    Lambda = Vertex_1.GetDecay(0);
    pion_plus = Vertex_1.GetDecay(1);

    // Getting total energy for second decay from Lambda
    TLorentzVector V2 = (TLorentzVector)*Lambda + Target;

    // Setting the second decay
    if (!Vertex_2.SetDecay(V2, 2, Masses_2)) continue;
    // Defining the phasespace for this decay
    Double_t Phasespace_Weight_2 = Vertex_2.Generate();
    resc_Lambda = Vertex_2.GetDecay(0);//Lambda
    resc_Proton = Vertex_2.GetDecay(1);//proton
    //cout <<Lambda->E()+Target->E()- resc_Proton->E()-resc_Lambda->E()<< endl;
    TLorentzVector V3 = (TLorentzVector)*resc_Lambda;
    Vertex_3.SetDecay(V3, 2, Masses_3);

    Double_t Phasespace_Weight_3 = Vertex_3.Generate();
    //cout << Phasespace_Weight_1 << " " << Phasespace_Weight_2 << " " << Phasespace_Weight_3 << endl;
    decay_Proton = Vertex_3.GetDecay(0);
    decay_Pi_Minus = Vertex_3.GetDecay(1);

    Lambda_Theta = Lambda->Theta();
    Lambda_Beta = Lambda.Beta();
    Lambda_Gamma = Lambda.Gamma();
    Lambda_Phi = Lambda->Phi();
    Lambda_Phi_V = TMath::ATan2(yverte,xverte);
    vertex_d = TMath::Sqrt(xverte * xverte + yverte * yverte);
    alpha = TMath::Pi() - Lambda_Phi_V + Lambda_Phi;
    
    t_path = vertex_d * TMath::Cos(alpha) + TMath::Sqrt(vertex_d*vertex_d * (TMath::Cos(alpha) * TMath::Cos(alpha)) + radius*radius - vertex_d*vertex_d);
    
    l_p = t_path / TMath::Tan(Lambda_Theta);

    z_exit = l_p + zverte;
  
    decay_path->SetParameter(0, Lambda_Beta);
    decay_path->SetParameter(1, Lambda_Gamma);
   
    if (z_exit > 0.4) {
        l_p = T_L - zverte;
        L = l_p / TMath::Cos(Lambda_Theta);
        t_path = l_p * TMath::Tan(Lambda_Theta);
       // h_t_path1->Fill(Lambda_Theta, L);
    }
    else {
        L = t_path / TMath::Sin(Lambda_Theta);
        //h_t_path2->Fill(Lambda_Theta, L);
    }
    
    
    GDecay = decay_path->GetRandom() / 100;
    
    if (GDecay < L) {
        L_Choice = GDecay;
    }
    else if (GDecay > L) {
        L_Choice = L;
    }
    
    h_L->Fill(Lambda_Theta, L);
    h_L_D->Fill(Lambda_Theta, L_Choice);
    h_l_p->Fill(Lambda_Theta, l_p);
    h_t_path->Fill(Lambda_Phi, t_path);

    //1 vertex decay

    h_pion_plus->Fill(pion_plus->Rho(), Phasespace_Weight_1);
    h_pion_plus_Angular->Fill(pion_plus->Rho(), pion_plus->Theta(), Phasespace_Weight_1);
    h_pion_plus_Phi->Fill(pion_plus->Rho(), pion_plus->Phi(), Phasespace_Weight_1);

    h_Lambda->Fill(Lambda->Rho(), Phasespace_Weight_1);
    h_Lambda_Angular->Fill(Lambda->Rho(), Lambda_Theta, Phasespace_Weight_1);
    h_Lambda_Phi->Fill(Lambda->Rho(), Lambda->Phi(), Phasespace_Weight_1);    

    //2 vertex decay

    h_resc_lambda_p->Fill(resc_Lambda->Rho(), Phasespace_Weight_1 * Phasespace_Weight_2);
    h_resc_lambda_angular->Fill(resc_Lambda->Rho(), resc_Lambda->Theta(), Phasespace_Weight_1 * Phasespace_Weight_2);
    h_resc_lambda_Phi->Fill(resc_Lambda->Rho(), resc_Lambda->Phi(), Phasespace_Weight_1 * Phasespace_Weight_2);
    
    h_resc_Proton->Fill(resc_Proton->Rho(), Phasespace_Weight_1 * Phasespace_Weight_2);
    h_resc_Proton_Angular->Fill(resc_Proton->Rho(), resc_Proton->Theta(), Phasespace_Weight_1 * Phasespace_Weight_2);
    h_resc_Proton_Phi->Fill(resc_Proton->Rho(), resc_Proton->Phi(), Phasespace_Weight_1 * Phasespace_Weight_2);
    
    //3 vertex decay

    h_decay_Proton_Angular->Fill(decay_Proton->Rho(), decay_Proton->Theta(), Phasespace_Weight_2 * Phasespace_Weight_3);
    h_decay_Proton->Fill(decay_Proton->Rho(), Phasespace_Weight_2 * Phasespace_Weight_3);
    h_decay_Proton_Phi->Fill(decay_Proton->Rho(), decay_Proton->Phi(), Phasespace_Weight_2 * Phasespace_Weight_3);
    
    h_decay_Pi_Minus->Fill(decay_Pi_Minus->Rho(), Phasespace_Weight_2 * Phasespace_Weight_3);
    h_decay_Pi_Minus_Angular->Fill(decay_Pi_Minus->Rho(), decay_Pi_Minus->Theta(), Phasespace_Weight_2 * Phasespace_Weight_3);
    h_decay_Pi_Minus_Phi->Fill(decay_Pi_Minus->Rho(), decay_Pi_Minus->Phi(), Phasespace_Weight_2 * Phasespace_Weight_3);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    h_K_Long6->Fill(BeamRandom, Phasespace_Weight_2 * Phasespace_Weight_3);
    h_K_Long10->Fill(Lambda->Rho(), Phasespace_Weight_2 * Phasespace_Weight_3);

   
    //One big if statement
    if (pion_plus->Rho() > 0.2 && resc_Proton->Rho() > 0.2 && decay_Proton->Rho() > 0.2 && decay_Pi_Minus->Rho() > 0.2) {
        if (((pion_plus->Theta() > 3 * TMath::DegToRad() && pion_plus->Theta() < 15 * TMath::DegToRad()) || (pion_plus->Theta() > 20 * TMath::DegToRad() && pion_plus->Theta() < 165 * TMath::DegToRad()))){
            if (((resc_Proton->Theta() > 3 * TMath::DegToRad() && resc_Proton->Theta() < 15 * TMath::DegToRad()) || (resc_Proton->Theta() > 20 * TMath::DegToRad() && resc_Proton->Theta() < 165 * TMath::DegToRad()))) {
                if (((decay_Proton->Theta() > 3 * TMath::DegToRad() && decay_Proton->Theta() < 15 * TMath::DegToRad()) || (decay_Proton->Theta() > 20 * TMath::DegToRad() && decay_Proton->Theta() < 165 * TMath::DegToRad()))) {
                    if (((decay_Pi_Minus->Theta() > 3 * TMath::DegToRad() && decay_Pi_Minus->Theta() < 15 * TMath::DegToRad()) || (decay_Pi_Minus->Theta() > 20 * TMath::DegToRad() && decay_Pi_Minus->Theta() < 165 * TMath::DegToRad()))) {
                        h_pion_plus_Angular_C->Fill(pion_plus->Rho(), pion_plus->Theta(), Phasespace_Weight_1);
                        h_resc_Proton_Angular_C->Fill(resc_Proton->Rho(), resc_Proton->Theta(), Phasespace_Weight_1 * Phasespace_Weight_2);
                        h_decay_Proton_Angular_C->Fill(decay_Proton->Rho(), decay_Proton->Theta(), Phasespace_Weight_2 * Phasespace_Weight_3);
                        h_decay_Pi_Minus_Angular_C->Fill(decay_Pi_Minus->Rho(), decay_Pi_Minus->Theta(), Phasespace_Weight_2 * Phasespace_Weight_3);
                        h_K_Long5->Fill(BeamRandom, Phasespace_Weight_2 * Phasespace_Weight_3);
                        h_K_Long9->Fill(Lambda->Rho(), Phasespace_Weight_2 * Phasespace_Weight_3);
                    }
                }
            }
        }
    }
   
   
}




h_pion_plus_Ac->Divide(h_pion_plus_Angular_C, h_pion_plus_Angular);
h_resc_Proton_Ac->Divide(h_resc_Proton_Angular_C, h_resc_Proton_Angular);
h_decay_Pi_Minus_Ac->Divide(h_decay_Pi_Minus_Angular_C, h_decay_Pi_Minus_Angular);
h_decay_Proton_Ac->Divide(h_decay_Proton_Angular_C, h_decay_Proton_Angular);

/*
decay_path_1->SetParameter(1, 0.1);
decay_path_2->SetParameter(1, 0.3);
decay_path_3->SetParameter(1, 0.5);
decay_path_4->SetParameter(1, 0.7);
decay_path_5->SetParameter(1, 0.9);

decay_path_1->SetLineColorAlpha(1, 0.5);
decay_path_2->SetLineColorAlpha(2, 0.5);
decay_path_3->SetLineColorAlpha(3, 0.5);
decay_path_4->SetLineColorAlpha(4, 0.5);
decay_path_5->SetLineColorAlpha(5, 0.5);

decay_path_1->Draw();
decay_path_2->Draw("SAME");
decay_path_3->Draw("SAME");
decay_path_4->Draw("SAME");
decay_path_5->Draw("SAME");
*/

h_K_Long7->Divide(h_K_Long5, h_K_Long6); // acceptanee vs beam momentum

h_K_Long11->Divide(h_K_Long9, h_K_Long10); // acceptance vs lambda momentum

MeanPL = h_L_D->GetMean(2);

 // lambda pi + events for the 100 days
//h_K_Long->Draw();

h_K_Long3->Scale(MeanPL * 100 * (0.602 / 1.007947) * 0.071 / 1000); // lambda proton events for the 100 days // new flux for lambda proton events

h_K_Long8->Multiply(h_K_Long7, h_K_Long); // accepted lambda pi + events for the 100 days
//h_K_Long8->Draw();

h_K_Long12->Multiply(h_K_Long11, h_K_Long3); //accepted lambda proton events for the 100 days


break;

TCanvas* c1 = new TCanvas("c1", "Accepted Events", 800, 600);
c1->Divide(3, 1);
c1->cd(1);
h_K_Long->Draw();
c1->cd(2);
h_K_Long3->Draw();
c1->cd(3);
h_K_Long12->Draw();

TCanvas* c6 = new TCanvas("c6", "Acceptance", 1600, 800);
c6->Divide(2, 1);
c6->cd(1);
gPad->SetLogz();
h_Lambda_Angular->Draw("colz");
c6->cd(2);
h_resc_lambda_angular->Draw("colz");


TCanvas* c6 = new TCanvas("c6", "Acceptance", 1600, 800);
c6->Divide(2, 2);
c6->cd(1);
h_pion_plus_Ac->Draw("colz");
c6->cd(2);
h_resc_Proton_Ac->Draw("colz");
c6->cd(3);
h_decay_Proton_Ac->Draw("colz");
c6->cd(4);
h_decay_Pi_Minus_Ac->Draw("colz");

TCanvas* c6 = new TCanvas("c6", "Acceptance", 1600, 800);
c6->Divide(2, 1);
c6->cd(1);
h_K_Long10->Draw();
c6->cd(2);
h_K_Long9->Draw();

TCanvas* c1 = new TCanvas("c1", "#Lambda #Rho unconstrained vs constrained", 800, 600);
c1->Divide(2, 2);
c1->cd(1);
h_pion_plus_Angular->Draw("colz");
c1->cd(2);
h_resc_Proton_Angular->Draw("colz");
c1->cd(3);
h_decay_Proton_Angular->Draw("colz");
c1->cd(4);
h_decay_Pi_Minus_Angular->Draw("colz");

/*
TCanvas* c1 = new TCanvas("c1", "Accepted Events", 800, 600);
c1->Divide(2, 1);
c1->cd(1);
h_K_Long->Draw();
c1->cd(2);
h_K_Long8->Draw();
break;
*/
//# of particles calculation
/*
std::cout << "Mean N From Histogram = " << h_N_total->GetMean(2) << endl;
std::cout << "Mean Rescattered N From Histogram = " << h_N_resc_total->GetMean(2) << endl;

 
//h_Lambda_N->Draw();
 
lambda_events = MeanPL * 100 * 10000 * (6.02E23 / 1.007947) * 0.071 * 2.2E-28; //luminosity of lambda * cross section in cm squared, flux is 10000

resc_events_lumi = lambda_events * MeanPL * 100 * (6.02E23 / 1.007947) * 0.071;
resc_events = resc_events_lumi * 1E-26;

std::cout << MeanPL << endl;

std::cout << "Number of Lambda Particles Created = " << lambda_events << endl;
std::cout << "Number of Lambda Events over 100 Days = " << lambda_events * 100 * 24 * 60 * 60 << endl;

std::cout << "Number of Lambda Rescattering Events = " << resc_events << endl;
std::cout << "Number of Lambda Rescattering Events over 100 Days = " << resc_events * 100 * 24 * 60 * 60 << endl;
*/
//Getting Acceptance



TCanvas* c1 = new TCanvas("c1", "Accepted Events", 800, 600);
c1->Divide(3, 1);
c1->cd(1);
h_K_Long->Draw();
c1->cd(2);
h_K_Long3->Draw();
c1->cd(3);
h_K_Long12->Draw();
break;

//h_pion_plus_AcN->Multiply(h_pion_plus_Ac, h_luminosity);
/*
TCanvas* c1 = new TCanvas("c1", "Vertex 1", 800, 600);
c1->Divide(3, 1);
c1->cd(1);
h_t_path->Draw("colz");
c1->cd(2);
h_t_path1->Draw("colz");
c1->cd(3);
h_t_path2->Draw("colz");
break;

TCanvas* c1 = new TCanvas("c1", "Vertex 1", 800, 600);
c1->Divide(2, 1);
c1->cd(1);
gPad->SetLogz();
h_L->Draw("colz");
c1->cd(2);
h_L_D->Draw("colz");


TCanvas* c1 = new TCanvas("c1", "Vertex 1", 800, 600);
//c1->Divide(1, 1);
c1->cd(1);
h_Lambda->Draw();
c1->cd(2);
gPad->SetLogz();
h_Lambda_Angular->Draw("colz");
c1->cd(3);
h_Lambda_Phi->Draw("colz");
c1->cd(2);
h_pion_plus->Draw();
c1->cd(4);

h_pion_plus_Angular->Draw("colz");

c1->cd(6);
h_pion_plus_Phi->Draw("colz");

TCanvas* c2 = new TCanvas("c2", "Vertex 2", 800, 600);
c2->Divide(1, 1);
c2->cd(1);
h_resc_lambda_p->Draw();
c2->cd(3);
h_resc_lambda_angular->Draw("colz");
c2->cd(5);
h_resc_lambda_Phi->Draw("colz");
c2->cd(2);
h_resc_Proton->Draw();
c2->cd(4);
h_resc_Proton_Angular->Draw("colz");
c2->cd(6);
h_resc_Proton_Phi->Draw("colz");

TCanvas* c3 = new TCanvas("c3", "Vertex 3", 800, 600);
c3->Divide(2, 3);
c3->cd(1);
h_decay_Proton->Draw();
c3->cd(3);
h_decay_Proton_Angular->Draw("colz");
c3->cd(5);
h_decay_Proton_Phi->Draw("colz");
c3->cd(2);
h_decay_Pi_Minus->Draw();
c3->cd(4);
h_decay_Pi_Minus_Angular->Draw("colz");
c3->cd(6);
h_decay_Pi_Minus_Phi->Draw("colz");

TCanvas* c4 = new TCanvas("c4", "Path Length Calculations", 800, 600);
c4->Divide(2, 1);
c4->cd(1);
h_L->Draw("colz");
c4->cd(2);
h_L_D->Draw("colz");

TCanvas* c5 = new TCanvas("c5", "Comparison of Theta and Constrained Theta", 1920, 1080);
c5->Divide(2,4);
c5->cd(1);
h_decay_Pi_Minus_Angular->Draw("colz");
c5->cd(2);
h_decay_Pi_Minus_Angular_C->Draw("colz");
c5->cd(4);
h_decay_Proton_Angular_C->Draw("colz");
c5->cd(3);
h_decay_Proton_Angular->Draw("colz");
c5->cd(5);
h_pion_plus_Angular->Draw("colz");
c5->cd(6);
h_pion_plus_Angular_C->Draw("colz");
c5->cd(7);
h_resc_Proton_Angular->Draw("colz");
c5->cd(8);
h_resc_Proton_Angular_C->Draw("colz");

TCanvas* c6 = new TCanvas("c6", "Acceptance", 1600, 800);
c6->Divide(2, 2);
c6->cd(1);
h_pion_plus_Ac->Draw("colz");
c6->cd(2);
h_resc_Proton_Ac->Draw("colz");
c6->cd(3);
h_decay_Pi_Minus_Ac->Draw("colz");
c6->cd(4);
h_decay_Proton_Ac->Draw("colz");

TCanvas* c7 = new TCanvas("c7", "N Events Per second", 800, 600);
c7->Divide(2, 1);
c7->cd(1);
h_K_Long->Draw();//h_N_total->Draw("colz");
c7->cd(2);
h_N_resc_total->Draw("colz");
c7->cd(3);
//xsection->Draw("colz");

*/
TCanvas* c8 = new TCanvas("c8", "Events Over 100 Days vs accepted", 800, 600);
c8->Divide(3, 1);
c8->cd(1);
h_K_Long->Draw();
c8->cd(2);
h_K_Long3->Draw();
c8->cd(3);
h_K_Long12->Draw();



}
  
