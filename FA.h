#ifndef FA_H
#define FA_H

//#include "TParameter.h"
#include <string>

//#include "include/GeneralIncludes.h"


namespace BBBA07{

//================== Parameters ==================
// q2 = -Q2
MnvH1D* h_q2qe_base;
MnvH1D* h_flux;
MnvH1D* h_enu;

MnvH2D* h_enu_q2qe_mc;
MnvH2D* h_enu_q2qe_data;
MnvH2D* h_enu_q2qe;

vector<vector<double>> par(20,vector<double>(14,0));
vector<vector<double>> par2(20,vector<double>(5,0));
vector<double> EnuMin({1.50296, 1.50629, 1.51129, 1.51795, 1.52461, 1.54126, 1.56791, 
    1.59455, 1.63452, 1.82687, 2.19686, 2.62256, 2.99946, 3.34313, 
    4.10667, 5.86014, 7.91494, 10.5563});

vector<double> EnuMax({20.003, 20.0063, 20.0113, 20.0179, 20.0246, 20.0413, 20.0679, 
    20.0946, 20.1345, 20.1878, 20.2677, 20.3743, 20.4809, 20.5875, 
    20.8539, 21.6, 22.6658, 24.2644});


//vector<double> EnuMax({
//    10.003, 10.0063, 10.0113, 10.0179, 10.0246, 10.0413, 10.0679, 
//    10.0946, 10.1345, 10.1878, 10.2677, 10.3743, 10.4809, 10.5875, 
//    10.8539, 11.6, 12.6658, 14.2644});

double scale=1.81348205392007398/2.0636;
bool doScale=false;

double Mp = 0.9382720881;
double Mn = 0.93956542052;

//double e = 1.602176634*pow(10,-19);
//double e = 1;
double alpha = 1/137.;
double e = 0.30282212;
double hbar= 6.58212*pow(10,-25); //GeVs

//double uN = e*hbar/(2*Mp);
//double uN = 1/(2*Mp);
double uN = 1;
double up = 2.7930*uN;
double un = -1.913042*uN;

double gA = -1.267;
//double Gf = 1.1803*pow(10,-5);
double Gf = 1.1663787*pow(10,-5);
//double cosThetaC = 0.9742;
double cosThetaC = 0.974000001;
double XI = 3.706*uN;
double Mv = 0.860;
double MVsq = Mv*Mv;//GeV^2

double Mmu = .105658372;
double Mpion = .1395702;



vector<double> aE({1,-.24});
vector<double> aM({1,0.1717});
vector<double> bE({10.98,12.82,21.97});
vector<double> bM({11.26,19.32,8.33});

vector<double>XiK({0,1/6., 1/3., 1/2., 2/3., 5/6., 1 } );
vector<double>AEp({1.,0.9927,0.9898,0.9975,0.9812,0.9340,1.});
vector<double>AMp({1.,1.0011,0.9992,0.9974,1.0010,1.0003,1.});
vector<double>AEpd({1.,0.9839,0.9632,0.9748,0.9136,0.5447,-0.2682});
vector<double>AMpd({1.,0.9916,0.9771,0.9801,1.0321,1.0429,0.5084});
vector<double>AMn25({1.,0.9958,0.9877,1.0193,1.0350,0.9164,0.7300});
vector<double>AMn43({1.,0.9958,0.9851,1.0187,1.0307,0.9080,0.9557});
vector<double>AEn25({1.,1.1011,1.1392,1.0203,1.1093,1.5429,0.9706});
vector<double>AEn43({1.,1.1019,1.1387,1.0234,1.1046,1.5395,1.2708});
vector<double>AFA25D({1.,0.9207,0.9795,1.0480,1.0516,1.2874,0.7707});

map<string, vector<double>*> Aparameters;


void DeclareFittingPars()
{
  par[1]=vector<double>({0.80054138,-0.00649664,-0.00323802,-0.00141958,0.00001152,5.74e-6,0.00070979,-5.76e-6,-7.4e-6,-0.00283917,0.00001152,0.79983159,-0.00649088,-0.00017028});
  par[2]=vector<double>({0.80196097,-0.00935488,-0.00466261,-0.00425875,0.00004968,0.00002476,0.00212938,-0.00002484,-2.33e-6,-0.0085175,0.00004968,0.79983159,-0.00933005,0.00062957});
  par[3]=vector<double>({0.80409034,-0.01366115,-0.00680892,-0.0085175,0.00014471,0.00007212,0.00425875,-0.00007235,0.00001475,-0.01703501,0.00014471,0.79983159,-0.0135888,0.00282707});
  par[4]=vector<double>({0.80692951,-0.01943811,-0.00968825,-0.01419584,0.00034196,0.00017044,0.00709792,-0.00017098,0.00005522,-0.02839168,0.00034196,0.79983159,-0.01926713,0.0076195});
  par[5]=vector<double>({0.80976868,-0.02525539,-0.01258766,-0.01987417,0.00061984,0.00030894,0.00993709,-0.00030992,0.00011591,-0.03974835,0.00061984,0.79983159,-0.02494547,0.01454043});
  par[6]=vector<double>({0.8168666,-0.03997495,-0.01992411,-0.03407001,0.00166728,0.000831,0.01703501,-0.00083364,0.0003561,-0.06814002,0.00166728,0.79983159,-0.0391413,0.04115491});
  par[7]=vector<double>({0.82822327,-0.0640503,-0.03192362,-0.05678335,0.00439132,0.0021887,0.02839168,-0.00219566,0.00100326,-0.1135667,0.00439132,0.79983159,-0.06185464,0.11140849});
  par[8]=vector<double>({0.83957994,-0.08877067,-0.04424462,-0.07949669,0.00840536,0.00418936,0.03974835,-0.00420268,0.00197395,-0.15899338,0.00840536,0.79983159,-0.08456798,0.21571796});
  par[9]=vector<double>({0.85661494,-0.1270606,-0.06332889,-0.1135667,0.0168452,0.0083959,0.05678335,-0.0084226,0.0040366,-0.2271334,0.0168452,0.79983159,-0.118638,0.43603697});
  par[10]=vector<double>({0.87932828,-0.18037136,-0.08989977,-0.15899338,0.03261336,0.01625499,0.07949669,-0.01630668,0.00791913,-0.31798676,0.03261336,0.79983159,-0.16406468,0.84899129});
  par[11]=vector<double>({0.91339829,-0.26517504,-0.13216719,-0.2271334,0.06594068,0.03286582,0.1135667,-0.03297034,0.01616936,-0.4542668,0.06594068,0.79983159,-0.2322047,1.72384197});
  par[12]=vector<double>({0.95882497,-0.38727669,-0.19302447,-0.31798676,0.12843727,0.06401505,0.15899338,-0.06421863,0.03169903,-0.63597352,0.12843727,0.79983159,-0.32305806,3.36709208});
  par[13]=vector<double>({1.00425165,-0.51969843,-0.25902544,-0.40884012,0.21157403,0.10545165,0.20442006,-0.10578702,0.0524051,-0.81768024,0.21157403,0.79983159,-0.41391142,5.55523651});
  par[14]=vector<double>({1.04967833,-0.66244026,-0.33017009,-0.49969348,0.31535097,0.15717562,0.24984674,-0.15767549,0.07828757,-0.99938697,0.31535097,0.79983159,-0.50476478,8.28827525});
  par[15]=vector<double>({1.16324504,-1.06444523,-0.53053535,-0.72682688,0.6650941,0.3314928,0.36341344,-0.33254705,0.16564051,-1.45365377,0.6650941,0.79983159,-0.73189818,17.50478474});
  par[16]=vector<double>({1.4812318,-2.53320209,-1.26258563,-1.36280041,2.33066077,1.16163602,0.6814002,-1.16533038,0.58234414,-2.72560082,2.33066077,0.79983159,-1.3678717,61.42874735});
  par[17]=vector<double>({1.9354986,-5.50863372,-2.74558505,-2.27133401,6.46445683,3.22198152,1.13566701,-3.23222841,1.61762929,-4.54266803,6.46445683,0.79983159,-2.27640531,170.49328232});
  par[18]=vector<double>({2.6168988,-11.90679784,-5.9345253,-3.63413442,16.53518424,8.24138198,1.81706721,-8.26759212,4.14113258,-7.26826884,16.53518424,0.79983159,-3.63920572,436.25776923});

  //par2[1] = vector<double>({ 0.799604 ,+0.800119 ,+0.000649899 ,-0.00141883 ,+0.00070959 });  
  //par2[2] = vector<double>({ 0.79915  ,+0.80069  ,+0.00194762  ,-0.004252   ,+0.00212756 });  
  //par2[3] = vector<double>({ 0.798471 ,+0.801542 ,+0.00388902  ,-0.00849051 ,+0.00425148 });  
  //par2[4] = vector<double>({ 0.797571 ,+0.802667 ,+0.00646788  ,-0.0141208  ,+0.00707771 });  
  //par2[5] = vector<double>({ 0.796677 ,+0.80378  ,+0.00903567  ,-0.0197272  ,+0.00989747 });  
  //par2[6] = vector<double>({ 0.794465 ,+0.806511 ,+0.0154068   ,-0.0336381  ,+0.0169186  });  
  //par2[7] = vector<double>({ 0.790998 ,+0.810723 ,+0.0254567   ,-0.0555835  ,+0.0280683  });  
  //par2[8] = vector<double>({ 0.78762  ,+0.814744 ,+0.0353297   ,-0.0771451  ,+0.0391145  });  
  //par2[9] = vector<double>({ 0.782718 ,+0.820415 ,+0.0498073   ,-0.108767   ,+0.0554899  });  
  //par2[10]= vector<double>({ 0.776492 ,+0.827305 ,+0.0684913   ,-0.149587   ,+0.0769615  });  
  //par2[11]= vector<double>({ 0.767817 ,+0.8362   ,+0.09519     ,-0.207937   ,+0.108393   });  
  //par2[12]= vector<double>({ 0.757489 ,+0.845372 ,+0.128311    ,-0.280361   ,+0.148853   });  
  //par2[13]= vector<double>({ 0.748577 ,+0.851473 ,+0.158599    ,-0.346642   ,+0.187657   });  
  //par2[14]= vector<double>({ 0.741081 ,+0.854502 ,+0.186057    ,-0.406781   ,+0.224805   });  
  //par2[15]= vector<double>({ 0.728534 ,+0.848637 ,+0.242311    ,-0.530251   ,+0.310433   });  
  //par2[16]= vector<double>({ 0.740481 ,+0.730088 ,+0.305669    ,-0.671714   ,+0.49514    });  
  //par2[17]= vector<double>({ 0.877894 ,+0.299655 ,+0.155487    ,-0.351649   ,+0.618277   });  
  //par2[18]= vector<double>({ 1.34949  ,-0.9219   ,-0.60073     ,+1.28026    ,+0.492548   });

 par2[1] = vector<double>({  0.780293,  +0.78092,  +0.000555456,  -0.00138479,  +0.000692453 });
 par2[2] = vector<double>({  0.779804,  +0.781682,  +0.00166567,  -0.00415106,  +0.00207606});
 par2[3] = vector<double>({  0.779071,  +0.78282,  +0.00332923,  -0.00829219,  +0.0041482});
 par2[4] = vector<double>({  0.778095,  +0.78433,  +0.00554404,  -0.0137983,  +0.00690497});
 par2[5] = vector<double>({  0.777121,  +0.785832,  +0.00775512,  -0.0192867,  +0.00965479});
 par2[6] = vector<double>({  0.774695,  +0.789546,  +0.0132664,  -0.0329305,  +0.0164989});
 par2[7] = vector<double>({  0.769419,  +0.793949,  +0.0219089,  -0.0544335,  +0.0273089});
 par2[8] = vector<double>({  0.765608,  +0.799653,  +0.0305693,  -0.0757162,  +0.0380386});
 par2[9] = vector<double>({  0.759946,  +0.807946,  +0.0434494,  -0.107114,  +0.0539256});
 par2[10]= vector<double>({  0.748503,  +0.814436,  +0.0595153,  -0.14726,  +0.0743259});
 par2[11]= vector<double>({  0.725159,  +0.816193,  +0.0807522,  -0.202962,  +0.10278});
 par2[12]= vector<double>({  0.688976,  +0.809742,  +0.104316,  -0.268545,  +0.136512});
 par2[13]= vector<double>({  0.644346,  +0.789531,  +0.120926,  -0.321425,  +0.163873});
 par2[14]= vector<double>({  0.605127,  +0.771238,  +0.135061,  -0.367143,  +0.187728});
 par2[15]= vector<double>({  0.498374,  +0.696251,  +0.150425,  -0.435036,  +0.223654});
 par2[16]= vector<double>({  0.255804,  +0.441098,  +0.124292,  -0.40583,  +0.211021});
 par2[17]= vector<double>({  0.0571808,  +0.125821,  +0.0356157,  -0.147653,  +0.0768695});
 par2[18]= vector<double>({  0.0130002,  +0.0383026,  +0.00914536,  -0.0531916,  +0.0274938});



  Aparameters["AEp"]=&AEp;
  Aparameters["AMp"]=&AMp;
  Aparameters["AEpd"]=&AEpd;
  Aparameters["AMpd"]=&AMpd;
  Aparameters["AMn25"]=&AMn25;
  Aparameters["AMn43"]=&AMn43;
  Aparameters["AEn25"]=&AEn25;
  Aparameters["AEn43"]=&AEn43;
  Aparameters["AFA25D"]=&AFA25D;

}

void SetPars(TFile* f, TFile* f_xs2D, MnvH1D* h_flux_external)
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();

  h_q2qe_base = (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen");
  h_q2qe_base->SetName( "h_q2qe_base");
  h_q2qe_base->Reset();
  h_q2qe_base->SetLineWidth(2);

  h_flux = (MnvH1D*) h_flux_external->Clone("h_flux");

  h_enu_q2qe_data = (MnvH2D*) f_xs2D->Get("h_enu_q2qe_data_nobck_hydrogen_unfold_effcor");
  h_enu_q2qe_mc = (MnvH2D*) f_xs2D->Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen_unfold_effcor");
  h_enu_q2qe_mc->ClearAllErrorBands();


  DeclareFittingPars();

};


vector<double> GetPar( double Q2 )
{
  int bin = h_q2qe_base->FindBin( Q2 );
  if (bin < 1 || bin>18 ) return vector<double>(14,0);
  return par[ bin ];
}
vector<double> GetPar2( double Q2 )
{
  int bin = h_q2qe_base->FindBin( Q2 );
  if (bin < 1 || bin>18 ) return vector<double>(5,0);
  return par2[ bin ];
}

double tau( double Q2, bool proton=true) 
{ 
  if(proton)return Q2/4/Mp/Mp; 
  else return Q2/4/Mn/Mn;
}
double Xi( double Q2, bool proton=true ) {return 2/(1+TMath::Sqrt(1+1/tau(Q2, proton)));}

double GKelly(double Q2, bool Ep = true)
{
  vector<double>*a = (Ep)? &aE : &aM;
  vector<double>*b = (Ep)? &bE : &bM;
  double t = tau(Q2,true);
  double num = 0;
  double den = 0;
  for( int i = 0; i<= 1; i++ ) num+=(a->at(i) )*pow(t,i);
  for( int i = 1; i<= 3; i++ ) den+=(b->at(i-1)) *pow(t,i);
  return num/(1+den);
}

double AN( string name, double xi )
{
  double A = 0;
  vector<double> *parA = Aparameters[name];
  int nPar = parA->size();

  for( int j = 1; j<=nPar; j++ )
  {
    double Pj = parA->at(j-1);
    for( int k = 1; k<=nPar; k++ ) 
      Pj*=(j==k)? 1:(xi - XiK[k-1])/(XiK[j-1] - XiK[k-1] ) ;
    A+=Pj;
  }
  return A;
}

double GD(double Q2, double MA, bool axial = true)
{
  double C = (axial)? gA:1;
  return C/pow(1+Q2/MA/MA,2);
}

double G(string name, double Q2)
{
  if(name == "Mp") 
  {
    double xi = Xi(Q2,true);
    return up*( AN("AMp",xi)*GKelly(Q2,false) );
  }
  if(name == "Ep") 
  {
    double xi = Xi(Q2,true);
    return AN("AEp", xi)*GKelly(Q2, true );
  }

  if(name == "Mn") 
  {
    double xi = Xi(Q2,false);
    double tn= tau(Q2,false);
    return un*AN("AMn25",xi)*( AN("AMp",Xi(Q2,false))*GKelly(Q2,false) );
  }

  if(name == "En") 
  {
    double tn= tau(Q2,false);
    double a = 1.7, b=3.3;
    double xi = Xi(Q2,false);
    return AN("AEn25", xi)*AN("AEp", xi)*GKelly(Q2, true )*(a*tn/(1+b*tn));
  }
  if(name== "FA")
  {
    double xi = Xi(Q2,false);
    double gd = GD(Q2,1.015, true );
    return AN("AFA25D",xi)*gd;

  }
  return -99999;
}
double GMp(double Q2) { return G("Mp",Q2); }
double GMn(double Q2) { return G("Mn",Q2); }
double GEp(double Q2) { return G("Ep",Q2); }
double GEn(double Q2) { return G("En",Q2); }

vector<double> FV12( double Q2, bool isProton = true )
{
  double t = tau(Q2,isProton);
  double GEV = GEp(Q2) - GEn(Q2);
  double GMV = GMp(Q2) - GMn(Q2);
  double FV1 = (GEV+t*GMV)/(1+t);
  double XIFV2 = (GMV-GEV)/(1+t);
  return vector<double>({FV1, XIFV2});
}

double FA(double Q2, double MA)
{
  return gA/pow( 1+Q2/MA/MA, 2 );
}


double FP(double fa, double Q2)
{
  double M = (Mp+Mn)/2;
  return fa*(2*M*M)/(Mpion*Mpion+Q2);
}





double GEp2(double Q2)
{
  double MV = TMath::Sqrt(0.71); 
  return AN("AEpd", Xi(Q2,true) )*GD(Q2,MV,false );
}
double GMp2(double Q2)
{
  double MV = TMath::Sqrt(0.71); 
  return AN("AMpd", Xi(Q2,true) )*GD(Q2,MV,false )*up;
}


#include "Xsec.h"


};

#endif
