


double A(double Q2, double M, double fa, double fv1, double xifv2 )
{
  double t = Q2/4/M/M;
  double f = Q2/M/M;
  double m = Mmu;//muon mass (GeV)
  double fp =FP(fa,Q2);

  double prefactor = (m*m+Q2)/(4*M*M);
  double p1 = (4+f)*TMath::Abs(fa*fa);
  double p2 = (4-f)*TMath::Abs(fv1*fv1);
  double p3 = f*TMath::Abs(xifv2*xifv2)*(1-t);
  double p4 = 4*f*(fv1*xifv2);
 
  double P1 = p1-p2+p3+p4;

  double pref2 = m*m/M/M;
  double pp1 = pow( fv1+xifv2, 2 );
  double pp2 = pow( fa+2*fp, 2 );
  double pp3 = (-f-4)*fp*fp;
  double P2 = pref2*( pp1+pp2+pp3 );

  return prefactor*( P1-P2 );
}


double B(double Q2, double M, double fa, double fv1, double xifv2 )
{
  return Q2/M/M*fa*(fv1+xifv2);
}

double C(double Q2, double M, double fa, double fv1, double xifv2 )
{
  return 0.25*( fa*fa + TMath::Abs(fv1*fv1) + Q2/M/M*TMath::Abs( pow( (xifv2/2), 2 ) ) );

}

double XsPrefactor( double Ev )
{
  double M = (Mp+Mn)/2;
  double fact = 1.973269788*pow(10,-14);
  return fact*fact*TMath::Power( M*Gf*cosThetaC, 2 )/(8*TMath::Pi() *Ev*Ev );
}

double XsMult( double Ev, double Q2, double fa, bool isNeutino )
{
  bool isProton = !isNeutino;
  //double fa=FA( Q2, MA);
  vector<double> fv12 = FV12( Q2, isProton );
  double fv1 = fv12[0];
  double xifv2=fv12[1];
  double M = isProton? Mp : Mn;
  double m = Mmu;
  double a = A(Q2,M, fa, fv1,xifv2 );
  double b = B(Q2,M, fa, fv1,xifv2 );
  double c = C(Q2,M, fa, fv1,xifv2 );

  double suM = (4*M*Ev - Q2 - m*m)/(M*M);
  //double s = M*M+2*M*Ev;
  //double u = M*M+m*m - 2M*Ev
  int sign = (isNeutino)? -1:1;
  return a+sign*b*suM+c*suM*suM;
}


double dAdFA(double Q2,double M,double FA, double f1v, double xsif2v)
{
  double m = Mmu;
  double fP=FP(1,Q2);
  double m2=m*m;
  double M2=M*M;
  //(FA (m^2 + Q2) (M^2 (4 M^2 + Q2) + 
  //   m^2 (-(1 + 4 fP) M^2 + fP^2 Q2)))/(2 M^6)

  return (FA*(m2 + Q2)*(M2*(4*M2 + Q2) + 
     m2*(-(1 + 4*fP)*M2 + fP*fP*Q2)))/(2*pow(M2,3));
}

double dBdFA(double Q2, double M, double FA, double f1v, double xsif2v)
{
  return (Q2*(f1v + xsif2v))/(M*M);
}

double dCdFA(double Q2, double M, double FA, double f1v, double xsif2v)
{
  return FA/2;
}



double dXsMult( double Ev, double Q2, double fa, bool isNeutino )
{
  bool isProton = !isNeutino;
  //double fa=FA( Q2, MA);
  vector<double> fv12 = FV12( Q2, isProton );
  double fv1 = fv12[0];
  double xifv2=fv12[1];
  double M = isProton? Mp : Mn;
  double m = Mmu;
  double a = dAdFA(Q2,M, fa, fv1,xifv2 );
  double b = dBdFA(Q2,M, fa, fv1,xifv2 );
  double c = dCdFA(Q2,M, fa, fv1,xifv2 );

  double suM = (4*M*Ev - Q2 - m*m)/(M*M);
  //double s = M*M+2*M*Ev;
  //double u = M*M+m*m - 2M*Ev
  int sign = (isNeutino)? -1:1;
  return a+sign*b*suM+c*suM*suM;
}

double XsLS( double Ev, double Q2, double fa, bool isNeutino )
{
  double s = (doScale)? scale:1;
  return XsPrefactor(Ev) * XsMult( Ev, Q2, fa, isNeutino )*s;
}

double dXsLS( double Ev, double Q2, double fa, bool isNeutino )
{
  double s = (doScale)? scale:1;
  return XsPrefactor(Ev) * dXsMult( Ev, Q2, fa, isNeutino )*s;
}

//double XsKevin( double Ev, double Q2, double MA,bool isNeutino )
//{
//  bool isProton = !isNeutino;
//
//  vector<double> coeff = GetPar( Q2 );
//
//  double M = isProton? Mp : Mn;
//  double m = Mmu;
//
//  double GVE = GEp(Q2) - GEn(Q2);
//  double GVM = GMp(Q2) - GMn(Q2);
//  double Fa=FA( Q2, MA);
//
//  double a= M/Ev;
//  double a2= a*a;
//  double Fa2 = Fa*Fa;
//  double GVE2 = GVE*GVE;
//  double GVM2 = GVM*GVM;
//  double FaGVM = Fa*GVM;
//  double GVEGVM = GVE*GVM;
//
//  vector<double> var({ Fa2, Fa2*a, Fa2*a2, GVE2, GVE2*a, GVE2*a2, FaGVM*a, FaGVM*a2, GVEGVM, GVEGVM*a, GVEGVM*a2, GVM2, GVM2*a, GVM2*a2 } );
//
//  double xs = 0;
//  for( uint i = 0; i<coeff.size(); ++i ) xs+=var[i]*coeff[i];
//  return xs*pow(10,-38);
//}

double XsKevin(double Ev, double Q2, double Fa,bool isNeutrino )
{

  vector<double> coeff = GetPar2( Q2 );

  double M = isNeutrino? Mn : Mp;
  double m = Mmu;

  double GVE = GEp(Q2) - GEn(Q2);
  double GVM = GMp(Q2) - GMn(Q2);

  double Fa2 = Fa*Fa;
  double GVE2 = GVE*GVE;
  double GVM2 = GVM*GVM;
  double FaGVM = Fa*GVM;
  double GVEGVM = GVE*GVM;

  vector<double> var({ Fa2,  GVE2,FaGVM, GVEGVM,  GVM2} );

  double xs = 0;
  for( uint i = 0; i<coeff.size(); ++i ) xs+=var[i]*coeff[i];
  return xs*pow(10,-38);
}


double dXsKevin(double Ev, double Q2, double Fa, bool isNeutrino)
{
  bool isProton = !isNeutrino;

  vector<double> coeff = GetPar2( Q2 );

  double M = isProton? Mp : Mn;
  double m = Mmu;

  double GVE = GEp(Q2) - GEn(Q2);
  double GVM = GMp(Q2) - GMn(Q2);

  double a= M/Ev;
  double a2= a*a;
  double Fa2 = Fa*Fa;
  double GVE2 = GVE*GVE;
  double GVM2 = GVM*GVM;
  double FaGVM = Fa*GVM;
  double GVEGVM = GVE*GVM;

  //vector<double> var({ Fa2, Fa2*a, Fa2*a2, GVE2, GVE2*a, GVE2*a2, FaGVM*a, FaGVM*a2, GVEGVM, GVEGVM*a, GVEGVM*a2, GVM2, GVM2*a, GVM2*a2 } );
  vector<double> dvar({ Fa,  0,GVM,0,  0} );

  double dxs = 0;
  for( uint i = 0; i<coeff.size(); ++i ) dxs+=dvar[i]*coeff[i];
  return dxs*pow(10,-38);
}

double XsIntegral( double Q2, double fa, bool isNeutino, bool useKevin = true )
{

  if(useKevin)
  {
    return XsKevin(0, Q2, fa, isNeutino );
  }
  else
  {
    double xs=0;
    int bin = h_q2qe_base->FindBin(Q2);
    h_enu = h_flux;
    double integral = h_enu->Integral("width");
    double minE = EnuMin[bin-1];
    double maxE = EnuMax[bin-1];
    if( integral == 0 ) return xs;
    for( int i = 1; i<= h_enu->GetNbinsX(); i++ )
    {
      double E = h_enu->GetBinCenter( i );
      if (E<minE || E>maxE) continue;
      double w = h_enu->GetBinContent( i ) * h_enu->GetBinWidth(i)/ integral;
      //cout<<integral<<", "<<w<<", "<<E<<endl;
      //double xsE = (useKevin)? XsKevin(E, Q2, fa, isNeutino ):XsLS(E, Q2,fa,isNeutino);
      double xsE = XsLS(E, Q2,fa,isNeutino);
      xs+= w*xsE;
    }
    return xs;
  }
}

double dXsIntegral( double Q2, double fa, bool isNeutino, bool useKevin = true )
{

  if(useKevin)
  {
    return dXsKevin(0, Q2, fa, isNeutino );
  }
  else
  {
    double dxs=0;
    int bin = h_q2qe_base->FindBin(Q2);
    h_enu = h_flux;
    double integral = h_enu->Integral("width");
    double minE = EnuMin[bin-1];
    double maxE = EnuMax[bin-1];
    if( integral == 0 ) return dxs;
    for( int i = 1; i<= h_enu->GetNbinsX(); i++ )
    {
      double E = h_enu->GetBinCenter( i );
      if (E<minE || E>maxE) continue;
      double w = h_enu->GetBinContent( i ) * h_enu->GetBinWidth(i)/ integral;
      //cout<<integral<<", "<<w<<", "<<E<<endl;
      //double xsE = (useKevin)? XsKevin(E, Q2, fa, isNeutino ):XsLS(E, Q2,fa,isNeutino);
      double dxsE = dXsLS(E, Q2,fa,isNeutino);
      dxs+= w*dxsE;
    }
    return dxs;
  }
}


