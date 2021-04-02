#include "fitZexp.h"

double runFitter(double v)
{
  TFitter fitter(4);
  lambda = v;
  fitter.SetFCN( minuitFunction );
  double v=1,verr=1e-8,vmin=-200, vmax=200;
  //for( int i = 0; i< 4; i++ ) fitter.SetParameter(i, Form("%d", i), v, verr, vmin, vmax );
  fitter.SetParameter(0, Form("%d", 0), 2.56521, verr, vmin, vmax );
  fitter.SetParameter(1, Form("%d", 1), 1.14, verr, vmin, vmax );
  fitter.SetParameter(2, Form("%d", 2), -17.7, verr, vmin, vmax );
  fitter.SetParameter(3, Form("%d", 3), 6.4, verr, vmin, vmax );

  //fitter.ExecuteCommand("SET PRINTOUT",0,1);
  fitter.ExecuteCommand("MIGRAD",0,0);

  double par[4]={ fitter.GetParameter(0), fitter.GetParameter(1), fitter.GetParameter(2), fitter.GetParameter(3) };

  cout<<par<<endl;
  return chiSqFunc(0, par);
}

TFitter *fitter2;
double parameters[4];
double runFitter2(double v)
{
  fitter2 = new TFitter(4);
  lambda = v;
  fitter2->SetFCN( minuitFunction );
  double v=1,verr=1e-8,vmin=-200, vmax=200;
  fitter2->SetParameter(0, Form("%d", 0), 2.56521, verr, vmin, vmax );
  fitter2->SetParameter(1, Form("%d", 1), 1.14, verr, vmin, vmax );
  fitter2->SetParameter(2, Form("%d", 2), -17.7, verr, vmin, vmax );
  fitter2->SetParameter(3, Form("%d", 3), 6.4, verr, vmin, vmax );

  //fitter2.ExecuteCommand("SET PRINTOUT",0,1);
  fitter2->ExecuteCommand("SIMPLEX",0,0);
  fitter2->ExecuteCommand("MIGRAD",0,0);
  parameters[0]=fitter2->GetParameter(0);
  parameters[1]=fitter2->GetParameter(1);
  parameters[2]=fitter2->GetParameter(2);
  parameters[3]=fitter2->GetParameter(3);



  return chiSqFunc(lambda,parameters);
}

TGraph *g = new TGraph();
double ProcessResults()
{
  double minX = 99999;
  double minY=999999;
  vector<double> x;
  for( int i=-10; i<10;i++ ) x.push_back( TMath::Power(10,i*1.0) );
  vector<double> y;

  for( int i = 0; i< x.size(); i++ ) 
  {
    double l = x[i];
    double yl = NumericalMinimization("Minuit2","MIGRAD", l);
    y.push_back(yl);
    g->SetPoint(i, l, yl );
    if( yl < minY )
    {
      minY = yl;
      minX = l;
    }
  }

  g->GetXaxis()->SetRangeUser(x[0], 1000);
  g->GetXaxis()->SetTitle("#lambda");
  g->GetYaxis()->SetTitle("#chi^{2}");
  g->SetTitle("L Curve");

  NumericalMinimization("Minuit2","SteepestDescent", minX);


  return minX;

}
