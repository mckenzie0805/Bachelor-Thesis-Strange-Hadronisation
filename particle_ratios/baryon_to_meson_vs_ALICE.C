// -------------------------------------------------------------
// Plot_K0sXi_Ratio_Monash_vs_Junctions_vs_ALICE_LineALICE.C
//
// Draw X/K0s vs pT from both Monash and Junctions PYTHIA files
// and overlay ALICE data (HEPData-ins1471838-v1-root.root),
// merged into 20%-wide multiplicity classes: 0-20, 20-40, 40-60,
// 60-80, 80-100.
//
// ALICE handling (from HEPData ins1471838, which has 52 total tables):
//  - This macro uses Tables 1-10 (K0s) and Tables 21-30 (Xi spectra)
//  - For each multiplicity class, merge the two V0M tables by inverse-variance
//    weighting at the *ratio* level.
//  - X uses published bins
//  - K0s is averaged over the *overlapping* K0s fine bins inside each
//    X bin using length-weighted mean; total K0s error combined in
//    quadrature from stat (e1) and sys,total (e2), length-weighted.
//  - Ratio R = (X yield) / (K0s averaged yield), error from standard
//    propagation. No additional processing.
//
// PYTHIA handling:
//  - Same percentile finder on fHistMultiplicity.
//  - Project fHistPtK0s and fHistPtXis onto pT.
//  - Subsampling error: split the Y-bin range for a class into S equal
//    contiguous slices (default S=20, reduced automatically if too few
//    bins). Compute the X/K0s ratio for each slice, then use the mean as
//    the point and the standard error of the mean as the vertical error.
//    This keeps your look but is robust to sparse tails.
//  - Horizontal errors = bin half-width.
//
// ALICE data drawn as blue line with no markers or error bars.
// Output:
//  - One canvas per class, PNGs saved for both tunes.
// 
// Usage:
//   .L Plot_K0sXi_Ratio_Monash_vs_Junctions_vs_ALICE_LineALICE.C
//   Plot_K0sXi_Ratio_Monash_vs_Junctions_vs_ALICE_LineALICE(
//       "/path/to/pythia_monash.root",
//       "/path/to/pythia_junctions.root",
//       "/path/to/HEPData-ins1471838-v1-root.root",
//       2, 10.0, 20)
// -------------------------------------------------------------

#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TString.h"

// ---------- cosmetics ----------
namespace {
void SetPubStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetLabelSize(0.045,"XY");
  gStyle->SetTitleOffset(1.15,"Y");
  gStyle->SetTitleOffset(1.05,"X");
}

template<class T>
T* GetObj(TDirectory* d, const char* name) {
  return d ? dynamic_cast<T*>(d->Get(name)) : nullptr;
}

std::pair<double,double> CombineIVW(double r1,double e1,double r2,double e2) {
  if (e1<=0 && e2<=0) return { (r1+r2)/2.0, 0.0 };
  if (e1<=0) return { r2, e2 };
  if (e2<=0) return { r1, e1 };
  double w1=1.0/(e1*e1), w2=1.0/(e2*e2);
  double r=(w1*r1 + w2*r2)/(w1+w2);
  double e=std::sqrt(1.0/(w1+w2));
  return {r,e};
}

std::pair<int,int> YRangeForPercentile(TH1* hMult, double pTop, double pBot) {
  const int nb = hMult->GetNbinsX();
  const double total = hMult->Integral(1, nb);
  const double thrTop = total*(pTop/100.0);
  const double thrBot = total*(pBot/100.0);
  double cum = 0.0;
  int yH = nb, yL = nb;
  for (int b = nb; b >= 1; --b) {
    cum += hMult->GetBinContent(b);
    if (cum >= thrTop) { yH = b; break; }
  }
  for (int b = yH; b >= 1; --b) {
    if (cum >= thrBot) { yL = b; break; }
    cum += hMult->GetBinContent(b-1);
  }
  if (yL > yH) std::swap(yL,yH);
  return {yL, yH};
}

TH1D* RatioWithSubsampling(TH2* hK2D, TH2* hXi2D, const std::pair<int,int>& yRange,
                           const char* name, int nSplits = 20) {
  TH1D* hK = hK2D->ProjectionX(TString::Format("%s_K_full",name),
                               yRange.first, yRange.second, "e");
  TH1D* hXi = hXi2D->ProjectionX(TString::Format("%s_Xi_full",name),
                                 yRange.first, yRange.second, "e");

  TH1D* hR = new TH1D(name, "", hK->GetNbinsX(),
                      hK->GetXaxis()->GetXmin(), hK->GetXaxis()->GetXmax());
  if (!hR->GetSumw2N()) hR->Sumw2();
  hR->GetYaxis()->SetTitle("#Xi / K^{0}_{S}");

  int y1=yRange.first, y2=yRange.second;
  int nY = y2-y1+1;
  int S  = std::min(std::max(1,nY), nSplits);

  std::vector<std::pair<int,int>> slices; slices.reserve(S);
  int base=nY/S, rem=nY%S, start=y1;
  for (int s=0; s<S; ++s){
    int len=base+(s<rem?1:0);
    int end=start+len-1;
    slices.push_back({start,end});
    start=end+1;
  }

  for (int bx=1; bx<=hR->GetNbinsX(); ++bx){
    double K = hK->GetBinContent(bx);
    double Xi = hXi->GetBinContent(bx);
    double R = (K>0? Xi/K : 0.0);
    hR->SetBinContent(bx,R);

    std::vector<double> rs; rs.reserve(S);
    for (auto sl : slices){
      TH1D* kS = hK2D->ProjectionX("k_tmp", sl.first, sl.second, "e");
      TH1D* xiS = hXi2D->ProjectionX("xi_tmp", sl.first, sl.second, "e");
      double Ks=kS->GetBinContent(bx), XiS=xiS->GetBinContent(bx);
      rs.push_back(Ks>0? XiS/Ks : 0.0);
      delete kS; delete xiS;
    }
    double mean=0; for(double v:rs) mean+=v; mean/=rs.size();
    double var=0; for(double v:rs) var+=(v-mean)*(v-mean);
    // Divide variance by number of subsamples squared to normalize error properly
    double sem = (rs.size()>1? std::sqrt(var/(rs.size()*rs.size())) : 0.0);
    hR->SetBinError(bx, sem);
  }
  return hR;
}

double Quad(double a,double b){ return std::sqrt(a*a+b*b); }

void K0sAvgOverInterval(TH1* hK, TH1* hK_e1, TH1* hK_e2,
                        double xl, double xh,
                        double& y, double& ey) {
  double W = xh - xl;
  double sum=0.0, var=0.0;
  int nb=hK->GetNbinsX();
  for (int b=1;b<=nb;++b){
    double a=hK->GetXaxis()->GetBinLowEdge(b);
    double c=hK->GetXaxis()->GetBinUpEdge(b);
    double ov=std::max(0.0, std::min(xh,c)-std::max(xl,a));
    if (ov<=0) continue;
    double w=ov/W;
    double yb=hK->GetBinContent(b);
    double eb=Quad(hK_e1->GetBinContent(b), hK_e2->GetBinContent(b));
    sum += w*yb;
    var += w*w*eb*eb;
  }
  y=sum; ey=std::sqrt(var);
}

TGraphErrors* AliceRatioForTables(TDirectory* dK, TDirectory* dXi, const char* gname){
  auto hK    = GetObj<TH1>(dK,"Hist1D_y1");
  auto hK_e1 = GetObj<TH1>(dK,"Hist1D_y1_e1");
  auto hK_e2 = GetObj<TH1>(dK,"Hist1D_y1_e2");
  auto hXi    = GetObj<TH1>(dXi,"Hist1D_y1");
  auto hXi_e1 = GetObj<TH1>(dXi,"Hist1D_y1_e1");
  auto hXi_e2 = GetObj<TH1>(dXi,"Hist1D_y1_e2");
  if (!hK || !hK_e1 || !hK_e2 || !hXi || !hXi_e1 || !hXi_e2) return nullptr;

  std::vector<double> vx,vy,vex,vey;
  const int nbXi=hXi->GetNbinsX();
  for (int b=1;b<=nbXi;++b){
    double xl=hXi->GetXaxis()->GetBinLowEdge(b);
    double xh=hXi->GetXaxis()->GetBinUpEdge(b);
    double x=0.5*(xl+xh), ex=0.5*(xh-xl);

    double yXi=hXi->GetBinContent(b);
    double eXi=Quad(hXi_e1->GetBinContent(b), hXi_e2->GetBinContent(b));
    if (yXi<=0) continue;

    double yK=0, eK=0; K0sAvgOverInterval(hK,hK_e1,hK_e2,xl,xh,yK,eK);
    if (yK<=0) continue;

    double R=yXi/yK;
    double eR=R*std::sqrt((eXi/yXi)*(eXi/yXi) + (eK/yK)*(eK/yK));

    vx.push_back(x); vex.push_back(ex);
    vy.push_back(R); vey.push_back(eR);
  }
  auto* g = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
  g->SetName(gname);
  return g;
}

TGraphErrors* MergeAliceGraphsIVW(const TGraphErrors* g1, const TGraphErrors* g2, const char* name){
  if (!g1 && !g2) return nullptr;
  if (g1 && !g2) { auto* o=(TGraphErrors*)g1->Clone(name); o->SetName(name); return o; }
  if (!g1 && g2) { auto* o=(TGraphErrors*)g2->Clone(name); o->SetName(name); return o; }

  std::vector<double> vx,vy,vex,vey;
  const double tol=1e-9;
  for (int i=0;i<g1->GetN();++i){
    double x1,y1; g1->GetPoint(i,x1,y1);
    double e1=g1->GetErrorY(i);
    double ex=g1->GetErrorX(i);
    int jm=-1;
    for (int j=0;j<g2->GetN();++j){
      double x2,y2; g2->GetPoint(j,x2,y2);
      if (std::abs(x1-x2)<tol){ jm=j; break; }
    }
    if (jm<0) continue;
    double x2,y2; g2->GetPoint(jm,x2,y2);
    double e2=g2->GetErrorY(jm);
    auto comb=CombineIVW(y1,e1,y2,e2);
    vx.push_back(x1); vex.push_back(ex);
    vy.push_back(comb.first); vey.push_back(comb.second);
  }
  auto* g=new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
  g->SetName(name);
  return g;
}

TGraphErrors* ClipGraphX(const TGraphErrors* gin, double xmax, const char* name){
  if (!gin) return nullptr;
  std::vector<double> x,y,ex,ey;
  for (int i=0;i<gin->GetN();++i){
    double xi,yi; gin->GetPoint(i,xi,yi);
    if (xi<=xmax){
      x.push_back(xi);
      y.push_back(yi);
      ex.push_back(gin->GetErrorX(i));
      ey.push_back(gin->GetErrorY(i));
    }
  }
  auto* g=new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetName(name);
  return g;
}

double MaxYInRange(const TH1* h1, const TH1* h2, const TGraphErrors* g, double xmax){
  double m=0.0;
  if (h1){
    int bmax = h1->GetXaxis()->FindBin(xmax-1e-6);
    for (int b=1;b<=std::min(bmax,h1->GetNbinsX());++b){
      m = std::max(m, h1->GetBinContent(b)+h1->GetBinError(b));
    }
  }
  if (h2){
    int bmax = h2->GetXaxis()->FindBin(xmax-1e-6);
    for (int b=1;b<=std::min(bmax,h2->GetNbinsX());++b){
      m = std::max(m, h2->GetBinContent(b)+h2->GetBinError(b));
    }
  }
  if (g){
    for (int i=0;i<g->GetN();++i){
      double xi,yi; g->GetPoint(i,xi,yi);
      if (xi<=xmax) m = std::max(m, yi+g->GetErrorY(i));
    }
  }
  return m;
}
} // anon

// ------------------------------------------------------------------
// Main entry: Xi/K0s ratio with Monash and Junctions
// ------------------------------------------------------------------
void Plot_K0sXi_Ratio_Monash_vs_Junctions_vs_ALICE_LineALICE(
    const char* pythiaMonashFile,
    const char* pythiaJunctionsFile,
    const char* aliceFile,
    int rebinPt = 2,
    double ptMax = 10.0,
    int nSubSplits = 20)
{
  SetPubStyle();

  const double XMAX = 11.0;

  // ========== PYTHIA MONASH ==========
  TFile* fpMon = TFile::Open(pythiaMonashFile,"READ");
  if (!fpMon) { printf("Could not open Monash file: %s\n", pythiaMonashFile); return; }
  TH1* hMultMon = GetObj<TH1>(fpMon,"fHistMultiplicity");
  TH2* hK2DMon  = GetObj<TH2>(fpMon,"fHistPtK0s");
  TH2* hXi2DMon  = GetObj<TH2>(fpMon,"fHistPtXis");
  if (!hMultMon || !hK2DMon || !hXi2DMon) { 
    printf("Missing Monash histos.\n"); return; 
  }

  // ========== PYTHIA JUNCTIONS ==========
  TFile* fpJun = TFile::Open(pythiaJunctionsFile,"READ");
  if (!fpJun) { printf("Could not open Junctions file: %s\n", pythiaJunctionsFile); return; }
  TH1* hMultJun = GetObj<TH1>(fpJun,"fHistMultiplicity");
  TH2* hK2DJun  = GetObj<TH2>(fpJun,"fHistPtK0s");
  TH2* hXi2DJun  = GetObj<TH2>(fpJun,"fHistPtXis");
  if (!hMultJun || !hK2DJun || !hXi2DJun) { 
    printf("Missing Junctions histos.\n"); return; 
  }

  // ========== ALICE ==========
  TFile* fa = TFile::Open(aliceFile,"READ");
  if (!fa) { printf("Could not open ALICE file: %s\n", aliceFile); return; }

  struct ClassDef { int k1,k2,xi1,xi2; double pTop, pBot; const char* tag; const char* label; };
  std::vector<ClassDef> classes = {
    { 1,  2, 21, 22,  0, 20, "0-20",   "[0-20%]" },
    { 3,  4, 23, 24, 20, 40, "20-40",  "[20-40%]" },
    { 5,  6, 25, 26, 40, 60, "40-60",  "[40-60%]" },
    { 7,  8, 27, 28, 60, 80, "60-80",  "[60-80%]" },
    { 9, 10, 29, 30, 80,100, "80-100", "[80-100%]"}
  };

  // ========== Build ALICE graphs ==========
  std::vector<TGraphErrors*> gAlice(classes.size(), nullptr);
  for (size_t ic=0; ic<classes.size(); ++ic){
    const auto& C = classes[ic];
    TDirectory* dK1 = (TDirectory*)fa->Get(TString::Format("Table %d", C.k1));
    TDirectory* dK2 = (TDirectory*)fa->Get(TString::Format("Table %d", C.k2));
    TDirectory* dXi1 = (TDirectory*)fa->Get(TString::Format("Table %d", C.xi1));
    TDirectory* dXi2 = (TDirectory*)fa->Get(TString::Format("Table %d", C.xi2));

    TGraphErrors* gR1 = AliceRatioForTables(dK1,dXi1,TString::Format("gR_%s_a",C.tag));
    TGraphErrors* gR2 = AliceRatioForTables(dK2,dXi2,TString::Format("gR_%s_b",C.tag));
    TGraphErrors* gM  = MergeAliceGraphsIVW(gR1,gR2,TString::Format("gR_%s",C.tag));
    gAlice[ic] = ClipGraphX(gM, XMAX, TString::Format("gR_%s_clip",C.tag));

    if (gAlice[ic]) {
      gAlice[ic]->SetLineColor(kBlue+1);
      gAlice[ic]->SetLineWidth(2);
      gAlice[ic]->SetLineStyle(2);     // dashed line
      gAlice[ic]->SetMarkerStyle(0);   // no markers
      // Set zero error bars to hide them
      for (int i=0; i<gAlice[ic]->GetN(); ++i) {
        gAlice[ic]->SetPointError(i, 0, 0);
      }
    }
    delete gM;
  }

  // ========== Build PYTHIA Monash ratios ==========
  std::vector<TH1D*> hMonash(classes.size(), nullptr);
  for (size_t ic=0; ic<classes.size(); ++ic){
    const auto& C = classes[ic];
    auto yr = YRangeForPercentile(hMultMon, C.pTop, C.pBot);
    TH1D* hR = RatioWithSubsampling(hK2DMon,hXi2DMon,yr,TString::Format("hR_mon_%s",C.tag), nSubSplits);
    if (rebinPt>1) hR->Rebin(rebinPt);
    hR->GetXaxis()->SetRangeUser(0.0, XMAX);
    int bxMax = hR->GetXaxis()->FindBin(XMAX-1e-6);
    for (int b=bxMax+1; b<=hR->GetNbinsX(); ++b){ hR->SetBinContent(b,0); hR->SetBinError(b,0); }

    // Monash: green open circles, solid line
    hR->SetLineColor(kGreen+2);
    hR->SetMarkerColor(kGreen+2);
    hR->SetMarkerStyle(24);   // open circle
    hR->SetMarkerSize(1.1);
    hR->SetLineWidth(2);
    hR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hR->GetYaxis()->SetTitle("#Xi / K^{0}_{S}");
    hMonash[ic]=hR;
  }

  // ========== Build PYTHIA Junctions ratios ==========
  std::vector<TH1D*> hJunctions(classes.size(), nullptr);
  for (size_t ic=0; ic<classes.size(); ++ic){
    const auto& C = classes[ic];
    auto yr = YRangeForPercentile(hMultJun, C.pTop, C.pBot);
    TH1D* hR = RatioWithSubsampling(hK2DJun,hXi2DJun,yr,TString::Format("hR_jun_%s",C.tag), nSubSplits);
    if (rebinPt>1) hR->Rebin(rebinPt);
    hR->GetXaxis()->SetRangeUser(0.0, XMAX);
    int bxMax = hR->GetXaxis()->FindBin(XMAX-1e-6);
    for (int b=bxMax+1; b<=hR->GetNbinsX(); ++b){ hR->SetBinContent(b,0); hR->SetBinError(b,0); }

    // Junctions: red filled circles, solid line
    hR->SetLineColor(kRed+1);
    hR->SetMarkerColor(kRed+1);
    hR->SetMarkerStyle(20);   // filled circle
    hR->SetMarkerSize(1.1);
    hR->SetLineWidth(2);
    hR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hR->GetYaxis()->SetTitle("#Xi / K^{0}_{S}");
    hJunctions[ic]=hR;
  }

  // ========== Draw canvases ==========
  for (size_t ic=0; ic<classes.size(); ++ic){
    const auto& C = classes[ic];

    // Dynamic Y range
    double ymax = MaxYInRange(hMonash[ic], hJunctions[ic], gAlice[ic], XMAX);
    if (ymax<=0) ymax = 0.15;
    ymax *= 1.15;

    TCanvas* c = new TCanvas(TString::Format("c_ratio_%s",C.tag),
                             TString::Format("Xi/K0s vs pT  %s", C.label),
                             900,650);
    
    // Draw Junctions first (red)
    hJunctions[ic]->SetTitle(TString::Format("#Xi / K^{0}_{S} vs p_{T}   %s", C.label));
    hJunctions[ic]->SetMinimum(0.0);
    hJunctions[ic]->SetMaximum(ymax);
    hJunctions[ic]->Draw("E P");

    // Overlay Monash (green)
    hMonash[ic]->Draw("E P SAME");
    
    // Overlay ALICE (blue line)
    if (gAlice[ic]) gAlice[ic]->Draw("L SAME");

    // Legend
    TLegend* leg = new TLegend(0.62, 0.72, 0.95, 0.92);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->AddEntry(hJunctions[ic],"Junctions","l");
    leg->AddEntry(hMonash[ic],"Monash","l");
    if (gAlice[ic]) leg->AddEntry(gAlice[ic],"ALICE","l");
    leg->Draw();

    TString out = TString::Format("Ratio_XiOverK0s_Monash_vs_Junctions_%s.png", C.tag);
    c->SaveAs(out);
    printf("Saved: %s\n", out.Data());
  }

  // ========== QUANTITATIVE ANALYSIS OUTPUT ==========
  std::ofstream txtFile("Xi_K0s_Ratio_Analysis.txt");
  txtFile << "======================================================================\n";
  txtFile << "QUANTITATIVE ANALYSIS: Xi/K0s Ratio\n";
  txtFile << "Monash vs Junctions vs ALICE\n";
  txtFile << "======================================================================\n\n";
  txtFile << "Chi-squared formula used:\n";
  txtFile << "  chi^2 = SUM_i [(PYTHIA_i - ALICE_i)^2 / (sigma_PYTHIA_i^2 + sigma_ALICE_i^2)]\n";
  txtFile << "  where sigma = combined error from both PYTHIA statistics and ALICE uncertainties\n";
  txtFile << "  Lower chi^2/NDF indicates better agreement with experimental data.\n";
  txtFile << "  chi^2/NDF ~ 1 is ideal (within statistical uncertainties).\n\n";

  // Define pT ranges for averaging
  const double ptRanges[][2] = {{0.5, 2.0}, {2.0, 5.0}, {5.0, 8.0}};
  const char* ptLabels[] = {"Low pT (0.5-2.0 GeV/c)", "Mid pT (2.0-5.0 GeV/c)", "High pT (5.0-8.0 GeV/c)"};

  for (size_t ic=0; ic<classes.size(); ++ic){
    const auto& C = classes[ic];
    txtFile << "\n" << std::string(70, '=') << "\n";
    txtFile << "MULTIPLICITY CLASS: " << C.label << "\n";
    txtFile << std::string(70, '=') << "\n\n";

    // --- Mean ratios in pT ranges ---
    txtFile << "Mean Xi/K0s Ratios:\n";
    txtFile << std::string(70, '-') << "\n";
    txtFile << std::setw(25) << "pT Range"
            << std::setw(15) << "Monash"
            << std::setw(15) << "Junctions"
            << std::setw(15) << "ALICE\n";
    txtFile << std::string(70, '-') << "\n";

    for (int ipt=0; ipt<3; ++ipt){
      double ptLow = ptRanges[ipt][0];
      double ptHigh = ptRanges[ipt][1];

      // Calculate mean for Monash
      double sumMon=0, countMon=0;
      if (hMonash[ic]) {
        for (int b=1; b<=hMonash[ic]->GetNbinsX(); ++b){
          double pt = hMonash[ic]->GetBinCenter(b);
          if (pt>=ptLow && pt<=ptHigh){
            double val = hMonash[ic]->GetBinContent(b);
            if (val>0) { sumMon+=val; countMon++; }
          }
        }
      }
      double meanMon = (countMon>0 ? sumMon/countMon : 0.0);

      // Calculate mean for Junctions
      double sumJun=0, countJun=0;
      if (hJunctions[ic]) {
        for (int b=1; b<=hJunctions[ic]->GetNbinsX(); ++b){
          double pt = hJunctions[ic]->GetBinCenter(b);
          if (pt>=ptLow && pt<=ptHigh){
            double val = hJunctions[ic]->GetBinContent(b);
            if (val>0) { sumJun+=val; countJun++; }
          }
        }
      }
      double meanJun = (countJun>0 ? sumJun/countJun : 0.0);

      // Calculate mean for ALICE
      double sumAlice=0, countAlice=0;
      if (gAlice[ic]) {
        for (int i=0; i<gAlice[ic]->GetN(); ++i){
          double pt, val;
          gAlice[ic]->GetPoint(i, pt, val);
          if (pt>=ptLow && pt<=ptHigh && val>0){
            sumAlice+=val; countAlice++;
          }
        }
      }
      double meanAlice = (countAlice>0 ? sumAlice/countAlice : 0.0);

      txtFile << std::setw(25) << ptLabels[ipt]
              << std::setw(15) << std::fixed << std::setprecision(4) << meanMon
              << std::setw(15) << meanJun
              << std::setw(15) << (countAlice>0 ? meanAlice : 0.0) << "\n";
    }

    // --- Chi-squared analysis ---
    txtFile << "\nChi-squared Goodness-of-Fit (vs ALICE):\n";
    txtFile << std::string(70, '-') << "\n";

    // Need to restore ALICE errors for chi-squared calculation
    // Re-read the ALICE data with errors intact
    TDirectory* dK1 = (TDirectory*)fa->Get(TString::Format("Table %d", C.k1));
    TDirectory* dK2 = (TDirectory*)fa->Get(TString::Format("Table %d", C.k2));
    TDirectory* dXi1 = (TDirectory*)fa->Get(TString::Format("Table %d", C.xi1));
    TDirectory* dXi2 = (TDirectory*)fa->Get(TString::Format("Table %d", C.xi2));

    TGraphErrors* gR1_err = AliceRatioForTables(dK1,dXi1,TString::Format("gR_%s_a_err",C.tag));
    TGraphErrors* gR2_err = AliceRatioForTables(dK2,dXi2,TString::Format("gR_%s_b_err",C.tag));
    TGraphErrors* gAlice_err = MergeAliceGraphsIVW(gR1_err,gR2_err,TString::Format("gR_%s_err",C.tag));

    if (gAlice_err) {
      double chi2Mon=0, chi2Jun=0;
      int ndfMon=0, ndfJun=0;

      for (int i=0; i<gAlice_err->GetN(); ++i){
        double ptAlice, valAlice;
        gAlice_err->GetPoint(i, ptAlice, valAlice);
        double errAlice = gAlice_err->GetErrorY(i);
        if (valAlice<=0 || errAlice<=0) continue;

        // Find Monash bin
        if (hMonash[ic]) {
          int bin = hMonash[ic]->FindBin(ptAlice);
          double valMon = hMonash[ic]->GetBinContent(bin);
          double errMon = hMonash[ic]->GetBinError(bin);
          if (valMon>0){
            // Chi-squared using combined errors: sqrt(errPYTHIA^2 + errALICE^2)
            double totalErr = std::sqrt(errMon*errMon + errAlice*errAlice);
            chi2Mon += std::pow((valMon - valAlice)/totalErr, 2);
            ndfMon++;
          }
        }

        // Find Junctions bin
        if (hJunctions[ic]) {
          int bin = hJunctions[ic]->FindBin(ptAlice);
          double valJun = hJunctions[ic]->GetBinContent(bin);
          double errJun = hJunctions[ic]->GetBinError(bin);
          if (valJun>0){
            // Chi-squared using combined errors: sqrt(errPYTHIA^2 + errALICE^2)
            double totalErr = std::sqrt(errJun*errJun + errAlice*errAlice);
            chi2Jun += std::pow((valJun - valAlice)/totalErr, 2);
            ndfJun++;
          }
        }
      }

      txtFile << "Monash vs ALICE:    Chi2/NDF = " << std::fixed << std::setprecision(2)
              << (ndfMon>0 ? chi2Mon/ndfMon : 0.0) << " (" << chi2Mon << "/" << ndfMon << ")\n";
      txtFile << "Junctions vs ALICE: Chi2/NDF = " << (ndfJun>0 ? chi2Jun/ndfJun : 0.0)
              << " (" << chi2Jun << "/" << ndfJun << ")\n";

      if (ndfMon>0 && ndfJun>0) {
        if (chi2Mon/ndfMon < chi2Jun/ndfJun) {
          txtFile << ">>> Monash has BETTER agreement with ALICE\n";
        } else {
          txtFile << ">>> Junctions has BETTER agreement with ALICE\n";
        }
      }

      delete gAlice_err;
    } else {
      txtFile << "ALICE data not available for this class.\n";
    }

    // --- Percent difference: Junctions vs Monash ---
    txtFile << "\nJunctions Enhancement (relative to Monash):\n";
    txtFile << std::string(70, '-') << "\n";

    double sumDiff=0, countDiff=0;
    if (hMonash[ic] && hJunctions[ic]) {
      for (int b=1; b<=hMonash[ic]->GetNbinsX(); ++b){
        double valMon = hMonash[ic]->GetBinContent(b);
        double valJun = hJunctions[ic]->GetBinContent(b);
        if (valMon>0 && valJun>0){
          double pctDiff = 100.0 * (valJun - valMon) / valMon;
          sumDiff += pctDiff;
          countDiff++;
        }
      }
    }
    double avgPctDiff = (countDiff>0 ? sumDiff/countDiff : 0.0);

    txtFile << "Average percent difference: " << std::fixed << std::setprecision(2)
            << avgPctDiff << "%\n";
    if (avgPctDiff > 0) {
      txtFile << ">>> Junctions produces " << avgPctDiff << "% MORE Xi/K0s than Monash\n";
    } else {
      txtFile << ">>> Junctions produces " << std::abs(avgPctDiff) << "% LESS Xi/K0s than Monash\n";
    }
  }

  txtFile << "\n======================================================================\n";
  txtFile << "END OF ANALYSIS\n";
  txtFile << "======================================================================\n";
  txtFile.close();

  printf("Saved quantitative analysis: Xi_K0s_Ratio_Analysis.txt\n");
  printf("\nDone! Generated all Xi/K0s comparison plots.\n");
}
