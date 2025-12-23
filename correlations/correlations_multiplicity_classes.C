// ============================================================================
//  plot_correlations_MultClasses_3panel.C - 3-panel side-by-side comparison
// ============================================================================
// Usage: root -l -q 'plot_correlations_MultClasses_3panel.C("corr_Xi_mult_monash_100M.root", "corr_Xi_mult_junctions_100M.root")'

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <iostream>

#define PI 3.14159265

void plot_correlations_MultClasses_3panel(const char* monashFile, const char* junctionsFile)
{
    gStyle->SetOptStat(0);

    // Open files
    TFile* fMon = TFile::Open(monashFile);
    TFile* fJunc = TFile::Open(junctionsFile);

    if (!fMon || fMon->IsZombie()) {
        std::cerr << "Cannot open Monash file: " << monashFile << std::endl;
        return;
    }
    if (!fJunc || fJunc->IsZombie()) {
        std::cerr << "Cannot open Junctions file: " << junctionsFile << std::endl;
        return;
    }

    // Select 3 classes to display: 60-80%, 20-40%, 0-20% (left to right)
    // This puts lowest multiplicity on left, highest on right
    const int nPanels = 3;
    const char* classNames[3] = {"60_80", "20_40", "0_20"};
    const char* classLabels[3] = {"60-80%", "20-40%", "0-20%"};

    // Create canvas with 3 panels side-by-side
    TCanvas* c1 = new TCanvas("c_3panel", "Xi Correlations - Multiplicity Classes", 1800, 600);
    // Don't use Divide - we'll create custom pads with no gaps

    // Global maximum for consistent y-axis scaling
    double globalMax = 0;

    // First find global maximum across all panels
    for(int i = 0; i < nPanels; i++){
        TH1D* hDPhiMon = (TH1D*)fMon->Get(Form("hDPhi_%s", classNames[i]));
        TH1D* hTrPtMon = (TH1D*)fMon->Get(Form("hTrPt_%s", classNames[i]));
        TH1D* hDPhiJunc = (TH1D*)fJunc->Get(Form("hDPhi_%s", classNames[i]));
        TH1D* hTrPtJunc = (TH1D*)fJunc->Get(Form("hTrPt_%s", classNames[i]));

        if(!hDPhiMon || !hTrPtMon || !hDPhiJunc || !hTrPtJunc) continue;

        // Normalize
        double nTrigMon = hTrPtMon->Integral();
        double nTrigJunc = hTrPtJunc->Integral();
        double binWidth = hDPhiMon->GetBinWidth(1);

        TH1D* hMon = (TH1D*)hDPhiMon->Clone(Form("hMon_temp_%s", classNames[i]));
        TH1D* hJunc = (TH1D*)hDPhiJunc->Clone(Form("hJunc_temp_%s", classNames[i]));
        hMon->Scale(1.0 / (nTrigMon * binWidth));
        hJunc->Scale(1.0 / (nTrigJunc * binWidth));

        double maxMon = hMon->GetMaximum();
        double maxJunc = hJunc->GetMaximum();
        if(maxMon > globalMax) globalMax = maxMon;
        if(maxJunc > globalMax) globalMax = maxJunc;

        delete hMon;
        delete hJunc;
    }

    // Create custom pads with no gaps
    // Each panel is 1/3 width
    double padWidth = 1.0 / 3.0;
    TPad* pads[3];

    for(int i = 0; i < nPanels; i++){
        double xlow = i * padWidth;
        double xup = (i + 1) * padWidth;
        pads[i] = new TPad(Form("pad%d", i), "", xlow, 0.0, xup, 1.0);

        // Set margins: only leftmost panel has left margin for y-axis
        if(i == 0){
            pads[i]->SetLeftMargin(0.18);   // More space for y-axis to prevent cutoff
            pads[i]->SetRightMargin(0.0);   // No gap between panels
        } else if(i == nPanels - 1){
            pads[i]->SetLeftMargin(0.0);    // No gap
            pads[i]->SetRightMargin(0.05);  // Small margin on far right
        } else {
            pads[i]->SetLeftMargin(0.0);    // No gap
            pads[i]->SetRightMargin(0.0);   // No gap
        }

        pads[i]->SetBottomMargin(0.15);
        pads[i]->SetTopMargin(0.08);
        pads[i]->Draw();
    }

    // Now draw each panel
    for(int i = 0; i < nPanels; i++){
        pads[i]->cd();
        gPad->SetGrid();
        gStyle->SetGridColor(kGray);
        gStyle->SetGridStyle(3);  // Dotted lines
        gStyle->SetGridWidth(1);

        TH1D* hDPhiMon = (TH1D*)fMon->Get(Form("hDPhi_%s", classNames[i]));
        TH1D* hTrPtMon = (TH1D*)fMon->Get(Form("hTrPt_%s", classNames[i]));
        TH1D* hDPhiJunc = (TH1D*)fJunc->Get(Form("hDPhi_%s", classNames[i]));
        TH1D* hTrPtJunc = (TH1D*)fJunc->Get(Form("hTrPt_%s", classNames[i]));

        if(!hDPhiMon || !hTrPtMon || !hDPhiJunc || !hTrPtJunc) continue;

        // Normalize
        double nTrigMon = hTrPtMon->Integral();
        double nTrigJunc = hTrPtJunc->Integral();
        double binWidth = hDPhiMon->GetBinWidth(1);

        TH1D* hMon = (TH1D*)hDPhiMon->Clone(Form("hMon_%s", classNames[i]));
        TH1D* hJunc = (TH1D*)hDPhiJunc->Clone(Form("hJunc_%s", classNames[i]));
        hMon->Scale(1.0 / (nTrigMon * binWidth));
        hJunc->Scale(1.0 / (nTrigJunc * binWidth));

        // Styling
        hMon->SetLineColor(kBlue+1);
        hMon->SetLineWidth(2);
        hJunc->SetLineColor(kRed+1);
        hJunc->SetLineWidth(2);

        // Set consistent y-axis range
        hMon->SetMaximum(globalMax * 1.15);
        hMon->SetMinimum(0);

        // Axis titles and labels
        hMon->SetTitle("");
        hMon->GetXaxis()->SetTitle("#Delta#phi (rad)");
        hMon->GetXaxis()->SetTitleSize(0.045);
        hMon->GetXaxis()->SetLabelSize(0.035);

        // Y-axis: only show labels and title on leftmost panel
        if(i == 0){
            hMon->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN}{d#Delta#phi}");
            hMon->GetYaxis()->SetTitleSize(0.045);  // Reduced from 0.055
            hMon->GetYaxis()->SetLabelSize(0.035);  // Reduced from 0.045
            hMon->GetYaxis()->SetTitleOffset(1.4);
        } else {
            hMon->GetYaxis()->SetTitle("");
            hMon->GetYaxis()->SetLabelSize(0);  // Hide y-axis labels on other panels
        }

        // Draw
        hMon->Draw("HIST E");
        hJunc->Draw("HIST E SAME");

        // Add panel label (smaller size)
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextFont(42);

        // Position labels at same distance from left edge of each panel
        double xPos;
        if(i == 0) {
            xPos = 0.35;  // First panel has left margin, so position further right
        } else {
            xPos = 0.15;  // Other panels have no left margin
        }
        latex.DrawLatex(xPos, 0.88, Form("(%c) %s", 'a'+i, classLabels[i]));

        // Legend (only on rightmost panel - 0-20%)
        if(i == nPanels - 1){
            TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.82);
            leg->SetBorderSize(1);       // Add box around legend
            leg->SetFillStyle(1001);     // Solid fill
            leg->SetFillColor(kWhite);   // White background
            leg->SetTextSize(0.038);     // Reduced from 0.045
            leg->AddEntry(hMon, "Monash", "l");
            leg->AddEntry(hJunc, "Junctions", "l");
            leg->Draw();
        }
    }

    c1->SaveAs("correlations_Xi_MultClasses_3panel.png");
    std::cout << "Created: correlations_Xi_MultClasses_3panel.png" << std::endl;

    fMon->Close();
    fJunc->Close();
    delete c1;
}
