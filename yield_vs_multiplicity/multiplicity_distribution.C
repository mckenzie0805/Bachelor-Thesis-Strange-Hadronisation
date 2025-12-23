void drawMultiplicity() {
    // Open your spectra ROOT file
    TFile *f = TFile::Open("/data/alice/mckenziebtr/spectra_monash_100M.root");
    
    // Access histogram
    TH1D *hMult = (TH1D*)f->Get("fHistMultiplicity");
    if (!hMult) {
        cout << "Error: could not find histogram 'fHistMultiplicity'!" << endl;
        return;
    }
    
    // Canvas setup
    TCanvas *c1 = new TCanvas("c1", "Multiplicity Distribution", 1000, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);
    
    // Style
    hMult->SetLineColor(kRed+1);
    hMult->SetLineWidth(3);
    hMult->SetTitle("Multiplicity Distribution"); // Remove title
    hMult->SetStats(0);  // Remove stats box
    hMult->GetXaxis()->SetTitle("Multiplicity N_{ch}");
    hMult->GetYaxis()->SetTitle("Counts");
    hMult->GetXaxis()->SetTitleSize(0.045);
    hMult->GetYaxis()->SetTitleSize(0.045);
    hMult->GetXaxis()->SetLabelSize(0.04);
    hMult->GetYaxis()->SetLabelSize(0.04);
    hMult->GetXaxis()->SetTitleOffset(1.1);
    hMult->GetYaxis()->SetTitleOffset(1.2);
    
    // Draw histogram first
    hMult->Draw("HIST");
    
    // ---- PERCENTILE SHADING ----
    Int_t nbins = hMult->GetNbinsX();
    Double_t total = hMult->Integral();
    Double_t cumulative = 0.0;
    
    // Define percentile boundaries (from high to low multiplicity)
    vector<double> percentiles = {0, 1, 5, 10, 20, 40, 60, 80, 100};
    vector<double> edges;
    
    // Calculate bin edges for percentiles
    for (int p = 0; p < percentiles.size(); ++p) {
        Double_t target = total * (percentiles[p] / 100.0);
        cumulative = 0.0;
        
        for (int i = nbins; i >= 1; --i) {  // Start from high multiplicity
            cumulative += hMult->GetBinContent(i);
            if (cumulative >= target) {
                edges.push_back(hMult->GetBinLowEdge(i));
                break;
            }
        }
    }
    
    // Color scheme - alternating: transparent, light red tint
    vector<int> fillStyles = {0, 3003, 0, 3003, 0, 3003, 0, 3003}; // 0 = no fill, 3003 = diagonal lines
    vector<int> fillColors = {0, kRed-10, 0, kRed-10, 0, kRed-10, 0, kRed-10};
    
    for (int i = 1; i < (int)edges.size(); ++i) {
	  if (fillStyles[i-1] == 0) continue;        // skip transparent ones
	
  // find max bin content in this x-interval
	  int b1 = hMult->FindFixBin(edges[i-1] + 1e-6);
 	 int b2 = hMult->FindFixBin(edges[i]   - 1e-6);
 	 double ymax = 0.0;
 	 for (int b = b1; b <= b2; ++b) {
  	  double v = hMult->GetBinContent(b);
 	   if (v > ymax) ymax = v;
  	  }
	  if (ymax <= 0) ymax = 0.5;                 // guard for log y

// Draw shaded bands between percentile regions
	for (int i = 1; i < (int)edges.size(); ++i) {
  		if (fillStyles[i-1] != 0) { // Only draw if not transparent
    		 TBox *box = new TBox(edges[i-1], 0.5,                 // Extend down to 0.5 (covers full y-range on log)
                         edges[i],   hMult->GetMaximum()*2); // Extend to top
    		 box->SetFillStyle(fillStyles[i-1]);
    		 box->SetFillColor(fillColors[i-1]);
    		 box->Draw("same");
  		}
	}	
    // Redraw histogram on top
    hMult->Draw("HIST same");
    
    // Add percentile text labels (vertical orientation, with cascade effect)
    TLatex latex;
    latex.SetTextSize(0.028);
    latex.SetTextAlign(22); // centered
    latex.SetTextColor(kBlack);
    latex.SetTextAngle(90); // Vertical text
    
    // Define y-positions for each label to follow the data trend (cascade DOWN from left to right)
    vector<double> yPositions = {2e1, 1e2, 4e2, 7e2, 3e3, 1e4, 1e4, 1e4}; // 0-1 lower, 1-5 lower, 5-10 lower
    
    for (int i = 1; i < edges.size(); ++i) {
        double xpos;
        
        // Special handling for specific bins
        if (i == edges.size() - 1) {
            // 80-100% bin: position at multiplicity ~2 (between 1 and 10)
            xpos = 2.0;
        } else if (i == 1) {
            // 0-1% bin: move MUCH further left (about 2cm on the plot)
            xpos = edges[i] * 1.3;  // Shift significantly left
        } else {
            xpos = sqrt(edges[i-1] * edges[i]); // geometric mean for log scale
        }
        
        double ypos = yPositions[i-1];
        
        TString label = Form("%.0f-%.0f%%", percentiles[i-1], percentiles[i]);
        latex.DrawLatex(xpos, ypos, label);
    }
    
    // Save the output
    c1->SaveAs("Multiplicity_Monash_Final4.png");
    cout << "✓ Saved: Multiplicity_Monash_Final4.png" << endl;
}
}
