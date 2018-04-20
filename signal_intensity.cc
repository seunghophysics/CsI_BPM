// still in progress.


#include "pco2root2x2.h"


using namespace std;


// Sets error for each bin in input 2D histogram
void set_error(TH2S *hist, double exp_time_min)
{
  double dc_err = 0.01 * exp_time_min * 60;
  double ro_err = 11.0;
  double sig_err;

  for(int i = 1; i <= 1600; i++){
    for(int j = 1; j <= 1200; j++){
      sig_err = sqrt(hist->GetBinContent(i, j));
      hist->SetBinError(i, j, sqrt(dc_err * dc_err + ro_err * ro_err + sig_err * sig_err));
    }
  }  
 }

// Takes the "t1" histogram out of a root file.
TH2S *get_hist(const char *filename, const char *set_histname)
{
  TFile *f = TFile::Open(filename);
  TH2S *hist = (TH2S*)f->Get("t1");
  hist->SetName(set_histname);

  return hist;
}

// Gets the threshold ADC value for determining whether a pixel is hot or not.
// Draws the amplitude spectrum of the dark (or data) picture if draw_opt is TRUE.
double get_threshold(TH2S *dark, double significance, bool draw_opt)
{
  int max_content = (int)dark->GetMaximum();
   
  if(gROOT->FindObject("amp_spectrum") != NULL)
    gROOT->FindObject("amp_spectrum")->Delete();
  
  TH1D *amp_spectrum = new TH1D("amp_spectrum", "Intensity Spectrum; Intensity (ADC counts); Frequency;", max_content, 0, max_content);
  
  for(int i = 1; i <= dark->GetNbinsX(); i++)
    for(int j = 1; j <= dark->GetNbinsY(); j++)
      amp_spectrum->Fill(dark->GetBinContent(i, j));
  
  TF1 *gf = new TF1("gf", "gaus", 0, max_content);
  amp_spectrum->Fit(gf, "Q0");
  
  double mean = gf->GetParameter(1);
  double sigma = gf->GetParameter(2);
  
  if(draw_opt) 
    amp_spectrum->Draw();

  return mean + significance * sigma;
}


TH1D *amp_spectrum(TH2S *dark, const char *spec_name)
{
  int max_content = (int)dark->GetMaximum();
 
  TH1D *amp_spectrum = new TH1D(spec_name, "Intensity Spectrum; Intensity (ADC counts); Frequency;", max_content, 0, max_content);
  
  for(int i = 1; i <= dark->GetNbinsX(); i++)
    for(int j = 1; j <= dark->GetNbinsY(); j++)
      amp_spectrum->Fill(dark->GetBinContent(i, j));

	return amp_spectrum;	
}

// Returns a 2D map of hot pixels, given a dark (or data) picture.
TH2S *hotpixel_map(TH2S *dark, double significance)
{ 
  double threshold = get_threshold(dark, significance, 0);
	cout << "Threshold : " << threshold << endl;

  /*
  if(gROOT->FindObject("hotpixel_map") != NULL)
    gROOT->FindObject("hotpixel_map")->Delete();
  */

	double xmin = dark->GetXaxis()->GetXmin();
	double xmax = dark->GetXaxis()->GetXmax();
	double ymin = dark->GetYaxis()->GetXmin();
	double ymax = dark->GetYaxis()->GetXmax();
	
  TH2S *map = new TH2S(Form("%s_%s", dark->GetName(), "hp"), "Hot Pixel Coordinates; x coordinate; y coordinate", dark->GetNbinsX(), xmin, xmax, dark->GetNbinsY(), ymin, ymax);
  
  for(int i = 1; i <= dark->GetNbinsX(); i++){
    for(int j = 1; j <= dark->GetNbinsY(); j++){
      if(dark->GetBinContent(i, j) < threshold)
        map->SetBinContent(i, j, 0);
      else
        map->SetBinContent(i, j, 1);
    }
  }      
  
  return map;
}

// Gives an average of neighboring non-hot pixels.
double get_average(vector<int> vec, vector<int> map)
{
  double sum = 0.0;
  double denominator = 0.0;

  for(int i = 0; i <= 8; i++){
    if(!map[i]){
      sum += vec[i];
      denominator++;
    }
  }
  
  denominator += 1e-5;

  return sum / denominator;
}

// Returns a moderated data, given a hot pixel map.
TH2S *moderate(TH2S *data, TH2S *hotpixel_map)
{
  if(gROOT->FindObject("moderated_data") != NULL)
    gROOT->FindObject("moderated_data")->Delete();
  
  TH2S *moderated_data = (TH2S*)data->Clone("moderated_data");

  vector<int> neighbors;
  vector<int> neighbor_map;
  vector<vector<double>> sp_hot;

  double avg;

  for(int i = 2; i <= 1599; i++){
    for(int j = 2; j <= 1199; j++){

      neighbors.clear();
      neighbor_map.clear();

      if(hotpixel_map->GetBinContent(i, j)){

        for(int k = -1; k <= 1; k++){
          for(int l = -1; l <= 1; l++){
            neighbors.push_back(data->GetBinContent(i + k, j + l));
            neighbor_map.push_back(hotpixel_map->GetBinContent(i + k, j + l));
          }
        }      
    avg = get_average(neighbors, neighbor_map);  
    
    if(avg == 0.0)
      sp_hot.push_back({(double)i, (double)j, data->GetBinContent(i, j)});
    else
      moderated_data->SetBinContent(i, j, avg);
      }
    }
  }
  
  int it = 0;
  double sp_sum = 0.0;
  for(vector<vector<double>>::iterator iter = sp_hot.begin(); iter != sp_hot.end(); ++iter){
    sp_sum += sp_hot[it][2];
    it++;
  }

  it = 0; 
  for(vector<vector<double>>::iterator iter = sp_hot.begin(); iter != sp_hot.end(); ++iter){
    moderated_data->SetBinContent(sp_hot[it][0], sp_hot[it][1], (moderated_data->Integral() - sp_sum) / (1600 * 1200));
    it++;
  }
      
  return moderated_data;
}

// Returns a cut 2D histogram (a region of interest)
TH2S *cut_hist(TH2S *org_histo, int start_binx, int end_binx, int start_biny, int end_biny)
{  
  TH2S *cut_histo = new TH2S(org_histo->GetName(), "Cut Histogram; x coordinate; y coordinate;", end_binx - start_binx + 1, start_binx - 1, end_binx, end_biny - start_biny + 1, start_biny - 1, end_biny);

  int k, l;
  k = 1;
  for(int i = start_binx; i <= end_binx; i++){
    l = 1;
    for(int j = start_biny; j <= end_biny; j++){
      cut_histo->SetBinContent(k, l, org_histo->GetBinContent(i, j));
      l++;
    }
    k++;
  }

  return cut_histo;
}


double *int_err(TH1D *prj_x)
{
  double start_x = prj_x->GetBinCenter(1) - 0.5;
  double end_x = prj_x->GetBinCenter(prj_x->GetNbinsX()) - 0.5;
 
  TF1 *sig_back_f = new TF1("sig_back_f", "gaus(0) + pol2(3)", start_x - 1, end_x);
  TF1 *sig_f = new TF1("sig_f", "gaus", start_x - 1, end_x);
  TF1 *back_f = new TF1("back_f", "pol2", start_x - 1, end_x);
  
  sig_back_f->SetParameters(prj_x->GetMaximum() - prj_x->GetMinimum(), start_x + prj_x->GetMaximumBin(), 50, prj_x->GetMinimum(), 0, -0.05);

  sig_back_f->SetParName(0, "Signal Amplitude");
  sig_back_f->SetParName(1, "Signal Max Position");
  sig_back_f->SetParName(2, "Signal Width (sigma)");
  sig_back_f->SetParName(3, "Background Offset");
  sig_back_f->SetParName(4, "Background 1st-or. Coeff.");
  sig_back_f->SetParName(5, "Background 2st-or. Coeff.");

  sig_back_f->SetParLimits(0, 0, 5e4);
  sig_back_f->SetParLimits(1, start_x + prj_x->GetMaximumBin() - 50, start_x + prj_x->GetMaximumBin() + 50);
  sig_back_f->SetParLimits(2, 10, 100);
  sig_back_f->SetParLimits(3, 0, prj_x->GetMinimum());
  sig_back_f->SetParLimits(4, 0, 20);
  sig_back_f->SetParLimits(5, -5e-2, 0);
 
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TFitResultPtr fit_result = prj_x->Fit("sig_back_f", "Sbe");

  sig_f->SetParameters(sig_back_f->GetParameter(0), sig_back_f->GetParameter(1), sig_back_f->GetParameter(2));
  sig_f->SetParError(0, sig_back_f->GetParError(0));
  sig_f->SetParError(1, sig_back_f->GetParError(1));
  sig_f->SetParError(2, sig_back_f->GetParError(2));

  back_f->SetParameters(sig_back_f->GetParameter(3), sig_back_f->GetParameter(4), sig_back_f->GetParameter(5));
  back_f->Draw("same");

  double *int_and_err = new double[2];
  int_and_err[0] = sig_f->Integral(start_x - 1, end_x);
  int_and_err[1] = sig_f->IntegralError(start_x - 1, end_x, 0, fit_result->GetCovarianceMatrix().GetSub(0, 2, 0, 2).GetMatrixArray()) + back_f->IntegralError(sig_f->GetParameter(1) - 3 * sig_f->GetParameter(2), sig_f->GetParameter(1) + sig_f->GetParameter(2), 0, fit_result->GetCovarianceMatrix().GetSub(3, 5, 3, 5).GetMatrixArray());

  return int_and_err;
}


// Calculates the signal intensity of a region of interest, given dark, data, and significance for hot pixels.
// Example usage: int_err = intensity_and_error(get_hist("folder/file_data.root", "data"), get_hist("folder/file_dark.root", "dark"), 2, 10, 400, 1200, 200, 600);
// In this case, int_err[0] stores the total intensity, and int_err[1] stores the error of the total intensity.
double *intensity_and_error(TH2S *data, TH2S *dark, double significance, double exp_time_min, int start_x, int end_x, int start_y, int end_y)
{
  if(gROOT->FindObject("back_f") != NULL)
    gROOT->FindObject("back_f")->Delete();
 
  TH2S *reg_interest = cut_hist(moderate(data, hotpixel_map(dark, significance)), start_x, end_x, start_y, end_y);
  set_error(reg_interest, exp_time_min);
  TH1D *prj_x = reg_interest->ProjectionX("_px", 0, -1, "e");
  prj_x->GetYaxis()->SetTitle("Intensity (ADC counts)");
 
  return int_err(prj_x); 
}

TH2S *hot_pixel_coincidence(vector<const char *> dark_name, double significance)
{
	const int n = dark_name.size();

	TH2S *dark_data[n];
	TH2S *dark_map[n];		
	int *n_fold_coincidence = new int[n + 1] {}; 

	for(int i = 0; i < n; i++){
		TString map_name = Form("dark_%d", i);
		dark_data[i] = cut_hist(get_hist(dark_name[i], map_name.Data()), 101, 1500, 101, 1100);
		dark_map[i] = hotpixel_map(dark_data[i], significance);
		//dark_map[i] = cut_hist(hotpixel_map(get_hist(dark_name[i], map_name.Data()), significance), 101, 1500, 101, 1100);
	}

	TH2S *coincidence = new TH2S("coincidence", "Hot Pixel Coincidence Map; Pixel x coordinate; Pixel y coordinate", 1400, 100, 1500, 1000, 100, 1100);
	TH2S *incoincidence = new TH2S("incoincidence", "Hot Pixel Incoincidence Map; Pixel x coordinate; Pixel y coordinate", 1400, 100, 1500, 1000, 100, 1100);
	TH2S *coin_data = new TH2S("coin_data", "Hot Pixel Coincidence Data; Pixel x coordinate; Pixel y coordinate", 1400, 100, 1500, 1000, 100, 1100);
	
	int sum, ii;
	double sq_sum, mean, stdev;

	vector<int> on_off(n,0);
	vector<int> dark_val(n, 0);

	for(int i = 1; i <= 1500; i++){
		for(int j = 1; j <= 1100; j++){

			fill(on_off.begin(), on_off.end(), 0);
			sum = 0;

			for(int k = 0; k < n; k++){
				on_off[k] += dark_map[k]->GetBinContent(i, j);
				sum = accumulate(on_off.begin(), on_off.end(), 0);
			}

			if(sum == 1){
				ii = distance(on_off.begin(), find(on_off.begin(), on_off.end(), 1));
				incoincidence->SetBinContent(i, j, dark_data[ii]->GetBinContent(i, j));
				//cout << dark_data[ii]->GetBinContent(i, j) << endl;
			}

			if(sum == n){
				coincidence->SetBinContent(i, j, 1);

				for(int k = 0; k < n; k++)
					dark_val[k] = dark_data[k]->GetBinContent(i, j);

				sq_sum = inner_product(dark_val.begin(), dark_val.end(), dark_val.begin(), 0.0);
				mean = accumulate(dark_val.begin(), dark_val.end(), 0.0) / n;
				coin_data->SetBinContent(i, j, mean);
				stdev = sqrt(sq_sum / n - mean * mean);

				if(stdev / mean > 0.2)
					cout << "location: (" << i << ", " << j << ")\t" << "mean : " << (int)mean << "\tstdev : " << (int)stdev << endl;
				}

			for(int k = 1; k < n; k++){
				if(sum == k)
					n_fold_coincidence[k]++;	
			}
		}
	}
	
	TH1D *inc_amp = amp_spectrum(incoincidence, "inc_amp");
	TH1D *coi_amp = amp_spectrum(coin_data, "coi_amp");
	inc_amp->SetBinContent(1, 0);
	coi_amp->SetBinContent(1, 0);

	coi_amp->SetLineColor(kRed);
	coi_amp->SetLineWidth(2);
	inc_amp->Draw();
	coi_amp->Draw("same");
	
	for(int k = 0; k < n; k++)
		cout << Form("Percentage %d : ", k + 1) << 100 * coincidence->Integral() / dark_map[k]->Integral() << endl;

	cout << endl;

	for(int k = 1; k < n; k++)
		cout << Form("%d", k) << "-fold coincidence counts : " << n_fold_coincidence[k] << endl;
	
	cout << Form("%d", n) << "-fold coincidence counts : " << coincidence->Integral() << endl;

	delete [] n_fold_coincidence;

	return coincidence;		
}


void show_coincidence(const char *dark1_name, const char *dark2_name, double significance)
{
	TH2S *dark1 = cut_hist(get_hist(dark1_name, "dark1"), 101, 1500, 101, 1100);
	TH2S *dark2 = cut_hist(get_hist(dark2_name, "dark2"), 101, 1500, 101, 1100);	

	vector<const char*> dark_name = {dark1_name, dark2_name};

	TH2S *coincidence_map = hot_pixel_coincidence(dark_name, significance);
	TH2S *dark1_map = hotpixel_map(dark1, significance);
	TH2S *dark2_map = hotpixel_map(dark2, significance);

	TH1D *dark1_amp = amp_spectrum(dark1, "dark1_amp");
	TH1D *dark2_amp = amp_spectrum(dark2, "dark2_amp");

	//TH2S *incoincidence = (TH2S*)coincidence_map->Clone();
	//incoincidence->Reset();

	TH1D *empty_hist = (TH1D*)dark2_amp->Clone();
	empty_hist->Reset();

  TF1 *gf1 = new TF1("gf1", "gaus(0) + landau(3)", 0, dark1->GetMaximum());
	TF1 *gf2 = new TF1("gf2", "gaus(0) + landau(3)", 0, dark2->GetMaximum());

	gf1->SetParameters(1e5, 99, 5.4);
	gf2->SetParameters(1e5, 99, 5.4);
	gf1->SetParLimits(1, 85, 105);
	gf1->SetParLimits(2, 5, 6);
	
	gf2->SetParLimits(1, 85, 105);
	gf2->SetParLimits(2, 5, 6);

  dark1_amp->Fit(gf1);
  dark2_amp->Fit(gf2);

  double mean1 = gf1->GetParameter(1);
  double sigma1 = gf1->GetParameter(2);
	double mean2 = gf2->GetParameter(1);
	double sigma2 = gf2->GetParameter(2);

	double threshold1 = mean1 + significance * sigma1;
	double threshold2 = mean2 + significance * sigma2;

	cout << endl;
	cout << "Threshold 1 : " << threshold1 << " Threshold 2 : " << threshold2 << endl;

	double max_range = max(threshold1, threshold2) + 100;

	if(threshold1 < 150 && threshold2 < 150){
		dark1_amp->GetXaxis()->SetRangeUser(50, 250);
		dark2_amp->GetXaxis()->SetRangeUser(50, 250);	
  }

	else{
		dark1_amp->GetXaxis()->SetRangeUser(50, max_range);
		dark2_amp->GetXaxis()->SetRangeUser(50, max_range);
	}

	TH1D *dark1_th = (TH1D*)dark1_amp->Clone();
	TH1D *dark2_th = (TH1D*)dark2_amp->Clone();
	dark1_th->Reset();
	dark2_th->Reset();
	dark1_th->SetLineColor(kBlue);
	dark2_th->SetLineColor(kGreen + 2);

	dark1_th->Fill(threshold1, 1e4);
	dark2_th->Fill(threshold2, 1e4);

	dark1_th->GetYaxis()->SetRangeUser(0, dark1_amp->GetMaximum() + 1e4);
	dark2_th->GetYaxis()->SetRangeUser(0, dark2_amp->GetMaximum() + 1e4);

	dark1_amp->SetLineColor(kBlue);
	dark2_amp->SetLineColor(kGreen + 2);
	
	coincidence_map->SetMarkerStyle(3);
	dark1_map->SetMarkerStyle(3);
	dark2_map->SetMarkerStyle(3);
	coincidence_map->SetMarkerColor(kRed);
	//coincidence_map->SetMarkerSize(0.1);
	dark1_map->SetMarkerColor(kBlue);
	dark2_map->SetMarkerColor(kGreen + 2);

	dark1_th->SetLineStyle(2);
	dark2_th->SetLineStyle(2);

	double dark1_pixel_counts = dark1_map->Integral();
	double dark2_pixel_counts = dark2_map->Integral();
	double coincidence_counts = coincidence_map->Integral();

	auto legend1 = new TLegend(0.5, 0.7, 0.88, 0.88);
	legend1->AddEntry(dark1_map, "Outlying Pixels in Dark Image 1", "p");
	legend1->AddEntry(dark2_map, "Outlying Pixels in Dark Image 2", "p");
	legend1->AddEntry(coincidence_map, "Outlying Pixels of Coincidence", "p");
	legend1->SetTextAlign(22);

	auto legend2 = new TLegend(0.6, 0.7, 0.88, 0.88);
	legend2->AddEntry(dark1_amp, "Dark Image 1", "l");
	legend2->AddEntry(dark2_amp, "Dark Image 2", "l"); 
	legend2->SetTextAlign(22);

	TCanvas *c1 = new TCanvas("hot_pix_stability", "Hot Pixel Stability", 1800, 600);
	c1->Divide(2, 0, 0.00001, 0.001);  

	gStyle->SetOptStat(0);

	c1->cd(1);
	gPad->SetGrid();
	dark1_map->Draw();
	dark2_map->Draw("same");
	coincidence_map->Draw("same");
	legend1->Draw("same");

	c1->cd(2);
	gPad->SetGridx();
	dark1_th->Draw("hist same");	
	dark2_th->Draw("hist same");
	dark1_amp->Draw("hist same");
	dark2_amp->Draw("hist same");
	empty_hist->Draw("same");
	legend2->Draw("same");
	
	cout << endl;
	cout << "Number of Hot Pixels (Dark 1) : " << dark1_pixel_counts << endl;
	cout << "Number of Hot Pixels (Dark 2) : " << dark2_pixel_counts << endl;
	cout << "Number of Coincidence Counts  : " << coincidence_counts << endl;
	cout << endl;
	cout << "Coincidence Percentage (Dark 1) : " << 100 * coincidence_counts / dark1_pixel_counts << " %" << endl;
	cout << "Coincidence Percentage (Dark 2) : " << 100 * coincidence_counts / dark2_pixel_counts << " %" << endl;
	cout << endl;
}

void dir_moderate(TH2S *data, double threshold)
{	
	TH2S *mod_data = (TH2S*)data->Clone();
	
	vector<int> neighbors;
	vector<int> neighbor_map;

	int value;
	double average;
	//double max;

	for(int i = 1; i < data->GetNbinsX(); i++){
		for(int j = 1; j < data->GetNbinsY(); j++){
		
		neighbors.clear();
		neighbor_map.clear();

		//cout << i << " " << j << endl;	

		value = data->GetBinContent(i, j);

		if(value > threshold){
			
			for(int k = -1; k <= 1; k++){
				for(int l = -1; l <= 1; l++){
					neighbors.push_back(data->GetBinContent(i + k, j + l));
					if(data->GetBinContent(i + k, j + l) > threshold){
						neighbor_map.push_back(1);
					}
					else neighbor_map.push_back(0);
				}
			}
			
			average = get_average(neighbors, neighbor_map);
			//average = (accumulate(neighbors.begin(), neighbors.end(), 0.0) - value) / 8.0;
		
			//if(value > 1.2 * average){			
			cout << "Coordinates : (" << i << ", " << j << ")\tValue : " << value << " Average : " << average << endl;
	
			mod_data->SetBinContent(i, j, average);
			//}
		}
/*
			for(int k = -1; k <= 1; k++)
				for(int l = -1; l <= 1; l++)
					neighbors.push_back(data->GetBinContent(i + k, j + l));

			value = neighbors[4];
			max = *max_element(neighbors.begin(), neighbors.end());	
			average = (accumulate(neighbors.begin(), neighbors.end(), 0.0) - value) / 8.0;

			if(max == value && value > 1.2 * average){
				cout << "Coordinates : (" << i << ", " << j << ")\tValue : " << max << endl;
				data->SetBinContent(i, j, average);
			}	
*/
		}
	}

	mod_data->SetName("mod_data");
	TCanvas *c1 = new TCanvas("c1", "Comparison", 1800, 600);
	c1->Divide(2, 0, 0.00001, 0.001);
	gStyle->SetOptStat(0);

	c1->cd(1);
	data->Draw("colz");
	c1->cd(2);
	mod_data->Draw("colz");
}

void locate(TH2S *data, double threshold)
{
	int tot = 0;

	for(int i = 1; i <= data->GetNbinsX(); i++){
		for(int j = 1; j <= data->GetNbinsY(); j++){
			if(data->GetBinContent(i, j) > threshold){
				tot++;
				cout << "over threshold: " << endl;
				cout << "location: (" << i << ", " << j << ")\tValue: " << data->GetBinContent(i, j) << endl;
			}
			

		}
	}
cout << "total number: " << tot << endl;
}

TH2S *th_map(TH2S *data, int adc)
{
	TH2S *map = (TH2S*)data->Clone();
	map->Reset();
	map->SetName(Form("%d_map", adc));
	int sum = 0;

	for(int i = 1; i <= data->GetNbinsX(); i++){
		for(int j = 1; j <= data->GetNbinsY(); j++){
			if(data->GetBinContent(i, j) == adc){
				map->SetBinContent(i, j, 1);
				sum++;
	     // cout << "(" << i << ", " << j << ")\t" << data->GetBinContent(i, j) << endl;
			}
		}
	}
	cout << "Number of Pixels : " << sum << endl;

	return map;
}

TH2S *coincidence_map(TH2S *map1, TH2S *map2)
{
	TH2S *map = (TH2S*)map1->Clone();
	map->Reset();

	for(int i = 1; i <= map1->GetNbinsX(); i++){
		for(int j = 1; j <= map1->GetNbinsY(); j++){
			if(map1->GetBinContent(i, j) + map2->GetBinContent(i, j) == 2)
				map->SetBinContent(i, j, 1);
		}
	}
	
	return map;
}


