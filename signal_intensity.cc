// still in progress.

#include "pco2root2x2.h"


using namespace std;

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
  
  TH1D *amp_spectrum = new TH1D("amp_spectrum", "Intensity Spectrum; Intensity; Counts;", max_content, 0, max_content);
  
  for(int i = 1; i <= 1600; i++)
    for(int j = 1; j <= 1200; j++)
      amp_spectrum->Fill(dark->GetBinContent(i, j));
  
  TF1 *gf = new TF1("gf", "gaus", 0, max_content);
  amp_spectrum->Fit(gf, "Q0");
  
  double mean = gf->GetParameter(1);
  double sigma = gf->GetParameter(2);
  
  if(draw_opt) 
    amp_spectrum->Draw();

  return mean + significance * sigma;
}

// Returns a 2D map of hot pixels, given a dark (or data) picture.
TH2S *hotpixel_map(TH2S *dark, double significance)
{ 
  double threshold = get_threshold(dark, significance, 0);
  
  if(gROOT->FindObject("hotpixel_map") != NULL)
    gROOT->FindObject("hotpixel_map")->Delete();
  
  TH2S *map = new TH2S("hotpixel_map", "Hot Pixel Coordinates; x coordinate; y coordinate", 1600, 0, 1600, 1200, 0, 1200);
  
  for(int i = 1; i <= 1600; i++){
    for(int j = 1; j <= 1200; j++){
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
  if(gROOT->FindObject("cut_histogram") != NULL)
    gROOT->FindObject("cut_histogram")->Delete();
 
  TH2S *cut_histo = new TH2S("cut_histogram", "Cut Histogram; x coordinate; y coordinate;", end_binx - start_binx + 1, start_binx - 1, end_binx, end_biny - start_biny + 1, start_biny - 1, end_biny);

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

// Calculates the signal intensity of a region of interest, given dark, data, and significance for hot pixels.
// Example usage: intensity(get_hist("folder/file.root"), get_hist("folder/file.root"), 5, 400, 1200, 200, 600);
double intensity(TH2S *data, TH2S *dark, double significance, int start_x, int end_x, int start_y, int end_y)
{
  if(gROOT->FindObject("back_f") != NULL)
    gROOT->FindObject("back_f")->Delete();
 
  TH2S *reg_interest = cut_hist(moderate(data, hotpixel_map(dark, significance)), start_x, end_x, start_y, end_y);
  TH1D *prj_x = reg_interest->ProjectionX("_px", 0, -1, "e");
  prj_x->GetYaxis()->SetTitle("Intensity (ADC counts)");
 
  TF1 *sig_back_f = new TF1("sig_back_f", "gaus(0) + pol2(3)", 0, 1600);
  TF1 *back_f = new TF1("back_f", "pol2", start_x - 1, end_x);

  sig_back_f->SetParameters(prj_x->GetMaximum() - prj_x->GetMinimum(), start_x + prj_x->GetMaximumBin(), 50, prj_x->GetMinimum(), -10, -0.05);
  sig_back_f->SetParLimits(0, 0, 5e4);
  sig_back_f->SetParLimits(1, start_x + prj_x->GetMaximumBin() - 10, start_x + prj_x->GetMaximumBin() + 10);
  sig_back_f->SetParLimits(2, 10, 1e2);
  sig_back_f->SetParLimits(4, 0, 5e1);
  sig_back_f->SetParLimits(5, -5e-1, 0);

  prj_x->Fit("sig_back_f", "be");
  back_f->SetParameters(sig_back_f->GetParameter(3), sig_back_f->GetParameter(4), sig_back_f->GetParameter(5));
  back_f->Draw("same");

  return reg_interest->Integral(1, end_x - start_x + 1) - back_f->Integral(start_x - 1, end_x);
}
