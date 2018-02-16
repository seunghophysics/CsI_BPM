using namespace std;

double get_threshold(TH2S *dark, double significance)
{
  int max_content = (int)dark->GetMaximum();
  
  TH1D *amp_spectrum = new TH1D("amp_spectrum", "Intensity Spectrum; Intensity; Counts;", max_content, 0, max_content);
  
  for(int i = 1; i <= 1600; i++)
    for(int j = 1; j <= 1200; j++)
      amp_spectrum->Fill(dark->GetBinContent(i, j));
  
  TF1 *gf = new TF1("gf", "gaus", 0, max_content);
  amp_spectrum->Fit(gf, "Q0");
  
  double mean = gf->GetParameter(1);
  double sigma = gf->GetParameter(2);
  
  amp_spectrum->Delete();

  return mean + significance * sigma;
}

TH2S *hotpixel_map(TH2S *dark, double significance)
{ 
  double threshold = get_threshold(dark, significance);
  
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

TH2S *moderate(TH2S *data, TH2S *hotpixel_map)
{
  TH2S *moderated_data = (TH2S*)data->Clone();

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

