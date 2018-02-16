double get_threshold(TH2S *dark, double significance)
{
  int max_content = (int)dark->GetMaximum();
  
  TH1D *amp_spectrum = new TH1D("amp_spectrum", "Pixel Intensity Spectrum", max_content, 0, max_content);
  
  for(int i = 1; i <= 1600; i++)
    for(int j = 1; j <= 1200; j++)
      amp_spectrum->Fill(dark->GetBinContent(i, j));
  
  TF1 *gf = new TF1("gf", "gaus", 0, max_content);
  amp_spectrum->Fit(gf, "Q0");
  
  double mean = gf->GetParameter(1);
  double sigma = gf->GetParameter(2);

  return mean + significance * sigma;
}

TH2S *hotpixel_map(TH2S *dark, double significance)
{ 
  double threshold = get_threshold(dark, significance);
 
  TH2S *map = new TH2S("hotpixel_map", "Hot Pixel Coordinates", 1600, 0, 1600, 1200, 0, 1200);
  
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

double get_average(std::vector<int> vec, std::vector<int> map)
{
  double sum = 0.0;
  double denominator = 0.0;

  for(int i = 0; i <= 8; i++){
    if(map[i]){
      sum += vec[i];
      std::cout << "denom: " << denominator << std::endl;
      denominator++;
    }
  }
  
  denominator += 1e-3;
  std::cout << "denom: " << denominator << std::endl;
  return sum / denominator;
}

TH2S *moderate(TH2S *data, TH2S *hotpixel_map)
{
  TH2S *moderated_data = (TH2S*)data->Clone();
  TH2S *map = (TH2S*)hotpixel_map->Clone();
  
  std::vector<int> neighbors;
  std::vector<int> neighbor_map;
  
  double average;

  for(int i = 2; i <= 1599; i++){
    for(int j = 2; j <= 1199; j++){

      std::cout << "i: " << i << " j: " << j << std::endl;

      if(map->GetBinContent(i, j)){

        for(int k = -1; k <= 1; k++){
          for(int l = -1; l <= 1; l++){
            printf("%d %d \n", k, l);
            std::cout << data->GetBinContent(i + k, j + l) << std::endl;
            std::cout << map->GetBinContent(i + k, j + l) << std::endl;
            neighbors.push_back(data->GetBinContent(i + k, j + l));
            neighbor_map.push_back(map->GetBinContent(i + k, j + l));
          }
        }      
      
      average = get_average(neighbors, neighbor_map);
      std::cout << "average: " << average << std::endl;
      moderated_data->SetBinContent(i, j, average);
      
      }
    }
  }
      
  return moderated_data;
}
