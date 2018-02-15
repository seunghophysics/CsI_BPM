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
