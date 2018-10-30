//
TGraphErrors* g[4];

void read_table(TString infile)
{
  float ptmin,ptmax;
  float scale;
  float xs[4];
  float stat[4];
  float sys[4];

  ifstream ifs(infile);
  for (int i=0;i<4;i++)
    {
      g[i] = new TGraphErrors();
      g[i]->SetName("g_cent_"+TString::Itoa(i,10));
    }

  int index=0;
  while (!ifs.eof())
    {
      ifs >> ptmin >> ptmax;
      ifs >> scale;
      for (int i=0;i<4;i++)
	{
	  ifs >> xs[i] >> stat[i] >> sys[i];
	  g[i]->SetPoint(index,(ptmax+ptmin)/2,scale*xs[i]);
	  g[i]->SetPointError(index,(ptmax-ptmin)/2.,scale*sqrt(stat[i]*stat[i]+sys[i]*sys[i]));
	}
      index++;
    }
  
}
