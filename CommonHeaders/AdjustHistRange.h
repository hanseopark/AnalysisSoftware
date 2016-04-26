//This function is used to automatically adjust the y-axis range of an histogram
//It retrieves the maximum and minimum value(ignoring empty bins) of the histogram
//and adds a range on both sides which is defined by the variable factor

//used to calculate min and max values of a set of histograms
Bool_t AdjustHistRange(TH1D* histogram, Double_t factorLow, Double_t factorHigh, Bool_t useBinError, Double_t* A_tmp, Double_t* B_tmp){

	if(histogram->GetEntries()<=0) return kFALSE;

	Int_t xfirst = 1;
	Double_t max = histogram->GetBinContent(xfirst);
	Double_t min = histogram->GetBinContent(xfirst);

	if(histogram->GetBinContent(xfirst)==0){
		do{
			xfirst++;
			max = histogram->GetBinContent(xfirst);
			min = histogram->GetBinContent(xfirst);
		}
		while(histogram->GetBinContent(xfirst)==0 && xfirst < histogram->GetNbinsX());
	}
	if(xfirst==histogram->GetNbinsX()){
		cout << "Warning in AdjustHistRange: Histogram with name: " << histogram->GetName() << " contains 0 in all bins, return..." << endl;
		return kFALSE;
	}

	Double_t valueMin = 0;
	Double_t valueMax = 0;

	Int_t xlast  = histogram->GetXaxis()->GetNbins();

	for(Int_t binx = xfirst; binx <= xlast; binx++){
		valueMax = histogram->GetBinContent(binx);
		if(useBinError) valueMax += histogram->GetBinError(binx);
		valueMin = histogram->GetBinContent(binx);
		if(useBinError) valueMin -= histogram->GetBinError(binx);
		if ((valueMin < min) && (valueMin!=0)) min = valueMin;
		if ((valueMax > max) && (valueMax!=0)) max = valueMax;
	}
	if(max<0) max/=factorHigh;
	else max*=factorHigh;
	if(min<=0) min*=factorLow;
	else min/=factorLow;

	*A_tmp = min;
	*B_tmp = max;
	return kTRUE;
}

void AdjustHistRange(TH1D* histogram, Double_t factorLow, Double_t factorHigh, Bool_t useBinError, Int_t fixRange = 0, Double_t range = 0){

	if(histogram->GetEntries()<=0){
		cout << "Warning in AdjustHistRange: No Entries in given histogram with name: " << histogram->GetName() << ", return..." << endl;
		return;
	}

	Double_t A, B;
	if(!AdjustHistRange(histogram,factorLow,factorHigh,useBinError,&A,&B)) {cout << "ERROR: AdjustHistRange(TH1D* histogram, Double_t factorLow, Double_t factorHigh, Bool_t useBinError, Double_t* A_tmp, Double_t* B_tmp) returned kFALSE! Returning..." << endl; return;}

	if(fixRange==0) histogram->GetYaxis()->SetRangeUser(A, B);
	else if(fixRange==1) histogram->GetYaxis()->SetRangeUser(range, B);
	else if(fixRange==2) histogram->GetYaxis()->SetRangeUser(A, range);
	else cout << "ERROR: AdjustRange.h wrong value of parameter fixRange given (valid 0,1 or 2): " << fixRange << endl;

    return;
}

void AdjustHistRange(TH2D* histogram, Double_t factorLow, Double_t factorHigh, Bool_t useBinError, Int_t fixRange = 0, Double_t range = 0){
	Double_t min = 0;
	Double_t max = 0;
	TH1D* hist = (TH1D*) histogram->ProjectionY("proj",1,histogram->GetNbinsX());
	AdjustHistRange(hist,factorLow,factorHigh,useBinError,min,max);

	if(fixRange==0) histogram->GetYaxis()->SetRangeUser(min, max);
	else if(fixRange==1) histogram->GetYaxis()->SetRangeUser(range, max);
	else if(fixRange==2) histogram->GetYaxis()->SetRangeUser(min, range);
	else cout << "ERROR: AdjustRange.h wrong value of parameter fixRange given (valid 0,1 or 2): " << fixRange << endl;

	return;
}



//overloaded function for AdjustHistRange
//This function calls AdjustHistRange for all histograms in the same vector
void AdjustHistRange(std::vector<TH1D*> vectorhist, Double_t factorLow, Double_t factorHigh, Bool_t useBinError, Int_t fixRange = 0, Double_t range = 0, Bool_t isYLog = kFALSE){
  
  if(!vectorhist.empty()){
    Double_t Max_global,Min_global;
    Double_t A,B;
    
	Bool_t successRange = AdjustHistRange(vectorhist.at(0),factorLow,factorHigh, useBinError, &A, &B);
    if(isYLog){
        if(A<0){
          successRange = kFALSE;
          A = 0;
        }
    }

	Int_t iRange = 1;
    if(!successRange && ((Int_t) vectorhist.size() > 1) ){
		while(!successRange){
			successRange = AdjustHistRange(vectorhist.at(iRange++),factorLow,factorHigh, useBinError, &A, &B);
			if(isYLog){
				if(A<0) successRange = kFALSE;
			}
            if(iRange>=(Int_t) vectorhist.size()){
                cout << "ERROR in AdjustHistRange, iRange>vectorhist.size(), hist: "<< ((TH1D*)vectorhist.at(0))->GetName() <<", returning..." << endl;
                return;
            }
		}
	}

	Min_global = A;
    Max_global = B;
    
	for(Int_t j=iRange; j<(Int_t) vectorhist.size(); j++){
	  if(!AdjustHistRange(vectorhist.at(j),factorLow,factorHigh, useBinError, &A, &B)) continue;
	  if(A<Min_global){
		  if(isYLog){
			  if(A>0) Min_global=A;
		  }else  Min_global=A;
	  }
      if(B>Max_global) Max_global=B;
	}

	if(fixRange==1) Min_global=range;
	else if(fixRange==2) Max_global=range;
	else if(fixRange!=0) cout << "ERROR: AdjustRange.h wrong value of parameter fixRange given (valid 0,1 or 2): " << fixRange << endl;

	for(Int_t i=0; i<(Int_t) vectorhist.size(); i++) vectorhist.at(i)->GetYaxis()->SetRangeUser(Min_global, Max_global);
  }
  
  return;
}

//overloaded function for AdjustHistRange
//This function calls AdjustHistRange for all histograms in the same vector
void AdjustHistRange(std::vector<TH1D*> vectorhist[], Double_t factorLow, Double_t factorHigh, Int_t h, Int_t nSets, Bool_t useBinError){

  Double_t Max_global,Min_global;
  Double_t A,B;

  Bool_t successRange = AdjustHistRange(vectorhist[0].at(h),factorLow,factorHigh, useBinError, &A, &B);

  Int_t iRange = 1;
  if(!successRange){
	  while(!successRange && iRange<nSets) successRange = AdjustHistRange(vectorhist[iRange++].at(h),factorLow,factorHigh, useBinError, &A, &B);
  }

  if(iRange>=nSets){
	  cout << "ERROR in AdjustHistRange, iRange>nSets, returning..." << endl;
	  return;
  }

  Min_global = A;
  Max_global = B;

  for(Int_t i=iRange; i<nSets; i++){
	  if(!AdjustHistRange(vectorhist[i].at(h),factorLow,factorHigh, useBinError, &A, &B)) continue;
	  if(A<Min_global) Min_global=A;
	  if(B>Max_global) Max_global=B;
  }

  for(Int_t i=0; i<nSets; i++) vectorhist[i].at(h)->GetYaxis()->SetRangeUser(Min_global, Max_global);

  return;
}





