#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <sys/types.h>
#include <dirent.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>

using namespace std;using namespace ROOT::Math;

struct runlist
{
	string runno;
	string b[3];
	int t[3];
	string date;
	string time;
};
runlist rl;
vector <runlist> RL;

struct finallist
{
	string runno;
	string b;
	int t;
	int s;
};

struct calib
{
	string pmt;
	double g;
	double e;
};

struct edata
{
	int ieta[144];
	int iphi[144];
	int depth[144];
	int pulse[144][10];
	int histo[144][60];
};

struct gain
{
	string *BoxName;
	string *BoxBarcode;
	string *PMTName;
	string *PMTSerial;
	Float_t l2014[2];//gain : error
	Float_t il[2];//gain : error
	Float_t ilV2[2];//gain : error
	Float_t ih[2];//gain : error
	Float_t ihV2[2];//gain : error
	Float_t fl[2][2];//[1+2 : 3+4] [gain : error]
	Float_t fl2[2][2];//[1+2 : 3+4] [gain : error]
	Float_t flV2[2][2];//[1+2 : 3+4] [gain : error]
	Float_t fh[2][2];//[1+2 : 3+4] [gain : error]
	Float_t fhV2[2][2];//[1+2 : 3+4] [gain : error]
};
gain g;

struct npelist
{
	string b;
	string runno[2];
	int t[2];
	int s[2];
};

//from JM 070516

double adc2fC_QIE10[256]={

  // =========== RANGE 0 ===========

  // --------- subrange 1 ---------
  -14.45,-11.35,-8.25,-5.15,-2.05,1.05,4.15,7.25,10.35,13.45,16.55,19.65,22.75,25.85,28.95,32.05,
  // --------- subrange 2 ---------
  36.7,42.9,49.1,55.3,61.5,67.7,73.9,80.1,86.3,92.5,98.7,104.9,111.1,117.3,123.5,129.7,135.9,142.1,148.3,154.5,
  // --------- subrange 3 ---------
  163.8,176.2,188.6,201.0,213.4,225.8,238.2,250.6,263.0,275.4,287.8,300.2,312.6,325.0,337.4,349.8,362.2,374.6,387.0,399.4,411.8,
  // --------- subrange 4 ---------
  430.4,455.2,480.0,504.8,529.6,554.4,579.2,
  // =========== RANGE 1 ===========

  // --------- subrange 1 ---------
  529.4,554.2,579.0,603.8,628.6,653.4,678.2,703.0,727.8,752.6,777.4,802.2,827.0,851.8,876.6,901.4,
  // --------- subrange 2 ---------
  938.6,988.2,1037.8,1087.4,1137.0,1186.6,1236.2,1285.8,1335.4,1385.0,1434.6,1484.2,1533.8,1583.4,1633.0,1682.6,1732.2,1781.8,1831.4,1881.0,
  // --------- subrange 3 ---------
  1955.4,2054.6,2153.8,2253.0,2352.2,2451.4,2550.6,2649.8,2749.0,2848.2,2947.4,3046.6,3145.8,3245.0,3344.2,3443.4,3542.6,3641.8,3741.0,3840.2,3939.4,
  // --------- subrange 4 ---------
  4088.2,4286.6,4485.0,4683.4,4881.8,5080.2,5278.6,
  // =========== RANGE 2 ===========

  // --------- subrange 1 ---------
  4879.2,5077.6,5276.0,5474.4,5672.8,5871.2,6069.6,6268.0,6466.4,6664.8,6863.2,7061.6,7260.0,7458.4,7656.8,7855.2,
  // --------- subrange 2 ---------
  8152.8,8549.6,8946.4,9343.2,9740.0,10136.8,10533.6,10930.4,11327.2,11724.0,12120.8,12517.6,12914.4,13311.2,13708.0,14104.8,14501.6,14898.4,15295.2,15692.0,
  // --------- subrange 3 ---------
  16287.2,17080.8,17874.4,18668.0,19461.6,20255.2,21048.8,21842.4,22636.0,23429.6,24223.2,25016.8,25810.4,26604.0,27397.6,28191.2,28984.8,29778.4,30572.0,31365.6,32159.2,
  // --------- subrange 4 ---------
  33349.6,34936.8,36524.0,38111.2,39698.4,41285.6,42872.8,
  // =========== RANGE 3 ===========

  // --------- subrange 1 ---------
  39693.5,41280.5,42867.5,44454.5,46041.5,47628.5,49215.5,50802.5,52389.5,53976.5,55563.5,57150.5,58737.5,60324.5,61911.5,63498.5,
  // --------- subrange 2 ---------
  65879.0,69053.0,72227.0,75401.0,78575.0,81749.0,84923.0,88097.0,91271.0,94445.0,97619.0,100793.0,103967.0,107141.0,110315.0,113489.0,116663.0,119837.0,123011.0,126185.0,
  // --------- subrange 3 ---------
  130946.0,137294.0,143642.0,149990.0,156338.0,162686.0,169034.0,175382.0,181730.0,188078.0,194426.0,200774.0,207122.0,213470.0,219818.0,226166.0,232514.0,238862.0,245210.0,251558.0,257906.0,
  // --------- subrange 4 ---------
  267428.0,280124.0,292820.0,305516.0,318212.0,330908.0,343604.0

};

//station - PMT - 1+2|3+4 - ieta|iphi|depth (or fiber|channel|uHTR)
int MAP2Ch[3][24][2][3]={{{{7,2,1},{10,3,1}},{{10,2,1},{7,3,1}},{{6,0,1},{9,1,1}},{{9,0,1},{6,1,1}},{{8,0,1},{11,1,1}},{{11,0,1},{8,1,1}},{{6,2,1},{9,3,1}},{{9,2,1},{6,3,1}},{{8,2,1},{11,3,1}},{{11,2,1},{8,3,1}},{{7,0,1},{10,1,1}},{{10,0,1},{7,1,1}},{{1,2,1},{4,3,1}},{{4,2,1},{1,3,1}},{{0,0,1},{3,1,1}},{{3,0,1},{0,1,1}},{{2,0,1},{5,1,1}},{{5,0,1},{2,1,1}},{{0,2,1},{3,3,1}},{{3,2,1},{0,3,1}},{{2,2,1},{5,3,1}},{{5,2,1},{2,3,1}},{{1,0,1},{4,1,1}},{{4,0,1},{1,1,1}}},{{{1,2,2},{4,3,2}},{{4,2,2},{1,3,2}},{{0,0,2},{3,1,2}},{{3,0,2},{0,1,2}},{{2,0,2},{5,1,2}},{{5,0,2},{2,1,2}},{{0,2,2},{3,3,2}},{{3,2,2},{0,3,2}},{{2,2,2},{5,3,2}},{{5,2,2},{2,3,2}},{{1,0,2},{4,1,2}},{{4,0,2},{1,1,2}},{{7,2,2},{10,3,2}},{{10,2,2},{7,3,2}},{{6,0,2},{9,1,2}},{{9,0,2},{6,1,2}},{{8,0,2},{11,1,2}},{{11,0,2},{8,1,2}},{{6,2,2},{9,3,2}},{{9,2,2},{6,3,2}},{{8,2,2},{11,3,2}},{{11,2,2},{8,3,2}},{{7,0,2},{10,1,2}},{{10,0,2},{7,1,2}}},{{{19,2,2},{22,3,2}},{{22,2,2},{19,3,2}},{{18,0,2},{21,1,2}},{{21,0,2},{18,1,2}},{{20,0,2},{23,1,2}},{{23,0,2},{20,1,2}},{{18,2,2},{21,3,2}},{{21,2,2},{18,3,2}},{{20,2,2},{23,3,2}},{{23,2,2},{20,3,2}},{{19,0,2},{22,1,2}},{{22,0,2},{19,1,2}},{{13,2,2},{16,3,2}},{{16,2,2},{13,3,2}},{{12,0,2},{15,1,2}},{{15,0,2},{12,1,2}},{{14,0,2},{17,1,2}},{{17,0,2},{14,1,2}},{{12,2,2},{15,3,2}},{{15,2,2},{12,3,2}},{{14,2,2},{17,3,2}},{{17,2,2},{14,3,2}},{{13,0,2},{16,1,2}},{{16,0,2},{13,1,2}}}};

string ChNames[24][3]={{"A2","A2_1+2","A2_3+4"},{"A4","A4_1+2","A4_3+4"},{"A6","A6_1+2","A6_3+4"},{"A8","A8_1+2","A8_3+4"},{"A1","A1_1+2","A1_3+4"},{"A3","A3_1+2","A3_3+4"},{"A5","A5_1+2","A5_3+4"},{"A7","A7_1+2","A7_3+4"},{"B2","B2_1+2","B2_3+4"},{"B4","B4_1+2","B4_3+4"},{"B6","B6_1+2","B6_3+4"},{"B8","B8_1+2","B8_3+4"},{"B1","B1_1+2","B1_3+4"},{"B3","B3_1+2","B3_3+4"},{"B5","B5_1+2","B5_3+4"},{"B7","B7_1+2","B7_3+4"},{"C2","C2_1+2","C2_3+4"},{"C4","C4_1+2","C4_3+4"},{"C6","C6_1+2","C6_3+4"},{"C8","C8_1+2","C8_3+4"},{"C1","C1_1+2","C1_3+4"},{"C3","C3_1+2","C3_3+4"},{"C5","C5_1+2","C5_3+4"},{"C7","C7_1+2","C7_3+4"}};

inline bool SPEexists (const string& name)
{
    ifstream f(name.c_str());
    return f.good();
}

int FindMode(string b,string r)
{
	vector <runlist> RL1;
	ifstream infile("RunList.txt");
	while(!infile.eof())
	{
		infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
		RL1.push_back(rl);
	}
	RL1.pop_back();
	infile.close();
	int rt=-1;
	for(int i1=0;i1<RL1.size();i1++)
	{
		if(RL1[i1].runno==r)
		{
			for(int i2=0;i2<3;i2++)
			{
				if(RL1[i1].b[i2]==b && (RL1[i1].t[i2]==1 || RL1[i1].t[i2]==2 || RL1[i1].t[i2]==12 || RL1[i1].t[i2]==13 || RL1[i1].t[i2]==14)){rt=1;break;}
				if(RL1[i1].b[i2]==b && (RL1[i1].t[i2]==3 || RL1[i1].t[i2]==4 || RL1[i1].t[i2]==5 || RL1[i1].t[i2]==6 || RL1[i1].t[i2]==7 || RL1[i1].t[i2]==15 || RL1[i1].t[i2]==16 || RL1[i1].t[i2]==17)){rt=2;break;}
			}
		}
		if(rt!=-1) break;
	}
	if(rt==-1)
	{
		for(int i1=0;i1<RL1.size();i1++)
		{
			if(RL1[i1].runno<=r)
			{
				for(int i2=0;i2<3;i2++)
				{
					if(RL1[i1].b[i2]==b && (RL1[i1].t[i2]==3 || RL1[i1].t[i2]==4 || RL1[i1].t[i2]==5 || RL1[i1].t[i2]==6 || RL1[i1].t[i2]==7 || RL1[i1].t[i2]==15 || RL1[i1].t[i2]==16 || RL1[i1].t[i2]==17)){rt=2;break;}
				}
				if(rt!=-1) break;
			}
		}
		if(rt==-1) rt=1;
	}
	return rt;
}

int plotpeds(string r, string b, int id, int s)
{
	char hname[500];
	edata ed;
	sprintf(hname,"../NTuples/N_%s.root",r.c_str());
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
// 	tree->SetBranchAddress("histo",&ed.histo);
	
	sprintf(hname,"../Histos/%s/Pedestals_%d_%s.root",b.c_str(),id,r.c_str());
	TFile* outroot=new TFile(hname,"recreate");
	
	TF1* tf=new TF1("gaus","gaus",1.,10.);
	TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
	TCanvas* cc2=new TCanvas("cc2","cc2",900,1200);
	gStyle->SetOptStat(0);
	
	int m=FindMode(b,r);
	if(m==1) cc1->Divide(4,6,0,0);
	else cc1->Divide(8,6,0,0);
	cc2->SetGridy();cc2->Divide(1,2);
	
	TGraphErrors* tg1[2];
	tg1[0]=new TGraphErrors();tg1[0]->SetName("Pedestals_Mean_1+2");
	tg1[0]->SetMarkerStyle(24);tg1[0]->SetMarkerColor(1);
	tg1[1]=new TGraphErrors();tg1[1]->SetName("Pedestals_Mean_3+4");
	tg1[1]->SetMarkerStyle(25);tg1[1]->SetMarkerColor(2);
	
	TGraphErrors* tg2[2];
	tg2[0]=new TGraphErrors();tg2[0]->SetName("Pedestals_Sigma_1+2");
	tg2[0]->SetMarkerStyle(24);tg2[0]->SetMarkerColor(1);
	tg2[1]=new TGraphErrors();tg2[1]->SetName("Pedestals_Sigma_3+4");
	tg2[1]->SetMarkerStyle(25);tg2[1]->SetMarkerColor(2);
	
	TH1F* ADCH[24][2];
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			sprintf(hname,"ADCperTS %s",ChNames[i1][m+i2-2].c_str());
			ADCH[i1][i2-1]=new TH1F(hname,hname,10,-0.5,9.5);
		}
	}
	
	vector <int> clist[3];
	tree->GetEntry(0);
	for(int iz1=0;iz1<144;iz1++)
	{
		for(int i1=0;i1<24;i1++)
		{
			for(int i2=1;i2<=m;i2++)
			{
				if(ed.ieta[iz1]==MAP2Ch[s][i1][i2-1][0] && ed.iphi[iz1]==MAP2Ch[s][i1][i2-1][1] && ed.depth[iz1]==MAP2Ch[s][i1][i2-1][2])
				{
					clist[0].push_back(iz1);
					clist[1].push_back(i1);
					clist[2].push_back(i2-1);
				}
			}
		}
	}
	
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		for(int i1=0;i1<clist[0].size();i1++)
		{
			for(int i2=0;i2<10;i2++)
			{
				ADCH[clist[1][i1]][clist[2][i1]]->Fill(ed.pulse[clist[0][i1]][i2]);
			}
		}
	}
	
	double ymax=0.;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			if(ADCH[i1][i2-1]->GetBinContent(ADCH[i1][i2-1]->GetMaximumBin())>ymax){ymax=ADCH[i1][i2-1]->GetBinContent(ADCH[i1][i2-1]->GetMaximumBin());}
		}
	}
	
	outroot->cd();
	int tci1=1;
	double tg1maxmin[2]={0.,10000.};
	double tg2maxmin[2]={0.,10000.};
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			cc1->cd(tci1);
			ADCH[i1][i2-1]->Draw();cc1->Update();
			ADCH[i1][i2-1]->SetLineWidth(2);ADCH[i1][i2-1]->SetLineColor(1);ADCH[i1][i2-1]->GetXaxis()->SetRangeUser(0.,10.);ADCH[i1][i2-1]->GetYaxis()->SetRangeUser(0.,ymax+20000.);
			ADCH[i1][i2-1]->Fit(tf,"q","q",1.,10.);
			int nn=tg1[i2-1]->GetN();
			tg1[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(1));
			tg1[i2-1]->SetPointError(nn,0,tf->GetParError(1));
			tg2[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(2));
			tg2[i2-1]->SetPointError(nn,0,tf->GetParError(2));
			if(tf->GetParameter(1)>tg1maxmin[0]){tg1maxmin[0]=tf->GetParameter(1);}if(tf->GetParameter(1)<tg1maxmin[1]){tg1maxmin[1]=tf->GetParameter(1);}
			if(tf->GetParameter(2)>tg2maxmin[0]){tg2maxmin[0]=tf->GetParameter(2);}if(tf->GetParameter(2)<tg2maxmin[1]){tg2maxmin[1]=tf->GetParameter(2);}
			tci1++;
			ADCH[i1][i2-1]->Write();
		}
	}
	
	sprintf(hname,"Pedestals_%d.pdf(",id);
	cc1->Print(hname);
	cc2->cd(1);
	tg1[0]->Draw("AP");tg1[0]->GetYaxis()->SetTitle("Mean Pedestal ADC");tg1[0]->GetYaxis()->CenterTitle();tg1[0]->GetXaxis()->SetTitle("PMT ID");tg1[0]->GetXaxis()->CenterTitle();
	tg1[0]->GetYaxis()->SetRangeUser(tg1maxmin[1]-0.5,tg1maxmin[0]+0.5);gPad->SetGridy(1);gPad->SetGridx(1);cc2->Update();
	if(m==2) tg1[1]->Draw("same P");
	cc2->cd(2);
	tg2[0]->Draw("AP");tg2[0]->GetYaxis()->SetTitle("Pedestal ADC Sigma");tg2[0]->GetYaxis()->CenterTitle();tg2[0]->GetXaxis()->SetTitle("PMT ID");tg2[0]->GetXaxis()->CenterTitle();
	tg2[0]->GetYaxis()->SetRangeUser(tg2maxmin[1]-0.5,tg2maxmin[0]+0.5);gPad->SetGridy(1);gPad->SetGridx(1);cc2->Update();
	if(m==2) tg2[1]->Draw("same P");
	cc2->SetGridy();cc2->RedrawAxis();cc2->Update();
	sprintf(hname,"Pedestals_%d.pdf)",id);
	cc2->Print(hname);
	if(id==9){sprintf(hname,"mv Pedestals_%d.pdf Pedestals_%d_%s.pdf",id,id,r.c_str());system(hname);}
	sprintf(hname,"mv Pedestals_*.pdf ../Plots/%s",b.c_str());system(hname);
	
	outroot->cd();
	for(int i2=1;i2<=m;i2++)
	{
		tg1[i2-1]->Write();
		tg2[i2-1]->Write();
	}
	
	inroot->Close();
	outroot->Close();
	delete cc1;delete cc2;
}

int plotleds(string r, string b, int id, int s)
{
	char hname[500];
	edata ed;
	sprintf(hname,"../NTuples/N_%s.root",r.c_str());
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
// 	tree->SetBranchAddress("histo",&ed.histo);
	
	sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",b.c_str(),id,r.c_str());
	TFile* outroot=new TFile(hname,"recreate");
	
	sprintf(hname,"LEDGains_%s_%d.txt",b.c_str(),id);
	ofstream outfile(hname);
	
	TF1* tf=new TF1("gaus","gaus",0.,100000.);
	
	int m=FindMode(b,r);
	
	TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
	gStyle->SetOptStat(0);
	if(m==1) {cc1->Divide(4,6,0,0);}
	else {cc1->Divide(8,6,0,0);}
	TCanvas* cc3=new TCanvas("cc3","cc3",600,600);
	
	TGraphErrors* tg1[2];
	tg1[0]=new TGraphErrors();tg1[0]->SetName("LED_Mean_1+2");
	tg1[0]->SetMarkerStyle(24);tg1[0]->SetMarkerColor(1);
	tg1[1]=new TGraphErrors();tg1[1]->SetName("LED_Mean_3+4");
	tg1[1]->SetMarkerStyle(25);tg1[1]->SetMarkerColor(2);
	
	TGraphErrors* tg2[2];
	tg2[0]=new TGraphErrors();tg2[0]->SetName("LED_Sigma_1+2");
	tg2[0]->SetMarkerStyle(24);tg2[0]->SetMarkerColor(1);
	tg2[1]=new TGraphErrors();tg2[1]->SetName("LED_Sigma_3+4");
	tg2[1]->SetMarkerStyle(25);tg2[1]->SetMarkerColor(2);
	
	TGraphErrors* tg3[2];
	tg3[0]=new TGraphErrors();tg3[0]->SetName("LEDNPE_Mean_1+2");
	tg3[0]->SetMarkerStyle(24);tg3[0]->SetMarkerColor(1);
	tg3[1]=new TGraphErrors();tg3[1]->SetName("LEDNPE_Mean_3+4");
	tg3[1]->SetMarkerStyle(25);tg3[1]->SetMarkerColor(2);
	
	TGraphErrors* tg4[2];
	tg4[0]=new TGraphErrors();tg4[0]->SetName("LEDNPE_Sigma_1+2");
	tg4[0]->SetMarkerStyle(24);tg4[0]->SetMarkerColor(1);
	tg4[1]=new TGraphErrors();tg4[1]->SetName("LEDNPE_Sigma_3+4");
	tg4[1]->SetMarkerStyle(25);tg4[1]->SetMarkerColor(2);
	
	TH1F* QPulse[24][2][2];TH1F* QTot[24][2];TH1F* NPETot[24][2];TH1F* NPETotC[24];TH1F* Frac[24];TH1F* ADCPH[24][2];
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			sprintf(hname,"QPulse %s",ChNames[i1][m+i2-2].c_str());
			QPulse[i1][i2-1][0]=new TH1F(hname,hname,10,-0.5,9.5);
			QPulse[i1][i2-1][0]->GetXaxis()->SetTitle("TS (x25 ns)");
			QPulse[i1][i2-1][0]->GetXaxis()->CenterTitle();
			QPulse[i1][i2-1][0]->GetYaxis()->SetTitle("Mean Charge per TS (fC)");
			QPulse[i1][i2-1][0]->GetYaxis()->CenterTitle();
			sprintf(hname,"QPulse %s norm",ChNames[i1][m+i2-2].c_str());
			QPulse[i1][i2-1][1]=new TH1F(hname,hname,10,-0.5,9.5);
			sprintf(hname,"QTot %s",ChNames[i1][m+i2-2].c_str());
			QTot[i1][i2-1]=new TH1F(hname,hname,800,0.,20000.);
			QTot[i1][i2-1]->GetXaxis()->SetTitle("Charge (fC)");
			QTot[i1][i2-1]->GetXaxis()->CenterTitle();
			QTot[i1][i2-1]->GetYaxis()->SetTitle("Entries / 25 fC");
			QTot[i1][i2-1]->GetYaxis()->CenterTitle();
			sprintf(hname,"NPETot %s",ChNames[i1][m+i2-2].c_str());
// 			NPETot[i1][i2-1]=new TH1F(hname,hname,250,0.,500.);
			NPETot[i1][i2-1]=new TH1F(hname,hname,500,0.,500.);
			NPETot[i1][i2-1]->GetXaxis()->SetTitle("Npe");
			NPETot[i1][i2-1]->GetXaxis()->CenterTitle();
			NPETot[i1][i2-1]->GetYaxis()->SetTitle("Entries / pe");
			NPETot[i1][i2-1]->GetYaxis()->CenterTitle();
			sprintf(hname,"ADCPH %s",ChNames[i1][m+i2-2].c_str());
			ADCPH[i1][i2-1]=new TH1F(hname,hname,256,-0.5,255.5);
			ADCPH[i1][i2-1]->GetXaxis()->SetTitle("ADC Pulse Height");
			ADCPH[i1][i2-1]->GetXaxis()->CenterTitle();
			ADCPH[i1][i2-1]->GetYaxis()->SetTitle("Entries / pule height");
			ADCPH[i1][i2-1]->GetYaxis()->CenterTitle();
		}
	}
	if(m==2)
	{
		for(int i1=0;i1<24;i1++)
		{
			sprintf(hname,"NPETot %s",ChNames[i1][0].c_str());
			NPETotC[i1]=new TH1F(hname,hname,500,0.,500.);
			NPETotC[i1]->GetXaxis()->SetTitle("Npe");
			NPETotC[i1]->GetXaxis()->CenterTitle();
			NPETotC[i1]->GetYaxis()->SetTitle("Entries / pe");
			NPETotC[i1]->GetYaxis()->CenterTitle();
			
			sprintf(hname,"Frac %s",ChNames[i1][0].c_str());
			Frac[i1]=new TH1F(hname,hname,60,0.,3.);
			Frac[i1]->GetXaxis()->SetTitle("3+4/1+2");
			Frac[i1]->GetXaxis()->CenterTitle();
			Frac[i1]->GetYaxis()->SetTitle("Entries / 0.05");
			Frac[i1]->GetYaxis()->CenterTitle();
		}
	}
	
	vector <int> clist[3];
	tree->GetEntry(0);
	for(int iz1=0;iz1<144;iz1++)
	{
		for(int i1=0;i1<24;i1++)
		{
			for(int i2=1;i2<=m;i2++)
			{
				if(ed.ieta[iz1]==MAP2Ch[s][i1][i2-1][0] && ed.iphi[iz1]==MAP2Ch[s][i1][i2-1][1] && ed.depth[iz1]==MAP2Ch[s][i1][i2-1][2])
				{
					clist[0].push_back(iz1);
					clist[1].push_back(i1);
					clist[2].push_back(i2-1);
				}
			}
		}
	}
	double Q[10]={0.};double Qtot=0.;double ped=0.;
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		for(int i1=0;i1<clist[0].size();i1++)
		{
			Qtot=0.;
			for(int i2=0;i2<10;i2++)
			{
				Q[i2]=adc2fC_QIE10[ed.pulse[clist[0][i1]][i2]];
				if(i2>=2 && i2<=6){Qtot+=Q[i2];}
				ADCPH[clist[1][i1]][clist[2][i1]]->Fill(ed.pulse[clist[0][i1]][i2]);
			}
			ped=(Q[0]+Q[1])/2.;
			Qtot-=(5.*ped);
			for(int i2=0;i2<10;i2++)
			{
				QPulse[clist[1][i1]][clist[2][i1]][0]->Fill(i2,Q[i2]);
				QPulse[clist[1][i1]][clist[2][i1]][1]->Fill(i2);
			}
			QTot[clist[1][i1]][clist[2][i1]]->Fill(Qtot);
		}
	}
	
	double ymaxpulse=0.;
	double ymaxQtot=0.;double xmaxQtot=0.;
	double tg1max=0.;double tg1min=10000.;
	double tg2max=0.;double tg2min=10000.;
	double G[24][2]={{0.}};double Ge[24][2]={{0.}};
// 	double tfmean=0.;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			QPulse[i1][i2-1][0]->Divide(QPulse[i1][i2-1][1]);
			if(QPulse[i1][i2-1][0]->GetBinContent(QPulse[i1][i2-1][0]->GetMaximumBin())>ymaxpulse){ymaxpulse=QPulse[i1][i2-1][0]->GetBinContent(QPulse[i1][i2-1][0]->GetMaximumBin());}
			if(QTot[i1][i2-1]->GetBinContent(QTot[i1][i2-1]->GetMaximumBin())>ymaxQtot){ymaxQtot=QTot[i1][i2-1]->GetBinContent(QTot[i1][i2-1]->GetMaximumBin());}
			if(QTot[i1][i2-1]->GetBinCenter(QTot[i1][i2-1]->FindLastBinAbove(0))>xmaxQtot){xmaxQtot=QTot[i1][i2-1]->GetBinCenter(QTot[i1][i2-1]->FindLastBinAbove(0));}
			
			tf->SetLineStyle(2);
			QTot[i1][i2-1]->GetXaxis()->SetRange(1,5);
			if(QTot[i1][i2-1]->GetMaximumBin()==1)
			{
				QTot[i1][i2-1]->GetXaxis()->SetRange(1,800);
				QTot[i1][i2-1]->Fit(tf,"q","q",50.,20000.);
				QTot[i1][i2-1]->Fit(tf,"q","q",(tf->GetParameter(1)-1.5*tf->GetParameter(2))>50.?(tf->GetParameter(1)-1.5*tf->GetParameter(2)):50.,tf->GetParameter(1)+1.5*tf->GetParameter(2));
			}
			else
			{
				QTot[i1][i2-1]->GetXaxis()->SetRange(1,800);
				QTot[i1][i2-1]->Fit(tf,"q","q",0.,20000.);
				QTot[i1][i2-1]->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
			}
// 			tfmean=QTot[i1][i2-1]->GetXaxis()->GetBinCenter(QTot[i1][i2-1]->GetMaximumBin());
// 			QTot[i1][i2-1]->GetXaxis()->SetRange(1,800);
// 			tf->SetParameter(1,tfmean);tf->SetParameter(2,70.);
// 			if(tfmean<100.) QTot[i1][i2-1]->Fit(tf,"q","q",0.,20000.);
// 			else QTot[i1][i2-1]->Fit(tf,"q","q",25.,20000.);
			
			int nn=tg1[i2-1]->GetN();
			tg1[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(1));
			tg1[i2-1]->SetPointError(nn,0,tf->GetParError(1));
			tg2[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(2));
			tg2[i2-1]->SetPointError(nn,0,tf->GetParError(2));
			
			if(tf->GetParameter(1)>tg1max){tg1max=tf->GetParameter(1);}if(tf->GetParameter(1)<tg1min){tg1min=tf->GetParameter(1);}
			if(tf->GetParameter(2)>tg2max){tg2max=tf->GetParameter(2);}if(tf->GetParameter(2)<tg2min){tg2min=tf->GetParameter(2);}
			
			G[i1][i2-1]=((pow(tf->GetParameter(2),2.)/(tf->GetParameter(1)*1.15))*10000./1.6);
			Ge[i1][i2-1]=(sqrt(pow(tf->GetParError(1),2.)+pow(sqrt(2.)*tf->GetParError(2),2.))*10000./(1.16*1.15));
			
			outfile<<ChNames[i1][m+i2-2]<<" "<<G[i1][i2-1]<<" "<<Ge[i1][i2-1]<<endl;
		}
	}
	
	double npe=0.;double Npe[24][2]={{0.}};
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		for(int i1=0;i1<clist[0].size();i1++)
		{
			Qtot=0.;
			for(int i2=0;i2<10;i2++)
			{
				Q[i2]=adc2fC_QIE10[ed.pulse[clist[0][i1]][i2]];
				if(i2>=2 && i2<=6){Qtot+=Q[i2];}
			}
			ped=(Q[0]+Q[1])/2.;
			Qtot-=(5.*ped);
			npe=Qtot/(G[clist[1][i1]][clist[2][i1]]*1.6/10000.);
			NPETot[clist[1][i1]][clist[2][i1]]->Fill(npe);
			Npe[clist[1][i1]][clist[2][i1]]=npe;
		}
		if(m==2)
		{
			for(int i1=0;i1<24;i1++)
			{
				NPETotC[i1]->Fill(Npe[i1][0]+Npe[i1][1]);
				if(Npe[i1][0]>0) Frac[i1]->Fill(Npe[i1][1]/Npe[i1][0]);
			}
		}
	}
	
	double ymaxNPE=0.;double xmaxNPE=0.;
	double tg3max=0.;double tg3min=10000.;
	double tg4max=0.;double tg4min=10000.;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			if(NPETot[i1][i2-1]->GetBinContent(NPETot[i1][i2-1]->GetMaximumBin())>ymaxNPE){ymaxNPE=NPETot[i1][i2-1]->GetBinContent(NPETot[i1][i2-1]->GetMaximumBin());}
			if(NPETot[i1][i2-1]->GetBinCenter(NPETot[i1][i2-1]->FindLastBinAbove(0))>xmaxNPE){xmaxNPE=NPETot[i1][i2-1]->GetBinCenter(NPETot[i1][i2-1]->FindLastBinAbove(0));}
			
			tf->SetLineStyle(2);
			NPETot[i1][i2-1]->GetXaxis()->SetRange(1,5);
			if(NPETot[i1][i2-1]->GetMaximumBin()==1)
			{
				NPETot[i1][i2-1]->GetXaxis()->SetRange(1,500);
				if((NPETot[i1][i2-1]->Integral(1,3)/NPETot[i1][i2-1]->Integral())>0.8)
				{
					NPETot[i1][i2-1]->Fit(tf,"q","q",0.,500.);
					NPETot[i1][i2-1]->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
				}
				else
				{
					NPETot[i1][i2-1]->Fit(tf,"q","q",2.,500.);
					NPETot[i1][i2-1]->Fit(tf,"q","q",(tf->GetParameter(1)-1.5*tf->GetParameter(2))>2.?(tf->GetParameter(1)-1.5*tf->GetParameter(2)):2.,tf->GetParameter(1)+1.5*tf->GetParameter(2));
				}
			}
			else
			{
				NPETot[i1][i2-1]->GetXaxis()->SetRange(1,800);
				NPETot[i1][i2-1]->Fit(tf,"q","q",0.,500.);
				NPETot[i1][i2-1]->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
			}
// 			tfmean=NPETot[i1][i2-1]->GetXaxis()->GetBinCenter(NPETot[i1][i2-1]->GetMaximumBin());
// 			NPETot[i1][i2-1]->GetXaxis()->SetRange(1,500);
// 			tf->SetParameter(1,tfmean);tf->SetParameter(2,3.);
// 			if(tfmean<5.) NPETot[i1][i2-1]->Fit(tf,"q","q",0.,500.);
// 			else NPETot[i1][i2-1]->Fit(tf,"q","q",1.,500.);
			
			int nn=tg3[i2-1]->GetN();
			tg3[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(1));
			tg3[i2-1]->SetPointError(nn,0,tf->GetParError(1));
			tg4[i2-1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(2));
			tg4[i2-1]->SetPointError(nn,0,tf->GetParError(2));
			
			if(tf->GetParameter(1)>tg3max){tg3max=tf->GetParameter(1);}if(tf->GetParameter(1)<tg3min){tg3min=tf->GetParameter(1);}
			if(tf->GetParameter(2)>tg4max){tg4max=tf->GetParameter(2);}if(tf->GetParameter(2)<tg4min){tg4min=tf->GetParameter(2);}
		}
	}
	
	int tci1=1;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			cc1->cd(tci1);
			QPulse[i1][i2-1][0]->Draw("hist");
			QPulse[i1][i2-1][0]->GetYaxis()->SetRangeUser(-10.,ymaxpulse+50.);
			QPulse[i1][i2-1][0]->SetLineColor(4);QPulse[i1][i2-1][0]->SetFillColor(4);
// 			cc1->Update();
			tci1++;
		}
	}
	sprintf(hname,"LEDs_%d.pdf(",id);
	cc1->Print(hname);
	
	tci1=1;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			cc1->cd(tci1);
			QTot[i1][i2-1]->Draw();
			QTot[i1][i2-1]->GetXaxis()->SetRangeUser(0.,xmaxQtot+50.);
			QTot[i1][i2-1]->GetYaxis()->SetRangeUser(0.,ymaxQtot+150.);
// 			cc1->Update();
			tci1++;
		}
	}
	sprintf(hname,"LEDs_%d.pdf",id);
	cc1->Print(hname);
// 	cc3->Update();
	tci1=1;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			cc1->cd(tci1);
			NPETot[i1][i2-1]->Draw();
			NPETot[i1][i2-1]->GetXaxis()->SetRangeUser(0.,xmaxNPE+10.);
			NPETot[i1][i2-1]->GetYaxis()->SetRangeUser(0.,ymaxNPE+150.);
// 			cc1->Update();
			tci1++;
		}
	}
// 	cc1->Update();
	sprintf(hname,"LEDs_%d.pdf",id);
	cc1->Print(hname);
	delete cc1;
	
	TCanvas* cc2=new TCanvas("cc2","cc2",900,1200);
	cc2->Divide(1,4);
	cc2->cd(1);
	tg1[0]->Draw("AP");tg1[0]->GetYaxis()->SetTitle("LED Charge Mean (fC)");tg1[0]->GetYaxis()->CenterTitle();tg1[0]->GetXaxis()->SetTitle("PMT ID");tg1[0]->GetXaxis()->CenterTitle();gPad->SetGridy(1);gPad->SetGridx(1);tg1[0]->GetYaxis()->SetRangeUser(tg1min-100.,tg1max+100.);
	if(m==2){tg1[1]->Draw("same P");}
	cc2->Update();
	cc2->cd(2);
	tg2[0]->Draw("AP");tg2[0]->GetYaxis()->SetTitle("LED Charge Sigma (fC)");tg2[0]->GetYaxis()->CenterTitle();tg2[0]->GetXaxis()->SetTitle("PMT ID");tg2[0]->GetXaxis()->CenterTitle();gPad->SetGridy(1);gPad->SetGridx(1);tg2[0]->GetYaxis()->SetRangeUser(tg2min-50.,tg2max+50.);
	if(m==2){tg2[1]->Draw("same P");}
	cc2->Update();
	
	cc2->cd(3);
	tg3[0]->Draw("AP");tg3[0]->GetYaxis()->SetTitle("Average Number of Photoelectrons");tg3[0]->GetYaxis()->CenterTitle();tg3[0]->GetXaxis()->SetTitle("PMT ID");tg3[0]->GetXaxis()->CenterTitle();gPad->SetGridy(1);gPad->SetGridx(1);tg3[0]->GetYaxis()->SetRangeUser(tg3min-10.,tg3max+10.);
	if(m==2){tg3[1]->Draw("same P");}
	cc2->Update();
	cc2->cd(4);
	tg4[0]->Draw("AP");tg4[0]->GetYaxis()->SetTitle("<NPE> Sigma");tg4[0]->GetYaxis()->CenterTitle();tg4[0]->GetXaxis()->SetTitle("PMT ID");tg4[0]->GetXaxis()->CenterTitle();gPad->SetGridy(1);gPad->SetGridx(1);tg4[0]->GetYaxis()->SetRangeUser(tg4min-5.,tg4max+5.);
	if(m==2){tg4[1]->Draw("same P");}
	cc2->Update();
	
	sprintf(hname,"LEDs_%d.pdf)",id);
	cc2->Print(hname);
	
	if(id==10){sprintf(hname,"mv LEDs_%d.pdf LEDs_%d_%s.pdf",id,id,r.c_str());system(hname);}
	sprintf(hname,"mv LEDs_*.pdf ../Plots/%s",b.c_str());
	system(hname);
	
	TGraphErrors* tg5[2];
	tg5[0]=new TGraphErrors();tg5[0]->SetName("TotalNPE");
	tg5[0]->SetMarkerStyle(24);tg5[0]->SetMarkerColor(1);
	tg5[1]=new TGraphErrors();tg5[1]->SetName("TotalNPE_Sigma");
	tg5[1]->SetMarkerStyle(25);tg5[1]->SetMarkerColor(2);
	if(m==2)
	{
		for(int i1=0;i1<24;i1++)
		{
			tf->SetLineStyle(2);
			NPETotC[i1]->Fit(tf,"q","q",2.,500.);
			int nn=tg5[0]->GetN();
			tg5[0]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(1));
			tg5[0]->SetPointError(nn,0,tf->GetParError(1));
			tg5[1]->SetPoint(nn,((double)(i1+1)),tf->GetParameter(2));
			tg5[1]->SetPointError(nn,0,tf->GetParError(2));
		}
	}
	
	outroot->cd();
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			QPulse[i1][i2-1][0]->Write();
			QTot[i1][i2-1]->Write();
			NPETot[i1][i2-1]->Write();
			ADCPH[i1][i2-1]->Write();
			if(i1==0)
			{
				tg1[i2-1]->Write();
				tg2[i2-1]->Write();
				tg3[i2-1]->Write();
				tg4[i2-1]->Write();
			}
		}
		if(m==2){NPETotC[i1]->Write();Frac[i1]->Write();}
	}
	if(m==2){tg5[0]->Write();tg5[1]->Write();}
	outroot->Close();
	
	inroot->Close();
	outfile.close();
	if(id==10){sprintf(hname,"mv LEDGains_%s_%d.txt LEDGains_%s_%d_%s.txt",b.c_str(),id,b.c_str(),id,r.c_str());system(hname);}
	sprintf(hname,"mv LEDGains_*.txt ../Plots/%s",b.c_str());
	system(hname);
	
// 	gain g;
	TFile* gainroot=new TFile("../Gains/Gains.root");
	TTree *gtree = (TTree*)gainroot->Get("Gains");
	gtree->SetBranchAddress("BoxName",&g.BoxName);
	gtree->SetBranchAddress("BoxBarcode",&g.BoxBarcode);
	gtree->SetBranchAddress("PMTName",&g.PMTName);
	gtree->SetBranchAddress("PMTSerial",&g.PMTSerial);
	gtree->SetBranchAddress("l2014",&g.l2014);
	gtree->SetBranchAddress("il",&g.il);
	gtree->SetBranchAddress("ilV2",&g.ilV2);
	gtree->SetBranchAddress("ih",&g.ih);
	gtree->SetBranchAddress("ihV2",&g.ihV2);
	gtree->SetBranchAddress("fl",&g.fl);
	gtree->SetBranchAddress("fl2",&g.fl2);
	gtree->SetBranchAddress("flV2",&g.flV2);
	gtree->SetBranchAddress("fh",&g.fh);
	gtree->SetBranchAddress("fhV2",&g.fhV2);
	
	TFile* goutroot=new TFile("../Gains/Gains_out.root","recreate");
	TTree* treeout = new TTree("Gains", "Gains");
	treeout->Branch("BoxName", &g.BoxName);
	treeout->Branch("BoxBarcode", &g.BoxBarcode);
	treeout->Branch("PMTName", &g.PMTName);
	treeout->Branch("PMTSerial", &g.PMTSerial);
	treeout->Branch("l2014", g.l2014, "l2014[2]/F");
	treeout->Branch("il", g.il, "il[2]/F");
	treeout->Branch("ilV2", g.ilV2, "ilV2[2]/F");
	treeout->Branch("ih", g.ih, "ih[2]/F");
	treeout->Branch("ihV2", g.ihV2, "ihV2[2]/F");
	treeout->Branch("fl", g.fl, "fl[2][2]/F");
	treeout->Branch("fl2", g.fl2, "fl2[2][2]/F");
	treeout->Branch("flV2", g.flV2, "flV2[2][2]/F");
	treeout->Branch("fh", g.fh, "fh[2][2]/F");
	treeout->Branch("fhV2", g.fhV2, "fhV2[2][2]/F");
	
	for(int i=0;i<gtree->GetEntries();i++)
	{
		gtree->GetEntry(i);
		if(*g.BoxBarcode==b)
		{
			for(int i1=0;i1<24;i1++)
			{
				if(*g.PMTName==ChNames[i1][0])
				{
					if(id==2){g.il[0]=G[i1][0];g.il[1]=Ge[i1][0];}
					else if(id==14){g.ilV2[0]=G[i1][0];g.ilV2[1]=Ge[i1][0];}
					else if(id==4){g.fl[0][0]=G[i1][0];g.fl[0][1]=Ge[i1][0];g.fl[1][0]=G[i1][1];g.fl[1][1]=Ge[i1][1];}
					else if(id==6){g.fl2[0][0]=G[i1][0];g.fl2[0][1]=Ge[i1][0];g.fl2[1][0]=G[i1][1];g.fl2[1][1]=Ge[i1][1];}
					else if(id==17){g.flV2[0][0]=G[i1][0];g.flV2[0][1]=Ge[i1][0];g.flV2[1][0]=G[i1][1];g.flV2[1][1]=Ge[i1][1];}
				}
			}
		}
		treeout->Fill();
	}
	goutroot->cd();
// 	treeout->Write();
	goutroot->Write();
	goutroot->Close();
	gainroot->Close();
	system("mv ../Gains/Gains_out.root ../Gains/Gains.root");
	
// 	g.BoxName=0;g.BoxBarcode=0;g.PMTName=0;g.PMTSerial=0;
// 	delete g.BoxName;delete g.BoxBarcode; delete g.PMTName;delete g.PMTSerial;
	
// 	delete gtree;delete treeout;
	
// 	delete inroot;delete outroot;delete gainroot;delete goutroot;
	
// 	delete cc1;
	delete cc2;delete cc3;
// 	delete cc3;delete cc5;
}

string ledcomp(vector <string> r, vector <int> s)//beware when comparing runs of different modes
{
	char hname[300];
	TF1* tf=new TF1("gaus","gaus",0.,100000.);
	TGraphErrors* tg1[3][2];
	TGraphErrors* tg2[3][2];
	TGraphErrors* tg3[3][2];
	TGraphErrors* tg4[3][2];
	TGraphErrors* tgg[3][2];
	
	vector <string> b;vector <int> id;vector <int> m;
	ifstream infile("RunList.txt");
	vector <runlist> RLc;
	while(!infile.eof())
	{
		infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
		RLc.push_back(rl);
	}
	RLc.pop_back();
	infile.close();
	for(int i2=0;i2<RLc.size();i2++)
	{
		for(int ii1=0;ii1<r.size();ii1++)
		{
			if(RLc[i2].runno==r[ii1])
			{
				b.push_back(RLc[i2].b[s[ii1]-1]);
				id.push_back(RLc[i2].t[s[ii1]-1]);
				m.push_back(FindMode(RLc[i2].b[s[ii1]-1],r[ii1]));
			}
		}
	}
	
	int colors[3]={1,2,4};
	TH1F* QTot[3][24][2];TH1F* NPE[3][24][2];
	double gains[3][24][2][2]={{{{0.}}}};
	TF1 *fit;
	
	for(int ii1=0;ii1<r.size();ii1++)
	{
		for(int i2=1;i2<=m[ii1];i2++)
		{
			tgg[ii1][i2-1]=new TGraphErrors();
			tgg[ii1][i2-1]->SetMarkerStyle(24+(i2-1));
			tgg[ii1][i2-1]->SetMarkerColor(colors[i2]);
		}
	}
	
	TFile* inroot[3];
	string gnames[2]={"1+2","3+4"};
	double qtotmax[2]={0.,0.};double npemax[2]={0.,0.};//x, y
	double tgsmaxmin[4][2]={{0.,10000.},{0.,10000.},{0.,10000.},{0.,10000.}};
	int dd=0;int de=0;
	for(int ii1=0;ii1<r.size();ii1++)
	{
		sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",b[ii1].c_str(),id[ii1],r[ii1].c_str());
		inroot[ii1]=new TFile(hname);
		for(int i1=0;i1<24;i1++)
		{
			for(int i2=1;i2<=m[ii1];i2++)
			{
				sprintf(hname,"QTot %s",ChNames[i1][m[ii1]+i2-2].c_str());
				inroot[ii1]->GetObject(hname,QTot[ii1][i1][i2-1]);
				QTot[ii1][i1][i2-1]->SetLineColor(colors[ii1]);
				if(QTot[ii1][i1][i2-1]->GetBinContent(QTot[ii1][i1][i2-1]->GetMaximumBin())>qtotmax[1]){qtotmax[1]=QTot[ii1][i1][i2-1]->GetBinContent(QTot[ii1][i1][i2-1]->GetMaximumBin());}
				if(QTot[ii1][i1][i2-1]->GetBinCenter(QTot[ii1][i1][i2-1]->FindLastBinAbove(0))>qtotmax[0]){qtotmax[0]=QTot[ii1][i1][i2-1]->GetBinCenter(QTot[ii1][i1][i2-1]->FindLastBinAbove(0));}
				fit=QTot[ii1][i1][i2-1]->GetFunction("gaus");
				gains[ii1][i1][i2-1][0]=pow(fit->GetParameter(2),2.)/(1.15*fit->GetParameter(1));
				gains[ii1][i1][i2-1][1]=gains[ii1][i1][i2-1][0]*sqrt(pow(2.*fit->GetParError(2)/fit->GetParameter(2),2.)+pow(fit->GetParError(1)/fit->GetParameter(1),2.));
				
				sprintf(hname,"NPETot %s",ChNames[i1][m[ii1]+i2-2].c_str());
				inroot[ii1]->GetObject(hname,NPE[ii1][i1][i2-1]);
				NPE[ii1][i1][i2-1]->SetLineColor(colors[ii1]);
				if(NPE[ii1][i1][i2-1]->GetBinContent(NPE[ii1][i1][i2-1]->GetMaximumBin())>npemax[1]){npemax[1]=NPE[ii1][i1][i2-1]->GetBinContent(NPE[ii1][i1][i2-1]->GetMaximumBin());}
				if(NPE[ii1][i1][i2-1]->GetBinCenter(NPE[ii1][i1][i2-1]->FindLastBinAbove(0))>npemax[0]){npemax[0]=NPE[ii1][i1][i2-1]->GetBinCenter(NPE[ii1][i1][i2-1]->FindLastBinAbove(0));}
				
				if(i1==0)
				{
					sprintf(hname,"LED_Mean_%s",gnames[i2-1].c_str());
					inroot[ii1]->GetObject(hname,tg1[ii1][i2-1]);
					tg1[ii1][i2-1]->SetMarkerStyle(24+i2-1);tg1[ii1][i2-1]->SetMarkerColor(colors[i2-1]);
					dd=TMath::LocMax(tg1[ii1][i2-1]->GetN(),tg1[ii1][i2-1]->GetY());
					if(tg1[ii1][i2-1]->GetY()[dd]>tgsmaxmin[0][0]){tgsmaxmin[0][0]=tg1[ii1][i2-1]->GetY()[dd];}
					de=TMath::LocMin(tg1[ii1][i2-1]->GetN(),tg1[ii1][i2-1]->GetY());
					if(tg1[ii1][i2-1]->GetY()[de]<tgsmaxmin[0][1]){tgsmaxmin[0][1]=tg1[ii1][i2-1]->GetY()[de];}
					
					sprintf(hname,"LED_Sigma_%s",gnames[i2-1].c_str());
					inroot[ii1]->GetObject(hname,tg2[ii1][i2-1]);
					tg2[ii1][i2-1]->SetMarkerStyle(24+i2-1);tg2[ii1][i2-1]->SetMarkerColor(colors[i2-1]);
					dd=TMath::LocMax(tg2[ii1][i2-1]->GetN(),tg2[ii1][i2-1]->GetY());
					if(tg2[ii1][i2-1]->GetY()[dd]>tgsmaxmin[1][0]){tgsmaxmin[1][0]=tg2[ii1][i2-1]->GetY()[dd];}
					de=TMath::LocMin(tg2[ii1][i2-1]->GetN(),tg2[ii1][i2-1]->GetY());
					if(tg2[ii1][i2-1]->GetY()[de]<tgsmaxmin[1][1]){tgsmaxmin[1][1]=tg2[ii1][i2-1]->GetY()[de];}
					
					sprintf(hname,"LEDNPE_Mean_%s",gnames[i2-1].c_str());
					inroot[ii1]->GetObject(hname,tg3[ii1][i2-1]);
					tg3[ii1][i2-1]->SetMarkerStyle(24+i2-1);tg3[ii1][i2-1]->SetMarkerColor(colors[i2-1]);
					dd=TMath::LocMax(tg3[ii1][i2-1]->GetN(),tg3[ii1][i2-1]->GetY());
					if(tg3[ii1][i2-1]->GetY()[dd]>tgsmaxmin[2][0]){tgsmaxmin[2][0]=tg3[ii1][i2-1]->GetY()[dd];}
					de=TMath::LocMin(tg3[ii1][i2-1]->GetN(),tg3[ii1][i2-1]->GetY());
					if(tg3[ii1][i2-1]->GetY()[de]<tgsmaxmin[2][1]){tgsmaxmin[2][1]=tg3[ii1][i2-1]->GetY()[de];}
					
					sprintf(hname,"LEDNPE_Sigma_%s",gnames[i2-1].c_str());
					inroot[ii1]->GetObject(hname,tg4[ii1][i2-1]);
					tg4[ii1][i2-1]->SetMarkerStyle(24+i2-1);tg4[ii1][i2-1]->SetMarkerColor(colors[i2-1]);
					dd=TMath::LocMax(tg4[ii1][i2-1]->GetN(),tg4[ii1][i2-1]->GetY());
					if(tg4[ii1][i2-1]->GetY()[dd]>tgsmaxmin[3][0]){tgsmaxmin[3][0]=tg4[ii1][i2-1]->GetY()[dd];}
					de=TMath::LocMin(tg4[ii1][i2-1]->GetN(),tg4[ii1][i2-1]->GetY());
					if(tg4[ii1][i2-1]->GetY()[de]<tgsmaxmin[3][1]){tgsmaxmin[3][1]=tg4[ii1][i2-1]->GetY()[de];}
				}
			}
		}
	}
	double tggmaxmin[2]={0.,1000.};
	for(int ii1=0;ii1<r.size();ii1++)
	{
		for(int i2=1;i2<=m[ii1];i2++)
		{
			for(int i1=0;i1<24;i1++)
			{
				tgg[ii1][i2-1]->SetPoint(i1,((double)(i1+1)),gains[ii1][i1][i2-1][0]);
				tgg[ii1][i2-1]->SetPointError(i1,0.,gains[ii1][i1][i2-1][1]);
				if(gains[ii1][i1][i2-1][0]>tggmaxmin[0]){tggmaxmin[0]=gains[ii1][i1][i2-1][0];}
				if(gains[ii1][i1][i2-1][0]<tggmaxmin[1]){tggmaxmin[1]=gains[ii1][i1][i2-1][0];}
			}
		}
	}
	
	TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
	gStyle->SetOptStat(0);
	if(m[0]==1) {cc1->Divide(4,6,0,0);}
	else {cc1->Divide(8,6,0,0);}
	
	char pdfname[500];
	if(r.size()==2) sprintf(pdfname,"LEDs_%s_%d_%s_%d.pdf",r[0].c_str(),s[0],r[1].c_str(),s[1]);
	else sprintf(pdfname,"LEDs_%s_%d_%s_%d_%s_%d.pdf",r[0].c_str(),s[0],r[1].c_str(),s[1],r[2].c_str(),s[2]);
	
	int tci1=1;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m[0];i2++)
		{
			cc1->cd(tci1);
			for(int ii1=0;ii1<r.size();ii1++)
			{
				if(ii1==0) QTot[ii1][i1][i2-1]->Draw();
				else QTot[ii1][i1][i2-1]->Draw("same");
				QTot[ii1][i1][i2-1]->GetXaxis()->SetRangeUser(0.,qtotmax[0]+100.);
				QTot[ii1][i1][i2-1]->GetYaxis()->SetRangeUser(0.,qtotmax[1]+150.);
			}
			tci1++;
		}
	}
	sprintf(hname,"%s(",pdfname);
	cc1->Print(hname);
	
	tci1=1;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m[0];i2++)
		{
			cc1->cd(tci1);
			for(int ii1=0;ii1<r.size();ii1++)
			{
				if(ii1==0) NPE[ii1][i1][i2-1]->Draw();
				else NPE[ii1][i1][i2-1]->Draw("same");
				NPE[ii1][i1][i2-1]->GetXaxis()->SetRangeUser(0.,npemax[0]+5.);
				NPE[ii1][i1][i2-1]->GetYaxis()->SetRangeUser(0.,npemax[1]+150.);
			}
			tci1++;
		}
	}
	sprintf(hname,"%s",pdfname);
	cc1->Print(hname);
	
	TCanvas* cc2=new TCanvas("cc2","cc2",900,1200);
	gStyle->SetOptStat(0);
	cc2->Divide(r.size(),5,0,0);
	tci1=1;
	for(int ii1=0;ii1<r.size();ii1++)
	{
		cc2->cd(tci1);
		tg1[ii1][0]->Draw("AP");
		tg1[ii1][0]->GetYaxis()->SetRangeUser(tgsmaxmin[0][1]-100.,tgsmaxmin[0][0]+100.);
		tg1[ii1][0]->GetXaxis()->SetTitle("PMT ID");tg1[ii1][0]->GetXaxis()->CenterTitle();
		tg1[ii1][0]->GetYaxis()->SetTitle("Mean LED Charge (fC)");tg1[ii1][0]->GetYaxis()->CenterTitle();
		gPad->SetGridy(1);gPad->SetGridx(1);
		if(m[ii1]==2) tg1[ii1][1]->Draw("same P");
		tci1++;
	}
	for(int ii1=0;ii1<r.size();ii1++)
	{
		cc2->cd(tci1);
		tg2[ii1][0]->Draw("AP");tg2[ii1][0]->GetYaxis()->SetRangeUser(tgsmaxmin[1][1]-50.,tgsmaxmin[1][0]+50.);
		tg2[ii1][0]->GetXaxis()->SetTitle("PMT ID");tg2[ii1][0]->GetXaxis()->CenterTitle();
		tg2[ii1][0]->GetYaxis()->SetTitle("LED Charge Sigma (fC)");tg2[ii1][0]->GetYaxis()->CenterTitle();
		gPad->SetGridy(1);gPad->SetGridx(1);
		if(m[ii1]==2) tg2[ii1][1]->Draw("same P");
		tci1++;
	}
	for(int ii1=0;ii1<r.size();ii1++)
	{
		cc2->cd(tci1);
		tg3[ii1][0]->Draw("AP");tg3[ii1][0]->GetYaxis()->SetRangeUser(tgsmaxmin[2][1]-10.,tgsmaxmin[2][0]+10.);
		tg3[ii1][0]->GetXaxis()->SetTitle("PMT ID");tg3[ii1][0]->GetXaxis()->CenterTitle();
		tg3[ii1][0]->GetYaxis()->SetTitle("Mean NPE");tg3[ii1][0]->GetYaxis()->CenterTitle();
		gPad->SetGridy(1);gPad->SetGridx(1);
		if(m[ii1]==2) tg3[ii1][1]->Draw("same P");
		tci1++;
	}
	for(int ii1=0;ii1<r.size();ii1++)
	{
		cc2->cd(tci1);
		tg4[ii1][0]->Draw("AP");tg4[ii1][0]->GetYaxis()->SetRangeUser(tgsmaxmin[3][1]-5.,tgsmaxmin[3][0]+5.);
		tg4[ii1][0]->GetXaxis()->SetTitle("PMT ID");tg4[ii1][0]->GetXaxis()->CenterTitle();
		tg4[ii1][0]->GetYaxis()->SetTitle("NPE Sigma");tg4[ii1][0]->GetYaxis()->CenterTitle();
		gPad->SetGridy(1);gPad->SetGridx(1);
		if(m[ii1]==2) tg4[ii1][1]->Draw("same P");
		tci1++;
	}
	for(int ii1=0;ii1<r.size();ii1++)
	{
		cc2->cd(tci1);
		tgg[ii1][0]->Draw("AP");tgg[ii1][0]->GetYaxis()->SetRangeUser(tggmaxmin[1]-5.,tggmaxmin[0]+5.);
		tgg[ii1][0]->GetXaxis()->SetTitle("PMT ID");tgg[ii1][0]->GetXaxis()->CenterTitle();
		tgg[ii1][0]->GetYaxis()->SetTitle("Gain (fC)");tgg[ii1][0]->GetYaxis()->CenterTitle();
		gPad->SetGridy(1);gPad->SetGridx(1);
		if(m[ii1]==2) tgg[ii1][1]->Draw("same P");
		tci1++;
	}
	sprintf(hname,"%s)",pdfname);
	cc2->Print(hname);
	
	for(int ii1=0;ii1<r.size();ii1++){inroot[ii1]->Close();}
	sprintf(hname,"mv LEDs_*pdf ../Plots/Others");
	system(hname);
	
	delete cc1;delete cc2;
	return b[0];
}

int ifcomp(vector <string> r, vector <int> s)
{
	char hname[300];
	TF1* tf=new TF1("gaus","gaus",0.,1000.);
	
	vector <string> b;vector <int> id;
	ifstream infile("RunList.txt");
	vector <runlist> RLc;
	while(!infile.eof())
	{
		infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
		RLc.push_back(rl);
	}
	RLc.pop_back();
	infile.close();
	for(int i2=0;i2<RLc.size();i2++)
	{
		for(int ii1=0;ii1<r.size();ii1++)
		{
			if(RLc[i2].runno==r[ii1])
			{
				b.push_back(RLc[i2].b[s[ii1]-1]);
				id.push_back(RLc[i2].t[s[ii1]-1]);
			}
		}
	}
	TFile* inroot[2];
	TGraphErrors* tg1[2];
	TGraphErrors* tg2;
	TGraphErrors* tg3[2];
	TGraphErrors* tg4;
	TGraphErrors* tg5[3];
	
	sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",b[0].c_str(),id[0],r[0].c_str());inroot[0]=new TFile(hname);
	sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",b[1].c_str(),id[1],r[1].c_str());inroot[1]=new TFile(hname);
	
	int dd=0;int de=0;double tg1maxmin[2]={0.,10000.};
	inroot[0]->GetObject("LEDNPE_Mean_1+2",tg1[0]);
	tg1[0]->SetMarkerStyle(24);tg1[0]->SetMarkerColor(1);
	dd=TMath::LocMax(tg1[0]->GetN(),tg1[0]->GetY());
	if(tg1[0]->GetY()[dd]>tg1maxmin[0]){tg1maxmin[0]=tg1[0]->GetY()[dd];}
	de=TMath::LocMin(tg1[0]->GetN(),tg1[0]->GetY());
	if(tg1[0]->GetY()[de]<tg1maxmin[1]){tg1maxmin[1]=tg1[0]->GetY()[de];}
	
	inroot[1]->GetObject("TotalNPE",tg1[1]);
	tg1[1]->SetMarkerStyle(25);tg1[1]->SetMarkerColor(2);
	dd=TMath::LocMax(tg1[1]->GetN(),tg1[1]->GetY());
	if(tg1[1]->GetY()[dd]>tg1maxmin[0]){tg1maxmin[0]=tg1[1]->GetY()[dd];}
	de=TMath::LocMin(tg1[1]->GetN(),tg1[1]->GetY());
	if(tg1[1]->GetY()[de]<tg1maxmin[1]){tg1maxmin[1]=tg1[1]->GetY()[de];}
	
	tg2=new TGraphErrors();//tg2->SetName("Dev");tg2->SetTitle("Dev");
	tg2->SetMarkerStyle(20);tg2->SetMarkerColor(1);
	for(int i1=0;i1<tg1[0]->GetN();i1++)
	{
		tg2->SetPoint(i1,((double)(i1+1)),(tg1[1]->GetY()[i1]-tg1[0]->GetY()[i1])/tg1[0]->GetY()[i1]);
		tg2->SetPointError(i1,0.,(tg1[1]->GetY()[i1]/tg1[0]->GetY()[i1])*sqrt(pow(tg1[1]->GetEY()[i1]/tg1[1]->GetY()[i1],2.)+pow(tg1[0]->GetEY()[i1]/tg1[0]->GetY()[i1],2.)));
	}
	
	inroot[1]->GetObject("LEDNPE_Mean_1+2",tg3[0]);
	inroot[1]->GetObject("LEDNPE_Mean_3+4",tg3[1]);
	tg4=new TGraphErrors();//tg4->SetName("Rat");tg4->SetTitle("Rat");
	tg4->SetMarkerStyle(21);tg4->SetMarkerColor(2);
	for(int i1=0;i1<tg3[0]->GetN();i1++)
	{
// 		tg4->SetPoint(i1,((double)(i1+1)),(tg3[1]->GetY()[i1]/tg3[0]->GetY()[i1]));
// 		tg4->SetPointError(i1,0.,(tg3[1]->GetY()[i1]/tg3[0]->GetY()[i1])*sqrt(pow(tg3[1]->GetEY()[i1]/tg3[1]->GetY()[i1],2.)+pow(tg3[0]->GetEY()[i1]/tg3[0]->GetY()[i1],2.)));
		
		tg4->SetPoint(i1,((double)(i1+1)),(tg3[0]->GetY()[i1]-tg3[1]->GetY()[i1])/(tg3[0]->GetY()[i1]+tg3[1]->GetY()[i1]));
		tg4->SetPointError(i1,0.,sqrt(2.*(pow(tg3[1]->GetY()[i1],2.)+pow(tg3[0]->GetY()[i1],2.))*(pow(tg3[1]->GetEY()[i1],2.)+pow(tg3[0]->GetEY()[i1],2.)))/(pow(tg3[1]->GetY()[i1]+tg3[0]->GetY()[i1],2.)));
	}
	
	TFile* gainroot=new TFile("../Gains/Gains.root");
	TTree *gtree = (TTree*)gainroot->Get("Gains");
	gtree->SetBranchAddress("BoxName",&g.BoxName);
	gtree->SetBranchAddress("BoxBarcode",&g.BoxBarcode);
	gtree->SetBranchAddress("PMTName",&g.PMTName);
	gtree->SetBranchAddress("PMTSerial",&g.PMTSerial);
	gtree->SetBranchAddress("l2014",&g.l2014);
	gtree->SetBranchAddress("il",&g.il);
	gtree->SetBranchAddress("ilV2",&g.ilV2);
	gtree->SetBranchAddress("ih",&g.ih);
	gtree->SetBranchAddress("ihV2",&g.ihV2);
	gtree->SetBranchAddress("fl",&g.fl);
	gtree->SetBranchAddress("fl2",&g.fl2);
	gtree->SetBranchAddress("flV2",&g.flV2);
	gtree->SetBranchAddress("fh",&g.fh);
	gtree->SetBranchAddress("fhV2",&g.fhV2);
	
	tg5[0]=new TGraphErrors();tg5[1]=new TGraphErrors();tg5[2]=new TGraphErrors();
	
// 	double XY5[3][2][24]={{{0.}}};
// 	double XYx[2][24]={{0.}};
	for(int i=0;i<gtree->GetEntries();i++)
	{
		gtree->GetEntry(i);
		if(*g.BoxBarcode==b[0])
		{
			for(int i1=0;i1<24;i1++)
			{
				if(*g.PMTName==ChNames[i1][0])
				{
// 					XY5[0][0][i1]=g.il[0]*1.6/10000.;XY5[0][1][i1]=g.il[1]*1.6/10000.;
// 					XY5[1][0][i1]=g.fl2[0][0]*1.6/10000.;XY5[1][1][i1]=g.fl2[0][1]*1.6/10000.;
// 					XY5[2][0][i1]=g.fl2[1][0]*1.6/10000.;XY5[2][1][i1]=g.fl2[1][1]*1.6/10000.;
// 					
// 					XYx[0][i1]=((double)(i1+1));
// 					XYx[1][i1]=0.;
					
					tg5[0]->SetPoint(i1,((double)(i1+1)),g.il[0]*1.6/10000.);
					tg5[0]->SetPointError(i1,0.,g.il[1]*1.6/10000.);
					tg5[1]->SetPoint(i1,((double)(i1+1)),g.fl2[0][0]*1.6/10000.);
					tg5[1]->SetPointError(i1,0.,g.fl2[0][1]*1.6/10000.);
					tg5[2]->SetPoint(i1,((double)(i1+1)),g.fl2[1][0]*1.6/10000.);
					tg5[2]->SetPointError(i1,0.,g.fl2[1][1]*1.6/10000.);
				}
			}
		}
	}
	gainroot->Close();
	
// 	for(int i1=0;i1<tg5[0]->GetN();i1++)
// 	{
// 		cout<<i1<<" "<<tg5[0]->GetX()[i1]<<" "<<tg5[0]->GetY()[i1]<<endl;
// 	}
	
// 	tg5[0]=new TGraphErrors(24,XYx[0],XY5[0][0],XYx[1],XY5[0][1]);
// 	tg5[1]=new TGraphErrors(24,XYx[0],XY5[1][0],XYx[1],XY5[1][1]);
// 	tg5[2]=new TGraphErrors(24,XYx[0],XY5[2][0],XYx[1],XY5[2][1]);
	
	TCanvas* cc2=new TCanvas("cc2","cc2",900,1200);
	gStyle->SetOptStat(0);
	cc2->Divide(1,4,0,0);
	cc2->cd(1);
	tg1[0]->Draw("AP");
	tg1[0]->GetYaxis()->SetRangeUser(tg1maxmin[1]-10.,tg1maxmin[0]+10.);
	tg1[0]->GetXaxis()->SetTitle("PMT ID");tg1[0]->GetXaxis()->CenterTitle();
	tg1[0]->GetYaxis()->SetTitle("<NPE>");tg1[0]->GetYaxis()->CenterTitle();
	gPad->SetGridy(1);gPad->SetGridx(1);
	tg1[0]->GetYaxis()->SetLabelSize(0.06);
	tg1[0]->GetYaxis()->SetTitleSize(0.05);
	tg1[1]->Draw("same P");
	cc2->cd(2);
	tg2->Draw("AP");
	tg2->GetXaxis()->SetTitle("PMT ID");tg2->GetXaxis()->CenterTitle();
	tg2->GetYaxis()->SetTitle("(Post-Pre)/Pre");tg2->GetYaxis()->CenterTitle();
	gPad->SetGridy(1);gPad->SetGridx(1);
	tg2->GetYaxis()->SetLabelSize(0.06);
	tg2->GetYaxis()->SetTitleSize(0.05);
	double tg2maxmin[2]={1.,-1.};
	dd=TMath::LocMax(tg2->GetN(),tg2->GetY());de=TMath::LocMin(tg2->GetN(),tg2->GetY());
	if(tg2->GetY()[dd]>tg2maxmin[0]){tg2maxmin[0]=tg2->GetY()[dd];}
	if(tg2->GetY()[de]<tg2maxmin[1]){tg2maxmin[1]=tg2->GetY()[de];}
	tg2->GetYaxis()->SetRangeUser(tg2maxmin[1],tg2maxmin[0]);
	cc2->cd(3);
	tg4->Draw("AP");
	tg4->GetXaxis()->SetTitle("PMT ID");tg4->GetXaxis()->CenterTitle();
// 	tg4->GetYaxis()->SetTitle("3+4/1+2");tg4->GetYaxis()->CenterTitle();
	tg4->GetYaxis()->SetTitle("((1+2)-(3+4))/((1+2)+(3+4))");tg4->GetYaxis()->CenterTitle();
	gPad->SetGridy(1);gPad->SetGridx(1);
	tg4->GetYaxis()->SetLabelSize(0.06);
	tg4->GetYaxis()->SetTitleSize(0.05);
	double tg4maxmin[2]={1.,-1.};
	dd=TMath::LocMax(tg4->GetN(),tg4->GetY());de=TMath::LocMin(tg4->GetN(),tg4->GetY());
	if(tg4->GetY()[dd]>tg4maxmin[0]){tg4maxmin[0]=tg4->GetY()[dd];}
	if(tg4->GetY()[de]<tg4maxmin[1]){tg4maxmin[1]=tg4->GetY()[de];}
	tg4->GetYaxis()->SetRangeUser(tg4maxmin[1],tg4maxmin[0]);
	cc2->cd(4);
	tg5[0]->SetMarkerStyle(24);tg5[0]->SetMarkerColor(1);tg5[0]->SetLineColor(1);
	tg5[1]->SetMarkerStyle(25);tg5[1]->SetMarkerColor(2);tg5[1]->SetLineColor(2);
	tg5[2]->SetMarkerStyle(26);tg5[2]->SetMarkerColor(4);tg5[2]->SetLineColor(4);
	tg5[0]->Draw("AP");tg5[1]->Draw("same P");tg5[2]->Draw("same P");
	tg5[0]->GetXaxis()->SetTitle("PMT ID");tg5[0]->GetXaxis()->CenterTitle();
	tg5[0]->GetYaxis()->SetTitle("SPE Charge (fC)");tg5[0]->GetYaxis()->CenterTitle();
	gPad->SetGridy(1);gPad->SetGridx(1);
	tg5[0]->GetYaxis()->SetLabelSize(0.06);
	tg5[0]->GetYaxis()->SetTitleSize(0.05);
	tg5[0]->GetXaxis()->SetLabelSize(0.06);
	tg5[0]->GetXaxis()->SetTitleSize(0.05);
	double tg5maxmin[2]={10.,50.};
	dd=TMath::LocMax(tg5[0]->GetN(),tg5[0]->GetY());de=TMath::LocMin(tg5[0]->GetN(),tg5[0]->GetY());
	if(tg5[0]->GetY()[dd]>tg5maxmin[0]){tg5maxmin[0]=tg5[0]->GetY()[dd];}
	if(tg5[0]->GetY()[de]<tg5maxmin[1]){tg5maxmin[1]=tg5[0]->GetY()[de];}
	dd=TMath::LocMax(tg5[1]->GetN(),tg5[1]->GetY());de=TMath::LocMin(tg5[1]->GetN(),tg5[1]->GetY());
	if(tg5[1]->GetY()[dd]>tg5maxmin[0]){tg5maxmin[0]=tg5[1]->GetY()[dd];}
	if(tg5[1]->GetY()[de]<tg5maxmin[1]){tg5maxmin[1]=tg5[1]->GetY()[de];}
	dd=TMath::LocMax(tg5[2]->GetN(),tg5[2]->GetY());de=TMath::LocMin(tg5[2]->GetN(),tg5[2]->GetY());
	if(tg5[2]->GetY()[dd]>tg5maxmin[0]){tg5maxmin[0]=tg5[2]->GetY()[dd];}
	if(tg5[2]->GetY()[de]<tg5maxmin[1]){tg5maxmin[1]=tg5[2]->GetY()[de];}
	tg5[0]->GetYaxis()->SetRangeUser(tg5maxmin[1]-5.,tg5maxmin[0]+5.);
	
	sprintf(hname,"IFComp_%s_%s_%d_%s_%d.pdf",b[0].c_str(),r[0].c_str(),s[0],r[1].c_str(),s[1]);
	cc2->Print(hname);
	
	inroot[0]->Close();inroot[1]->Close();
	sprintf(hname,"mv IFComp_*pdf ../Plots/Others");
	system(hname);
	
	delete cc2;
}

int plothistos(string r, string b, int id, int s)
{
	char hname[500];
	edata ed;
	sprintf(hname,"../NTuples/N_%s.root",r.c_str());
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
// 	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("histo",&ed.histo);
	
	sprintf(hname,"../Histos/%s/PedestalHistos_%d_%s.root",b.c_str(),id,r.c_str());
	TFile* outroot=new TFile(hname,"recreate");
	
	sprintf(hname,"SPEGains_%s_%d.txt",b.c_str(),id);
	ofstream outfile(hname);
// 	TF1* t2=new TF1("t2","gaus(0)+gaus(3)",-50.,36.);
	TF1* t2=new TF1("t2","gaus(0)+landau(3)",-50.,36.);
// 	TF1* t2=new TF1("t2","gaus(0)+gaus(3)+gaus(6)",-50.,50.);
// 	TF1* t2=new TF1("t2","gaus(0)+gaus(3)+gaus(6)+gaus(9)",-50.,50.);
	TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
	TCanvas* cc3=new TCanvas("cc3","cc3",4500,6000);
	TCanvas* cc2=new TCanvas("cc2","cc2",900,1200);
	gStyle->SetOptStat(0);
	
	int m=FindMode(b,r);
	
	if(m==1) {cc1->Divide(4,6,0,0);cc3->Divide(4,6,0,0);}
	else {cc1->Divide(8,6,0,0);cc3->Divide(8,6,0,0);}
	cc2->Divide(1,2);
	
	TGraphErrors* tg1[2];
	tg1[0]=new TGraphErrors();tg1[0]->SetName("SPE_Mean_1+2");
	tg1[0]->SetMarkerStyle(24);tg1[0]->SetMarkerColor(1);
	tg1[1]=new TGraphErrors();tg1[1]->SetName("SPE_Mean_3+4");
	tg1[1]->SetMarkerStyle(25);tg1[1]->SetMarkerColor(2);
	
	TGraphErrors* tg2[2];
	tg2[0]=new TGraphErrors();tg2[0]->SetName("SPE_Sigma_1+2");
	tg2[0]->SetMarkerStyle(24);tg2[0]->SetMarkerColor(1);
	tg2[1]=new TGraphErrors();tg2[1]->SetName("SPE_Sigma_3+4");
	tg2[1]->SetMarkerStyle(24);tg2[1]->SetMarkerColor(2);
	
	TH1F* QH[24][2];
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			sprintf(hname,"QH %s",ChNames[i1][m+i2-2].c_str());
			QH[i1][i2-1]=new TH1F(hname,hname,160,-40,600);
			QH[i1][i2-1]->GetXaxis()->SetTitle("Charge (fC)");
			QH[i1][i2-1]->GetXaxis()->CenterTitle();
			QH[i1][i2-1]->GetYaxis()->SetTitle("Entries / 4 fC");
			QH[i1][i2-1]->GetYaxis()->CenterTitle();
		}
	}
	
	vector <int> clist[3];
	tree->GetEntry(0);
	for(int iz1=0;iz1<144;iz1++)
	{
		for(int i1=0;i1<24;i1++)
		{
			for(int i2=1;i2<=m;i2++)
			{
				if(ed.ieta[iz1]==MAP2Ch[s][i1][i2-1][0] && ed.iphi[iz1]==MAP2Ch[s][i1][i2-1][1] && ed.depth[iz1]==MAP2Ch[s][i1][i2-1][2])
				{
					clist[0].push_back(iz1);
					clist[1].push_back(i1);
					clist[2].push_back(i2-1);
				}
			}
		}
	}
	
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		for(int i1=0;i1<clist[0].size();i1++)
		{
			for(int i2=0;i2<60;i2++)
			{
				QH[clist[1][i1]][clist[2][i1]]->Fill(adc2fC_QIE10[i2],ed.histo[clist[0][i1]][i2]);
// 				cout<<i2<<" "<<adc2fC_QIE10[i2]<<" "<<ed.histo[clist[0][i1]][i2]<<endl;
			}
		}
	}
	
	double ymax=0.;
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			if(QH[i1][i2-1]->GetBinContent(QH[i1][i2-1]->GetMaximumBin())>ymax){ymax=QH[i1][i2-1]->GetBinContent(QH[i1][i2-1]->GetMaximumBin());}
		}
	}
	
	outroot->cd();
	int tci1=1;
	double tg1maxmin[2]={0.,10000.};
	double tg2maxmin[2]={0.,10000.};
	double G[24][2]={{0.}};double Ge[24][2]={{0.}};
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=1;i2<=m;i2++)
		{
			cc1->cd(tci1);gPad-> SetLogy();
			QH[i1][i2-1]->Draw();
			QH[i1][i2-1]->GetXaxis()->SetRangeUser(-20.,100.);QH[i1][i2-1]->GetYaxis()->SetRangeUser(1.,ymax+5000.);
			cc1->Update();
			
			t2->SetParLimits(0,1e7,1e10);
			t2->SetParLimits(1,-10,10);
			t2->SetParLimits(2,2,10);
			t2->SetParLimits(3,1000,10000);
			t2->SetParLimits(4,10,30);
			t2->SetParLimits(5,2,10);
			
			t2->SetParameter(0,1e9);
			t2->SetParameter(1,0.);
			t2->SetParameter(2,2.5);
			t2->SetParameter(3,4000.);
			t2->SetParameter(4,20.);
			t2->SetParameter(5,4.);
			
			double df=28.;double dfopt=0.;double eerr=1000.;
			while(df<=36.)
			{
				QH[i1][i2-1]->Fit(t2,"q","q",-50.,df);
				if(t2->GetParError(4)<eerr){eerr=t2->GetParError(4);dfopt=df;}
				df+=4.;
			}
			
			t2->SetLineStyle(2);
// 				h[i1][i2]->Fit(t2,"q","q",-50.,35.);
// 					h[i1][i2]->Fit(t2,"q","q",-50.,24.);
			QH[i1][i2-1]->Fit(t2,"q","q",-50.,dfopt);
// 					if(t2->GetParError(4)<fp[5])
			
			int nn=tg1[i2-1]->GetN();
			tg1[i2-1]->SetPoint(nn,((double)(i1+1)),t2->GetParameter(4)-t2->GetParameter(1));
			tg1[i2-1]->SetPointError(nn,0,sqrt(pow(t2->GetParError(1),2.)+pow(t2->GetParError(4),2.)));
			tg2[i2-1]->SetPoint(nn,((double)(i1+1)),t2->GetParameter(5));
			tg2[i2-1]->SetPointError(nn,0,sqrt(pow(t2->GetParError(2),2.)+pow(t2->GetParError(5),2.)));
			
			if((t2->GetParameter(4)-t2->GetParameter(1))>tg1maxmin[0]){tg1maxmin[0]=(t2->GetParameter(4)-t2->GetParameter(1));}if((t2->GetParameter(4)-t2->GetParameter(1))<tg1maxmin[1]){tg1maxmin[1]=(t2->GetParameter(4)-t2->GetParameter(1));}
			if(t2->GetParameter(5)>tg2maxmin[0]){tg2maxmin[0]=t2->GetParameter(5);}if(t2->GetParameter(5)<tg2maxmin[1]){tg2maxmin[1]=t2->GetParameter(5);}
			
			G[i1][i2-1]=((t2->GetParameter(4)-t2->GetParameter(1))*10000./1.6);
			Ge[i1][i2-1]=(sqrt(pow(t2->GetParError(1),2.)+pow(t2->GetParError(4),2.))*10000./1.6);
			
			outfile<<ChNames[i1][m+i2-2]<<" "<<G[i1][i2-1]<<" "<<Ge[i1][i2-1]<<endl;
			
			tci1++;
			
			QH[i1][i2-1]->Write();
		}
	}
	
	sprintf(hname,"Histograms_%d.pdf(",id);
	cc1->Print(hname);
	cc2->cd(1);
	tg1[0]->Draw("AP");tg1[0]->GetYaxis()->SetTitle("Mean SPE Charge (fC)");tg1[0]->GetYaxis()->CenterTitle();tg1[0]->GetXaxis()->SetTitle("Channel");tg1[0]->GetXaxis()->CenterTitle();
	tg1[0]->GetYaxis()->SetRangeUser(tg1maxmin[1]-5.,tg1maxmin[0]+5.);gPad->SetGridy(1);gPad->SetGridx(1);
	cc2->Update();
	if(m==2) tg1[1]->Draw("same P");
	cc2->cd(2);
	tg2[0]->Draw("AP");tg2[0]->GetYaxis()->SetTitle("SPE Charge Sigma (fC)");tg2[0]->GetYaxis()->CenterTitle();tg2[0]->GetXaxis()->SetTitle("Channel");tg2[0]->GetXaxis()->CenterTitle();
	tg2[0]->GetYaxis()->SetRangeUser(tg2maxmin[1]-2.,tg2maxmin[0]+2.);gPad->SetGridy(1);gPad->SetGridx(1);
	cc2->Update();
	if(m==2) tg2[1]->Draw("same P");
	sprintf(hname,"Histograms_%d.pdf)",id);
	cc2->Print(hname);
	
	if(id==11){sprintf(hname,"mv Histograms_%d.pdf Histograms_%d_%s.pdf",id,id,r.c_str());system(hname);}
	sprintf(hname,"mv Histograms_*.pdf ../Plots/%s",b.c_str());
	system(hname);
	
	outfile.close();
	if(id==11){sprintf(hname,"mv SPEGains_%s_%d.txt SPEGains_%s_%d_%s.txt",b.c_str(),id,b.c_str(),id,r.c_str());system(hname);}
	sprintf(hname,"mv SPEGains_*.txt ../Plots/%s",b.c_str());
	system(hname);
	
	inroot->Close();
	
	outroot->cd();
	tg1[0]->Write();if(m==2)tg1[1]->Write();
	tg2[0]->Write();if(m==2)tg2[1]->Write();
	outroot->Close();
	
// 	gain g;
	TFile* gainroot=new TFile("../Gains/Gains.root");
	TTree *gtree = (TTree*)gainroot->Get("Gains");
	gtree->SetBranchAddress("BoxName",&g.BoxName);
	gtree->SetBranchAddress("BoxBarcode",&g.BoxBarcode);
	gtree->SetBranchAddress("PMTName",&g.PMTName);
	gtree->SetBranchAddress("PMTSerial",&g.PMTSerial);
	gtree->SetBranchAddress("l2014",&g.l2014);
	gtree->SetBranchAddress("il",&g.il);
	gtree->SetBranchAddress("ilV2",&g.ilV2);
	gtree->SetBranchAddress("ih",&g.ih);
	gtree->SetBranchAddress("ihV2",&g.ihV2);
	gtree->SetBranchAddress("fl",&g.fl);
	gtree->SetBranchAddress("fl2",&g.fl2);
	gtree->SetBranchAddress("flV2",&g.flV2);
	gtree->SetBranchAddress("fh",&g.fh);
	gtree->SetBranchAddress("fhV2",&g.fhV2);
	
	TFile* goutroot=new TFile("../Gains/Gains_out.root","recreate");
	TTree* treeout = new TTree("Gains", "Gains");
	treeout->Branch("BoxName", &g.BoxName);
	treeout->Branch("BoxBarcode", &g.BoxBarcode);
	treeout->Branch("PMTName", &g.PMTName);
	treeout->Branch("PMTSerial", &g.PMTSerial);
	treeout->Branch("l2014", g.l2014, "l2014[2]/F");
	treeout->Branch("il", g.il, "il[2]/F");
	treeout->Branch("ilV2", g.ilV2, "ilV2[2]/F");
	treeout->Branch("ih", g.ih, "ih[2]/F");
	treeout->Branch("ihV2", g.ihV2, "ihV2[2]/F");
	treeout->Branch("fl", g.fl, "fl[2][2]/F");
	treeout->Branch("fl2", g.fl2, "fl2[2][2]/F");
	treeout->Branch("flV2", g.flV2, "flV2[2][2]/F");
	treeout->Branch("fh", g.fh, "fh[2][2]/F");
	treeout->Branch("fhV2", g.fhV2, "fhV2[2][2]/F");
	
	for(int i=0;i<gtree->GetEntries();i++)
	{
		gtree->GetEntry(i);
		if(*g.BoxBarcode==b)
		{
			for(int i1=0;i1<24;i1++)
			{
				if(*g.PMTName==ChNames[i1][0])
				{
					if(id==12){g.ih[0]=G[i1][0];g.ih[1]=Ge[i1][0];}
					else if(id==13){g.ihV2[0]=G[i1][0];g.ihV2[1]=Ge[i1][0];}
					else if(id==15){g.fh[0][0]=G[i1][0];g.fh[0][1]=Ge[i1][0];g.fh[1][0]=G[i1][1];g.fh[1][1]=Ge[i1][1];}
					else if(id==16){g.fhV2[0][0]=G[i1][0];g.fhV2[0][1]=Ge[i1][0];g.fhV2[1][0]=G[i1][1];g.fhV2[1][1]=Ge[i1][1];}
				}
			}
		}
		treeout->Fill();
	}
	goutroot->cd();
	goutroot->Write();
	goutroot->Close();
	gainroot->Close();
	system("mv ../Gains/Gains_out.root ../Gains/Gains.root");
	
	delete cc1;delete cc2;delete cc3;
}

int finalplots(string b)
{
	vector <finallist> FL;
	finallist fl;
	char hname[300];
	for(int i2=0;i2<RL.size();i2++)
	{
		for(int i3=0;i3<3;i3++)
		{
			if(RL[i2].b[i3]==b)
			{
				fl.runno=RL[i2].runno;
				fl.b=RL[i2].b[i3];
				fl.t=RL[i2].t[i3];
				fl.s=i3;
				FL.push_back(fl);
			}
		}
	}
	sprintf(hname,"../Histos/%s/Finals_%s.root",b.c_str(),b.c_str());
	TFile* outroot=new TFile(hname,"recreate");
	TFile* ff[3];int ss[3]={0};
	for(int i1=0;i1<FL.size();i1++)
	{
		if(FL[i1].t==2)//initial LED
		{
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",FL[i1].b.c_str(),FL[i1].t,FL[i1].runno.c_str());
			ff[0]=new TFile(hname);
			ss[0]=FL[i1].s;
			cout<<hname<<" Initial LED"<<endl; 
		}
		if(FL[i1].t==4)//final LED 1
		{
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",FL[i1].b.c_str(),FL[i1].t,FL[i1].runno.c_str());
			ff[1]=new TFile(hname);
			ss[1]=FL[i1].s;
			cout<<hname<<" Final LED 1"<<endl;
		}
		if(FL[i1].t==6)//final LED 2
		{
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",FL[i1].b.c_str(),FL[i1].t,FL[i1].runno.c_str());
			ff[2]=new TFile(hname);
			ss[2]=FL[i1].s;
			cout<<hname<<" Final LED 2"<<endl;
		}
	}
	TCanvas* cc=new TCanvas("cc","cc",500,500);
	gStyle->SetOptStat(0);
	TH1F* hNPE[24][3];TH1F* hFrac[24][2];
	int colors[3]={1,2,4};
	double npexymax[2]={0.};
	for(int i1=0;i1<24;i1++)
	{
		for(int i2=0;i2<3;i2++)
		{
			sprintf(hname,"NPETot %s",ChNames[i1][0].c_str());
			ff[i2]->GetObject(hname,hNPE[i1][i2]);hNPE[i1][i2]->SetLineColor(colors[i2]);
			if(i2!=1)
			{
				if(hNPE[i1][i2]->GetBinCenter(hNPE[i1][i2]->FindLastBinAbove(0.))>npexymax[0]){npexymax[0]=hNPE[i1][i2]->GetBinCenter(hNPE[i1][i2]->FindLastBinAbove(0.));}
				if(hNPE[i1][i2]->GetBinContent(hNPE[i1][i2]->GetMaximumBin())>npexymax[1]){npexymax[1]=hNPE[i1][i2]->GetBinContent(hNPE[i1][i2]->GetMaximumBin());}
			}
			if(i2>0)
			{
				sprintf(hname,"Frac %s",ChNames[i1][0].c_str());
				ff[i2]->GetObject(hname,hFrac[i1][i2-1]);hFrac[i1][i2-1]->SetLineColor(colors[i2]);
			}
		}
	}
	outroot->cd();
	for(int i1=0;i1<24;i1++)
	{
		hNPE[i1][0]->Draw();
// 		hNPE[i1][1]->Draw("same");
		hNPE[i1][2]->Draw("same");
		hNPE[i1][0]->GetXaxis()->SetRangeUser(0.,npexymax[0]+10.);
		hNPE[i1][0]->GetYaxis()->SetRangeUser(0.,npexymax[1]+200.);
		hNPE[i1][0]->GetYaxis()->SetTitleOffset(1.5);
		sprintf(hname,"LEDComp_%d.png",i1);
		cc->SaveAs(hname);
// 		hFrac[i1][0]->Draw();hFrac[i1][1]->Draw("same");
		hFrac[i1][1]->Draw();
// 		hFrac[i1][0]->GetYaxis()->SetTitleOffset(1.5);
		hFrac[i1][1]->GetYaxis()->SetTitleOffset(1.5);
		sprintf(hname,"LEDfrac_%d.png",i1);
		cc->SaveAs(hname);
		sprintf(hname,"%s Initial LED",hNPE[i1][0]->GetName());hNPE[i1][0]->SetName(hname);
		sprintf(hname,"%s Final LED 1",hNPE[i1][1]->GetName());hNPE[i1][1]->SetName(hname);
		sprintf(hname,"%s Final LED 2",hNPE[i1][2]->GetName());hNPE[i1][2]->SetName(hname);
		hNPE[i1][0]->Write();hNPE[i1][1]->Write();hNPE[i1][2]->Write();
		sprintf(hname,"%s Final LED 1",hFrac[i1][0]->GetName());hFrac[i1][0]->SetName(hname);
		sprintf(hname,"%s Final LED 2",hFrac[i1][1]->GetName());hFrac[i1][1]->SetName(hname);
		hFrac[i1][0]->Write();hFrac[i1][1]->Write();
	}
	ff[0]->Close();ff[1]->Close();ff[2]->Close();
	
	TFile* ov[50];int ssov[50]={0};int nov=0;
	for(int i1=0;i1<FL.size();i1++)
	{
		if(FL[i1].t==7)//overnight LEDs
		{
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",FL[i1].b.c_str(),FL[i1].t,FL[i1].runno.c_str());
			ov[nov]=new TFile(hname);
			ss[nov]=FL[i1].s;
			nov++;
			cout<<hname<<" Overnight LED "<<nov<<endl;
		}
	}
	TH2F* ovH[24];TH2F* ovP[24];
	TH1F* hovH;
	for(int i1=0;i1<24;i1++)
	{
		sprintf(hname,"ADCperTS %s",ChNames[i1][0].c_str());
		ovP[i1]=new TH2F(hname,hname,nov,-0.5,nov-0.5,556,-0.5,555.5);
		ovP[i1]->GetXaxis()->SetTitle("Run ID");ovP[i1]->GetXaxis()->CenterTitle();
		ovP[i1]->GetYaxis()->SetTitle("1+2                           3+4");ovP[i1]->GetYaxis()->CenterTitle();ovP[i1]->GetYaxis()->SetTitleOffset(1.5);
		for(int i3=0;i3<nov;i3++)
		{
			sprintf(hname,"ADCPH %s",ChNames[i1][1].c_str());
			ov[i3]->GetObject(hname,hovH);
			for(int i4=1;i4<=hovH->GetNbinsX();i4++)
			{
				ovP[i1]->Fill(i3,hovH->GetBinCenter(i4),hovH->GetBinContent(i4));
			}
			sprintf(hname,"ADCPH %s",ChNames[i1][2].c_str());
			ov[i3]->GetObject(hname,hovH);
			for(int i4=1;i4<=hovH->GetNbinsX();i4++)
			{
				ovP[i1]->Fill(i3,hovH->GetBinCenter(i4)+300,hovH->GetBinContent(i4));
			}
		}
		gStyle->SetPalette(kRainBow);
		ovP[i1]->Draw("colz");
		sprintf(hname,"OvernightPulses_%d.png",i1);
		cc->SaveAs(hname);
		
		sprintf(hname,"OvernightNPEs %s",ChNames[i1][0].c_str());
		ovH[i1]=new TH2F(hname,hname,nov,-0.5,nov-0.5,1100,0,1100);
		ovH[i1]->GetXaxis()->SetTitle("Run ID");ovH[i1]->GetXaxis()->CenterTitle();
		ovH[i1]->GetYaxis()->SetTitle("1+2                           3+4");ovH[i1]->GetYaxis()->CenterTitle();ovH[i1]->GetYaxis()->SetTitleOffset(1.5);
		for(int i3=0;i3<nov;i3++)
		{
			sprintf(hname,"NPETot %s",ChNames[i1][1].c_str());
			ov[i3]->GetObject(hname,hovH);
			for(int i4=1;i4<=hovH->GetNbinsX();i4++)
			{
				ovH[i1]->Fill(i3,hovH->GetBinCenter(i4),hovH->GetBinContent(i4));
			}
			sprintf(hname,"NPETot %s",ChNames[i1][2].c_str());
			ov[i3]->GetObject(hname,hovH);
			for(int i4=1;i4<=hovH->GetNbinsX();i4++)
			{
				ovH[i1]->Fill(i3,hovH->GetBinCenter(i4)+600,hovH->GetBinContent(i4));
			}
		}
		gStyle->SetPalette(kRainBow);
		ovH[i1]->Draw("colz");
		sprintf(hname,"OvernightNPEs_%d.png",i1);
		cc->SaveAs(hname);
		outroot->cd();
		ovH[i1]->Write();
		ovP[i1]->Write();
	}
	
	for(int i3=0;i3<nov;i3++)
	{
		ov[i3]->Close();
	}
	sprintf(hname,"mv *.png ../Plots/%s",FL[0].b.c_str());
	system(hname);
	
	vector <string> rif;vector <int> sif;
	int sif1[2]={0};string rif1[2]={"",""};
	string bif="";
	for(int i1=0;i1<FL.size();i1++)
	{
		if(FL[i1].t==2)//initial LED
		{
			sif1[0]=FL[i1].s+1;rif1[0]=FL[i1].runno;
			bif=FL[i1].b;
		}
		if(FL[i1].t==6)//final LED 2
		{
			sif1[1]=FL[i1].s+1;rif1[1]=FL[i1].runno;
		}
	}
	rif.push_back(rif1[0]);rif.push_back(rif1[1]);
	sif.push_back(sif1[0]);sif.push_back(sif1[1]);
	ifcomp(rif,sif);
	char ifname[2000];
	sprintf(ifname,"mv ../Plots/Others/IFComp_%s_%s_%d_%s_%d.pdf ../Plots/%s/IFCompare.pdf",bif.c_str(),rif[0].c_str(),sif[0],rif[1].c_str(),sif[1],bif.c_str());
	system(ifname);
}

int npediff()
{
	vector <npelist> NPEL[2];
	vector <npelist> NPEL2ch;
	npelist npel;
	int dl[4]={2,14,6,17};
	char hname[300];
	for(int i2=0;i2<RL.size();i2++)
	{
		for(int i3=0;i3<3;i3++)
		{
			for(int i5=0;i5<2;i5++)
			{
				if(RL[i2].t[i3]==dl[2*i5] || RL[i2].t[i3]==dl[2*i5+1])
				{
					int ned=-1;
					for(int i4=0;i4<NPEL[i5].size();i4++)
					{
						if(NPEL[i5][i4].b==RL[i2].b[i3]){ned=i4;}
					}
					if(ned!=-1)
					{
						NPEL[i5][ned].runno[1]=RL[i2].runno;
						NPEL[i5][ned].t[1]=RL[i2].t[i3];
						NPEL[i5][ned].s[1]=i3+1;
					}
					else
					{
						npel.b=RL[i2].b[i3];
						npel.runno[0]=RL[i2].runno;
						npel.t[0]=RL[i2].t[i3];
						npel.s[0]=i3+1;
						NPEL[i5].push_back(npel);
					}
				}
			}
		}
	}
	TFile* inroot[2];TGraphErrors* tg[2];
	TFile* outroot=new TFile("../Histos/Others/NPEDiff.root","recreate");
	TH1F* NPEDiff[3];
	NPEDiff[0]=new TH1F("NPEDiff1Ch","NPEDiff1Ch",50,-10.,10.);
	NPEDiff[1]=new TH1F("NPEDiff2Ch","NPEDiff2Ch",50,-10.,10.);
	NPEDiff[2]=new TH1F("NPEDiffAll","NPEDiffAll",50,-10.,10.);
	for(int i1=0;i1<3;i1++)
	{
		NPEDiff[i1]->GetXaxis()->SetTitle("NPE_{OV2+100}-NPE_{OV2}");NPEDiff[i1]->GetXaxis()->CenterTitle();
		NPEDiff[i1]->GetYaxis()->SetTitle("Entries / 0.4");NPEDiff[i1]->GetYaxis()->CenterTitle();
	}
	TH2F* NPEDiff_vs_NPE=new TH2F("NPEDiff_vs_NPE","NPEDiff_vs_NPE",200,0.,200.,50,-10.,10.);
	NPEDiff_vs_NPE->GetXaxis()->SetTitle("NPE");NPEDiff_vs_NPE->GetXaxis()->CenterTitle();
	NPEDiff_vs_NPE->GetYaxis()->SetTitle("NPE_{OV2+100}-NPE_{OV2}");NPEDiff_vs_NPE->GetYaxis()->CenterTitle();
	for(int i2=0;i2<2;i2++)
	{
		for(int i1=0;i1<NPEL[i2].size();i1++)
		{
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",NPEL[i2][i1].b.c_str(),NPEL[i2][i1].t[0],NPEL[i2][i1].runno[0].c_str());
			inroot[0]=new TFile(hname);
			sprintf(hname,"../Histos/%s/LEDs_%d_%s.root",NPEL[i2][i1].b.c_str(),NPEL[i2][i1].t[1],NPEL[i2][i1].runno[1].c_str());
			inroot[1]=new TFile(hname);
			inroot[0]->GetObject("LEDNPE_Mean_1+2",tg[0]);
			inroot[1]->GetObject("LEDNPE_Mean_1+2",tg[1]);
			for(int i3=0;i3<tg[0]->GetN();i3++)
			{
				NPEDiff[i2]->Fill(tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
				NPEDiff[2]->Fill(tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
				NPEDiff_vs_NPE->Fill(tg[0]->GetY()[i3],tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
			}
			if(i2==1)
			{
				inroot[0]->GetObject("LEDNPE_Mean_3+4",tg[0]);
				inroot[1]->GetObject("LEDNPE_Mean_3+4",tg[1]);
				for(int i3=0;i3<tg[0]->GetN();i3++)
				{
					NPEDiff[i2]->Fill(tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
					NPEDiff[2]->Fill(tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
					NPEDiff_vs_NPE->Fill(tg[0]->GetY()[i3],tg[1]->GetY()[i3]-tg[0]->GetY()[i3]);
				}
			}
			inroot[0]->Close();
			inroot[1]->Close();
		}
	}
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		NPEDiff[i1]->Fit("gaus");
// 	}
	outroot->Write();
	outroot->Close();
	
// 	for(int i2=0;i2<2;i2++)
// 	{
// 		for(int i1=0;i1<NPEL[i2].size();i1++)
// 		{
// 			cout<<i1<<" "<<NPEL[i2][i1].b<<" "<<NPEL[i2][i1].runno[0]<<" "<<NPEL[i2][i1].t[0]<<" "<<NPEL[i2][i1].s[0]<<" "<<NPEL[i2][i1].runno[1]<<" "<<NPEL[i2][i1].t[1]<<" "<<NPEL[i2][i1].s[1]<<endl;
// 		}
// 	}
	
}

int insertPreSignoff(vector <string> r, vector <int> s)
{
	string b=ledcomp(r,s);
	cout<<b<<endl;
	char hname[1000];
	sprintf(hname,"mv ../Plots/Others/LEDs_%s_%d_%s_%d.pdf ../Plots/%s/PreSignOff.pdf",r[0].c_str(),s[0],r[1].c_str(),s[1],b.c_str());
	system(hname);
	sprintf(hname,"../Plots/%s/LEDGains_%s_10_%s.txt",b.c_str(),b.c_str(),r[1].c_str());
	ifstream infile(hname);
	calib cn;double G[24][3]={{0.}};double Ge[24][3]={{0.}};
	while(!infile.eof())
	{
		infile>>cn.pmt>>cn.g>>cn.e;
		for(int i1=0;i1<24;i1++)
		{
			for(int i2=0;i2<3;i2++)
			{
				if(cn.pmt==ChNames[i1][i2]){G[i1][i2]=cn.g;Ge[i1][i2]=cn.e;}
			}
		}
	}
	
	TFile* gainroot=new TFile("../Gains/Gains.root");
	TTree *gtree = (TTree*)gainroot->Get("Gains");
	gtree->SetBranchAddress("BoxName",&g.BoxName);
	gtree->SetBranchAddress("BoxBarcode",&g.BoxBarcode);
	gtree->SetBranchAddress("PMTName",&g.PMTName);
	gtree->SetBranchAddress("PMTSerial",&g.PMTSerial);
	gtree->SetBranchAddress("l2014",&g.l2014);
	gtree->SetBranchAddress("il",&g.il);
	gtree->SetBranchAddress("ilV2",&g.ilV2);
	gtree->SetBranchAddress("ih",&g.ih);
	gtree->SetBranchAddress("ihV2",&g.ihV2);
	gtree->SetBranchAddress("fl",&g.fl);
	gtree->SetBranchAddress("fl2",&g.fl2);
	gtree->SetBranchAddress("flV2",&g.flV2);
	gtree->SetBranchAddress("fh",&g.fh);
	gtree->SetBranchAddress("fhV2",&g.fhV2);
	
	TFile* goutroot=new TFile("../Gains/Gains_out.root","recreate");
	TTree* treeout = new TTree("Gains", "Gains");
	treeout->Branch("BoxName", &g.BoxName);
	treeout->Branch("BoxBarcode", &g.BoxBarcode);
	treeout->Branch("PMTName", &g.PMTName);
	treeout->Branch("PMTSerial", &g.PMTSerial);
	treeout->Branch("l2014", g.l2014, "l2014[2]/F");
	treeout->Branch("il", g.il, "il[2]/F");
	treeout->Branch("ilV2", g.ilV2, "ilV2[2]/F");
	treeout->Branch("ih", g.ih, "ih[2]/F");
	treeout->Branch("ihV2", g.ihV2, "ihV2[2]/F");
	treeout->Branch("fl", g.fl, "fl[2][2]/F");
	treeout->Branch("fl2", g.fl2, "fl2[2][2]/F");
	treeout->Branch("flV2", g.flV2, "flV2[2][2]/F");
	treeout->Branch("fh", g.fh, "fh[2][2]/F");
	treeout->Branch("fhV2", g.fhV2, "fhV2[2][2]/F");
	
	for(int i=0;i<gtree->GetEntries();i++)
	{
		gtree->GetEntry(i);
		if(*g.BoxBarcode==b)
		{
			for(int i1=0;i1<24;i1++)
			{
				if(*g.PMTName==ChNames[i1][0])
				{
					cout<<ChNames[i1][0]<<" "<<g.fl2[0][0]<<" "<<G[i1][1]<<" "<<g.fl2[1][0]<<" "<<G[i1][2]<<endl;
					g.fl2[0][0]=G[i1][1];g.fl2[0][1]=Ge[i1][1];g.fl2[1][0]=G[i1][2];g.fl2[1][1]=Ge[i1][2];
				}
			}
		}
		treeout->Fill();
	}
	goutroot->cd();
	goutroot->Write();
	goutroot->Close();
	gainroot->Close();
	system("mv ../Gains/Gains_out.root ../Gains/Gains.root");
}

int main(int argc, char *argv[])
{
	int opt=atoi(argv[1]);
	
	if(opt==1)//plot the newly added runs
	{
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		infile.close();
		infile.open("RunDiff.txt");
		vector <string> rd;string r;
		while(!infile.eof())
		{
			infile>>r;
			rd.push_back(r);
		}
		if(rd.size()>1) {rd.pop_back();}
		infile.close();
		for(int i1=0;i1<rd.size();i1++)
		{
			for(int i2=0;i2<RL.size();i2++)
			{
				if(rd[i1]==RL[i2].runno)
				{
					cout<<"Run: "<<RL[i2].runno<<endl;
					for(int i3=0;i3<3;i3++)
					{
						if(RL[i2].t[i3]==1 || RL[i2].t[i3]==3 || RL[i2].t[i3]==5 || RL[i2].t[i3]==9){plotpeds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
						if(RL[i2].t[i3]==2 || RL[i2].t[i3]==4 || RL[i2].t[i3]==6 || RL[i2].t[i3]==7 || RL[i2].t[i3]==8 ||RL[i2].t[i3]==14 || RL[i2].t[i3]==17 || RL[i2].t[i3]==10){plotleds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
						if(RL[i2].t[i3]==12 || RL[i2].t[i3]==13 || RL[i2].t[i3]==15 || RL[i2].t[i3]==16 || RL[i2].t[i3]==11){plothistos(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
						if(RL[i2].t[i3]==17){finalplots(RL[i2].b[i3]);}
					}
				}
			}
		}
	}
	else if(opt==2)//single run
	{
		string RunNo=argv[2];
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		for(int i2=0;i2<RL.size();i2++)
		{
			if(RunNo==RL[i2].runno)
			{
				cout<<"Run: "<<RL[i2].runno<<endl;
				for(int i3=0;i3<3;i3++)
				{
					if(RL[i2].t[i3]==1 || RL[i2].t[i3]==3 || RL[i2].t[i3]==5 || RL[i2].t[i3]==9){plotpeds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
					if(RL[i2].t[i3]==2 || RL[i2].t[i3]==4 || RL[i2].t[i3]==6 || RL[i2].t[i3]==7 || RL[i2].t[i3]==8 ||RL[i2].t[i3]==14 || RL[i2].t[i3]==17 || RL[i2].t[i3]==10){plotleds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
					if(RL[i2].t[i3]==12 || RL[i2].t[i3]==13 || RL[i2].t[i3]==15 || RL[i2].t[i3]==16 || RL[i2].t[i3]==11){plothistos(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
					if(RL[i2].t[i3]==17){finalplots(RL[i2].b[i3]);}
				}
			}
		}
	}
	else if(opt==3)//compare up to 3 runs
	{
		vector <string> RunNos;
		vector <int> stations;
		RunNos.push_back(argv[2]);stations.push_back(atoi(argv[3]));
		if(argc>4){RunNos.push_back(argv[4]);stations.push_back(atoi(argv[5]));}
		if(argc>6){RunNos.push_back(argv[6]);stations.push_back(atoi(argv[7]));}
		ledcomp(RunNos,stations);
	}
	else if(opt==4)//plot everything
	{
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		infile.close();
		for(int i2=0;i2<RL.size();i2++)
		{
			cout<<"Run: "<<RL[i2].runno<<endl;
			for(int i3=0;i3<3;i3++)
			{
				if(RL[i2].t[i3]==1 || RL[i2].t[i3]==3 || RL[i2].t[i3]==5 || RL[i2].t[i3]==9){plotpeds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(RL[i2].t[i3]==2 || RL[i2].t[i3]==4 || RL[i2].t[i3]==6 || RL[i2].t[i3]==7 || RL[i2].t[i3]==8 ||RL[i2].t[i3]==14 || RL[i2].t[i3]==17 || RL[i2].t[i3]==10){plotleds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(RL[i2].t[i3]==12 || RL[i2].t[i3]==13 || RL[i2].t[i3]==15 || RL[i2].t[i3]==16 || RL[i2].t[i3]==11){plothistos(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(RL[i2].t[i3]==17){finalplots(RL[i2].b[i3]);}
			}
		}
	}
	else if(opt==5)//compare initial and final
	{
		vector <string> RunNos;
		vector <int> stations;
		RunNos.push_back(argv[2]);stations.push_back(atoi(argv[3]));
		if(argc>4){RunNos.push_back(argv[4]);stations.push_back(atoi(argv[5]));}
		ifcomp(RunNos,stations);
	}
	else if(opt==6)//print modes
	{
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		infile.close();
		for(int i1=0;i1<RL.size();i1++)
		{
			cout<<RL[i1].runno<<" "<<RL[i1].b[0]<<" "<<FindMode(RL[i1].b[0],RL[i1].runno)<<" "<<RL[i1].b[1]<<" "<<FindMode(RL[i1].b[1],RL[i1].runno)<<" "<<RL[i1].b[2]<<" "<<FindMode(RL[i1].b[2],RL[i1].runno)<<endl;
		}
	}
	else if(opt==7)//NPE diff for OV2 and OV2+100
	{
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		infile.close();
		npediff();
	}
	else if(opt==8)//plot specific types // 0 ped, 1 led (not overnight), 2 histos, 3 overnight, 4 finals
	{
		int opt2=atoi(argv[2]);
		ifstream infile("RunList.txt");
		while(!infile.eof())
		{
			infile>>rl.runno>>rl.b[0]>>rl.t[0]>>rl.b[1]>>rl.t[1]>>rl.b[2]>>rl.t[2]>>rl.date>>rl.time;
			RL.push_back(rl);
		}
		RL.pop_back();
		infile.close();
		for(int i2=0;i2<RL.size();i2++)
		{
			cout<<"Run: "<<RL[i2].runno<<endl;
			for(int i3=0;i3<3;i3++)
			{
				if(opt2==0 && (RL[i2].t[i3]==1 || RL[i2].t[i3]==3 || RL[i2].t[i3]==5 || RL[i2].t[i3]==9)){plotpeds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(opt2==1 && (RL[i2].t[i3]==2 || RL[i2].t[i3]==4 || RL[i2].t[i3]==6 ||RL[i2].t[i3]==14 || RL[i2].t[i3]==17 || RL[i2].t[i3]==10)){plotleds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(opt2==2 && (RL[i2].t[i3]==12 || RL[i2].t[i3]==13 || RL[i2].t[i3]==15 || RL[i2].t[i3]==16 || RL[i2].t[i3]==11)){plothistos(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(opt2==3 && (RL[i2].t[i3]==7 || RL[i2].t[i3]==8)){plotleds(RL[i2].runno,RL[i2].b[i3],RL[i2].t[i3],i3);}
				if(opt2==4 && RL[i2].t[i3]==17){finalplots(RL[i2].b[i3]);}
			}
		}
	}
	else if(opt==9)//replace fl2 with pre sign-off
	{
		vector <string> RunNos;
		vector <int> stations;
		RunNos.push_back(argv[2]);stations.push_back(atoi(argv[3]));
		if(argc>4){RunNos.push_back(argv[4]);stations.push_back(atoi(argv[5]));}
		insertPreSignoff(RunNos,stations);
	}
}

















