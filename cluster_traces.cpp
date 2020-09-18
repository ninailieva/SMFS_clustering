// To compile the code:  module purge
//                       module load intel/18.0.3
// 			 icc -std=c++0x -O3 cluster_traces.cpp -fopenmp -o cluster_traces.x
// To run the code:      ./cluster_traces.x input.txt > output.txt 
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
using namespace std;
string filename;                        // Stores the name of the input file read from the command line.  
int    n_cont=20;			// Minimum number of points with force>0; used to detect the starting point in the contact part;
double sigma_noise=5.67;		// Standard deviation of the noise;
double sigma_cut=4.*sigma_noise; 	// Threshold used for non-contact part removal;
double sigma_tail=2.*sigma_noise; 	// Threshold used to detect spurious traces (traces having strange tails);
double len_min=50.; 			// Threshold on the trace length used to discard short traces [nm];
double abnormal_f=5000.;		// Force values above abnormal_f are considered abnormal;
double abnormal_x=5000.;		// Extension values above abnormal_x are considered abnormal;
double dx_interp=1.0; 			// Grid width used in trace interpolation [nm];
double l_p=0.4; 			// Persistence length [nm];
double F_min=30.;      			// Lower bound limit of the WLC force range;
double F_max=500.;    			// Upper bound limit of the WLC force range;
double lim_p_max=0.01;			// Peaks in the Lc histogram above this value are considered relevant; if below they get peak score 0;
double peak_width=75.;		        // Used in score assignation: score is assigned to points within this peak width;
double bin_width=8.; 			// Lc histogram bin width;
double score_thresh=0.5; 		// the threshold on the ratio between global score and trace length;
double F_scoring=4.*sigma_noise; 	// Match/mismatch threshold used to compute the similarity score in the distance;
double peaks_diff=10000.0; 			// Distances are computed only between traces that differ in peaks by no more than this number;
int Npeaks_min=1;                       //traces with less than Npeaks_min are discarded
int Npeaks_max=100000;                  //traces with more than Npeaks_maxn are discarded
double len_diff=0.2;                    // Distances are computed only between traces that differ in length by no more than this percentage; 
int    kNN=6; 				// Number of nearest neighbours in the density estimator;
double r_cut=0.3;                       // Cutoff distance used in the clustering to determine the cluster core;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int debug=0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double delta1 = 0.002*dx_interp;        // dynamic programming gap penalty in the first 10 nm;
double delta2 = 0.8*dx_interp;          // dynamic programming gap penalty after 10 nm;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double kbT=4.11433259; 			// Boltzmann constant*Temperature for T=298K [pN*nm]
int    NCEN=0;				// Total number of clusters;
string trace_input;
// Print output variables.
int print_f=0;				// Print the original F-x curves in format: "f" tr_no x F(x)
int print_f_interp=0;                   // Print the interpolated F-x curves in format: "fi" tr_no x F(x)    
int print_Lc=0;				// Print Lc values for each point in the F-x curves in format: "Lc" tr_no x F(x) Lc
int print_hist=0;			// Print Lc histograms in format: "h" tr_no bincenter probability
int print_w=0;				// Print WLC score for each point in the F-x curves in format: "w" tr_no x F(x) w
int print_globalscore=0;		// Print the global score and trace length for each trace in format: "gl" tr_no score length
int print_dist=0;			// Print the distance matrix in format: "dist" i j distance(i,j)
int print_trlen=0;			// Print the length of each trace in format: "l" tr_no length
int print_stdout=0;			// Print the reason for a trace to be discarded in format: "rejected" tr_no <reason explained>
int print_knndist=0;			// Prints the distances between each traces and its k nearest neighbors in format: "knndist" i k distance(i, k)
int print_align=0;			// Prints trace pairs aligned
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
struct clust_s{
  int size=0;
  int center=0;
  vector <int> members;
  vector <double> distcenter;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
struct bord_s{
  int isborder=0;
  int cl1;
  int cl2;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct dp_s{
  int i_ini;
  int cen;
  int numb;
  int nn;   //number of neighbors within rcut
  double rho;
  double delta;
  double gamma;
  vector <int> ind_nn;
  vector <double> dist_nn;
  double dist_center;
  double dist_other;
  int cl=-1;
  int core=0;
  int nn_hrho=-1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool mt(const dp_s &a, const dp_s  &b)
{
  return a.rho > b.rho;
};
bool ltgamma(const dp_s &a, const dp_s  &b)
{
  return a.gamma > b.gamma;
};
bool index(const dp_s &a, const dp_s  &b)
{
  return a.i_ini < b.i_ini;
};
bool ltsize(const clust_s &a, const clust_s  &b)
{
  return a.size> b.size;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct dist_s{
  double dist;
  int ind;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool lt(const dist_s &a, const dist_s  &b)
{
  return a.dist < b.dist;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The read_input function reads all input parameters from input.txt together with various options for printing more detailed output information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_input(string& filename, string& trace_input, int& n_cont, double& sigma_noise, double& sigma_cut, double& sigma_tail, double& len_min, double& abnormal_f, double& abnormal_x, double& dx_interp, double& l_p, double& F_min, double& F_max, double& lim_p_max, double& peak_width, double& bin_width, double& score_thresh, double& delta1, double& delta2, double& F_scoring, double& peaks_diff, double& len_diff, int& Npeaks_min, int& Npeaks_max, int& kNN, double& r_cut, int& print_f, int& print_f_interp, int& print_Lc, int& print_hist, int&  print_w, int& print_globalscore, int& print_dist, int& print_trlen, int& print_stdout,int& print_knndist, int& print_align, int& debug ){
  //string filename="input.txt";
  ifstream fileIN;
  fileIN.open(filename.c_str());

  string line;
  string keyword;
  string value;
  
  // Read the trace input filename from the first line.
  getline(fileIN, line);
  stringstream sline(line);
  sline >> trace_input;

  while(fileIN.good()){
    while(getline(fileIN, line)) {
      stringstream sline(line);
      sline >> keyword >> value;
      if (keyword=="n_cont" && atoi(value.c_str())>0){
         n_cont=atoi(value.c_str());
      }
      else if (keyword=="sigma_noise" && atof(value.c_str())>0.){
         sigma_noise=atof(value.c_str());
         sigma_cut=4.*sigma_noise;
         sigma_tail=2.*sigma_noise;
         F_scoring=4.*sigma_noise;
      }
      else if (keyword=="sigma_cut" && atof(value.c_str())>0.){
         sigma_cut=atof(value.c_str())*sigma_noise;
      }
      else if (keyword=="sigma_tail" && atof(value.c_str())>0.){
         sigma_tail=atof(value.c_str())*sigma_noise;
      }
      else if (keyword=="len_min" && atof(value.c_str())>0.){
         len_min=atof(value.c_str());
      }
      else if (keyword=="abnormal_f" && atof(value.c_str())>0.){
         abnormal_f=atof(value.c_str());
      }
      else if (keyword=="abnormal_x" && atof(value.c_str())>0.){
         abnormal_x=atof(value.c_str());
      }
      else if (keyword=="dx_interp" && atof(value.c_str())>0.){
         dx_interp=atof(value.c_str());
      }
      else if (keyword=="l_p" && atof(value.c_str())>0.){
         l_p=atof(value.c_str());
      }
      else if (keyword=="F_min" && atof(value.c_str())>0.){
         F_min=atof(value.c_str());
      }
      else if (keyword=="F_max" && atof(value.c_str())>0.){
         F_max=atof(value.c_str());
      }
      else if (keyword=="lim_p_max" && atof(value.c_str())>0.){
         lim_p_max=atof(value.c_str());
      }
      else if (keyword=="peak_width" && atof(value.c_str())>0.){
         peak_width=atof(value.c_str());
      }
      else if (keyword=="bin_width" && atof(value.c_str())>0.){
         bin_width=atof(value.c_str());
      }
      else if (keyword=="score_thresh" && atof(value.c_str())>0.){
         score_thresh=atof(value.c_str());
      }
      else if (keyword=="delta1" && atof(value.c_str())>0.){
         delta1=atof(value.c_str())*dx_interp;
      }
      else if (keyword=="delta2" && atof(value.c_str())>0.){
         delta2=atof(value.c_str())*dx_interp;
      }
      else if (keyword=="F_scoring" && atof(value.c_str())>0.){
         F_scoring=atof(value.c_str())*sigma_noise;
      }
      else if (keyword=="peaks_diff" && atof(value.c_str())>0.){
         peaks_diff=atof(value.c_str());
      }
      else if (keyword=="len_diff" && atof(value.c_str())>0.){
         len_diff=atof(value.c_str());
      }
      else if (keyword=="Npeaks_min" && atof(value.c_str())>0.){
         Npeaks_min=atof(value.c_str());
      }
      else if (keyword=="Npeaks_max" && atof(value.c_str())>0.){
         Npeaks_max=atof(value.c_str());
      }
      else if (keyword=="kNN" && atoi(value.c_str())>0){
         kNN=atoi(value.c_str());
      }
      else if (keyword=="r_cut" && atof(value.c_str())>0.){
         r_cut=atof(value.c_str());
      }
      else if (keyword=="print_f" && atoi(value.c_str())==1){
	print_f=1;
      }
      else if (keyword=="print_f_interp" && atoi(value.c_str())==1){
	print_f_interp=1;
      }
      else if (keyword=="print_Lc" && atoi(value.c_str())==1){
	print_Lc=1;
      }
      else if (keyword=="print_hist" && atoi(value.c_str())==1){
	print_hist=1;
      }
      else if (keyword=="print_w" && atoi(value.c_str())==1){
	print_w=1;
      }
      else if (keyword=="print_globalscore" && atoi(value.c_str())==1){
	print_globalscore=1;
      }
      else if (keyword=="print_dist" && atoi(value.c_str())==1){
	print_dist=1;
      }
      else if (keyword=="print_trlen" && atoi(value.c_str())==1){
	print_trlen=1;
      }
      else if (keyword=="print_stdout" && atoi(value.c_str())==1){
	print_stdout=1;
      }
      else if (keyword=="print_knndist" && atoi(value.c_str())==1){
	print_knndist=1;
      }
      else if (keyword=="print_align" && atoi(value.c_str())==1){
	print_align=1;
      }
      else if (keyword=="debug" && atoi(value.c_str())>0){
         debug=atoi(value.c_str());
      }
    } 
  } 
  fileIN.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The cut_trace function removes the negative contact part and the non-contact part from each trace.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cut_trace(vector<double> & data_x,vector<double> & data_f, int& start_f, int& tr_len_f, int& tr_no){
  double err=0.;
  int Ltr=data_x.size();
  // Detect traces containing extremely large force and/or extension values.  
  for (int i=0;i<Ltr;i++) {
    if (data_f[i]>abnormal_f){
      if(print_stdout) {cout << "rejected " << tr_no << " this traces contains abnormal force" << endl;}
      tr_len_f=0;
      break;
    }
    if (data_x[i]>abnormal_x){
      if(print_stdout) {cout << "rejected " <<tr_no << " this trace contains abnormal extension" << endl;}
      tr_len_f=0;
      break;
    }
  }
  // Find the starting point: where the initial negative contact part ends and the meaningful positive contact part begins. 
  int c=0;
  int last_i=0;
  for (int i=0;i<Ltr-1;i++) {
    if (data_f[i+1]>0 && data_f[i]>0 && data_x[i]>0) {
      c=c+1;
      last_i=i+1;
      if (c==n_cont) {break;}
    }
    else {
      c=0;
    }
  }
  if (c==0 || c<n_cont ) {
    if(print_stdout) {cout << "rejected " <<tr_no << " this is a completely negative trace" << endl;}
    start_f=0;
    tr_len_f=0;
    return;
  }
  start_f=last_i-c;
  int  l_window=50;
  int  l_test=10;
  int  l_min=20;
  // Assuming the last 50 measures are noise... 
  tr_len_f=Ltr-l_window;  
  double F_av=0.;
  for(int j=tr_len_f;j<Ltr;j++) F_av+=data_f[j];
  F_av/=double(l_window);
  // If the trace doesn't have a tail, we simply do no cut it: tr_len_f=Ltr.
  // If the trace doesn't have a tail, we simply discarde it
  if(F_av>F_min){	
    //tr_len_f=Ltr;
    tr_len_f=0;
    if(print_stdout) {cout << "rejected " <<tr_no << " this trace has no tail" << endl;}
  }
  else {
    while(sqrt(err*err)<sigma_cut && tr_len_f>l_min){ 
      double xa=0.; double ya=0.; double x2a=0.; double y2a=0.; double xya=0.; int n=0;
      for (int j=tr_len_f;j<=Ltr-1;j++){ 
	n++;
	xa+=data_x[j];
	ya+=data_f[j];
	x2a+=data_x[j]*data_x[j];
	y2a+=data_f[j]*data_f[j];
	xya+=data_x[j]*data_f[j];
      }
      x2a/=double(n); y2a/=double(n); xa/=double(n); ya/=double(n); xya/=double(n); 
      double bb=(xya-xa*ya)/(x2a-xa*xa);
      double aa=ya-bb*xa;
      err=0.;
      for (int j=tr_len_f-l_test;j<tr_len_f;j++){
	double fmod=aa+bb*data_x[j];
	err+=(fmod-data_f[j]);
      }
      err/=double(l_test);
      tr_len_f=tr_len_f-l_test;
    } 
  }
  int ntt=0;
  double sigma_tail0=0.;
  for(int i=tr_len_f+l_test;i<Ltr;i++){ 
     ntt++;
     sigma_tail0+=data_f[i]*data_f[i]; 
  }
  sigma_tail0=sqrt(sigma_tail0/ntt);
  if (sigma_tail0>sigma_tail) {
//      tr_len_f=0;
      if(print_stdout) {cout << "rejected " <<tr_no << " this trace has strange tail" << endl;}
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compute_w function computes the WLC score by calculating Lc values, estimating Lc histograms and finally computing and assigning the score to all points.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_w(vector<double> &  xx, vector<double> &  ff, vector<double> &  ww, int& tr_no, int& tr_len, double& Npeaks){
  // Compute Lc values.
  double alpha=kbT/l_p; 
  complex<double> beta(-0.5,0.5*sqrt(3.));
  vector<double> hist;
  vector<double> Lc_all;
  vector<int> ind_Lc;
  hist.resize(1000);
  int nt=0;
  double Lcmin=0.;
  double Lcmax=0.;
  for(int i=0;i<xx.size();i++){
    if(ff[i]>F_min && ff[i]<F_max){
      complex<double> omega=(4*ff[i]/alpha)-3.;
      complex<double> gamma=pow((pow(omega,3)-216.+12.*sqrt(complex<double>(324.)-3.*pow(omega,3))),(1./3.));
      complex<double>  p1=omega/12.*(1.+omega/(beta*gamma));
      complex<double> Lcc=(1./(omega/12.*(1.+omega/(beta*gamma))+(beta*gamma)/12.+1.))*xx[i];
      double Lc=real(Lcc);
      if (print_Lc==1) {cout << "Lc " << tr_no << " " << xx[i] << " " << ff[i] << " "<< Lc << endl;}
      Lc_all.push_back(Lc);
      ind_Lc.push_back(i);
      if(Lc>Lcmax)Lcmax=Lc;
      nt++;
    }
  }
  // Discard traces for which the histogram cannot be computed: traces with less than 10 Lc points.
  if(nt<10) return; 
  double renorm=double(nt)/double(xx.size());
  
  // Compute histograms.
  int n_bins=int((Lcmax-Lcmin)/bin_width)+10;
  int psize=n_bins-1;
  double bins[n_bins];
  vector<double> prob;
  prob.resize(psize);
  vector<double> bincenters;
  bincenters.resize(psize);
  bins[0]=Lcmin;
  for (int i=1;i<n_bins;i++){bins[i]=bins[i-1]+10.;}
  for (int i=0;i<psize;i++){ 
    bincenters[i]=0.5*(bins[i]+bins[i+1]);
    prob[i]=0.0;
  }
    
  double st=bins[0];
  double end=st+bin_width;
  for (int i=0;i<psize;i++){
    for (int j=0;j<nt;j++){
      if (Lc_all[j]>=st && Lc_all[j]<=end){
	prob[i]=prob[i]+1.;
      }
    }
    st=end;
    end=end+bin_width;
  }
  // Eliminate the first peak, if it starts from the first bin.
  double Lcstart=0.;
  int found_start=0;
  if(prob[0]==0)found_start=1;
  for (int i=1;i<psize-1;i++){ // find the first minimum of the histogram, start from there
     if((prob[i-1]>=prob[i] && prob[i+1]>=prob[i])  && found_start==0){
        Lcstart=bincenters[i];
        found_start=1;
     }
  }
  if(found_start==0)Lcstart=bincenters[psize-1]; //the trace is discarded

  double norm=0.0;
  for (int i=0;i<psize;i++){
    norm=norm+prob[i];
  }
  for (int i=0;i<psize;i++){
    prob[i]=prob[i]/norm;
  }
  if(print_hist==1){
    for (int i=0;i<psize;i++){
      cout << "h "<<  tr_no<< " "<< bincenters[i] << " " << prob[i]*renorm << endl;
    }
  }
  // Find minima and maxima.
  vector<int> max_indices;
  vector<int> min_indices;
    
  double df=prob[1]-prob[0];
  int is=0;
  
  // In case there are many 0 points in the beginning, find the first non-zero point and start minima and maxima detection from there.
  while(df==0){	
    is++;
    df=prob[is+1]-prob[is];
  }
    
  int increasing=1;
  min_indices.push_back(0);
  if(df<0){		
    prob[is]=prob[is+1]-1.e-8;
  }
    
  for (int i=is+1;i<psize-2;i++) {
    df=prob[i+1]-prob[i];
    if(increasing==1 && df<0){
      increasing=0;
      max_indices.push_back(i);
    }
    if(increasing==0 && df>0){
      increasing=1;
      min_indices.push_back(i);
    }
  }
  
  // Find the last minimum.
  df=prob[psize-1]-prob[psize-2];	
  if(df<=0){		
    min_indices.push_back(psize-1);
  }
  else {
    if(print_stdout) {cout << "rejected " <<tr_no << " " << tr_len<< " the last point of prob is a maximum"<<endl;}
    tr_len=0;
  }
  int n_mins=min_indices.size();
  int n_maxs=max_indices.size();
  
  if (print_hist==1){
    for (int i=0;i<n_maxs;i++) {
      cout <<"h " << tr_no << " m "<< bincenters[min_indices[i]] << " " << prob[min_indices[i]]*renorm << endl;
      cout <<"h " << tr_no << " M "<< bincenters[max_indices[i]] << " " << prob[max_indices[i]]*renorm << endl;
    }
    cout <<"h " << tr_no <<  " m "<< bincenters[min_indices[n_mins-1]] << " " << prob[min_indices[n_mins-1]]*renorm << endl;
  }

  if(n_mins-n_maxs != 1 && n_mins==2){
    if(print_stdout) {cout << "rejected " <<tr_no << " " << tr_len << " wrong min max assignation"<<endl;}
    tr_len=0;
  }

  int nmm=n_mins+n_maxs;
  int total_list[nmm];
  // In xm[n_maxs] we store the xmax value for each peak. It is used later for w score assignation.
  double xm[n_maxs];	
  for (int i=0;i<n_maxs;i++){
    total_list[i]=max_indices[i];
    xm[i]=0.0;
  }
  for (int i=0;i<n_mins;i++){
    total_list[i+n_maxs]=min_indices[i];
  }
  sort(total_list, total_list+nmm);
  
  // Compute the w score.
  int n_maxs0=0; 		// n_maxs0 counts the maxima different from zero
  double min_to_max=0.0;
  double wi=0.0;
  double max_w_score[n_maxs];
  int m=0;
  for (int i=0;i<n_maxs;i++){
    m=max_indices[i];
    for (int j=0;j<nmm;j++){
      if (m==total_list[j]){
       // A peak is considered as meaningful if defined by  more than 5 points.
       if (prob[total_list[j]]*norm>5){ 
	if(prob[total_list[j]]*renorm > lim_p_max){
	  min_to_max=0.5*(((prob[total_list[j-1]]/prob[total_list[j]]))+((prob[total_list[j+1]]/prob[total_list[j]])));
	  min_to_max=exp(-((min_to_max)*(min_to_max))/(2*0.5*0.5));
	  wi=min_to_max;
	  max_w_score[i]=wi;
          n_maxs0++; 
	}
        else {
	  max_w_score[i]=0.0;
	}
       } 
       else {	
          // Assign 0 score to a peak defined by 5 or less than 5 Lc points.
          max_w_score[i]=0.0;
       }	
      } 
    } 
  }  
  Npeaks=n_maxs0;
  // Assign the same w score to all points belonging to each peak. 
  for (int i=0;i<nt;++i){
    int kmax=-1;
    for(int k=1;k<n_mins;k++){
      if(Lc_all[i]>bin_width*min_indices[k-1] && Lc_all[i]<bin_width*min_indices[k]){ 
	kmax=k-1;
        // Find xmax for each peak.
	if (xx[ind_Lc[i]]>xm[k-1]) xm[k-1]=xx[ind_Lc[i]];
      }
    }
    if(kmax==-1){
      cout<< bin_width*min_indices[n_mins-1]<<endl;
      cout<< i <<" "<<Lc_all[i]<<endl;
      cout<< tr_no << " max not found"<<endl;
      exit(0);
    }
    else {
      ww[ind_Lc[i]]=max_w_score[kmax];
    }
    if(Lc_all[i]<Lcstart)ww[i]=200000.;
  }
  for (int i=0;i<xx.size();++i){
    if (ww[i]==0) {
      for (int j=0;j<n_maxs;j++){
	if(xx[i]<xm[j] and ff[i]>0 and xm[j]-xx[i]<peak_width){
	  ww[i]=max_w_score[j];
	  break; 
	}
      }  
    }
    if(ff[i]<0)ww[i]=0.;
    if(ww[i]==200000.)ww[i]=0.;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compute_global_score function sums the WLC point scores in each trace.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double compute_global_score(vector<double> &  ww, int& tr_no){
  double global_score=0.0;
  for (int k=0;k<ww.size();k++){
    global_score=global_score+ww[k];
  }
  return global_score;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The similarity_score function computes the match/mismatch score for the alignment.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double similarity_score(double a,double b){
  double result=0.;
  if (abs(a-b)<F_scoring){
     result=1-0.5*(abs(a-b)/(F_scoring));
  }
  else
  {
     result=(-0.5*abs(a-b)/(F_scoring));
  }
  return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The find_array_max function finds the maximum value in an array. It is used inside compute_dist_align.
double find_array_max(double array[],int length, int &ind){

  double max = array[0];            
  ind = 0;

  for(int i = 1; i<length; i++){
    if(array[i] > max){
      max = array[i];
      ind = i;
    }
  }
  return max;                    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compute_dist_align function performs global alignment and computes the similarity distance for a pair of traces a and b.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double compute_dist_align(double* forcea,double* wwa, int N_a, int a, double* forceb,double* wwb, int N_b, int b,int cl){

  int N_max=max(N_a,N_b);	
  double Fma=*max_element(forcea,forcea+N_a);
  double Fmb=*max_element(forceb,forceb+N_b);
  double F_av =0.5*(Fma+Fmb);
  vector <double> H((N_a+1)*(N_b+1),0.);	// Dynamic programming matrix
  vector <int> I_i((N_a+1)*(N_b+1),0);		// Index matrices used for backtracking
  vector <int> I_j((N_a+1)*(N_b+1),0); 		
			
  // Initialize first row and first column with -i*delta.
  int row0=0;
  double penalty=0;
  for (int col_n=0; col_n<=N_b;col_n++){
    int ngap=int(10./dx_interp);
    if(col_n < ngap) penalty=penalty-delta1;
    else penalty = penalty - delta2;
    H.at((N_b+1)*row0+col_n)=penalty;
  }

  int col0=0;
  penalty=0;
  for (int row_n=0; row_n<=N_a;row_n++){
    int ngap=int(10./dx_interp);
    if(row_n < ngap) penalty=penalty-delta1;
    else penalty = penalty - delta2;
    H.at(col0+row_n*(N_b+1))=penalty;
  }
			
  double temp[3];
  // Here comes the actual algorithm.
  for(int i=1;i<=N_a;i++){            
    for(int j=1;j<=N_b;j++){
      temp[0] = H.at((N_b+1)*(i-1)+(j-1))+similarity_score(forcea[i-1],forceb[j-1]); 
      if (i<10){
	temp[1] = H.at((N_b+1)*(i-1)+j)-delta1;}
      else {
	temp[1] = H.at((N_b+1)*(i-1)+j)-delta2;}
      if (j<10) {                 
	temp[2] = H.at((N_b+1)*(i)+(j-1))-delta1;}
      else {
	temp[2] = H.at((N_b+1)*(i)+(j-1))-delta2;}                
      int ind;
      H.at((N_b+1)*(i)+j) = find_array_max(temp,3,ind);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
	I_i.at((N_b+1)*i+j) = i-1;			
	I_j.at((N_b+1)*i+j) = j-1;
	break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence a

	I_i.at((N_b+1)*i+j) = i-1;                 
	I_j.at((N_b+1)*i+j) = j;
	break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence b
	I_i.at((N_b+1)*i+j) = i;
	I_j.at((N_b+1)*i+j) = j-1;
	break;
      }		   
    }			    
  } 		

  // By definition the highest score corresponding to the best alignment is in the last cell.	
  int i_max=N_a;
  int j_max=N_b;
  double H_max=H.at((N_b+1)*i_max+j_max);
  // Backtracking from H_max...
  int current_i=i_max,current_j=j_max;
  int next_i=I_i.at((N_b+1)*current_i+current_j);
  int next_j=I_j.at((N_b+1)*current_i+current_j);
  int tick=0;
  int consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2]; 
  int consensus_x_a[N_a+N_b+2],consensus_x_b[N_a+N_b+2]; 
  int indices_tb_a[N_a+N_b+2],indices_tb_b[N_a+N_b+2];


  while(((current_i!=next_i) || (current_j!=next_j)) && ((next_j!=0) || (next_i!=0))){

    if(next_i==current_i && next_i>10)  {
      consensus_a[tick] = -20000;		// deletion in a with gap penalty delta2
      consensus_x_a[tick] = -20000;
      indices_tb_a[tick]=-1;}  			
    else if(next_i==current_i && next_i<10)  {
      consensus_a[tick] = -10000;     		// deletion in a with gap penalty delta1
      consensus_x_a[tick] = -10000;
      indices_tb_a[tick]=-1;}  			
    else  {
      consensus_a[tick] = forcea[current_i-1];	// match/mismatch in a
      consensus_x_a[tick] =current_i-1;
      indices_tb_a[tick]=current_i-1;}		

    if(next_j==current_j && next_j>10)  {
      consensus_b[tick] = -20000;		// deletion in b with gap penalty delta2
      consensus_x_b[tick] =-20000;
      indices_tb_b[tick]=-1;} 	
    else if(next_j==current_j && next_j<10)  {
      consensus_b[tick] = -10000;	       // deletion in b with gap penalty delta1
      consensus_x_b[tick] = -10000;
      indices_tb_b[tick]=-1;} 
    else  {
      consensus_b[tick] = forceb[current_j-1];	// match/mismatch in b
      consensus_x_b[tick] =  current_j-1;
      indices_tb_b[tick]=current_j-1;}	
    current_i = next_i;
    current_j = next_j;
    next_i = I_i.at((N_b+1)*current_i+current_j);
    next_j = I_j.at((N_b+1)*current_i+current_j);

    tick++;
  }  

  double H_i[tick];
  double w_i_a[tick];
  double w_i_b[tick];
  int gg=0;

  for(int i=tick-1;i>=0;i--) {
    // Find local scores H_i[i].
    if(consensus_a[i]==-10000 || consensus_b[i]==-10000){
      H_i[i]=-delta1;} 
    else if (consensus_a[i]==-20000 || consensus_b[i]==-20000){
      H_i[i]=-delta2;} 
    else {
      H_i[i]=similarity_score(consensus_a[i],consensus_b[i]);}
    // Find local quality scores, w_i in sequences a and b. 
    if(consensus_a[i]==-10000 || consensus_a[i]==-20000){
      w_i_a[i]=0.0;}
    else {
      w_i_a[i]=wwa[indices_tb_a[i]];}
    if(consensus_b[i]==-10000 || consensus_b[i]==-20000){
      w_i_b[i]=0.0;}
    else {
      w_i_b[i]=wwb[indices_tb_b[i]];}
  }
  if(print_align){
    for(int i=tick-1;i>=0;i--) {
      cout << cl << " "<< a << " " << b << " " << indices_tb_a[i] << " "<< indices_tb_b[i] <<" ";
      if(consensus_a[i]==-20000 || consensus_a[i]==-10000 ) {
	cout << "- " << w_i_a[i] << " ";}
      else
	cout << consensus_a[i]<<" " << w_i_a[i]<< " ";
      if(consensus_b[i]==-20000 || consensus_b[i]==-10000) {
	cout << "- " << w_i_b[i] <<" " << H_i[i] << endl;}
      else
	cout << consensus_b[i]  << " " << w_i_b[i] << " " << H_i[i] << endl;
    }
  }
  // Compute distance.
  double distf=1-(H_max/N_max); 	
  return distf;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char** argv)
{
  if(argc>1)
  {
      filename=argv[1];      
  }
  else
  {
      cout << "ERROR: No input file has been passed in the command line." << endl;
      exit(0);
  }
  string trace_input;
  read_input(filename,trace_input,n_cont,sigma_noise,sigma_cut,sigma_tail,len_min,abnormal_f,abnormal_x,dx_interp,l_p,F_min,F_max,lim_p_max,peak_width,bin_width,score_thresh,delta1,delta2,F_scoring,peaks_diff,len_diff,Npeaks_min,Npeaks_max,kNN,r_cut,print_f,print_f_interp,print_Lc,print_hist,print_w,print_globalscore,print_dist,print_trlen,print_stdout,print_knndist,print_align,debug);
  // Print the parameter values used for this particular run.
  cout << "# List of parameters used for this run." << endl;
  cout << "n_cont= " << n_cont << endl;			
  cout << "sigma_noise= " << sigma_noise << endl;
  cout << "sigma_cut= " << sigma_cut << endl;	        
  cout << "sigma_tail= " << sigma_tail << endl;		
  cout << "len_min= " << len_min << endl;	
  cout << "abnormal_f= " << abnormal_f << endl;	
  cout << "abnormal_x= " << abnormal_x << endl;	
  cout << "dx_interp= " <<  dx_interp << endl; 		
  cout << "l_p= " <<  l_p << endl; 		
  cout << "F_min= " <<  F_min << endl; 		
  cout << "F_max= " <<  F_max << endl; 		
  cout << "lim_p_max= " <<  lim_p_max << endl; 		
  cout << "peak_width= " <<  peak_width << endl; 		
  cout << "bin_width= " << bin_width << endl; 		
  cout << "score_threshold= " << score_thresh << endl;	
  cout << "delta1= " <<  delta1 << endl; 		
  cout << "delta2= " <<  delta2 << endl; 		
  cout << "F_scoring= " <<  F_scoring << endl; 		
  cout << "peaks_diff= " <<  peaks_diff << endl; 		
  cout << "len_diff= " <<  len_diff << endl; 		
  cout << "Npeaks_min= " <<  Npeaks_min << endl; 		
  cout << "Npeaks_max= " <<  Npeaks_max << endl; 		
  cout << "kNN= " << kNN << endl; 			
  cout << "r_cut= " << r_cut << endl; 			
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //if (debug==1) return 0; //DEBUG
  int t=0;
  ifstream fileIN;
  fileIN.open(trace_input.c_str());
  int nl=0;
  int num_traces=0;
  string line;
  vector <double>  data_x;
  vector <double>  data_f;
  vector <double>  force;
  vector <double>  ww;
  vector <double>  global_score;
  vector <int>     start_tr;
  vector <int>     Length_tr;
  vector <double>  Length_tr_nm;
  vector <int>     orig_index;
  vector <double>  n_peaks;
  int Ltot=0;
  int start;
  int tr_len;
  int tr_no;
  double xold=0.;
  double Npeaks;					  // Total number of high quality peaks in the Lc histogram;
  while (getline(fileIN,line)) {      
    istringstream streamA(line);
    double xr;
    while(streamA >> xr) {
      if(nl%2==0)data_f.push_back(xr*1.e12);
      if(nl%2==1)data_x.push_back(xr*1.e9);
    }
    if(nl%2==0) {tr_no=int(nl/2)+1;}	      		 // tr_no is the original trace number;  
    if(nl%2==1){
      vector <double>  xt;
      vector <double>  ft;
      for(int i=0;i<data_x.size()-2;i++){ 
	if(data_x[i]>0 && data_f[i]>0){	
	  if(xt.size()==0 || data_x[i]>xt[xt.size()-1]){ 	
	    xt.push_back(data_x[i]);
	    ft.push_back(data_f[i]);
	    if (print_f==1) {cout << "f " << tr_no << " " << data_x[i] << " " << data_f[i] << endl;}
	  }
	}
      }
      // Remove the negative contact anf the non-contact part (the "tail") of the trace.
      cut_trace(xt,ft,start,tr_len,tr_no); 
      double ww_sum=0.0;
      double ratio=0.0;
      vector <double>  f_interp,x_interp,wwl;
      // Remove traces with gaps (negative force) longer that 10 nm.
      for(int i=1;i<tr_len;i++)if(xt[i]-xt[i-1]>10.)tr_len=0.;  
      if (tr_len>20){  
	int iO=int(xt[0]/dx_interp);
	for(int i=0;i<=start;i++) {ft[i]=0.0;}
	for(int i=1;i<tr_len;i++){ 
	  int iN=int(xt[i]/dx_interp);
	  if(iN>iO){
	    for(int k=iO+1;k<=iN;k++){
	      double f_int=ft[i-1]+(ft[i]-ft[i-1])*(double(k)*dx_interp-xt[i-1])/(xt[i]-xt[i-1]);
	      f_interp.push_back(f_int);
	      x_interp.push_back(double(k)*dx_interp);
	      wwl.push_back(0);
	    }
	  }
	  iO=iN;
	}
	tr_len=x_interp.size();
	double x_len=x_interp[x_interp.size()-1];
        // Find maximum force value in each trace.
        double Fm=-1000.;
        for (int k=50;k<tr_len;k++){
	    if (f_interp[k]>Fm) Fm=f_interp[k]; 
        }
        if (print_f_interp==1){
	   for (int i=0;i<tr_len;i++){
	     cout << "fi " << tr_no << " " << x_interp[i] << " " << f_interp[i] << endl;
	   }
        }
        if(print_stdout) {
           if (tr_len*dx_interp<=len_min){
              cout << "rejected " <<tr_no << " " << tr_len << " this trace is too short" << endl;
           }   
        } 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
	// Compute w score.
        if (tr_len*dx_interp>len_min){  
          compute_w(x_interp,f_interp,wwl,tr_no,tr_len,Npeaks);  
	  ww_sum=compute_global_score(wwl,tr_no); 
	  ratio=ww_sum/tr_len;
        }
        if (print_w==1){
          for(int i=0;i<tr_len;i++){
              cout<<"w " << tr_no <<" "<< x_interp[i] <<" "<< f_interp[i]<<" " << wwl[i]<<endl; 
          }
        }
        if(print_stdout) {
           if (tr_len*dx_interp>len_min && ratio<score_thresh){
              cout << "rejected " <<tr_no << " " << ratio << " this trace has ratio<0.5" << endl;
           }
        }
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(Npeaks>Npeaks_max){
           tr_len=0;
           if(print_stdout)cout << "rejected " <<tr_no << " this trace has more  peaks than " <<Npeaks_max<< endl;
        }
        if(Npeaks<Npeaks_min){
           tr_len=0;
           if(print_stdout)cout << "rejected " <<tr_no << " this trace has less peaks than " <<Npeaks_min<< endl;
        }
        if(tr_len>len_min && ratio>score_thresh){  
	  Length_tr.push_back(tr_len);
	  Length_tr_nm.push_back(x_len);
	  global_score.push_back(ww_sum);
	  start_tr.push_back(Ltot);
	  Ltot=Ltot+tr_len;
	  orig_index.push_back(tr_no);
	  n_peaks.push_back(Npeaks);
	  ++num_traces;
	  force.insert(force.end(),f_interp.begin(),f_interp.begin()+tr_len);
	  ww.insert(ww.end(),wwl.begin(),wwl.begin()+tr_len);
	}
      } 
      if (print_globalscore) {
	cout << "gl " << tr_no << " " << ww_sum << " " << tr_len << endl;
      }
      data_x.clear();
      data_f.clear();
      xt.clear();
      ft.clear();
      f_interp.clear();
      x_interp.clear();
    } 
    nl++;
  }
  fileIN.close();
  cout << "Total number of traces passing the quality filter= " << num_traces << endl;
  if (print_w==1){
    for(int i=0;i<num_traces;i++){
      for(int j=start_tr[i];j<start_tr[i]+Length_tr[i];j++){
	cout<<"w " << orig_index[i] <<" "<< (j-start_tr[i]+1)*dx_interp <<" "<< force[j]<<" " << ww[j]<<endl; 
      }
    }
  }
  if (print_trlen==1){
    for(int i=0;i<num_traces;i++){
      cout << "l " << orig_index[i] << " " << Length_tr[i] << " " << Length_tr_nm[i] <<endl;	
    }
  }
  //if (debug==1) return 0; //DEBUG
  // Compute distances and densities.
  vector <dist_s> fakedist(num_traces);
  vector <dist_s> truedist(num_traces); 
  vector <double> rho(num_traces);
  vector <dp_s> DP(num_traces);
  for (int i=0;i<num_traces;i++){
    DP[i].dist_nn.resize(num_traces); 
    DP[i].ind_nn.resize(num_traces); 
  }
  if(print_align==1){
    int i=0;
    for (int nn=1; nn<10; nn++){
      double dist=compute_dist_align(&force[start_tr[i]],&ww[start_tr[i]],Length_tr[i],i,&force[start_tr[nn]],&ww[start_tr[nn]],Length_tr[nn],nn,0);
    }
  }
  else {
    for(int i=0;i<num_traces;i++){
#pragma omp parallel 
      {
#pragma omp for
      for (int k=0;k<num_traces;++k){ 
          double len_a=Length_tr[i];
          double len_b=Length_tr[k];
          double len_max=max(len_a,len_b);
          double diff_len=abs(len_a-len_b);
          diff_len=diff_len/len_max;
          double diff_peaks=abs(n_peaks[i]-n_peaks[k]);
          if (i==k || diff_len>len_diff || diff_peaks>peaks_diff) {
	    double dist=5.;
            truedist[k].dist=dist;
	    truedist[k].ind=k;           
	    continue;
	  }
	  double dist=compute_dist_align(&force[start_tr[i]],&ww[start_tr[i]],Length_tr[i],i,&force[start_tr[k]],&ww[start_tr[k]],Length_tr[k],k,0);
	  truedist[k].dist=dist;         
	  truedist[k].ind=k;           
	}
     }
     sort(truedist.begin(), truedist.end(),lt);
     if (print_dist) {
      	for (int j=0;j<num_traces;++j){ 
        	  cout << "distance " << orig_index[i] << " "  << orig_index[truedist[j].ind] << " " << truedist[j].dist << endl;
     	 }
      } 
      if (print_knndist) {
	for (int l=0;l<kNN;l++){ 
	  cout << "knndist " << orig_index[i]  << " " << orig_index[truedist[l].ind] << " " <<  truedist[l].dist << endl;
	}
      }
      DP[i].nn=0;
      DP[i].i_ini=i;
      for (int k=0;k<num_traces;++k){ 
	DP[i].ind_nn[k]=truedist[k].ind;
	DP[i].dist_nn[k]=truedist[k].dist;
        if(truedist[k].dist<r_cut)DP[i].nn++;
      }
      DP[i].rho=log(1./truedist[kNN-1].dist);
//      DP[i].rho=DP[i].nn*(global_score[i]/Length_tr_nm[i]);  
    }
 }
  //////////////////////////////////////////////////////////////////////////
  // Compute delta.
  for(int i=0;i<num_traces;i++){
      DP[i].cl=-1;
      DP[i].cen=0;
      DP[i].delta=10.;
      for (int k=0;k<num_traces;++k){
         int ik=DP[i].ind_nn[k];
         double distk=DP[i].dist_nn[k];
         if(DP[ik].rho > DP[i].rho && distk<DP[i].delta)DP[i].delta=distk; 
      }
      DP[i].gamma=DP[i].rho*DP[i].delta;
      if(DP[i].nn<kNN)DP[i].gamma=0.;
  }
  // Find the cluster centers.
  sort(DP.begin(), DP.end(),ltgamma);

  int  nmaxcen=num_traces/10;
  if(nmaxcen>100)nmaxcen=100;
  int addcen=1;
  NCEN=1;
  DP[0].cl=1;
  DP[0].cen=1;
  for (int i=1;i<num_traces;++i){
     if(DP[i].gamma>0){
        addcen=1;
        int ii=DP[i].i_ini;
        for (int k=0;k<i;++k){
           if(DP[k].cen==1){
              int icen=DP[k].i_ini;
              double dist=compute_dist_align(&force[start_tr[ii]],&ww[start_tr[ii]],Length_tr[ii],ii,&force[start_tr[icen]],&ww[start_tr[icen]],Length_tr[icen],icen,0);
              if(dist<r_cut)addcen=0;
           }
        }
        if(addcen==1){
           NCEN++;
           DP[i].cl=NCEN;
           DP[i].cen=1;
        }
        if(NCEN==nmaxcen){
           cout<<"More clusters than the maximum allowed value (min(100,num_traces/10))"<<endl;
           addcen=0;
        }
     }
  }
/*
  while (addcen==1){
     addcen=1;
     for (int k=0;k<num_traces;++k){
         int ik=DP[NCEN].ind_nn[k];
         double distk=DP[NCEN].dist_nn[k];
         for(int i=0;i<NCEN;i++){
            int icen=DP[i].i_ini;
            if(icen==ik && distk<r_cut)addcen=0;
         }
     }
     if(addcen==1)NCEN++;
     if(NCEN==nmaxcen){
        cout<<"More clusters than the maximum allowed value (min(100,num_traces/10))"<<endl;
        addcen=0;
     }
  }

  for(int i=0;i<NCEN;i++){
      DP[i].cl=i+1;
      DP[i].cen=1;
  }
*/
  //////////////////////////////////////////////////////////////////////////
  sort(DP.begin(), DP.end(),index);
  for(int i=0;i<num_traces;i++){
     if(DP[i].i_ini!=i) cout<<"error in the order"<<endl;
  }

  // Find the nearest neighbour of higher density
  for(int i=0;i<num_traces;i++){
      DP[i].nn_hrho=-100;
      double distm=10000.;
      for (int k=0;k<num_traces;++k){
         int ik=DP[i].ind_nn[k];
         if(DP[ik].rho > DP[i].rho){
            double distk=DP[i].dist_nn[k];
            if(distk<distm){
               distm=distk;
               DP[i].nn_hrho=ik;
            }
         }
      }
  }
  // Perform cluster assignation.
  //  
  int nfound=NCEN;
  while (nfound<num_traces){
    for(int i=0;i<num_traces;i++){
      if(DP[i].cl==-1){
	int j=DP[i].nn_hrho;
	if(DP[j].cl>-1){
	  DP[i].cl=DP[j].cl;
	  nfound++;
	}
      }
    }
  }

  // Compute the distance of the trace from the center of its cluster, and from all the other cluster centers.
  cout << "#COL2 \t Trace number" << endl; 
  cout << "#COL3 \t Length [nm]" << endl; 
  cout << "#COL4 \t Quality score" << endl; 
  cout << "#COL5 \t Number of peaks" << endl; 
  cout << "#COL6 \t Cluster number (=0 if COL6>r_cut)" << endl;
  cout << "#COL7 \t Cluster number without filters" << endl;
  cout << "#COL8 \t Distance from cluster center" << endl;
  for(int i=0;i<num_traces;i++){
      double dist_center=100.;
      double dist_other=100.;
      for (int k=0;k<num_traces;++k){
         int ik=DP[i].ind_nn[k];
         if(DP[ik].cen == 1){
            double distk=DP[i].dist_nn[k];
            if(DP[ik].cl==DP[i].cl)DP[i].dist_center=distk;
            if(DP[ik].cl!=DP[i].cl){
                 if(distk<dist_other) {
                   dist_other=distk;
                   DP[i].dist_other=dist_other;
                 }
            }
         }
      }
      if(DP[i].cen == 1)DP[i].dist_center=0.;
      for (int j=1;j<=NCEN;j++){
          if (DP[i].cl==j){
               DP[i].core=0;
               if(DP[i].dist_center<=r_cut && DP[i].dist_center<DP[i].dist_other)DP[i].core=DP[i].cl;
               cout << "tr "<< orig_index[i]<<"\t"<< Length_tr[i] << "\t" << global_score[i] << "\t" ;
               cout <<  n_peaks[i] << "\t" << DP[i].core<< "\t" << DP[i].cl << "\t" << DP[i].dist_center <<  endl;
          }
      }
  }
  vector <clust_s> cluster(NCEN+1);
  for(int cl=0;cl<=NCEN;cl++){
    cluster[cl].members.resize(0);
  }
  for(int i=0;i<num_traces;i++){
     int cl=DP[i].core;
     if (DP[i].cen==1)cluster[cl].center=orig_index[i];
     if(cl!=0){
        vector<int>::iterator it=cluster[cl].members.begin();
        vector<double>::iterator jt=cluster[cl].distcenter.begin();
        double dcen=DP[i].dist_center;
	int ipos=0;
        int sz=cluster[cl].size;
        if(sz==0)ipos=0;
        if(sz==1)if(dcen>cluster[cl].distcenter[0])ipos=1;
        if(sz>1)for(int j=1;j<sz;j++)if(dcen>cluster[cl].distcenter[j-1] && dcen<=cluster[cl].distcenter[j])ipos=j;
        if(sz>1)if(dcen>cluster[cl].distcenter[sz-1])ipos=sz;
        cluster[cl].members.insert(it+ipos,orig_index[i]);
        cluster[cl].distcenter.insert(jt+ipos,DP[i].dist_center);
        cluster[cl].size++;
     }
  }
  // Print summary information about the clusters.
  cout << "# Clusters" << endl;
  cout << "Total number of clusters: " << NCEN << endl;
  sort(cluster.begin(), cluster.end(),ltsize);
  for(int cl=0;cl<NCEN;cl++){
     if(cluster[cl].size>=kNN){
        cout << " cluster " << cl+1 << " center: " << cluster[cl].center ;
        cout << " size= " << cluster[cl].size <<  " members: ";
        for(int i=0;i<cluster[cl].size;i++)cout << cluster[cl].members[i]<<" ";
        cout<< endl;
     }
  }
  return 0;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}



