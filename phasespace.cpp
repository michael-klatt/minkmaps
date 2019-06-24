#include "phasespace.h"

// C is measured in units \chi
// P is measured in units P/2
// n_bins_P (2*(N/2+1)*N+3), n_bins_C ((N/2+1)*(N/2+1)+N*(N/2+1)+4), abs_min_C (N*(N/2+1)+1)
// max_P (n_bins_P), max_C ((N/2+1)*(N/2+1)+1), min_C ((-1)*abs_min_C)

PhaseSpace::PhaseSpace(const unsigned &N, const unsigned &A,
		       const double &nmbr_of_microstates, const unsigned &NmbrMacrostates,
		       const double &f_initial, const double &f_mod, const double &f_epsilon, 
		       const double &hist_dev,
		       const std::string &prefix_if, const std::string &prefix_of,
		       const bool &verbose,
		       const unsigned &seed, const unsigned &seedwght) :
  N (N), N2 (N*N), A (A), W (N2-A), nmbr_of_microstates (nmbr_of_microstates), NmbrMacrostates (NmbrMacrostates), timescale (1.0/NmbrMacrostates), f_initial (f_initial), f_mod (f_mod), f (f_initial), f_epsilon (f_epsilon), hist_dev (hist_dev), sample (BinField<bool> (N+2, false)), black_pixels_x (std::vector<unsigned> (A, 0)), black_pixels_y (std::vector<unsigned> (A, 0)), white_pixels_x (std::vector<unsigned> (W, 0)), white_pixels_y (std::vector<unsigned> (W, 0)), xb (0), yb (0), xw (0), yw (0), P0 (0), C0 (0), bin0 (0), P1 (0), C1 (0), bin1 (0), n_bins_P (2*(N/2+1)*N+3), n_bins_C ((N/2+1)*(N/2+1)+N*(N/2+1)+4), abs_min_C (N*(N/2+1)+1), max_P (n_bins_P), max_C ((N/2+1)*(N/2+1)+1), min_C ((-1)*abs_min_C), log_DoS (std::vector<double> (n_bins_P * n_bins_C,0)), log_DoS_all (std::vector< std::vector<double> > (20,std::vector<double> (n_bins_P * n_bins_C,0))), DoS (std::vector<double> (n_bins_P * n_bins_C,0)), hist (std::vector<unsigned> (n_bins_P * n_bins_C,0)), non_zero_hist (0), hist_min(0), hist_mean (0), total_count (0), prefix_if (prefix_if), prefix_of (prefix_of), verbose (verbose), seed (seed+A), seedwght (seedwght+seed+A*2), original_seed (seed), original_seedwght (seedwght)
{
  std::stringstream intermediate_DoS_stst;
  intermediate_DoS_stst << prefix_of << "intermediate_log_DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
			<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string intermediate_DoS_st= intermediate_DoS_stst.str();
  std::ifstream intermediate_DoS( intermediate_DoS_st.c_str() );
  
  if( intermediate_DoS.good() ){
    
    // Read in: non_zero_hist; log_DoS; f; total_count; black_pixels_x, black_pixels_y, white_pixels_x, white_pixels_y;
    std::string dummy_line="";
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    intermediate_DoS >> f;
    std::cout << f << std::endl;
    getline(intermediate_DoS,dummy_line);
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    intermediate_DoS >> total_count;
    std::cout << total_count << std::endl;
    getline(intermediate_DoS,dummy_line);
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    for(unsigned ai = 0; ai < A; ai++){
      intermediate_DoS >> black_pixels_x[ai];
      std::cout << std::setw(5) << black_pixels_x[ai];
    }
    std::cout << std::endl;
    getline(intermediate_DoS,dummy_line);
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    for(unsigned ai = 0; ai < A; ai++){
      intermediate_DoS >> black_pixels_y[ai];
      std::cout << std::setw(5) << black_pixels_y[ai];
    }
    std::cout << std::endl;
    getline(intermediate_DoS,dummy_line);
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    for(unsigned wi = 0; wi < W; wi++){
      intermediate_DoS >> white_pixels_x[wi];
      std::cout << std::setw(5) << white_pixels_x[wi];
    }
    std::cout << std::endl;
    getline(intermediate_DoS,dummy_line);
    
    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    for(unsigned wi = 0; wi < W; wi++){
      intermediate_DoS >> white_pixels_y[wi];
      std::cout << std::setw(5) << white_pixels_y[wi];
    }
    std::cout << std::endl;
    getline(intermediate_DoS,dummy_line);

    getline(intermediate_DoS,dummy_line);
    std::cout << dummy_line << " ";
    intermediate_DoS >> non_zero_hist;
    std::cout << non_zero_hist << std::endl;
    getline(intermediate_DoS,dummy_line);
 
    getline(intermediate_DoS,dummy_line);
    // "i log_DoS[i]"
    unsigned i=0;
    for(unsigned ui = 0; ui < non_zero_hist; ui++){
      intermediate_DoS >> i;
      intermediate_DoS >> log_DoS[i];
    }
    
    for(unsigned ai = 0; ai < A; ai++)
      sample.assign(black_pixels_x[ai],black_pixels_y[ai],true);
    
  }
  else{

    std::vector<unsigned> all_pixels_x(N2, 0);
    std::vector<unsigned> all_pixels_y(N2, 0);
    for(unsigned xi = 0; xi < N; xi++)
      for(unsigned yi = 0; yi < N; yi++){
	all_pixels_x[xi*N+yi] = xi+1;
	all_pixels_y[xi*N+yi] = yi+1;
      }
  
    for(unsigned ai = 0; ai < A; ai++){
      unsigned ni = RandomGSLInteger(N2-ai, seed);
      black_pixels_x[ai] = all_pixels_x[ni];
      black_pixels_y[ai] = all_pixels_y[ni];
      all_pixels_x.erase(all_pixels_x.begin()+ni);
      all_pixels_y.erase(all_pixels_y.begin()+ni);
      sample.assign(black_pixels_x[ai],black_pixels_y[ai],true);
    }
  
    for(unsigned wi = 0; wi < W; wi++){
      white_pixels_x[wi] = all_pixels_x[wi];
      white_pixels_y[wi] = all_pixels_y[wi];
    }
    
  }
  intermediate_DoS.close();
    
  int A0 = area_wbc_pix(sample);
  P0 = perimeter_wbc_pix(sample)/16;
  C0 = euler_wbc_pix(sample)/8;
  P1 = P0;
  C1 = C0;
  
  bin0 = n_bin(P0,C0);
  bin1 = bin0;
  
  std::cout << "Initial sample:" << std::endl;
  display_sample();
  
  if( A0 != int(A*8) ){
    std::cerr << "ERROR: initialization or read in failed; area is " << A0 << ", compare given A = " << A << ";" << std::endl;
    exit(-1);
  }
}

void PhaseSpace::reset(){
  f=f_initial;
  black_pixels_x=std::vector<unsigned> (A, 0);
  black_pixels_y=std::vector<unsigned> (A, 0);
  white_pixels_x=std::vector<unsigned> (W, 0);
  white_pixels_y=std::vector<unsigned> (W, 0);
  P0=0;
  C0=0;
  bin0=0;
  P1=0;
  C1=0;
  bin1=0;
  log_DoS=std::vector<double> (n_bins_P * n_bins_C,0);
  DoS=std::vector<double> (n_bins_P * n_bins_C,0);
  hist=std::vector<unsigned> (n_bins_P * n_bins_C,0);
  non_zero_hist=0;
  hist_min=0;
  hist_mean=0;
    
  std::vector<unsigned> all_pixels_x(N2, 0);
  std::vector<unsigned> all_pixels_y(N2, 0);
  for(unsigned xi = 0; xi < N; xi++)
    for(unsigned yi = 0; yi < N; yi++){
      all_pixels_x[xi*N+yi] = xi+1;
      all_pixels_y[xi*N+yi] = yi+1;
      sample.assign(xi+1,yi+1,false);
    }
  
  for(unsigned ai = 0; ai < A; ai++){
    unsigned ni = RandomGSLInteger(N2-ai, seed);
    black_pixels_x[ai] = all_pixels_x[ni];
    black_pixels_y[ai] = all_pixels_y[ni];
    all_pixels_x.erase(all_pixels_x.begin()+ni);
    all_pixels_y.erase(all_pixels_y.begin()+ni);
    sample.assign(black_pixels_x[ai],black_pixels_y[ai],true);
  }
  
  for(unsigned wi = 0; wi < W; wi++){
    white_pixels_x[wi] = all_pixels_x[wi];
    white_pixels_y[wi] = all_pixels_y[wi];
  }
  
  int A0 = area_wbc_pix(sample);
  P0 = perimeter_wbc_pix(sample)/16;
  C0 = euler_wbc_pix(sample)/8;
  P1 = P0;
  C1 = C0;
  
  bin0 = n_bin(P0,C0);
  bin1 = bin0;
  
  std::cout << "Initial sample:" << std::endl;
  display_sample();
  
  if((A0 != int(A*8)) || (all_pixels_x.size() != W)){
    std::cerr << "ERROR: initialization failed, area is " << A0/8 << ", compare given A = " << A << ";" << std::endl
	      << "       and number of pixels left is " << all_pixels_x.size() << ", compare given W = " << W << ";" << std::endl;
    exit(-1);
  }  
}

void PhaseSpace::CHANGE_STATE(){
  int bi = RandomGSLInteger(A, seed);
  int wi = RandomGSLInteger(W, seed);
  
  //compute P1, C1;
  
  xb = black_pixels_x[bi];
  yb = black_pixels_y[bi];
  
  xw = white_pixels_x[wi];
  yw = white_pixels_y[wi];
  
  change_macrostate(xb, yb, xw, yw);
  bin1 = n_bin(P1,C1);
  
  if( log_DoS[bin1] < log_DoS[bin0] ){
    P0 = P1;
    C0 = C1;
    bin0 = bin1;
    
    sample.assign(xw,yw,true);
    black_pixels_x[bi] = xw;
    black_pixels_y[bi] = yw;
    white_pixels_x[wi] = xb;
    white_pixels_y[wi] = yb;
  }
  else if( RandomRational(seedwght) < exp(log_DoS[bin0] - log_DoS[bin1]) ){
    P0 = P1;
    C0 = C1;
    bin0 = bin1;
    
    sample.assign(xw,yw,true);
    black_pixels_x[bi] = xw;
    black_pixels_y[bi] = yw;
    white_pixels_x[wi] = xb;
    white_pixels_y[wi] = yb;
  }
  else{
    P1 = P0;
    C1 = C0;
    bin1 = bin0;
    sample.assign(xb,yb,true);
  }
  
  // hist_sum++;
  // if( hist[bin1] < 0.1 )
  //   non_zero_hist++;
  
  hist[bin1] += 1;
  log_DoS[bin1] += f;
}

void PhaseSpace::CHANGE_STATE_wo_hist(){
  int bi = RandomGSLInteger(A, seed);
  int wi = RandomGSLInteger(W, seed);
  
  //compute P1, C1;
  
  xb = black_pixels_x[bi];
  yb = black_pixels_y[bi];
  
  xw = white_pixels_x[wi];
  yw = white_pixels_y[wi];
  
  change_macrostate(xb, yb, xw, yw);
  bin1 = n_bin(P1,C1);
  
  if( log_DoS[bin1] < log_DoS[bin0] ){
    P0 = P1;
    C0 = C1;
    bin0 = bin1;
    
    sample.assign(xw,yw,true);
    black_pixels_x[bi] = xw;
    black_pixels_y[bi] = yw;
    white_pixels_x[wi] = xb;
    white_pixels_y[wi] = yb;
  }
  else if( RandomRational(seedwght) < exp(log_DoS[bin0] - log_DoS[bin1]) ){
    P0 = P1;
    C0 = C1;
    bin0 = bin1;
    
    sample.assign(xw,yw,true);
    black_pixels_x[bi] = xw;
    black_pixels_y[bi] = yw;
    white_pixels_x[wi] = xb;
    white_pixels_y[wi] = yb;
  }
  else{
    P1 = P0;
    C1 = C0;
    bin1 = bin0;
    sample.assign(xb,yb,true);
  }
  
  log_DoS[bin1] += f;
}

void PhaseSpace::EVOLVE_classic_20times(){ 
  for(unsigned i=0; i<20; i++){
    EVOLVE_classic();
    log_DoS_all[i]=log_DoS;
    reset();
  }
}

void PhaseSpace::EVOLVE_classic(){
  // output_runtime_stats_DoS();
  std::stringstream stats_filename_stst;
  stats_filename_stst << prefix_of << "DoS_runtime_stats_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
		      << "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string stats_filename_st = stats_filename_stst.str();
  std::ofstream stats( stats_filename_st.c_str(), std::ofstream::app );
  
  stats << "# Runtime statistics" << std::endl
	<< "# " << nmbr_of_microstates << " microstates" << std::endl
	<< "# " << 2*log_DoS.size() << " inital steps" << std::endl
	<< "# f_mod = " << f_mod << " ; f_epsilon " << f_epsilon << " ; Hdev = " << hist_dev << " ; f_initial = " << f_initial << " ;" << std::endl
	<< "#" << std::setw(34) << "f" 
	<< std::setw(35) << "Time" << std::setw(35) << "Walltime" << std::setw(35) << "CPU-Time"
	<< std::setw(35) << "count" << std::setw(35) << "total_count" << std::endl;
  
  unsigned N_mirco_check = 16777215;
  
  double evolve_time=0;
  time_t begin;
  time(&begin);
  time_t now;
  
  unsigned long int count = 0;
  while(f > f_epsilon){
    if(count == 0)
      for(unsigned i = 0; i < 2*log_DoS.size(); i++){
	CHANGE_STATE();
	count++;
      }
    
    CHANGE_STATE();
    count++;
    
    if( (count&N_mirco_check) == 0 ){
      
      for(std::vector<unsigned>::size_type i = 0; i < hist.size(); i++)
	if(hist[i] > 0.1){
	  hist_min = hist[i];
	  break;
	}
      non_zero_hist=0;
      for(std::vector<unsigned>::size_type i = 0; i < hist.size(); i++)
	if(hist[i] > 0.1){
	  non_zero_hist++;
	  if(hist[i] < hist_min)
	    hist_min = hist[i];
	}
      hist_mean = count*1.0/non_zero_hist;
      
      if( hist_min > (hist_dev*hist_mean) ){
	hist_min = 0;
	hist = std::vector<unsigned> (n_bins_P * n_bins_C,0);
	
	// Normalize log_DoS
	double log_DoS_min = std::numeric_limits<double>::max();
	for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
	  if( (log_DoS[i] > 0.1) && (log_DoS[i] < log_DoS_min) )
	    log_DoS_min = log_DoS[i];
	log_DoS_min--;
	
	for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
	  if( log_DoS[i] > 0.1 )
	    log_DoS[i] -= log_DoS_min;
	
	total_count += count;
	time(&now);
	evolve_time=difftime(now,begin);
	char rtime[18];
	strftime (rtime, 18, "%Y%m%d.%H%M%S", gmtime(&now));

	int who = RUSAGE_SELF;
	struct rusage usage;
	getrusage(who, &usage);

	stats << std::setw(35) << f
	      << std::setw(35) << rtime << std::setw(35) << evolve_time << std::setw(35) << usage.ru_utime.tv_sec + usage.ru_stime.tv_sec
	      << std::setw(35) << count << std::setw(35) << total_count << std::endl;
	
	if(f < 5e-5)
	  {
	    // output_ammendment_log_DoS();
	    std::stringstream filename_stst;
	    filename_stst << prefix_of << "log_DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
			  << "_seed_" << original_seed << "_seedwght_" << original_seedwght << "_f_" << f << ".dat";
	    std::string filename_st = filename_stst.str();
	    std::ofstream Data( filename_st.c_str() );
	  
	    for(unsigned i = 0; i < log_DoS.size(); i++)
	      if(log_DoS[i] > 0){
		int P = 2*(i%n_bins_P);
		int C = i/n_bins_P - abs_min_C;
	      
		Data << std::setw(34) << P << std::setw(34) << C 
		     << std::setw(34) << std::setprecision(20) << log_DoS[i] << std::endl;
	      
		if(P==0){
		  std::cerr << "ERROR in assignment of P and C:" << std::endl;
		  std::cerr << std::setw(10) << "n_bins_P" << std::setw(10) << "n_bins_C" << std::setw(10) <<  "abs_min_C" << std::endl;
		  std::cerr << std::setw(10) << n_bins_P << std::setw(10) << n_bins_C << std::setw(10) <<  abs_min_C << std::endl;
		  std::cerr << "i: " << i << std::endl;
		
		  stats.close();
		  Data.close();
		  exit(-1);
		}
	    
	      }
	    Data.close();
	  
	    std::stringstream intermediate_DoS_stst;
	    intermediate_DoS_stst << prefix_of << "intermediate_log_DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
				  << "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
	    std::string intermediate_DoS_st = intermediate_DoS_stst.str();
	    std::ofstream intermediate_DoS( intermediate_DoS_st.c_str() );
	  
	    // black_pixels_x, black_pixels_y, white_pixels_x, white_pixels_y;
	    intermediate_DoS << "# f:" << std::endl
			     << std::setprecision(20) << f << std::endl
			     << "# total_count:" << std::endl
			     << std::setprecision(20) << total_count << std::endl;
	  
	    intermediate_DoS << "# black_pixels_x:" << std::endl
			     << black_pixels_x[0];
	    for(unsigned ai = 1; ai < A; ai++)
	      intermediate_DoS << " " << black_pixels_x[ai];
	    intermediate_DoS << std::endl
			     << "# black_pixels_y:" << std::endl
			     << black_pixels_y[0];
	    for(unsigned ai = 1; ai < A; ai++)
	      intermediate_DoS << " " << black_pixels_y[ai];
	    intermediate_DoS << std::endl
			     << "# white_pixels_x:" << std::endl
			     << white_pixels_x[0];
	    for(unsigned wi = 1; wi < W; wi++)
	      intermediate_DoS << " " << white_pixels_x[wi];
	    intermediate_DoS << std::endl
			     << "# white_pixels_y:" << std::endl
			     << white_pixels_y[0];
	    for(unsigned wi = 1; wi < W; wi++)
	      intermediate_DoS << " " << white_pixels_y[wi];
	    intermediate_DoS << std::endl;

	    intermediate_DoS << "# non_zero_hist:" << std::endl
			     << std::setprecision(20) << non_zero_hist << std::endl
			     << "# i log_DoS[i]:" << std::endl;
	  
	    for(unsigned i = 0; i < log_DoS.size(); i++)
	      if(log_DoS[i] > 0){
		intermediate_DoS << i << " " << std::setprecision(20) << log_DoS[i] << std::endl;  
	      }
	  
	    intermediate_DoS.close();
	  
	  }
	
	f *= f_mod;
	
	if(verbose)
	  std::cout << "A flat distribution is reached at count = " 
		    << count*1e-6 << " * 10^6; setting f = " << f << std::endl;
	count = 0;
      }
      else{
	if(count > 30000000000){
	  hist_min = 0;
	  hist = std::vector<unsigned> (n_bins_P * n_bins_C,0);
	  
	  // Normalize log_DoS
	  double log_DoS_min = std::numeric_limits<double>::max();
	  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
	    if( (log_DoS[i] > 0.1) && (log_DoS[i] < log_DoS_min) )
	      log_DoS_min = log_DoS[i];
	  log_DoS_min--;
	
	  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
	    if( log_DoS[i] > 0.1 )
	      log_DoS[i] -= log_DoS_min;
	  
	  total_count += count;
	  time(&now);
	  evolve_time=difftime(now,begin);
	  char rtime[18];
	  strftime (rtime, 18, "%Y%m%d.%H%M%S", gmtime(&now));
	  
	  int who = RUSAGE_SELF;
	  struct rusage usage;
	  getrusage(who, &usage);
	  
	  stats << std::setw(35) << f
		<< std::setw(35) << rtime << std::setw(35) << evolve_time << std::setw(35) << usage.ru_utime.tv_sec + usage.ru_stime.tv_sec
		<< std::setw(35) << count << std::setw(35) << total_count << std::endl;
	  
	  std::stringstream intermediate_DoS_stst;
	  intermediate_DoS_stst << prefix_of << "intermediate_log_DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
				<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
	  std::string intermediate_DoS_st = intermediate_DoS_stst.str();
	  std::ofstream intermediate_DoS( intermediate_DoS_st.c_str() );
	  
	  // black_pixels_x, black_pixels_y, white_pixels_x, white_pixels_y;
	  intermediate_DoS << "# f:" << std::endl
			   << std::setprecision(20) << f << std::endl
			   << "# total_count:" << std::endl
			   << std::setprecision(20) << total_count << std::endl;
	  
	  intermediate_DoS << "# black_pixels_x:" << std::endl
			   << black_pixels_x[0];
	  for(unsigned ai = 1; ai < A; ai++)
	    intermediate_DoS << " " << black_pixels_x[ai];
	  intermediate_DoS << std::endl
			   << "# black_pixels_y:" << std::endl
			   << black_pixels_y[0];
	  for(unsigned ai = 1; ai < A; ai++)
	    intermediate_DoS << " " << black_pixels_y[ai];
	  intermediate_DoS << std::endl
			   << "# white_pixels_x:" << std::endl
			   << white_pixels_x[0];
	  for(unsigned wi = 1; wi < W; wi++)
	    intermediate_DoS << " " << white_pixels_x[wi];
	  intermediate_DoS << std::endl
			   << "# white_pixels_y:" << std::endl
			   << white_pixels_y[0];
	  for(unsigned wi = 1; wi < W; wi++)
	    intermediate_DoS << " " << white_pixels_y[wi];
	  intermediate_DoS << std::endl;

	  intermediate_DoS << "# non_zero_hist:" << std::endl
			   << std::setprecision(20) << non_zero_hist << std::endl
			   << "# i log_DoS[i]:" << std::endl;
	  
	  for(unsigned i = 0; i < log_DoS.size(); i++)
	    if(log_DoS[i] > 0){
	      intermediate_DoS << i << " " << std::setprecision(20) << log_DoS[i] << std::endl;  
	    }
	  
	  intermediate_DoS.close();
	  
	  if(verbose)
	    std::cout << "No flat distribution is reached at count = " 
		      << count*1e-9 << " * 10^9; restetting histogram at f =" << f << std::endl;
	  count = 0;
	}
      }
      
      if( verbose ){
	std::cout << "count = " << count*1e-6 << " * 10^6; f = " << f << std::endl;
	display_sample();
      }
      
    } // if( count%N_mirco_check == 0 )
    
  } // while
}

unsigned PhaseSpace::evaluate_NmbrMacrostates() const{
  unsigned NMacro=0;
  
  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
    if(log_DoS[i] > 0.1)
      NMacro++;
  
  return NMacro;
}

void PhaseSpace::EVOLVE(const unsigned &slowdown){
  for(unsigned i = 0; i < 10000*N*N; i++)
    CHANGE_STATE_wo_hist();
  unsigned long int count = 0;
  double time = (10000*N*N+1)*timescale;
  
  unsigned sweep=1000/timescale;
  do{
    CHANGE_STATE();
    count++;
    
    if(count%sweep == 0){
      if(verbose)
	display_sample();
      for(std::vector<unsigned>::size_type i = 0; i < hist.size(); i++)
	if(hist[i] > 0.1){
	  hist_min = hist[i];
	  break;
	}
      non_zero_hist=0;
      for(std::vector<unsigned>::size_type i = 0; i < hist.size(); i++)
	if(hist[i] > 0.1){
	  non_zero_hist++;
	  if(hist[i] < hist_min)
	    hist_min = hist[i];
	}
      hist_mean = count*1.0/non_zero_hist;
      
      if( hist_min > (hist_dev*hist_mean) ){
	hist_min = 0;
	hist = std::vector<unsigned> (n_bins_P * n_bins_C,0);
		
	f *= f_mod;
	
	if(verbose){
	  std::cout << "A flat distribution is reached at count = " << count*1e-6 << " * 10^6; set f = " << f << std::endl;
	  if(f*time > slowdown)
	    std::cout << "Time = " << time*1e-4 << " * 10^4; f*time = " << f*time << " > " << slowdown << "; " 
		      << "keeping to classic Wang-Landau algorithm" << std::endl;
	  else
	    std::cout << "Time = " << time*1e-4 << " * 10^4; f*time = " << f*time << " <= " << slowdown << "; " 
		      << "to avoid saturation of error: new scaling f=1/t" << std::endl;
	}
	count=0;
      }//if( hist_min > (hist_dev*hist_mean) )
      else{
	time+=timescale;
      }
      
    }//if(count%sweep == 0)
    else{
      time+=timescale;
    }
  }while(f*time > slowdown);
  
  while(f > f_epsilon){
    f=slowdown/time;
    CHANGE_STATE_wo_hist();
    time+=timescale;
    
    if( unsigned(time/timescale)%(2*sweep) == 0 && verbose ){
      std::cout << "Time = " << time*1e-4 << " * 10^4; f = " << f << ";" << std::endl;
      display_sample();
    }
  }
  // MISSING: decide about parameters
}

void PhaseSpace::EVOLVE_analytic(){
  if(A <= N2/2){
    // initialize sample again
    for(unsigned ni = 0; ni < (N+2)*(N+2); ni++)
      sample.assign(ni,false);
    
    std::vector<unsigned> black_pixels_pos(A,0);
    for(unsigned ni = 0; ni < A; ni++){
      black_pixels_pos[ni] = ni;
    
      sample.assign(ni/N+1,ni%N+1,true);
    }
  
    std::cout << "Analytic determination of DoS" << std::endl
	      << "Initial configuration:" << std::endl;
    display_sample();
    
    int A0 = area_wbc_pix(sample);
    if(A0 != int(A*8)){
      std::cerr << "ERROR: initialization failed, area is " << A0 << ", compare given A = " << A << ";" << std::endl;
      exit(-1);
    }
    
    P1 = perimeter_wbc_pix(sample)/16;
    C1 = euler_wbc_pix(sample)/8;
    bin1 = n_bin(P1,C1);
    
    DoS[bin1]++;
    
    // Change to all other microstates
    if(A > 0){
      const unsigned j_max = A-1;
      const unsigned N2minusA = N2-A;
      unsigned j = j_max;
      unsigned original_site = 0;
      while(black_pixels_pos.front() < N2-A){
	
	j=j_max;
	while(black_pixels_pos[j] >= N2minusA+j)
	  j--;
	
	original_site = black_pixels_pos[j];
	black_pixels_pos[j]++;
	xb = original_site/N+1;
	yb = original_site%N+1;
	
	xw = black_pixels_pos[j]/N+1;
	yw = black_pixels_pos[j]%N+1;
	
	change_macrostate(xb, yb, xw, yw);
	bin1 = n_bin(P1,C1);
	sample.assign(xw,yw,true);
	
	for(unsigned ni = j+1; ni < A; ni++){
	  original_site = black_pixels_pos[ni];
	  black_pixels_pos[ni] = black_pixels_pos[j] + ni - j;
	  if(original_site != black_pixels_pos[ni]){
	    xb = original_site/N+1;
	    yb = original_site%N+1;
    
	    xw = black_pixels_pos[ni]/N+1;
	    yw = black_pixels_pos[ni]%N+1;
	  
	    change_macrostate(xb, yb, xw, yw);
	    bin1 = n_bin(P1,C1);
	    sample.assign(xw,yw,true);
	  }
	}	
	
	DoS[bin1]++;
      }// all steps
    }
  }// equal or less black pixel than white pixel
  else{
    // initialize sample again
    for(unsigned xi = 1; xi < (N+1); xi++)
      for(unsigned yi = 1; yi < (N+1); yi++)
	sample.assign(xi,yi,true);
    
    std::vector<unsigned> white_pixels_pos(W,0);
    for(unsigned ni = 0; ni < W; ni++){
      white_pixels_pos[ni] = ni;
      
      sample.assign(ni/N+1,ni%N+1,false);
    }
  
    std::cout << "Analytic determination of DoS" << std::endl
	      << "Initial configuration:" << std::endl;
    display_sample();
    
    int A0 = area_wbc_pix(sample);
    if(A0 != int(A*8)){
      std::cerr << "ERROR: initialization failed, area is " << A0 << ", compare given A = " << A << ";" << std::endl;
      exit(-1);
    }
  
    P1 = perimeter_wbc_pix(sample)/16;
    C1 = euler_wbc_pix(sample)/8;
    bin1 = n_bin(P1,C1);
    
    DoS[bin1]++;
        
    // Change to all other microstates
    if(A < N2){
      const unsigned j_max = W-1;
      const unsigned N2minusW = N2-W;
      unsigned j = j_max;
      unsigned original_site = 0;
      while(white_pixels_pos.front() < N2-W){
    
	j=j_max;
	while(white_pixels_pos[j] >= N2minusW+j)
	  j--;
      
	original_site = white_pixels_pos[j];
	white_pixels_pos[j]++;
	xb = original_site/N+1;
	yb = original_site%N+1;
    
	xw = white_pixels_pos[j]/N+1;
	yw = white_pixels_pos[j]%N+1;
	
	change_macrostate(xw, yw, xb, yb);
	bin1 = n_bin(P1,C1);
	sample.assign(xb,yb,true);
	
	for(unsigned ni = j+1; ni < W; ni++){
	  original_site = white_pixels_pos[ni];
	  white_pixels_pos[ni] = white_pixels_pos[j] + ni - j;
	  if(original_site != white_pixels_pos[ni]){
	    xb = original_site/N+1;
	    yb = original_site%N+1;
    
	    xw = white_pixels_pos[ni]/N+1;
	    yw = white_pixels_pos[ni]%N+1;
	    
	    change_macrostate(xw, yw, xb, yb);
	    bin1 = n_bin(P1,C1);
	    sample.assign(xb,yb,true);
	  }
	}
	
	DoS[bin1]++;
      }// all steps
    }
  }

}

void PhaseSpace::NORMALIZE(){
  
  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
    if(log_DoS[i] > 0)
      DoS[i] = exp(log_DoS[i]);
  
}

void PhaseSpace::NORMALIZE_with_total_sum(){
  
  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
    if(log_DoS[i] > 0)
      DoS[i] = exp(log_DoS[i]);
  
  double sum = 0;
  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
    sum += DoS[i];
  
  for(std::vector<double>::size_type i = 0; i < log_DoS.size(); i++)
    DoS[i] *= nmbr_of_microstates/sum;
  
}

void PhaseSpace::output_DoS(const bool &numerical){
  std::stringstream filename_stst;
  filename_stst << prefix_of << "DoS_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
		<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "# Structure distribution" << std::endl
       << "# f_initial: " << f_initial << std::endl
       << "# f_mod: " << f_mod << std::endl
       << "# f_epsilon: " << f_epsilon << std::endl
       << "# hist_dev: " << hist_dev << std::endl;
  if(numerical) 
    Data << "# result derived via numerical simulation with Wang-Landau algorithm" << std::endl;
  else
    Data << "# result derived analytically via brute force approach" << std::endl;
  Data << "#" << std::setw(33) << "P" << std::setw(34) << "C" << std::setw(34) << "DoS" 
       << std::endl;
  
  for(unsigned i = 0; i < DoS.size(); i++){
    int P = 2*(i%n_bins_P);
    int C = i/n_bins_P - abs_min_C;
    
    Data << std::setw(34) << P << std::setw(34) << C 
	 << std::setw(34) << std::setprecision(20) << DoS[i] << std::endl;  
  }
  
  Data.close();
}

void PhaseSpace::output_ammendment_DoS(){
  std::stringstream filename_stst;
  filename_stst << prefix_of << "DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
		<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  for(unsigned i = 0; i < DoS.size(); i++)
    if(DoS[i] > 0){
      int P = 2*(i%n_bins_P);
      int C = i/n_bins_P - abs_min_C;
      
      Data << std::setw(34) << P << std::setw(34) << C
	   << std::setw(34) << std::setprecision(20) << log_DoS[i]
	   << std::setw(34) << std::setprecision(20) << DoS[i] << std::endl;  
    }
  
  Data.close();
}

void PhaseSpace::output_ammendment_log_DoS(){
  std::stringstream filename_stst;
  filename_stst << prefix_of << "log_DoS_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
		<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  for(unsigned i = 0; i < log_DoS.size(); i++)
    if(log_DoS[i] > 0){
      int P = 2*(i%n_bins_P);
      int C = i/n_bins_P - abs_min_C;
      
      Data << std::setw(34) << P << std::setw(34) << C 
	   << std::setw(34) << std::setprecision(20) << log_DoS[i] << std::endl;  
    }
  
  Data.close();
}

void PhaseSpace::output_ammendment_log_DoS_all(){
  std::stringstream filename_stst;
  filename_stst << prefix_of << "log_DoS_all_ammendment_A_" << A << "_PC_" << "pix" << "_" << "wbc" << "_" << N << "x" << N
		<< "_seed_" << original_seed << "_seedwght_" << original_seedwght << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  for(unsigned i = 0; i < log_DoS.size(); i++){
    bool non_zero=false;
    for(std::vector<unsigned>::size_type run = 0; run < log_DoS_all.size(); run++)
      if(log_DoS_all[run][i] > 0){
	non_zero=true;
	break;
      }
    if(non_zero){
      int P = 2*(i%n_bins_P);
      int C = i/n_bins_P - abs_min_C;
      
      Data << std::setw(34) << P << std::setw(34) << C;
      for(std::vector<unsigned>::size_type run = 0; run < log_DoS_all.size(); run++)
	Data << std::setw(34) << std::setprecision(20) << log_DoS_all[run][i];
      Data << std::endl;
    }
  }
  
  Data.close();
}

void PhaseSpace::display_sample() const{
  for(unsigned i = 0; i < N+4; i++)
    if(i == (N-1))
      std::cout << "==";
    else
      std::cout << "===";
  std::cout << std::endl;
  for(unsigned yi = 0; yi < N+2; yi++){
    std::cout << "|| ";
    for(unsigned xi = 0; xi < N+2; xi++)
      if(sample.call(xi,N+1-yi))
	if((xi < (N+1)) && (xi > 0) && (yi < (N+1)) && (yi > 0))
	  std::cout << "<> ";
	else
	  std::cout << ">< ";     
      else
	if((xi < (N+1)) && (xi > 0) && (yi < (N+1)) && (yi > 0))
	  std::cout << "__ ";
	else
	  std::cout << "   "; 
    std::cout << "||" << std::endl;
  }
  for(unsigned i = 0; i < N+4; i++)
    if(i == (N+1))
      std::cout << "==";
    else
      std::cout << "===";
  std::cout << std::endl;
  
}

void PhaseSpace::change_macrostate(const unsigned &xb_, const unsigned &yb_, const unsigned &xw_, const unsigned &yw_){
  
  if( sample.call(xb_,yb_+1) && sample.call(xb_,yb_-1) )
    P1 += 1;    
  else if( !sample.call(xb_,yb_+1) && !sample.call(xb_,yb_-1) )
    P1 -= 1;
  
  if( sample.call(xb_+1,yb_) && sample.call(xb_-1,yb_) )
    P1 += 1;    
  else if( !sample.call(xb_+1,yb_) && !sample.call(xb_-1,yb_) )
    P1 -= 1;
  
  C1 += change_in_C_b2w[sample.call(xb_-1,yb_+1)*pow_of_2[0] + sample.call(xb_,yb_+1)*pow_of_2[1] + sample.call(xb_+1,yb_+1)*pow_of_2[2] + sample.call(xb_-1,yb_)*pow_of_2[3] + sample.call(xb_+1,yb_)*pow_of_2[4] + sample.call(xb_-1,yb_-1)*pow_of_2[5] + sample.call(xb_,yb_-1)*pow_of_2[6] + sample.call(xb_+1,yb_-1)*pow_of_2[7]];
  
  sample.assign(xb_,yb_,false);
  
  if( sample.call(xw_,yw_+1) && sample.call(xw_,yw_-1) )
    P1 -= 1;    
  else if( !sample.call(xw_,yw_+1) && !sample.call(xw_,yw_-1) )
    P1 += 1;
  
  if( sample.call(xw_+1,yw_) && sample.call(xw_-1,yw_) )
    P1 -= 1;    
  else if( !sample.call(xw_+1,yw_) && !sample.call(xw_-1,yw_) )
    P1 += 1;
  
  C1 += change_in_C_w2b[sample.call(xw_-1,yw_+1)*pow_of_2[0] + sample.call(xw_,yw_+1)*pow_of_2[1] + sample.call(xw_+1,yw_+1)*pow_of_2[2] + sample.call(xw_-1,yw_)*pow_of_2[3] + sample.call(xw_+1,yw_)*pow_of_2[4] + sample.call(xw_-1,yw_-1)*pow_of_2[5] + sample.call(xw_,yw_-1)*pow_of_2[6] + sample.call(xw_+1,yw_-1)*pow_of_2[7]];
  
  assert( (P1 >= n_bins_P) || (C1 <= min_C) || (C1 >= max_C) );
  /*
    if( (P1 >= n_bins_P) || (C1 <= min_C) || (C1 >= max_C) ){
    std::cerr << "ERROR in choice of the ranges of P and C:" << std::endl;
    std::cerr << std::setw(10) << "P_max = " << std::setw(10) << n_bins_P << std::endl;
    std::cerr << std::setw(10) << "P1    = " << std::setw(10) << P1 << std::endl;
    std::cerr << std::setw(10) << "C_max = " << std::setw(10) << max_C << std::endl;
    std::cerr << std::setw(10) << "C_min = " << std::setw(10) << min_C << std::endl;
    std::cerr << std::setw(10) << "C1    = " << std::setw(10) << C1 << std::endl;
    
    exit(-1);
    }
  */
}

unsigned PhaseSpace::n_bin(const int &P, const int &C) const{
  return P + n_bins_P*(C + abs_min_C);
}

void PhaseSpace::evaluate_single_deviations(const unsigned &start_seed, const unsigned &start_seedwght,
					    const unsigned &end_seed, const unsigned &end_seedwght,
					    const unsigned &bins, const double &max_dev){
  // Choose PhaseSpace Numerics
  unsigned stdseed_tmp = 101;
  unsigned stdseedwgth_tmp = 1001;
  if(N == 7){
    stdseed_tmp = 12;
    stdseedwgth_tmp = 24;
  }
  else if(N == 8){
    stdseed_tmp = 17;
    stdseedwgth_tmp = 23;
  }  
  const unsigned stdseed = stdseed_tmp;
  const unsigned stdseedwgth = stdseedwgth_tmp;
  
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  const double max_dev_check = max_dev-0.001;
  const double dev_width=2*max_dev/bins;
  
  // Read in DoS
  {
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << N << "x" << N << "/DoS_ammendment_A_"
	      << A << "_PC_pix_wbc_" << N << "x" << N << "_seed_" << stdseed << "_seedwght_" << stdseedwgth << ".dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << N << "x" << N << "/DoS_ammendment_A_"
		<< A << "_PC_pix_wbc_" << N << "x" << N << "_seed_" << stdseed << "_seedwght_" << stdseedwgth << ".dat"
		<< std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
    
      double log_DoS = 0;
      double DoS = 0;
      structure_distribution >> log_DoS;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(A,P,C,DoS));
    }
    structure_distribution.close();
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  std::vector< std::vector<double> > dev(nmbr_of_macrostates,std::vector<double>(bins,0));
  double norm = 1.0/(end_seed - start_seed + 1)/(end_seedwght - start_seedwght + 1)/dev_width;
  std::vector<double> delta(bins,0);
  double delta_0=dev_width/2.0 - max_dev;
  for(unsigned i = 0; i < bins; i++)
    delta[i] = delta_0 + i*dev_width;
  
  for(unsigned sd = start_seed; sd <= end_seed; sd++)
    for(unsigned sdw = start_seedwght; sdw <= end_seedwght; sdw++){
      std::stringstream infile_st;
      infile_st << "input/DoS_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		<< "_seed_" << start_seed << "_" << end_seed << "_seedwght_" << start_seedwght << "_" << end_seedwght
		<< "/DoS_ammendment_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		<< "_seed_" << sd << "_seedwght_" << sdw << ".dat";
      std::string infile = infile_st.str();
      std::ifstream structure_distribution((infile).c_str ());
      if(structure_distribution.fail()){
	std::cerr << "ERROR: ifstream failed to read " << "input/DoS_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		  << "_seed_" << start_seed << "_" << end_seed << "_seedwght_" << start_seedwght << "_" << end_seedwght
		  << "/DoS_ammendment_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		  << "_seed_" << sd << "_seedwght_" << sdw << ".dat"
		  << std::endl;
	exit(-1);
      }
      unsigned i = 0;
      while(!structure_distribution.eof()){
	int P(0), C(0);
	structure_distribution >> P;
	structure_distribution >> C;
    
	double log_DoS = 0;
	double DoS = 0;
	structure_distribution >> log_DoS;
	structure_distribution >> DoS;
	if(  DoS && (P == DoS_APC.at(i).P()) && (C == DoS_APC.at(i).C()) ){
	  double single_dev = (DoS_APC[i].cl()-DoS)/DoS_APC.at(i).cl();
	  if(fabs(single_dev) > max_dev_check){
	    std::cerr << "ERROR: deviation exceeds maximum" << std::endl;
	    exit(-1);
	  }
	  unsigned bin_pos = (single_dev+max_dev)/dev_width;
	  (dev[i][bin_pos]) += norm;
	}
	else if( DoS ){
	  std::cerr << "ERROR: unknown structure: seed = " << sd << std::endl
		    << P << std::endl
		    << DoS_APC.at(i).P() << std::endl
		    << C << std::endl
		    << DoS_APC.at(i).C() << std::endl;
	  exit(-1);
	}
	i++;
      }
      structure_distribution.close();
    }
    
  std::stringstream filename_stst;
  filename_stst << prefix_of << "DoS_distribution_of_single_deviations_" << N << "x" << N << "_A_" << A
		<< "_seed_" << start_seed << "_" << end_seed << "_seedwght_" << start_seedwght << "_" << end_seedwght << "_" << bins << "_bins.dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "#" << std::setw(14) << "delta";
  for(unsigned m = 0; m < nmbr_of_macrostates; m++)
    Data << std::setw(15) << DoS_APC[m].cl();
  Data << std::endl;    
  for(unsigned i = 0; i < bins; i++){
    Data << std::setw(15) << delta[i];
    for(unsigned m = 0; m < nmbr_of_macrostates; m++)
      Data << std::setw(15) << dev[m][i];
    Data << std::endl;
  }
  
  Data.close();
  
}

void PhaseSpace::evaluate_all_maximum_deviations(const unsigned &start_seed, const unsigned &start_seedwght,
						 const unsigned &end_seed, const unsigned &end_seedwght,
						 const unsigned &end_actual_seed, const unsigned &end_actual_seedwght){
  
  // Choose PhaseSpace Numerics
  unsigned stdseed_tmp = 101;
  unsigned stdseedwgth_tmp = 1001;
  if(N == 7){
    stdseed_tmp = 12;
    stdseedwgth_tmp = 24;
  }
  else if(N == 8){
    stdseed_tmp = 17;
    stdseedwgth_tmp = 23;
  }  
  const unsigned stdseed = stdseed_tmp;
  const unsigned stdseedwgth = stdseedwgth_tmp;
  
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  std::vector<double> DoS_vec; DoS_vec.clear();
  
  // Read in DoS
  {
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << N << "x" << N << "/DoS_ammendment_A_"
	      << A << "_PC_pix_wbc_" << N << "x" << N << "_seed_" << stdseed << "_seedwght_" << stdseedwgth << ".dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << N << "x" << N << "/DoS_ammendment_A_"
		<< A << "_PC_pix_wbc_" << N << "x" << N << "_seed_" << stdseed << "_seedwght_" << stdseedwgth << ".dat"
		<< std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double log_DoS = 0;
      double DoS = 0;
      structure_distribution >> log_DoS;
      structure_distribution >> DoS;
      if(DoS){
	DoS_APC.push_back(stateclassifier(A,P,C,DoS));
	DoS_vec.push_back(DoS);
      }
    }
    structure_distribution.close();
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  std::vector<double> max_dev_vec(nmbr_of_macrostates,std::numeric_limits<double>::min());
  std::vector<double> min_dev_vec(nmbr_of_macrostates,std::numeric_limits<double>::max());
  
  for(unsigned sd = start_seed; sd <= end_actual_seed; sd++)
    for(unsigned sdw = start_seedwght; sdw <= end_actual_seedwght; sdw++){
      std::stringstream infile_st;
      infile_st << "input/DoS_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		<< "_seed_" << start_seed << "_" << end_seed << "_seedwght_" << start_seedwght << "_" << end_seedwght
		<< "/DoS_ammendment_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		<< "_seed_" << sd << "_seedwght_" << sdw << ".dat";
      std::string infile = infile_st.str();
      std::ifstream structure_distribution((infile).c_str ());
      if(structure_distribution.fail()){
	std::cerr << "ERROR: ifstream failed to read " << "input/DoS_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		  << "_seed_" << start_seed << "_" << end_seed << "_seedwght_" << start_seedwght << "_" << end_seedwght
		  << "/DoS_ammendment_A_" << A << "_PC_pix_wbc_" << N << "x" << N
		  << "_seed_" << sd << "_seedwght_" << sdw << ".dat"
		  << std::endl;
	exit(-1);
      }
      unsigned i = 0;
      while(!structure_distribution.eof()){
	int P(0), C(0);
	structure_distribution >> P;
	structure_distribution >> C;
    
	double log_DoS = 0;
	double DoS = 0;
	structure_distribution >> log_DoS;
	structure_distribution >> DoS;
	if(  DoS && (P == DoS_APC.at(i).P()) && (C == DoS_APC.at(i).C()) ){
	  double relative_error=(DoS_APC[i].cl()-DoS)/DoS_APC.at(i).cl();
	  if(max_dev_vec[i] < relative_error)
	    max_dev_vec[i] = relative_error;
	  if(min_dev_vec[i] > relative_error)
	    min_dev_vec[i] = relative_error;
	}
	else if( DoS ){
	  std::cerr << "ERROR: unknown structure: seed = " << sd << std::endl
		    << P << std::endl
		    << DoS_APC.at(i).P() << std::endl
		    << C << std::endl
		    << DoS_APC.at(i).C() << std::endl;
	  exit(-1);
	}
	i++;
      }
      structure_distribution.close();
    }
  
  std::stringstream filename_stst;
  filename_stst << prefix_of << "DoS_distribution_of_maximum_deviations_" << N << "x" << N << "_A_" << A
		<< "_seed_" << start_seed << "_" << end_actual_seed
		<< "_seedwght_" << start_seedwght << "_" << end_actual_seedwght << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "#" << std::setw(14) << "DoS"
       << std::setw(15) << "minimum deviation"
       << std::setw(15) << "maximum deviation"
       << std::endl;
  
  for(unsigned i = 0; i < nmbr_of_macrostates; i++){
    Data << std::setw(15) << DoS_vec[i]
	 << std::setw(15) << min_dev_vec[i]
	 << std::setw(15) << max_dev_vec[i]
	 << std::endl;
  }
  
  Data.close();
  
}

