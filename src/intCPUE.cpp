// Integrated CPUE standardization algorithm
#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#define _USE_MATH_DEFINES
#include <cmath>

// List Of Sparse Matrices (LOSM): read a List from R into vector<SparseMatrix<> >
// following sdmTMB, taken from kaskr, https://github.com/kaskr/adcomp/issues/96
using namespace Eigen;
using namespace tmbutils;
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){
    this->resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};

// Bias corrected log-normal
template<class Type>
Type dlnorm_bc(const Type& x, const Type& meanlog, const Type& sdlog, int give_log = 0) {
    Type adjusted_meanlog = meanlog - pow(sdlog, 2) / 2;
    Type logres = dnorm(log(x), adjusted_meanlog, sdlog, true) - log(x);
	
    if (give_log)
        return logres;
    else
        return exp(logres);
}

// PC prios on range and sigma
template <class Type>
Type pc_prior_matern(        
  Type range, 
  Type sigma,   
  Type matern_range,   // r₀       — the user‐supplied “Pr(range<r₀)=range_prob” cutoff  
  Type matern_sigma,   // σ₀       — the user‐supplied “Pr(σ>σ₀)=SD_prob” cutoff
  Type range_prob,     // Pr(range<r₀)  
  Type sigma_prob,     // Pr(σ>σ₀)  
  int  give_log=0,     // if 1, return log‐density; if 0, return density  
  int  share_range=0   // if 1, omit the range‐prior term (useful in hierarchical models) 
){
  // dimension and half‐dimension
  Type d     = 2.;
  Type dhalf = d / 2.;

  // compute λ parameters
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(sigma_prob) / matern_sigma;

  // log‐density of PC‐prior for range (unless share_range==1)
  Type range_ll = log(dhalf) 
                + log(lam1) 
                + (-1.0 - dhalf) * log(range)
                - lam1 * pow(range, -dhalf);

  // log‐density of PC‐prior for σ
  Type sigma_ll = log(lam2) 
                - lam2 * sigma;

  // combine
  Type penalty = sigma_ll;
  if(!share_range) penalty += range_ll;

  if(give_log)   return penalty;  
  else           return exp(penalty);
}


// Main TMB objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Dimensions
  DATA_INTEGER(n_i); // Number of observations (stacked across all years)
  DATA_INTEGER(n_t); // Number of time-indices
  DATA_INTEGER(n_v); // Number of vessels (i.e., levels for the factor explaining overdispersion)
  DATA_INTEGER(n_f); // Number of flags (e.g., 0 = the reference, 1 = other 1, 2 = other 2, etc.)
  DATA_INTEGER(n_g); // Number of extrapolation-grid cells
  
  // Data
  DATA_VECTOR(b_i); // Response (Positive CPUE, biomass / effort) for each observation
  DATA_IVECTOR(e_i); // Response (0 or 1) for each observation
  DATA_IVECTOR(t_i); // Time index for each observation
  DATA_IVECTOR(v_i); // Vessel for each observation
  DATA_IVECTOR(f_i); // Flag (gear) type (e.g., 0 = the reference, 1 = other 1, 2 = other 2, etc.) for each observation  
  
  DATA_VECTOR(area_g); // Area for each extrapolation-grid cell
  
  // Projection matrices from knots s to data i or extrapolation-grid cells x
  DATA_SPARSE_MATRIX(A_is); // Project vertices to samples 
  DATA_SPARSE_MATRIX(A_gs); // Project vertices to integration points
  DATA_IMATRIX(Ais_ij);
  DATA_VECTOR(Ais_x);
  
  // PC priors on SPDE hyperparamters
  DATA_SCALAR(matern_range);    // r₀: e.g. diff(range(mesh$loc[,1]))/5
  DATA_SCALAR(range_prob);      // P(range < r₀)
  DATA_SCALAR(matern_sigma_0);  // σ₀ for omega
  DATA_SCALAR(matern_sigma_t);  // σ₀ for epsilon
  DATA_SCALAR(matern_sigma_flag);  // σ₀ for epsilon
  DATA_SCALAR(sigma_prob);      // P(σ > σ₀)
  
  // Smooth objects
  DATA_INTEGER(has_smooths);
  DATA_MATRIX(Xs);                 // n_i x K_smooth (linear part of smooths)
  DATA_STRUCT(Zs, LOSM_t);   // list of n_i x k_s sparse basis matrices
  DATA_IVECTOR(b_smooth_start);    // length = n_smooth, start index in b_smooth for each smooth block

  // Aniso SPDE objects
  DATA_STRUCT(spde, spde_aniso_t);

  // Parameters  
  PARAMETER(ln_sd); // SD for Log-normal distribution
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  
  PARAMETER(ln_range_1); // SPDE hyper-parameter
  PARAMETER(ln_sigma_0_1); // Time-constant SPDE hyper-parameter
  PARAMETER(ln_sigma_t_1); // Time-varying SPDE hyper-parameter
  PARAMETER(ln_sigma_flag_1); // SPDE hyper-parameter for differences in gear-specific catchability over space
  PARAMETER(ln_range_2); // SPDE hyper-parameter
  PARAMETER(ln_sigma_0_2); // Time-constant SPDE hyper-parameter
  PARAMETER(ln_sigma_t_2); // Time-varying SPDE hyper-parameter
  PARAMETER(ln_sigma_flag_2); // SPDE hyper-parameter for differences in gear-specific catchability over space

  // Encounter probability (Poisson)
  PARAMETER_VECTOR(ves_v_1); // Random vessel effect
  PARAMETER_VECTOR(yq_t_1); // Year effect
  PARAMETER_VECTOR(flag_f_1); // Differences in gear-specific catchability (that is, systematic difference)
  PARAMETER_MATRIX(flag_s_1); // Differences in gear-specific catchability over space [s,n_f-1]
  PARAMETER_VECTOR(omega_s_1); // Time-constant SPDE spatial field (over mesh nodes)
  PARAMETER_MATRIX(epsilon_st_1); // Time-varying SPDE spatial field (over mesh nodes) across time [s,c]
  PARAMETER(ves_ln_std_dev_1); // SD for random vessel effect
  PARAMETER(flag_ln_std_dev_1); // SD for gear-specific catchability (systematic difference)
  
  // Positive catch (Lognormal)
  PARAMETER_VECTOR(ves_v_2); // Random vessel effect 
  PARAMETER_VECTOR(yq_t_2);  // Year effect
  PARAMETER_VECTOR(flag_f_2); // Differences in gear-specific catchability
  PARAMETER_MATRIX(flag_s_2); // Differences in gear-specific catchability over space [s,n_f-1]
  PARAMETER_VECTOR(omega_s_2); // Time-constant SPDE spatial field (over mesh nodes) 
  PARAMETER_MATRIX(epsilon_st_2); // Time-varying SPDE spatial field (over mesh nodes) across time [s,c]
  PARAMETER(ves_ln_std_dev_2); // SD for random vessel effect
  PARAMETER(flag_ln_std_dev_2); // SD for gear-specific catchability

  // Smoothers (two components: encounter=0, positive=1)
  // bs: coefficients for Xs (unpenalized part / linear part of smooth)
  PARAMETER_MATRIX(bs);              // [K_smooth, 2]

  // b_smooth: coefficients for each Zs block stacked together
  PARAMETER_MATRIX(b_smooth);        // [sum_k, 2]

  // ln_smooth_sigma: SD for each smooth block's random-effect coefficients
  PARAMETER_MATRIX(ln_smooth_sigma); // [n_smooth, 2]
  
  // Global variables
  Type nll = 0;
  Type nll_prior = 0;
  Type nll_penalty = 0;
  Type sd = exp(ln_sd); 
  Type range_1      = exp(ln_range_1);
  Type sigma_0_1    = exp(ln_sigma_0_1);
  Type sigma_t_1    = exp(ln_sigma_t_1);
  Type sigma_flag_1 = exp(ln_sigma_flag_1);  
  Type ves_std_dev_1 = exp(ves_ln_std_dev_1); 
  Type flag_std_dev_1 = exp(flag_ln_std_dev_1); 
  Type range_2      = exp(ln_range_2);
  Type sigma_0_2    = exp(ln_sigma_0_2);
  Type sigma_t_2    = exp(ln_sigma_t_2);
  Type sigma_flag_2 = exp(ln_sigma_flag_2); 
  Type ves_std_dev_2 = exp(ves_ln_std_dev_2);
  Type flag_std_dev_2 = exp(flag_ln_std_dev_2); 
  
  // SPDE hyper transforms
  Type kappa_1 = sqrt(Type(8.0)) / range_1;
  Type tau_0_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_0_1);
  Type tau_t_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_t_1);
  Type tau_flag_1 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_1*sigma_flag_1);
  Type kappa_2 = sqrt(Type(8.0)) / range_2;
  Type tau_0_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_0_2);
  Type tau_t_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_t_2);
  Type tau_flag_2 = Type(1.0) / (Type(2.0)*sqrt(M_PI)*kappa_2*sigma_flag_2);
  
  // PC‐priors
  nll_prior -= pc_prior_matern(
    range_1, sigma_0_1,
    matern_range, matern_sigma_0,
    range_prob, sigma_prob,
    1,             // give_log = 1  (return log‐density)
    0              // share_range = 0 (penalize range)
  );
  
  nll_prior -= pc_prior_matern(
    range_1, sigma_t_1,
    matern_range, matern_sigma_t,
    range_prob, sigma_prob,
    1,  // give_log
    1   // share_range = 1 (skip range penalty here)
  );
  
  nll_prior -= pc_prior_matern(
    range_1, sigma_flag_1,
    matern_range, matern_sigma_flag,
    range_prob, sigma_prob,
    1,  // give_log
    1   // share_range = 1 (skip range penalty here)
  );
  
  nll_prior -= pc_prior_matern(
    range_2, sigma_0_2,
    matern_range, matern_sigma_0,
    range_prob, sigma_prob,
    1,             // give_log = 1  (return log‐density)
    0              // share_range = 0 (penalize range)
  );
  
  nll_prior -= pc_prior_matern(
    range_2, sigma_t_2,
    matern_range, matern_sigma_t,
    range_prob, sigma_prob,
    1,  // give_log
    1   // share_range = 1 (skip range penalty here)
  );
  
  nll_prior -= pc_prior_matern(
    range_2, sigma_flag_2,
    matern_range, matern_sigma_flag,
    range_prob, sigma_prob,
    1,  // give_log
    1   // share_range = 1 (skip range penalty here)
  );
  
  
  // Need to parameterize H matrix such that det(H)=1 (preserving volume) 
  // Note that H appears in (20) in Lindgren et al 2011
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  
  SparseMatrix<Type> Q_1 = Q_spde(spde, kappa_1, H); // Precision matrix
  SparseMatrix<Type> Q_2 = Q_spde(spde, kappa_2, H); // Precision matrix
  
  
  // Encounter probability
  // Time-constant SPDE (Omega)
  nll += SCALE(GMRF(Q_1), 1./tau_0_1)(omega_s_1);

  // Time-varying SPDE (Epsilon)
  nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(0)); // t=0 
  for(int t=1; t<n_t; t++){
	nll += SCALE(GMRF(Q_1), 1./tau_t_1)(epsilon_st_1.col(t) - epsilon_st_1.col(t-1));
  }
   
  // Random vessel effect
  for(int i=0; i<n_v; i++){
    nll -= dnorm(ves_v_1(i), Type(0.0), ves_std_dev_1, true);
  }
  
  // Differences in gear-specific catchability
  for(int i=0; i<n_f-1; i++){
    nll -= dnorm(flag_f_1(i), Type(0.0), flag_std_dev_1, true);
  }

  // Differences in gear-specific catchability over space
  for(int i=0; i<n_f-1; i++){
    nll += SCALE(GMRF(Q_1), 1./tau_flag_1)(flag_s_1.col(i));
  }
  
  
  // Positive catch
  // Time-constant SPDE (Omega)
  nll += SCALE(GMRF(Q_2), 1./tau_0_2)(omega_s_2);
  
  // Time-varying SPDE (Epsilon)
  nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(0)); //t=0 
  for(int t=1; t<n_t; t++){
	nll += SCALE(GMRF(Q_2), 1./tau_t_2)(epsilon_st_2.col(t) - epsilon_st_2.col(t-1));
  }
   
  // Random vessel effect
  for(int i=0; i<n_v; i++){
    nll -= dnorm(ves_v_2(i), Type(0.0), ves_std_dev_2, true);
  }

  // Differences in gear-specific catchability
  for(int i=0; i<n_f-1; i++){
    nll -= dnorm(flag_f_2(i), Type(0.0), flag_std_dev_2, true);
  }
  
  // Differences in gear-specific catchability over space
  for(int i=0; i<n_f-1; i++){
    nll += SCALE(GMRF(Q_2), 1./tau_flag_2)(flag_s_2.col(i));
  }
             
	
  // ---- Smooth contributions (eta_smooth_1, eta_smooth_2) + priors ----
  vector<Type> eta_smooth_1(n_i); // encounter part
  vector<Type> eta_smooth_2(n_i); // positive part
  eta_smooth_1.setZero();
  eta_smooth_2.setZero();

  if (has_smooths) {
    int n_smooth = b_smooth_start.size();

    // 1) penalized spline RE part: sum_s Zs[s] * beta_s
    for (int s = 0; s < n_smooth; s++) {
      int k_s = Zs(s).cols();          // #coeffs in this smooth block
      int start = b_smooth_start(s);   // start index in stacked b_smooth

      // encounter component (col 0)
      vector<Type> beta0(k_s);
      for (int j = 0; j < k_s; j++) {
        beta0(j) = b_smooth(start + j, 0);
        nll -= dnorm(beta0(j), Type(0.0), exp(ln_smooth_sigma(s, 0)), true);
      }
      eta_smooth_1 += Zs(s) * beta0;

      // positive component (col 1)
      vector<Type> beta1(k_s);
      for (int j = 0; j < k_s; j++) {
        beta1(j) = b_smooth(start + j, 1);
        nll -= dnorm(beta1(j), Type(0.0), exp(ln_smooth_sigma(s, 1)), true);
      }
      eta_smooth_2 += Zs(s) * beta1;
    }

    // 2) unpenalized / linear part: Xs * bs
    if (Xs.cols() > 0) {
      eta_smooth_1 += Xs * vector<Type>(bs.col(0));
      eta_smooth_2 += Xs * vector<Type>(bs.col(1));
    }
  }

	
  // Response
  vector<Type> s_effect_1 = A_is * omega_s_1;
  vector<Type> s_effect_2 = A_is * omega_s_2;
  vector<Type> st_effect_1(n_i);
  vector<Type> st_effect_2(n_i);
  vector<Type> flag_s_effect_1(n_i);
  vector<Type> flag_s_effect_2(n_i);

  st_effect_1.setZero();
  st_effect_2.setZero();
  flag_s_effect_1.setZero();
  flag_s_effect_2.setZero();
  
  for(int r=0; r<Ais_ij.rows(); r++){
    int i = Ais_ij(r, 0);  // observation
    int s = Ais_ij(r, 1);  // knot
	int t_id = t_i(i);
	int f_id = f_i(i);
	  
	st_effect_1(i) += Ais_x(r) * epsilon_st_1(s, t_id);
	st_effect_2(i) += Ais_x(r) * epsilon_st_2(s, t_id);
	
	if(f_id > 0){
	  flag_s_effect_1(i) += Ais_x(r) * flag_s_1(s, f_id-1);
	  flag_s_effect_2(i) += Ais_x(r) * flag_s_2(s, f_id-1);
	}
  }
  
  for(int i=0; i<n_i; i++){	
	int tid = t_i(i);
	int vid = v_i(i);
	int fid = f_i(i);
	
    Type ves_effect_1 = ves_v_1(vid); 
	Type yq_effect_1 = yq_t_1(tid);
	Type flag_effect_1 = (fid == 0) ? Type(0) : flag_f_1(fid-1);  // if fid=0，flag_effect_1 = 0
	
    Type ves_effect_2 = ves_v_2(vid); 
	Type yq_effect_2 = yq_t_2(tid);
	Type flag_effect_2 = (fid == 0) ? Type(0) : flag_f_2(fid-1);  // if fid=0，flag_effect_2 = 0
	
	Type eta1 = ves_effect_1 + yq_effect_1 + flag_effect_1 + flag_s_effect_1(i) + s_effect_1(i) + st_effect_1(i) + eta_smooth_1(i); 					   
    Type eta2 = ves_effect_2 + yq_effect_2 + flag_effect_2 + flag_s_effect_2(i) + s_effect_2(i) + st_effect_2(i) + eta_smooth_2(i); 
	  
    // Poisson-link; i.e., complementary log–log (cloglog)
    Type log_one_minus_p = -1.0*exp(eta1);
    Type logp = logspace_sub( Type(0.0), log_one_minus_p );

    // Positive
	Type logcat = eta1 + eta2 - logp;
   
	//Likelihood
	// Probability
	nll -= e_i(i) * logp + (1 - e_i(i)) * log_one_minus_p;
				 
	//Positive catch
	if(e_i(i) == 1){
	  nll -= dlnorm_bc( b_i(i), logcat, sd, true ); 
	}
  }
  
  
  //Projection
  vector<Type> s_effect_proj_1 = A_gs * omega_s_1;
  vector<Type> s_effect_proj_2 = A_gs * omega_s_2;
  matrix<Type> st_effect_proj_1 = A_gs * epsilon_st_1; // [g, t]
  matrix<Type> st_effect_proj_2 = A_gs * epsilon_st_2;
  
  matrix<Type> cpue_density(n_g, n_t);
  vector<Type> mu_total(n_t); // Area-weighted CPUE in weight
  vector<Type> link_total(n_t);
  
  mu_total.setZero();

  for(int t=0; t<n_t; t++){ 
	Type yq_effect_proj_1 = yq_t_1(t); 
	Type yq_effect_proj_2 = yq_t_2(t); 
	
    for(int g=0; g<n_g; g++){	
	  Type eta1_proj = yq_effect_proj_1 + s_effect_proj_1(g) + st_effect_proj_1(g, t);
	  Type eta2_proj = yq_effect_proj_2 + s_effect_proj_2(g) + st_effect_proj_2(g, t);
  	      
	  // CPUE in weight
	  Type cpue = exp(eta1_proj + eta2_proj);
	  cpue_density(g, t) = cpue;
	  
	  Type adds = cpue * area_g(g);
		
	  // Area-weighted CPUE
      mu_total(t) += adds;
    }
	
	link_total(t) = log(mu_total(t));
  }
  
  PARAMETER_VECTOR(eps_index);
  
  if (eps_index.size() > 0) {
    Type S;
    for (int t=0; t < n_t; t++) {
      S = mu_total(t);
      //S = newton::Tag(S); // Set lowrank tag on S = sum(exp(x))
      nll_penalty += eps_index(t) * S;
    }
  }
  
  nll += nll_prior;   
  nll += nll_penalty; 
  
  // Reporting;
  REPORT(nll_prior);
  REPORT(nll_penalty);
  REPORT(cpue_density);
  ADREPORT(link_total);
  return nll;
}