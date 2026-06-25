#include <RcppArmadillo.h>
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
using namespace Rcpp;
using namespace arma;
using namespace Numer;
// [[Rcpp::export]]
IntegerVector find_cutoff(double eff_lower,
                          int n_stage1,
                          int n_stage2,
                          double Cs2,
                          int cutoff_step 
                          ){
  double a0 = 0.1;
  IntegerVector n_stage1_seq = seq(0, n_stage1);
  IntegerVector n_total_seq = seq(0, (n_stage1 + n_stage2));
  NumericVector stage1_vec = as<NumericVector>(n_stage1_seq) + a0;
  NumericVector total_vec = as<NumericVector>(n_total_seq) + a0;
  int response1 = -1;
  int response2 = -1;
  for (int i=0; i<n_stage1; i++) {
    double tmp_preff = R::pbeta(eff_lower, stage1_vec[i], stage1_vec[(n_stage1-i)], true, false);
    
    if(tmp_preff < Cs2){
      response1 = std::max(0, (i-1));
      break;
    }
  }
  
  for (int i=0; i<(n_stage1 + n_stage2); i++) {
    double tmp_preff = R::pbeta(eff_lower, total_vec[i], total_vec[(n_stage1+n_stage2-i)], true, false);
    if(tmp_preff < Cs2){
      response2 = std::max(0, (i-1));
      break;
    }
  }
  
  IntegerVector output;
  if(response1==-1 | response2==-1){
    Rcerr << "Please input a proper Cs2.\n";
    output = IntegerVector::create(-1,-1);
  }else{
    int stage1_cutoff = response1 + cutoff_step; 
    int total_cutoff = response2; 
    output = IntegerVector::create(stage1_cutoff, total_cutoff);
  }
  return output;
}



DataFrame DGP_comb_s1(int Ndose,
                         int n,
                         NumericVector prob,
                         int type = 1) {
  
  NumericVector CDF(n + 1, 1.0);
  if(type == 1){
    for (int i=0; i<Ndose; i++) {
      NumericVector prod = pbinom(seq(0, n), n, prob[i]);
      CDF = CDF * prod;
    }
  }else if(type == 2){
    for (int i=0; i<Ndose; i++) {
      NumericVector prod = pbinom(seq(0, n), n, prob[i]);
      CDF = CDF * (1 - prod);
    }
    CDF = 1 - CDF;
  }
  CDF.push_front(0);
  NumericVector PMF = diff(CDF);
  DataFrame df = DataFrame::create(Named("m") = seq(0, n), 
                                   Named("PMF") = PMF );
  return df;
}

DataFrame DGP_comb_total(int Ndose, 
                         int Nstage1, 
                         int Nstage2, 
                         NumericVector prob_s1,
                         double prob_s2,
                         double eff_lower,
                         double Cs2,
                         int cutoff_step,
                         int type) {
  
  int Ntotal = Nstage1 + Nstage2; 
  DataFrame PMFstage1 = DGP_comb_s1(Ndose, Nstage1, prob_s1, type);
  NumericVector m = PMFstage1["m"];
  NumericVector PMF = PMFstage1["PMF"];
  
  IntegerVector cutoff = find_cutoff(eff_lower, Nstage1, Nstage2, Cs2, cutoff_step);
  NumericVector PMF_adm = clone(PMF);
  PMF_adm[m <= cutoff[0]] = 0; 
  
  auto calc_PMF = [&](int r) {
    int idx1 = std::min(r, Nstage1);
    int idx2 = std::max(0, r - Nstage1);
    NumericVector PMF_temp = PMF_adm[m <= idx1];
    NumericVector prod; 
    prod = dbinom(rev(seq(idx2, r)), Nstage2, prob_s2); 
    double res =  sum(PMF_temp * prod);
    return res;
  };
  
  NumericVector PMF_total = sapply(seq(0, Ntotal), calc_PMF);
  IntegerVector m_total =  seq(0, Ntotal);
  
  PMF_total[m_total <= cutoff[1]] = 0;
  
  DataFrame df = DataFrame::create(Named("m") = m_total,
                                   Named("PMF") = PMF_total );
  return df;
}

// [[Rcpp::export]]
double prob_inadm(int Ndose, 
                  int Nstage1, 
                  int Nstage2, 
                  NumericVector prob_s1,
                  double prob_s2,
                  double eff_lower,
                  double Cs2, 
                  int cutoff_step){
  IntegerVector cutoff = find_cutoff(eff_lower, Nstage1, Nstage2, Cs2, cutoff_step);
  
  DataFrame PMFstage1 = DGP_comb_s1(Ndose, Nstage1, prob_s1, 1);
  NumericVector m = PMFstage1["m"];
  NumericVector PMF_s1=PMFstage1["PMF"];
  double prob_s1_inadm = sum(PMF_s1[seq(0, cutoff[0])]);

  double PMF_temp;
  double log_prod;
  double prob_s2_inadm = 0;
  int upper_idx = std::min((cutoff[1] + 1), (Nstage1 + 1));
  for (int ii=(cutoff[0]+1); ii < upper_idx; ii++) {
    PMF_temp = log(PMF_s1[ii]);
    log_prod = R::pbinom((cutoff[1]-ii), Nstage2, prob_s2, true, true); 
    prob_s2_inadm =  prob_s2_inadm + exp(PMF_temp + log_prod);
  }

  return prob_s1_inadm + prob_s2_inadm; 
  
}


class BetaIntegrate: public Func
{
private:
  double a;
  double b;
  double a_obd;
  double b_obd;
public:
  BetaIntegrate(double a_, double b_,
                double a_obd_, double b_obd_) :
  a(a_), b(b_), a_obd(a_obd_), b_obd(b_obd_) {}

  double operator()(const double& x) const
  {
    double log_val = R::dbeta(x, a_obd, b_obd, true) + R::pbeta(x, a, b, true, true);
    return exp(log_val);
  }
};


double integrate_prob(double a, double b,
                    double a_obd, double b_obd)
{
  const double lower = 0, upper = 1 - 1e-15;

  BetaIntegrate f(a, b, a_obd, b_obd);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  
  return res;
}

DataFrame expand_prob(int n1, int n2, double a0, double b0){
  IntegerVector a1 = seq(0,n1);
  IntegerVector a2 = seq(0,n2);
  IntegerVector b1 = rev(seq(0,n1));
  IntegerVector b2 = rev(seq(0,n2));
  Function expand("expand.grid");
  DataFrame a_value_comb = expand(a1,a2);
  DataFrame b_value_comb = expand(b1,b2);
  NumericVector a_grp1 = a_value_comb[0];
  NumericVector a_grp2 = a_value_comb[1];
  NumericVector b_grp1 = b_value_comb[0];
  NumericVector b_grp2 = b_value_comb[1];
  a_grp1 = a_grp1 + a0;
  a_grp2 = a_grp2 + a0;
  b_grp1 = b_grp1 + b0;
  b_grp2 = b_grp2 + b0;
  
  auto calc_prob = [&](int r) {
    return integrate_prob(a_grp1[r], b_grp1[r], a_grp2[r], b_grp2[r]);
  };
  
  NumericVector probability = sapply(seq(0, a_grp1.size()-1), calc_prob);
  
  DataFrame combinations = DataFrame::create(Named("a.value.grp1") = a_grp1,
                                             Named("b.value.grp1") = b_grp1,
                                             Named("a.value.grp2") = a_grp2,
                                             Named("b.value.grp2") = b_grp2,
                                             Named("probability") = probability);
  
  return combinations;
}

// [[Rcpp::export]]
DataFrame generate_df(int n_stage1, int n_stage2, int N_dose, int N_arm,
                         NumericVector eff_prob, double eff_lower, double Cs2,
                         int cutoff_step, int type, bool upfront){
  Function merge("merge");
  Function subset("[.data.frame");
  double a0 = 0.1;
  double b0 = 0.1;
  int n_total = n_stage1 + n_stage2;
  DataFrame all_cases;
  if(upfront){
    all_cases = expand_prob(n_total, n_total, a0, b0);
  }else{
    all_cases = expand_prob(n_stage2, n_total, a0, b0);
  }

  IntegerVector n_total_seq = seq(0, n_total);
  IntegerVector n_stage2_seq = seq(0, n_stage2);
  NumericVector a_value = as<NumericVector>(n_total_seq) + a0;
  NumericVector a_value_stage2 = as<NumericVector>(n_stage2_seq) + a0;
  NumericVector PMF_ctrl;
  DataFrame PMF_grp_ctrl;
  if(upfront){
    PMF_ctrl = dbinom(n_total_seq, n_total, eff_prob[0]);
    PMF_grp_ctrl = DataFrame::create(
      Named("a.value.grp1") = a_value,
      Named("PMF.grp1") = PMF_ctrl
    );
  }else{
    PMF_ctrl = dbinom(n_stage2_seq, n_stage2, eff_prob[0]);
    PMF_grp_ctrl = DataFrame::create(
      Named("a.value.grp1") = a_value_stage2,
      Named("PMF.grp1") = PMF_ctrl
    );
  }

  DataFrame PMF_total;
  if(type==1){
    PMF_total = DGP_comb_total(
      N_dose,n_stage1, n_stage2,
      eff_prob[Range((N_arm - 1), (N_arm + N_dose -2))],
      max(eff_prob[Range((N_arm - 1), (N_arm + N_dose -2))]),
      eff_lower, Cs2, cutoff_step,type);
  }else if(type==2){
    PMF_total = DGP_comb_total(
      N_dose,n_stage1, n_stage2,
      eff_prob[Range((N_arm - 1), (N_arm + N_dose -2))],
      min(eff_prob[Range((N_arm - 1), (N_arm + N_dose -2))]),
      eff_lower, Cs2, cutoff_step, type);
  }
  NumericVector PMF_comb = PMF_total["PMF"];
  DataFrame PMF_grp2 = DataFrame::create(
    Named("a.value.grp2") = a_value,
    Named("PMF.grp2") = PMF_comb
  );

  DataFrame all_cases_ctrl = merge(all_cases, PMF_grp_ctrl, "a.value.grp1");
  all_cases_ctrl = merge(all_cases_ctrl, PMF_grp2, "a.value.grp2");
  NumericVector pmf1 = all_cases_ctrl["PMF.grp1"];
  NumericVector pmf2 = all_cases_ctrl["PMF.grp2"];
  all_cases_ctrl["PMF"] = pmf1 * pmf2;
  all_cases_ctrl = DataFrame(all_cases_ctrl);

  auto calc_pmf = [&](int r) {
    double prob1 = R::dbinom(r, n_total, eff_prob[1], false) * R::dbinom(r, n_total, eff_prob[2], false);
    double prob2 = R::dbinom(r, n_total, eff_prob[1], false) * R::pbinom((r - 1), n_total, eff_prob[2], true, false);
    double prob3 = R::dbinom(r, n_total, eff_prob[2], false) * R::pbinom((r - 1), n_total, eff_prob[1], true, false);
    return prob1 + prob2 + prob3;
  };

  auto calc_pmf_stage2 = [&](int r) {
    double prob1 = R::dbinom(r, n_stage2, eff_prob[1], false) * R::dbinom(r, n_stage2, eff_prob[2], false);
    double prob2 = R::dbinom(r, n_stage2, eff_prob[1], false) * R::pbinom((r - 1), n_stage2, eff_prob[2], true, false);
    double prob3 = R::dbinom(r, n_stage2, eff_prob[2], false) * R::pbinom((r - 1), n_stage2, eff_prob[1], true, false);
    return prob1 + prob2 + prob3;
  };

  NumericVector PMF_single;
  DataFrame PMF_grp_single;
  if(upfront){
    if(N_arm == 4){
      PMF_single = sapply(seq(0,n_total), calc_pmf);
    }else if(N_arm == 3){
      PMF_single = dbinom(n_total_seq, n_total, eff_prob[1]);
    }else if(N_arm == 2){
      PMF_single = rep(0, (n_total+1));
    }
    PMF_grp_single = DataFrame::create(
      Named("a.value.grp1") = a_value,
      Named("PMF.grp1") = PMF_single
    );
  }else{
    if(N_arm == 4){
      PMF_single = sapply(seq(0,n_stage2), calc_pmf_stage2);
    }else if(N_arm == 3){
      PMF_single = dbinom(n_stage2_seq, n_stage2, eff_prob[1]);
    }else if(N_arm == 2){
      PMF_single = rep(0, (n_stage2+1));
    }
    PMF_grp_single = DataFrame::create(
      Named("a.value.grp1") = a_value_stage2,
      Named("PMF.grp1") = PMF_single
    );
  }


  DataFrame all_cases_single = merge(all_cases, PMF_grp_single, "a.value.grp1");

 return Rcpp::List::create(all_cases_ctrl, all_cases_single);
}



