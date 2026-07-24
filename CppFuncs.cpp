// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// General C++ Code functions                                        //
// in the paper:                                                     //
// Robust Best Subset Selection via Fast Approximate MM-Estimation   //
//                                                                   //
// Author: Martin Huang (martin.huang@sydney.edu.au)                 //          
//         School of Mathematics & Statistics, University of Sydney  //          
//         AUSTRALIA                                                 //          
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <random>
#include <chrono>
using namespace Rcpp;




// [[Rcpp::export]]
arma::vec wlsfit_cpp(const arma::mat& X,
                     const arma::vec& y,
                     const arma::vec& w,
                     bool add_int = false) {
  arma::mat Xnew = add_int
  ? arma::join_horiz(arma::ones<arma::vec>(X.n_rows), X)
    : X;
  arma::vec sqrtw = arma::sqrt(w);
  arma::mat Xw = Xnew.each_col() % sqrtw;
  arma::vec yw = y % sqrtw;
  return arma::solve(Xw.t() * Xw, Xw.t() * yw);
}

// [[Rcpp::export]]
Rcpp::NumericVector Rho_MM(const arma::vec& x, double cc = 4.685061) {
  arma::vec tmp  = x / cc;
  
  arma::vec tmp2 = std::pow(cc, 2.0) *
    arma::square(tmp) %
    (1.0 - arma::square(tmp) + arma::pow(tmp, 4) / 3.0) / 2.0;
  
  tmp2.elem(arma::find(arma::abs(tmp) > 1.0)).fill(std::pow(cc, 2.0) / 6.0);
  
  return Rcpp::NumericVector(tmp2.begin(), tmp2.end());
} 

// [[Rcpp::export]]
double mad(const arma::vec& resid) {
  const double constant = 1.4826;
  int n = resid.n_elem;
  
  if (n == 0) return NA_REAL;
  
  double med = arma::median(resid);
  arma::vec abs_dev = arma::abs(resid - med);
  return constant * arma::median(abs_dev);
} 


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List msFAMMcpp(const arma::mat& X, 
               const arma::vec& y, 
               const arma::imat& A, 
               int strata, 
               const std::string& params, 
               int B, 
               int m,
               bool crit = false,
               bool MMweights = false) { 
  int n = X.n_rows;
  int p = X.n_cols;
  int num_models = A.n_rows;
  
  // start timer
  Function sys_time = Environment::base_env()["Sys.time"];
  Function difft    = Environment::base_env()["difftime"];
  Function rep_f    = Environment::base_env()["rep"];
  NumericVector start = sys_time();
  
  // 1) Full-model MM fit
  arma::mat X_full(n, p+1);
  X_full.col(0).ones();
  X_full.cols(1, p) = X;
  Environment MASS = Environment::namespace_env("MASS");
  Function rlm_def = MASS["rlm.default"];
  List RLM = rlm_def(_["x"]=wrap(X_full), _["y"]=y,
                     _["method"]="MM", _["maxit"]=500);
  arma::vec weights = as<arma::vec>(RLM["w"]);
  arma::vec resid   = as<arma::vec>(RLM["residuals"]);
  
  arma::vec w_crit(n, arma::fill::ones);
  if (MMweights) {
    w_crit = weights;
  }
  
  
  // 2) residual ranks & strata assignment
  // ord[k] = original index of the k-th smallest residual (the sort permutation).
  // We need the inverse: rnk[i] = 0-based rank of observation i.
  arma::uvec ord = arma::sort_index(resid);
  arma::uvec rnk(n);
  rnk.elem(ord) = arma::regspace<arma::uvec>(0, n - 1);   // invert the permutation
  arma::ivec res_ind = arma::conv_to<arma::ivec>::from(
    arma::floor(arma::conv_to<arma::vec>::from(rnk) * strata / double(n))) + 1;
  
  // 3) split into strata
  std::vector<arma::mat> X_list(strata);
  std::vector<arma::vec> y_list(strata), w_list(strata);
  for(int h = 1; h <= strata; ++h) {
    arma::uvec idx = arma::find(res_ind == h);
    X_list[h-1] = X.rows(idx);
    y_list[h-1] = y.elem(idx);
    w_list[h-1] = weights.elem(idx);
  } 
  
  // 4) scale computation
  double s_scale = 0.0, mad_full = 0.0;
  int b_huber = 2;
  if(params == "MSB") {
    int best_i = arma::index_max(arma::sum(A,1));
    arma::uvec keep = arma::find(A.row(best_i) == 1);
    arma::mat X_sub(n, keep.n_elem+1);
    X_sub.col(0).ones();
    for(size_t k=0; k<keep.n_elem; ++k)
      X_sub.col(k+1) = X.col(keep[k]);
    if((int)keep.n_elem != p) {
      RLM = rlm_def(_["x"]=wrap(X_sub), _["y"]=y,
                    _["method"]="MM", _["maxit"]=500);
    } 
    Function summary_rlm = MASS["summary.rlm"];
    List sumL = summary_rlm(RLM);
    s_scale = as<double>(sumL["sigma"]);
  } else {
    mad_full = mad(resid); 
  }
  
  // 5) precompute model columns
  std::vector<arma::uvec> modelCols(num_models);
  for(int i=0; i<num_models; ++i)
    modelCols[i] = arma::find(A.row(i) == 1);
  
  // 6) bootstrap 
  std::vector<double> MSC2(num_models, 0.0);
  unsigned seed = static_cast<unsigned>(
    std::chrono::high_resolution_clock::now()
    .time_since_epoch().count());
    std::mt19937_64 rng(seed);
    
    for(int b = 0; b < B; ++b) {
      int nSample = m / strata;
      arma::mat Xboot(nSample*strata, p);
      arma::vec yboot(nSample*strata), wboot(nSample*strata);
      int row = 0;
      for(int h=0; h<strata; ++h) {
        int nh = X_list[h].n_rows;
        std::uniform_int_distribution<int> dist(0, nh-1);
        arma::uvec idx(nSample);
        for(int j=0; j<nSample; ++j) idx[j] = dist(rng);
        Xboot.rows(row, row+nSample-1) = X_list[h].rows(idx);
        yboot.subvec(row, row+nSample-1) = y_list[h].elem(idx);
        wboot.subvec(row, row+nSample-1) = w_list[h].elem(idx);
        row += nSample;
      } 
      for(int pp=0; pp<num_models; ++pp) {
        arma::mat Xb = Xboot.cols(modelCols[pp]);
        arma::vec coef = wlsfit_cpp(Xb, yboot, wboot, true);
        arma::mat Xf(n, modelCols[pp].n_elem+1);
        Xf.col(0).ones();
        for(size_t k=0; k<modelCols[pp].n_elem; ++k)
          Xf.col(k+1) = X.col(modelCols[pp][k]);
        arma::vec res = y - Xf * coef;
        if(params=="MSB") {
          NumericVector rho = Rho_MM(res / s_scale);
          arma::vec rho_a = as<arma::vec>(rho);
          MSC2[pp] += arma::dot(w_crit, rho_a);
        } else { 
          double cutoff = std::pow(b_huber * mad_full, 2);
          arma::uvec small = arma::find(arma::square(res) < cutoff);
          arma::uvec large = arma::find(arma::square(res) >= cutoff);
          MSC2[pp] += (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
            cutoff * arma::accu(w_crit.elem(large))) / n;   // (weighted) reduces to original when w_i = 1
        } 
      } 
    } 
    
    // 7) compute MSC5
    std::vector<double> MSC5(num_models);
    for(int pp=0; pp<num_models; ++pp) {
      MSC5[pp] = MSC2[pp] / B;
      arma::mat Xf(n, modelCols[pp].n_elem+1);
      Xf.col(0).ones();
      for(size_t k=0; k<modelCols[pp].n_elem; ++k)
        Xf.col(k+1) = X.col(modelCols[pp][k]);
      arma::vec coef = wlsfit_cpp(Xf, y, weights, false);
      arma::vec res = y - Xf * coef;
      if(params=="MSB") {
        NumericVector rho = Rho_MM(res / s_scale);
        arma::vec rho_a = as<arma::vec>(rho);                   // (added)
        double pen = (arma::dot(w_crit, rho_a) +                // (weighted) in-sample loss; penalty term unweighted
                      std::log(n) * modelCols[pp].n_elem * std::pow(s_scale,2)) / n;
        MSC5[pp] += pen;
      } else { 
        double cutoff = std::pow(b_huber * mad_full, 2);
        arma::uvec small = arma::find(arma::square(res) < cutoff);
        arma::uvec large = arma::find(arma::square(res) >= cutoff);           // (added) explicit large-residual index set
        double pen = (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
                      cutoff * arma::accu(w_crit.elem(large)) +
                      2 * std::log(n) * modelCols[pp].n_elem * std::pow(mad_full,2)) / n;  // (weighted) loss; penalty unweighted
        MSC5[pp] += pen;
      } 
    } 
    
    if (crit == false){
      // 8) select winners
      IntegerVector M05ind(num_models, 0), M02ind(num_models, 0);
      int best5 = std::distance(MSC5.begin(), std::min_element(MSC5.begin(), MSC5.end()));
      int best2 = std::distance(MSC2.begin(), std::min_element(MSC2.begin(), MSC2.end()));
      M05ind[best5] = 1;
      M02ind[best2] = 1;
      
      // end timer
      NumericVector end   = sys_time();
      NumericVector dt    = difft(_["time1"]=end, _["time2"]=start, _["units"]="secs");
      NumericVector Time  = rep_f(_["x"]=dt, _["times"]=2);
      
      // 9) output
      CharacterVector Criterion = {"M05.ind","M02.ind"};
      IntegerVector TrueCol   = {M05ind[0], M02ind[0]};
      List out;
      out["Criterion"] = Criterion;
      out["True"]      = TrueCol;
      for(int j=1; j<num_models; ++j) {
        std::string nm = "M" + std::to_string(j+1);
        out[nm] = IntegerVector::create(M05ind[j], M02ind[j]);
      }
      out["Time"]      = Time;
      out["Method"]    = CharacterVector::create("FAMM","FAMM");
      out.attr("class")     = CharacterVector::create("data.frame");
      out.attr("row.names") = IntegerVector::create(NA_INTEGER, -2);
      return out;
    }else{
      return List::create(_["MSC2"] = MSC2,
                          _["MSC5"] = MSC5);
    }
    
}


// [[Rcpp::export]]
List msMMcpp(const arma::mat& X, 
             const arma::vec& y, 
             const arma::imat& A, 
             int strata, 
             const std::string& params, 
             int B, 
             int m,
             bool crit = false,
             bool MMweights = false) { 
  int n = X.n_rows;
  int p = X.n_cols;
  int num_models = A.n_rows;
  
  // start timer
  Function sys_time = Environment::base_env()["Sys.time"];
  Function difft    = Environment::base_env()["difftime"];
  Function rep_f    = Environment::base_env()["rep"];
  NumericVector start = sys_time();
  
  // 1) Full-model MM fit
  arma::mat X_full(n, p+1);
  X_full.col(0).ones();
  X_full.cols(1, p) = X;
  Environment MASS = Environment::namespace_env("MASS");
  Function rlm_def = MASS["rlm.default"];
  List RLM = rlm_def(_["x"]=wrap(X_full), _["y"]=y,
                     _["method"]="MM", _["maxit"]=500);
  arma::vec weights = as<arma::vec>(RLM["w"]);
  arma::vec resid   = as<arma::vec>(RLM["residuals"]);
  
  arma::vec w_crit(n, arma::fill::ones);
  if (MMweights) {
    w_crit = weights;
  }
  
  
  // 2) residual ranks & strata assignment
  // ord[k] = original index of the k-th smallest residual (the sort permutation).
  // We need the inverse: rnk[i] = 0-based rank of observation i.
  arma::uvec ord = arma::sort_index(resid);
  arma::uvec rnk(n);
  rnk.elem(ord) = arma::regspace<arma::uvec>(0, n - 1);   // invert the permutation
  arma::ivec res_ind = arma::conv_to<arma::ivec>::from(
    arma::floor(arma::conv_to<arma::vec>::from(rnk) * strata / double(n))) + 1;
  
  // 3) split into strata
  std::vector<arma::mat> X_list(strata);
  std::vector<arma::vec> y_list(strata);
  for(int h = 1; h <= strata; ++h) {
    arma::uvec idx = arma::find(res_ind == h);
    X_list[h-1] = X.rows(idx);
    y_list[h-1] = y.elem(idx);
  } 
  
  // 4) scale computation
  double s_scale = 0.0, mad_full = 0.0;
  int b_huber = 2;
  if(params == "MSB") {
    int best_i = arma::index_max(arma::sum(A,1));
    arma::uvec keep = arma::find(A.row(best_i) == 1);
    arma::mat X_sub(n, keep.n_elem+1);
    X_sub.col(0).ones();
    for(size_t k=0; k<keep.n_elem; ++k)
      X_sub.col(k+1) = X.col(keep[k]);
    if((int)keep.n_elem != p) {
      RLM = rlm_def(_["x"]=wrap(X_sub), _["y"]=y,
                    _["method"]="MM", _["maxit"]=500);
    } 
    Function summary_rlm = MASS["summary.rlm"];
    List sumL = summary_rlm(RLM);
    s_scale = as<double>(sumL["sigma"]);
  } else {
    mad_full = mad(resid); 
  }
  
  // 5) precompute model columns
  std::vector<arma::uvec> modelCols(num_models);
  for(int i=0; i<num_models; ++i)
    modelCols[i] = arma::find(A.row(i) == 1);
  
  // 6) bootstrap MSC2
  std::vector<double> MSC2(num_models, 0.0);
  unsigned seed = static_cast<unsigned>(
    std::chrono::high_resolution_clock::now()
    .time_since_epoch().count());
    std::mt19937_64 rng(seed);
    
    for(int b = 0; b < B; ++b) {
      int nSample = m / strata;
      arma::mat Xboot(nSample*strata, p);
      arma::vec yboot(nSample*strata);
      int row = 0;
      for(int h=0; h<strata; ++h) {
        int nh = X_list[h].n_rows;
        std::uniform_int_distribution<int> dist(0, nh-1);
        arma::uvec idx(nSample);
        for(int j=0; j<nSample; ++j) idx[j] = dist(rng);
        Xboot.rows(row, row+nSample-1) = X_list[h].rows(idx);
        yboot.subvec(row, row+nSample-1) = y_list[h].elem(idx);
        row += nSample;
      } 
      for(int pp=0; pp<num_models; ++pp) {
        // arma::mat Xb = Xboot.cols(modelCols[pp]);
        arma::mat Xb(Xboot.n_rows, modelCols[pp].n_elem + 1);
        Xb.col(0).ones();
        for (size_t k = 0; k < modelCols[pp].n_elem; ++k)
          Xb.col(k + 1) = Xboot.col(modelCols[pp][k]);
        
        List r = rlm_def(_["x"]=wrap(Xb), _["y"]=yboot,
                         _["method"]="MM", _["maxit"]=500);
        arma::vec coef  = as<arma::vec>(r["coefficients"]);
        arma::mat Xf(n, modelCols[pp].n_elem+1);
        Xf.col(0).ones();
        for(size_t k=0; k<modelCols[pp].n_elem; ++k)
          Xf.col(k+1) = X.col(modelCols[pp][k]);
        arma::vec res = y - Xf * coef;
        if(params=="MSB") {
          NumericVector rho = Rho_MM(res / s_scale);
          arma::vec rho_a = as<arma::vec>(rho);                 // (added)
          MSC2[pp] += arma::dot(w_crit, rho_a); 
        } else { 
          double cutoff = std::pow(b_huber * mad_full, 2);
          arma::uvec small = arma::find(arma::square(res) < cutoff);
          arma::uvec large = arma::find(arma::square(res) >= cutoff);
          MSC2[pp] += (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
            cutoff * arma::accu(w_crit.elem(large))) / n;   // (weighted) reduces to original when w_i = 1
        } 
      } 
    } 
    
    // 7) compute MSC5
    std::vector<double> MSC5(num_models);
    for(int pp=0; pp<num_models; ++pp) {
      MSC5[pp] = MSC2[pp] / B;
      arma::mat Xf(n, modelCols[pp].n_elem+1);
      Xf.col(0).ones();
      for(size_t k=0; k<modelCols[pp].n_elem; ++k)
        Xf.col(k+1) = X.col(modelCols[pp][k]);
      
      List r5 = rlm_def(_["x"]=wrap(Xf), _["y"]=y,
                        _["method"]="MM", _["maxit"]=500);
      arma::vec coef  = as<arma::vec>(r5["coefficients"]);
      arma::vec res = y - Xf * coef;
      if(params=="MSB") {
        NumericVector rho = Rho_MM(res / s_scale);
        arma::vec rho_a = as<arma::vec>(rho);                   // (added)
        double pen = (arma::dot(w_crit, rho_a) +                // (weighted) in-sample loss; penalty term unweighted
                      std::log(n) * modelCols[pp].n_elem * std::pow(s_scale,2)) / n;
        MSC5[pp] += pen;
      } else { 
        double cutoff = std::pow(b_huber * mad_full, 2);
        arma::uvec small = arma::find(arma::square(res) < cutoff);
        arma::uvec large = arma::find(arma::square(res) >= cutoff); 
        double pen = (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
                      cutoff * arma::accu(w_crit.elem(large)) +
                      2 * std::log(n) * modelCols[pp].n_elem * std::pow(mad_full,2)) / n;  // (weighted) loss; penalty unweighted
        MSC5[pp] += pen;
      } 
    } 
    
    if (crit == false){
      // 8) select winners
      IntegerVector M05ind(num_models, 0), M02ind(num_models, 0);
      int best5 = std::distance(MSC5.begin(), std::min_element(MSC5.begin(), MSC5.end()));
      int best2 = std::distance(MSC2.begin(), std::min_element(MSC2.begin(), MSC2.end()));
      M05ind[best5] = 1;
      M02ind[best2] = 1;
      
      // end timer
      NumericVector end   = sys_time();
      NumericVector dt    = difft(_["time1"]=end, _["time2"]=start, _["units"]="secs");
      NumericVector Time  = rep_f(_["x"]=dt, _["times"]=2);
      
      // 9) output
      CharacterVector Criterion = {"M05.ind","M02.ind"};
      IntegerVector TrueCol   = {M05ind[0], M02ind[0]};
      List out;
      out["Criterion"] = Criterion;
      out["True"]      = TrueCol;
      for(int j=1; j<num_models; ++j) {
        std::string nm = "M" + std::to_string(j+1);
        out[nm] = IntegerVector::create(M05ind[j], M02ind[j]);
      }
      out["Time"]      = Time;
      out["Method"]    = CharacterVector::create("MM","MM");
      out.attr("class")     = CharacterVector::create("data.frame");
      out.attr("row.names") = IntegerVector::create(NA_INTEGER, -2);
      return out;
    }else{
      return List::create(_["MSC2"] = MSC2,
                          _["MSC5"] = MSC5);
    }
    
}


// [[Rcpp::export]]
double aic_cpp(arma::mat X, arma::vec y) {
  int n = X.n_rows;
  
  // Estimate coefficients
  arma::vec coef = arma::solve(X.t() * X, X.t() * y);
  arma::vec resid = y - X * coef;
  
  double rss = arma::dot(resid, resid);
  int p = X.n_cols; 
  
  double logLik = -0.5 * n * (std::log(2 * M_PI) + 1 - std::log(n) + std::log(rss));
  
  // df = p + 1  (includes sigma^2)
  double k = p + 1;
  
  double Aic = 2*k - 2 * logLik;
  return Aic;
}



// [[Rcpp::export]]
NumericVector computeAICcpp(const arma::mat& A, const arma::mat& X, const arma::vec& y, int B = 100) {
  int nA = A.n_rows;
  int nX = X.n_rows;
  NumericVector aic(nA);
  
  for (int i = 0; i < nA; i++) {
    NumericVector AICa(B);
    arma::rowvec a = A.row(i);
    arma::uvec idx_a = find(a == 1);
    
    for (int b = 0; b < B; b++) {
      IntegerVector idx = sample(nX, nX / 2, true) - 1;
      arma::uvec idx_u = as<arma::uvec>(idx);
      
      arma::mat Xa;
      if (idx_a.n_elem == 0) {
        Xa = arma::ones(idx_u.n_elem, 1);  // intercept-only model
      } else {
        Xa = X.rows(idx_u);
        Xa = Xa.cols(idx_a);
        Xa.insert_cols(0, arma::ones(idx_u.n_elem, 1));  // add intercept
      }
      arma::vec ya = y.elem(idx_u);
      
      AICa[b] = aic_cpp(Xa, ya);
    }
    
    aic[i] = mean(AICa);
  }
  
  return aic;
}

// [[Rcpp::export]]
List critFamm(const arma::mat& X, 
              const arma::vec& y, 
              const arma::imat& A, 
              int strata, 
              const std::string& params, 
              int B, 
              int m,
              const arma::vec& weights,
              const arma::vec& resid,
              bool MMweights = false) { 
  int n = X.n_rows;
  int p = X.n_cols;
  int num_models = A.n_rows;
  
  // start timer
  Function sys_time = Environment::base_env()["Sys.time"];
  Function difft    = Environment::base_env()["difftime"];
  Function rep_f    = Environment::base_env()["rep"];
  NumericVector start = sys_time();
  
  Environment MASS = Environment::namespace_env("MASS");
  
  arma::vec w_crit(n, arma::fill::ones);
  if (MMweights) {
    w_crit = weights;
  }
  
  // 2) residual ranks & strata assignment
  // ord[k] = original index of the k-th smallest residual (the sort permutation).
  // We need the inverse: rnk[i] = 0-based rank of observation i.
  arma::uvec ord = arma::sort_index(resid);
  arma::uvec rnk(n);
  rnk.elem(ord) = arma::regspace<arma::uvec>(0, n - 1);   // invert the permutation
  arma::ivec res_ind = arma::conv_to<arma::ivec>::from(
    arma::floor(arma::conv_to<arma::vec>::from(rnk) * strata / double(n))) + 1;
  
  // 3) split into strata
  std::vector<arma::mat> X_list(strata);
  std::vector<arma::vec> y_list(strata), w_list(strata);
  for(int h = 1; h <= strata; ++h) {
    arma::uvec idx = arma::find(res_ind == h);
    X_list[h-1] = X.rows(idx);
    y_list[h-1] = y.elem(idx);
    w_list[h-1] = weights.elem(idx);
  } 
  
  // 4) scale computation
  double s_scale = 0.0, mad_full = 0.0;
  int b_huber = 2;
  mad_full = mad(resid); 
  
  // 5) precompute model columns
  std::vector<arma::uvec> modelCols(num_models);
  for(int i=0; i<num_models; ++i)
    modelCols[i] = arma::find(A.row(i) == 1);
  
  // 6) bootstrap MSC2
  std::vector<double> MSC2(num_models, 0.0);
  unsigned seed = static_cast<unsigned>(
    std::chrono::high_resolution_clock::now()
    .time_since_epoch().count());
    std::mt19937_64 rng(seed);
    
    for(int b = 0; b < B; ++b) {
      int nSample = m / strata;
      arma::mat Xboot(nSample*strata, p);
      arma::vec yboot(nSample*strata), wboot(nSample*strata);
      int row = 0;
      for(int h=0; h<strata; ++h) {
        int nh = X_list[h].n_rows;
        std::uniform_int_distribution<int> dist(0, nh-1);
        arma::uvec idx(nSample);
        for(int j=0; j<nSample; ++j) idx[j] = dist(rng);
        Xboot.rows(row, row+nSample-1) = X_list[h].rows(idx);
        yboot.subvec(row, row+nSample-1) = y_list[h].elem(idx);
        wboot.subvec(row, row+nSample-1) = w_list[h].elem(idx);
        row += nSample;
      } 
      for(int pp=0; pp<num_models; ++pp) {
        arma::mat Xb = Xboot.cols(modelCols[pp]);
        arma::vec coef = wlsfit_cpp(Xb, yboot, wboot, true);
        arma::mat Xf(n, modelCols[pp].n_elem+1);
        Xf.col(0).ones();
        for(size_t k=0; k<modelCols[pp].n_elem; ++k)
          Xf.col(k+1) = X.col(modelCols[pp][k]);
        arma::vec res = y - Xf * coef;
        if(params=="MSB") {
          NumericVector rho = Rho_MM(res / s_scale);
          arma::vec rho_a = as<arma::vec>(rho);                 // (added)
          MSC2[pp] += arma::dot(w_crit, rho_a);   
        } else { 
          double cutoff = std::pow(b_huber * mad_full, 2);
          arma::uvec small = arma::find(arma::square(res) < cutoff); 
          arma::uvec large = arma::find(arma::square(res) >= cutoff);
          MSC2[pp] += (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
            cutoff * arma::accu(w_crit.elem(large))) / n;   // (weighted) reduces to original when w_i = 1
        } 
      } 
    } 
    
    // 7) compute MSC5
    std::vector<double> MSC5(num_models);
    for(int pp=0; pp<num_models; ++pp) {
      MSC5[pp] = MSC2[pp] / B;
      arma::mat Xf(n, modelCols[pp].n_elem+1);
      Xf.col(0).ones();
      for(size_t k=0; k<modelCols[pp].n_elem; ++k)
        Xf.col(k+1) = X.col(modelCols[pp][k]);
      arma::vec coef = wlsfit_cpp(Xf, y, weights, false);
      arma::vec res = y - Xf * coef;
      if(params=="MSB") {
        NumericVector rho = Rho_MM(res / s_scale);
        arma::vec rho_a = as<arma::vec>(rho);                   // (added)
        double pen = (arma::dot(w_crit, rho_a) +                // (weighted) in-sample loss; penalty term unweighted
                      std::log(n) * modelCols[pp].n_elem * std::pow(s_scale,2)) / n;
        MSC5[pp] += pen;
      } else { 
        double cutoff = std::pow(b_huber * mad_full, 2);
        arma::uvec small = arma::find(arma::square(res) < cutoff);
        arma::uvec large = arma::find(arma::square(res) >= cutoff);           // (added) explicit large-residual index set
        double pen = (arma::accu(w_crit.elem(small) % arma::square(res.elem(small))) +
                      cutoff * arma::accu(w_crit.elem(large)) +
                      2 * std::log(n) * modelCols[pp].n_elem * std::pow(mad_full,2)) / n;  // (weighted) loss; penalty unweighted
        MSC5[pp] += pen;
      } 
    } 
    
    // 8) select winners
    IntegerVector M05ind(num_models, 0), M02ind(num_models, 0);
    int best5 = std::distance(MSC5.begin(), std::min_element(MSC5.begin(), MSC5.end()));
    int best2 = std::distance(MSC2.begin(), std::min_element(MSC2.begin(), MSC2.end()));
    M05ind[best5] = 1;
    M02ind[best2] = 1;
    
    return List::create(
      _["MSC2"] = MSC2,
      _["MSC5"] = MSC5
    );
}