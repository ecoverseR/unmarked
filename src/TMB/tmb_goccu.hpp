#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// Adapted from Stan code by Maxwell B. Joseph,
// https://discourse.mc-stan.org/t/divergent-transition-every-iteration-multi-scale-occupancy-model/13739/5

// name of function below **MUST** match filename
template<class Type>
Type tmb_goccu(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  DATA_INTEGER(T);

  DATA_INTEGER(link);

  DATA_MATRIX(X_psi);
  DATA_MATRIX(X_phi);
  DATA_MATRIX(X_det);

  DATA_VECTOR(known_present);
  DATA_MATRIX(known_available);
  DATA_MATRIX(missing_session);

  PARAMETER_VECTOR(beta_psi);
  PARAMETER_VECTOR(beta_phi);
  PARAMETER_VECTOR(beta_det);

  vector<Type> psi = X_psi * beta_psi;
  vector<Type> phi = X_phi * beta_phi;
  vector<Type> p = X_det * beta_det;

  if(link == 1){
    psi = cloglog(psi);
  } else {
    psi = invlogit(psi);
  }

  phi = invlogit(phi);
  p = invlogit(p);

  Type loglik = 0.0;

  int M = y.rows();
  int J = y.cols() / T;

  Type obs_lp;
  Type poss_lp;

  int ystart;
  vector<Type> ysite;
  vector<Type> psite;
  vector<Type> ysub;
  vector<Type> psub;

  int tstart;
  int pstart;

  for (int i=0; i<M; i++){

    tstart = i*T;
    pstart = i*(T*J);

    ysite = y.row(i);
    psite = p.segment(pstart, T*J);

    if(known_present(i) == 1){
      for (int t=0; t<T; t++){

        if(missing_session(i, t) == 1) continue;

        ystart = t*J;

        psub = psite.segment(ystart, J);
        ysub = ysite.segment(ystart, J);

        obs_lp = log(phi(tstart+t));
        for (int j=0; j<J; j++){
            if(R_IsNA(asDouble(ysub(j)))) continue;
            obs_lp += dbinom(ysub(j), Type(1), psub(j), true);
        }
        if(known_available(i, t) == 1){
          loglik += obs_lp;
        } else {
          loglik += logspace_add(obs_lp, log(1-phi(tstart+t)));
        }
      }
      loglik += log(psi(i));
    } else {
      // Original math in Nichols: 
      // Pr(000 000) = (1−ψ) + ψ(1−θ₁)(1−θ₂) + ψθ₁q₁(1−θ₂) + ψ(1−θ₁)θ₂q₂ + ψθ₁θ₂q₁q₂
      //
      // Factorize and simplify with help from Claude:
      // = (1−ψ) + ψ · [ (1−θ₁)(1−θ₂) + θ₁q₁(1−θ₂) + (1−θ₁)θ₂q₂ + θ₁q₁·θ₂q₂ ]
      //
      //  θ_t is the probability the species is available (phi in the code).
      //  q_t = ∏_k (1 − p_tk) is the probability of missing it on 
      //  every secondary sample in period t, given that it’s available.
      //  a_t = 1 − θ_t is the probability that it’s unavailable.
      //  b_t = θ_t · q_t is the probability that it’s available but never detected.
      //
      //  Term	        Period 1	          Period 2	          In a/b form
      // (1−θ₁)(1−θ₂)	  unavailable	        unavailable	        a₁a₂
      //  θ₁q₁(1−θ₂)	  available, missed	  unavailable	        b₁a₂
      //  (1−θ₁)θ₂q₂	  unavailable	        available, missed	  a₁b₂
      //  θ₁q₁θ₂q₂	    available, missed	  available, missed	  b₁b₂
      //
      //  = (a₁ + b₁)·a₂ + (a₁ + b₁)·b₂
      //  = (a₁ + b₁)(a₂ + b₂)
      //
      //  Thus:
      //  Pr(000 000) = (1−ψ) + ψ · [(1−θ₁) + θ₁q₁] · [(1−θ₂) + θ₂q₂]
      // 
      //  In general: Pr(all zeros) = (1−ψ) + ψ · ∏ₜ [(1−θ_t) + θ_t q_t]
      //  log Pr = log[ exp( log ψ + Σₜ log(a_t + b_t) ) + (1−ψ) u]
      //  Missing periods contribute nothing to likelihood and are skipped
      //
      poss_lp = log(psi(i));
      for (int t=0; t<T; t++){
        if(missing_session(i, t) == 1) continue;
        ystart = t*J;
        psub = psite.segment(ystart, J);
        ysub = ysite.segment(ystart, J);
        obs_lp = log(phi(tstart+t));
        for (int j=0; j<J; j++){
          if(R_IsNA(asDouble(ysub(j)))) continue;
          obs_lp += dbinom(ysub(j), Type(1), psub(j), true);
        }
        poss_lp += logspace_add(obs_lp, log(1-phi(tstart+t)));
      }
      loglik += logspace_add(poss_lp, log(1-psi(i)));
    }
  }

  return -loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
