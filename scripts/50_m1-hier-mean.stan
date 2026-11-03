// LOCAL EXPERIMENT copy of m1.stan (49278c9): as hier-sd, plus a shared mean m_theta0 ~ N(0,1) on theta0.
// cytoJointR M1: donor-paired log-linear model of exact marker-pattern counts.
// eta_sk = sum_j B_kj a_sj + sum_{j<l} B_kj B_kl theta_sjl, p_s = softmax(eta_s).
// z and u_raw are standardised (noncentered, default) or raw (centered) effects.
// a_sj  = alpha_j + X_s gamma_j + (donor effect, loading W_s) + (donor x visit effect).
// overdispersed = 0: y_s ~ multinomial(mu_s), exactly M1.
// overdispersed = 1: y_s ~ dirichlet_multinomial(phi * mu_s), phi = 1 / inv_phi, so
//   p_s | y_s ~ Dirichlet(phi * mu_s + y_s) (conjugate) and inv_phi -> 0 recovers M1.
// theta_s = theta0 + Xs_s delta.   All design matrices are built in R (cj_design()).
functions {
  matrix DE_of(vector tau, matrix L_Omega, matrix z, int centered) {
    if (centered || num_elements(z) == 0) return z;
    return diag_pre_multiply(tau, L_Omega) * z;
  }
  vector u_of(vector tau_u, vector u_raw, int centered) {
    if (centered || num_elements(u_raw) == 0) return u_raw;
    return tau_u[1] * u_raw;
  }
  // DE: Dd x Ndon realised donor effects; u: realised donor x visit deviations.
  matrix eta_matrix(vector alpha, matrix gamma, vector theta0, matrix delta,
                    matrix DE, vector u, matrix X, matrix Xs, matrix W, array[] int donor_of,
                    array[] int dv_of, matrix Bt, matrix Bpt) {
    int S = rows(X);
    int J = rows(Bt);
    vector[S] shift = rep_vector(0, S);
    if (num_elements(DE) > 0) shift += columns_dot_product(W', DE[:, donor_of])';
    if (num_elements(u) > 0) shift += u[dv_of];
    matrix[S, J] A = rep_matrix(alpha', S) + X * gamma + shift * rep_row_vector(1, J);
    matrix[S, cols(delta)] TH = rep_matrix(theta0', S) + Xs * delta;
    return A * Bt + TH * Bpt;
  }
}
data {
  int<lower=1> S;                        // samples
  int<lower=2> K;                        // patterns, 2^J
  int<lower=2> J;                        // markers
  int<lower=1> M;                        // stimulations
  int<lower=1> P;                        // marker pairs
  int<lower=1> Pf;                       // fixed-effect columns (no intercept)
  int<lower=1> Ndon;
  int<lower=1> G;                        // donor x visit groups
  int<lower=0, upper=1> use_donor;
  int<lower=0, upper=1> use_dv;
  int<lower=0, upper=1> centered;        // 1: z and u_raw are the raw effects
  int<lower=0, upper=1> overdispersed;   // 1: Dirichlet-multinomial sample residual
  array[S, K] int<lower=0> y;
  matrix[S, Pf] X;
  matrix[S, M] Xs;                       // stimulation indicators
  matrix[S, 1 + M] W;                    // donor loadings: (1, Xs)
  array[S] int<lower=1, upper=Ndon> donor_of;
  array[S] int<lower=1, upper=G> dv_of;
  matrix[J, K] Bt;
  matrix[P, K] Bpt;
  // priors
  vector[J] alpha_mean;
  real<lower=0> alpha_sd;
  vector<lower=0>[Pf] gamma_sd;
  real<lower=0> theta0_sd;
  real<lower=0> delta_sd;
  real<lower=0> tau_sd;
  real<lower=0> tau_u_sd;
  real<lower=0> lkj_eta;
  real<lower=0> inv_phi_sd;              // half-normal scale of inv_phi
}
transformed data {
  int Dd = use_donor ? 1 + M : 0;
  int Nz = use_donor ? Ndon : 0;
  int Gu = use_dv ? G : 0;
}
parameters {
  vector[J] alpha;
  matrix[Pf, J] gamma;
  vector[P] theta0_raw;
  real m_theta0;
  matrix[M, P] delta_raw;
  real<lower=0> s_theta0;
  real<lower=0> s_delta;
  vector<lower=0>[Dd] tau;
  cholesky_factor_corr[Dd] L_Omega;
  matrix[Dd, Nz] z;
  vector<lower=0>[use_dv] tau_u;
  vector[Gu] u_raw;
  vector<lower=0>[overdispersed] inv_phi;
}
transformed parameters {
  vector[P] theta0 = m_theta0 + s_theta0 * theta0_raw;
  matrix[M, P] delta = s_delta * delta_raw;
}
model {
  alpha ~ normal(alpha_mean, alpha_sd);
  for (j in 1:J) gamma[, j] ~ normal(0, gamma_sd);
  m_theta0 ~ normal(0, 1);
  s_theta0 ~ normal(0, theta0_sd);
  s_delta ~ normal(0, delta_sd);
  theta0_raw ~ std_normal();
  to_vector(delta_raw) ~ std_normal();
  if (use_donor) {
    tau ~ normal(0, tau_sd);
    L_Omega ~ lkj_corr_cholesky(lkj_eta);
    if (centered) {
      matrix[Dd, Dd] LS = diag_pre_multiply(tau, L_Omega);
      for (i in 1:Nz) z[, i] ~ multi_normal_cholesky(rep_vector(0, Dd), LS);
    } else {
      to_vector(z) ~ std_normal();
    }
  }
  if (use_dv) {
    tau_u ~ normal(0, tau_u_sd);
    if (centered) u_raw ~ normal(0, tau_u[1]);
    else u_raw ~ std_normal();
  }
  matrix[S, K] eta = eta_matrix(alpha, gamma, theta0, delta, DE_of(tau, L_Omega, z, centered),
                                u_of(tau_u, u_raw, centered), X, Xs, W, donor_of, dv_of, Bt, Bpt);
  if (overdispersed) {
    inv_phi ~ normal(0, inv_phi_sd);
    for (s in 1:S) y[s] ~ dirichlet_multinomial(softmax(eta[s]') / inv_phi[1]);
  } else {
    for (s in 1:S) y[s] ~ multinomial_logit(eta[s]');
  }
}
generated quantities {
  array[S] vector[K] p;                  // sample-level probabilities
  array[S] vector[K] mu;                 // structured probabilities softmax(eta_s)
  vector[S] log_lik;
  vector[overdispersed] phi;
  vector[(Dd * (Dd - 1)) %/% 2] rho;
  {
    matrix[S, K] eta = eta_matrix(alpha, gamma, theta0, delta, DE_of(tau, L_Omega, z, centered),
                                  u_of(tau_u, u_raw, centered), X, Xs, W, donor_of, dv_of, Bt, Bpt);
    if (overdispersed) phi[1] = 1 / inv_phi[1];
    for (s in 1:S) {
      mu[s] = softmax(eta[s]');
      if (overdispersed) {
        // conjugate posterior of p_s; a draw that underflows to nothing falls back to mu_s
        p[s] = dirichlet_rng(phi[1] * mu[s] + to_vector(y[s]));
        if (!(sum(p[s]) > 0)) p[s] = mu[s];
        log_lik[s] = dirichlet_multinomial_lpmf(y[s] | phi[1] * mu[s]);
      } else {
        p[s] = mu[s];
        log_lik[s] = multinomial_logit_lpmf(y[s] | eta[s]');
      }
    }
    if (Dd > 1) {
      matrix[Dd, Dd] Omega = multiply_lower_tri_self_transpose(L_Omega);
      int c = 0;
      for (a in 1:(Dd - 1)) for (b in (a + 1):Dd) { c += 1; rho[c] = Omega[a, b]; }
    }
  }
}
