functions {
  // credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower,
                                     real point_scale_lower) {
    // works only for [2x2 matrix]
    real lpdf = lkj_corr_lpdf(rho | 1)
                + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }
  
  // pH and fi at a given time at column inlet
  vector gra_state(real t, vector hplcparam) {
    vector[2] sol;
    real tg = hplcparam[1];
    real td = hplcparam[2];
    real fio = hplcparam[5];
    real fik = hplcparam[6];
    real pHo = hplcparam[8];
    real alpha1 = hplcparam[9];
    real alpha2 = hplcparam[10];
    real fi;
    
    fi = fio + (fik - fio) / tg * (t - td);
    
    if (t < td) {
      fi = fio;
    } else if (t > tg + td) {
      fi = fik;
    }
    
    sol[1] = fi;
    sol[2] = pHo + alpha1 * fi + alpha2 * fi ^ 2;
    
    return sol;
  }
  
 real funlnki(vector logkw, vector apH, vector S1, real S2, vector pKaw, vector alpha,
                  int nDiss, real fi, real pH) {
    real lnki;
    vector[3] logkix;
    vector[2] pHmpKa;
    
    logkix = log(10) *(logkw - S1*(1+S2) * fi / (1 + S2 * fi) + apH * (pH - 7));
    pHmpKa = log(10) *(pH - (pKaw + alpha * fi));
    
    if (nDiss == 0) {
      lnki = logkix[1];
    } else if (nDiss == 1) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1] + logkix[2] - logkix[1]) -
             log1p_exp(pHmpKa[1]);
    } else if (nDiss == 2) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1] + log1p_exp(pHmpKa[2]+logkix[3]-logkix[2])) -
             log1p_exp(pHmpKa[1] + log1p_exp(pHmpKa[2]));
    }
    
    return lnki;
  }  
 
  vector areaandslope(real dt, real lnki1, real lnki2, real invki1, real invki2) {
    vector[2] cki_b;
    real bo;
    real cki;
    
    if (invki2 > 1.001 * invki1) {
      bo = (lnki1 - lnki2) / dt;
      cki = (invki2 - invki1) / bo;
    }
    
    else if (invki1 > 1.001 * invki2) {
      bo = (lnki2 - lnki1) / dt;
      cki = (invki1 - invki2) / bo;
    }
    else {
      bo = 0.001 / dt;
      cki = dt * (invki2 + invki1) / 2;
    }
    
    cki_b[1] = cki;
    cki_b[2] = bo;
    
    return cki_b;
  }
  
  real chromgratrapz(int steps, vector logkw, vector apH, vector S1,  real S2, vector pKaw,
                     vector alpha, int nDiss, vector hplcparam) {
                       
    real tg = hplcparam[1];
    real td = hplcparam[2];
    real to = hplcparam[3];
    real te = hplcparam[4];
    
    vector[1] sol;
    real time1;
    real time2;
    vector[2] fipH1;
    vector[2] fipH2;
    real lnki1;
    real lnki2;
    real invki1;
    real invki2;
    vector[2] cki_b;
    real cumki1;
    real cumki2;
    real bo;
    real tr;
    real dt;
    
    dt = tg / steps;
    
    time1 = 0;
    time2 = td;
    
    fipH1 = gra_state(time1, hplcparam);
    lnki1 = funlnki(logkw, apH, S1, S2, pKaw, alpha, nDiss, fipH1[1], fipH1[2]);
    lnki2=lnki1;
    
    invki1 = exp(-lnki1)/to;
    invki2 = invki1;
    
    cumki1 = 0;
    cumki2 = td*invki1; // cumulative area
    
    bo = 0.001 / td;    // slope
    
    for (x in 1 : steps) {
      if (cumki2 >= 1) {
        continue;
      }
      time1 = time2;
      time2 += dt;
      fipH2 = gra_state(time2, hplcparam);
      lnki1 = lnki2;
      lnki2 = funlnki(logkw, apH, S1, S2, pKaw, alpha, nDiss, fipH2[1], fipH2[2]);
      invki1 = invki2;
      invki2 = exp(-lnki2)/to;
      cki_b = areaandslope(dt, lnki1, lnki2, invki1, invki2);
      cumki1 = cumki2;
      cumki2 += cki_b[1]; // cumulative area
      bo = cki_b[2];     //slope
    }
    
    if (cumki2 >= 1 && cumki1==0) {
      tr = te+to+1/invki2;
    } else if (cumki2 >= 1) {
      tr = te+to+time1+log1p_exp(log((1-cumki1)*bo*to) + lnki1)/bo;
    } else if (cumki2 < 1) {
      tr = te+to+time2+(1-cumki2)/invki2;
    }
    
    return tr;
  }
  
  real partial_sum(array[] int ind, int start, int end, 
                   vector trobs,
                   array[] int steps,
                   array[] vector hplcparam,
                   array[] int analyte,
                   array[] int column,
                   array[] int modifier,
                   array[] int R,
                   array[,] vector logkw,
                   array[,] vector apH,
                   array[, ,] vector S1,
                   array[,] real S2,
                   array[] vector pKaw,
                   array[,] vector alpha,
                   vector dlogkT,
                   array[] vector sigma) {
                     
    real lp = 0;

    for (z in start : end) {
      
    real y_hat = chromgratrapz(steps[z], 
                                logkw[analyte[z], column[z],  : ] + dlogkT[analyte[z]]  * hplcparam[z, 11],
                                apH[analyte[z], column[z], : ],
                                S1[analyte[z], modifier[z], column[z],  : ],
                                S2[modifier[z], column[z]],
                                pKaw[analyte[z],  : ],
                                alpha[analyte[z],modifier[z],  : ],
                                R[analyte[z]],
                                hplcparam[z]);
                                                    
     lp = lp + student_t_lpdf(trobs[z] | 3, y_hat,  sigma[analyte[z],column[z]]);
      
    }
    return lp;
  }
}

data {
  int nAnalytes;           // number of analytes
  int nColumns;             // number of columns
  int nModifiers;           // number of org. modifiers
  int nObs;                // number of observations
  array[nObs] int analyte;   // analyte indexes
  array[nObs] int column;   // column indexes
  array[nObs] int modifier;  // modifier indexes
  array[nObs] int<lower=1> steps;        // steps for gradient retention time aproximation
  array[nObs] vector[12] hplcparam; // [tg, td, to, te, fio, fik, org modifier, pHo, alpha1, alpha2, (temp-25)/10, column]
  vector[nAnalytes] logPobs;
  int<lower=0, upper=2> maxR; //
  array[nAnalytes] int<lower=0, upper=2> R;

  int<lower=1> nGroupsA;
  int<lower=1> nGroupsB;
  vector[nGroupsA] pKaslitA;
  vector[nGroupsB] pKaslitB;
  array[nGroupsA,2] int idxGroupsA;
  array[nGroupsB,2] int idxGroupsB;
  
  vector[nObs] trobs;                   // observed retention factors
    
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
  
  int nexpid;                // number of individual experiments
  array[nObs] int expid;   // individual experiments indexes
  int nAnalytessim;
  array[nObs] int analytesim;
}

transformed data {
  int grainsize = 1;
  array[nObs] int ind = rep_array(1, nObs);
}

parameters {
  
  real logkwHat; // typical logkw [Neutral]
  real S1mHat;   // effect of MeOH on logkw [Neutral]
  real dS1Hat;   // effect of ACN on S1m [Neutral]
  array[2] real dlogkwHat; // effect of dissociation on logkw [Acids, Bases]
  array[2] real dS1mHat;   // effect of dissociation on S1m [Acids, Bases]
  array[2] real ddS1Hat;   // effect of dissociation on dS1 [Acids, Bases] 
  real logS2mHat; // typical value of S2m (log10 scale)
  real dlogS2Hat; // effect of ACN on logS2m
  vector[2] beta; // effect of logP on logkw and S1m
  real dlogkTHat; // effect of temperature on logkw
  vector[2] apH;  // effect of pH on logkw [Acids, Bases]
  array[2] real alphamHat;   // effect of MeOH on pKa [Acids, Bases]
  array[2] real dalphaHat;   // effect of ACN on alpham [Acids, Bases]
  vector<lower=0>[3] tau;   // sd for between analyte variability of pKa's
  vector<lower=0>[3] omega; // sd of BAV [logkw,S1m, dS1]
  corr_matrix[2] rho;      // correlation matrix [logkw vs. S1m] 
  real<lower=0> omegaT;     // sd of BAV [dlogkT]
  vector<lower=0>[3] kappa; // sd of BAV [dlogkw,dS1m,ddS1]
  
  // 2nd column
  vector[nColumns-1] clogkwHat;          // effect of column on logkw [Neutral]
  vector[nColumns-1] cS1mHat;            // effect of column on S1m [Neutral]
  vector[nColumns-1] cdS1Hat;            // effect of column on dS1 [Neutral]
  matrix[nColumns-1,2] cdlogkwHat;      // effect of column on logkw [Acids, Bases]
  matrix[nColumns-1,2] cdS1mHat;        // effect of column on dS1m [Acids, Bases]
  matrix[nColumns-1,2] cddS1Hat;        // effect of column on ddS1 [Acids, Bases] 
  matrix[nColumns-1,2]  cbeta;
  vector[nColumns-1]   cdlogkTHat; // ffect of column on dlogkTHat
  matrix[nColumns-1,2] capH;            // effect of column on apH  [Acids, Bases]
  matrix<lower=0>[nColumns-1,3] comega;  // sd of BAV [clogkw,cS1m,cdS1]
  vector<lower=0>[nColumns-1]  comegaT;  // sd of BAV [cdlogkT]
  matrix<lower=0>[nColumns-1,3] ckappa;  //sd of BAV [cdlogkw,cdS1m,cddS1]

  // 1st column
  array[nAnalytes] vector[2] paramN;
  vector[nAnalytes] dS1N;
  vector[nAnalytes] dlogkT;
  vector[nGroupsA] dlogkwA;
  vector[nGroupsB] dlogkwB;
  vector[nGroupsA] dS1mA;
  vector[nGroupsB] dS1mB;
  vector[nGroupsA] dS1A;
  vector[nGroupsB] dS1B;

  // 2nd column
  array[nColumns-1] vector[nAnalytes] clogkwN;
  array[nColumns-1] vector[nAnalytes] cS1mN;
  array[nColumns-1] vector[nAnalytes] cdS1N;
  array[nColumns-1] vector[nAnalytes] etacdlogkT;
  array[nColumns-1] vector[nGroupsA] etacdlogkwA;
  array[nColumns-1] vector[nGroupsB] etacdlogkwB;
  array[nColumns-1] vector[nGroupsA] etacdS1mA;
  array[nColumns-1] vector[nGroupsB] etacdS1mB;
  array[nColumns-1] vector[nGroupsA] etacdS1A;
  array[nColumns-1] vector[nGroupsB] etacdS1B;
   
  // Dissociation
  vector[nGroupsA] pKawA;
  vector[nGroupsB] pKawB;
  vector[nGroupsA] etaalphamA;
  vector[nGroupsB] etaalphamB;
  vector[nGroupsA] etadalphaA;
  vector[nGroupsB] etadalphaB;
  
  // residual variability for the 1st and 2nd column
  real<lower=0> msigma;   // typical sigma for the 1st column]
  real<lower=0> ssigma;
  vector[nAnalytes] logsigma; 
    
  vector[nColumns-1] clogmsigma; // effect of column on log(msigma)]
  vector<lower=0>[nColumns-1] cssigma; ; //sd of residual [1st,2nd column]
  array[nColumns-1] vector[nAnalytes] etaclogsigma;
}
transformed parameters {
  
  cov_matrix[2] Omega;
  array[nAnalytes] vector[3] miu;
  array[nColumns-1,nAnalytes] vector[3] cmiu;
  array[nAnalytes,nColumns]   vector[maxR + 1] logkwx;
  array[nAnalytes,nModifiers, nColumns] vector[maxR + 1] S1x;
  array[nModifiers,nColumns]  real S2x;
  array[nAnalytes,nColumns]   vector[maxR + 1] apHx;
  array[nAnalytes,nModifiers] vector[maxR] alphax;
  array[nAnalytes, nColumns]  real dlogkTx;
  array[nAnalytes] vector[maxR] pKawx;
  array[nAnalytes] vector[nColumns] sigmax;
  
  array[nColumns-1] vector[nGroupsA] cdlogkwA;
  array[nColumns-1] vector[nGroupsB] cdlogkwB;
  array[nColumns-1] vector[nGroupsA] cdS1mA;
  array[nColumns-1] vector[nGroupsB] cdS1mB;
  array[nColumns-1] vector[nGroupsA] cdS1A;
  array[nColumns-1] vector[nGroupsB] cdS1B;
   array[nColumns-1] vector[nAnalytes] cdlogkT;
   array[nColumns-1] vector[nAnalytes] clogsigma;
  vector[nGroupsA] alphamA;
  vector[nGroupsB] alphamB;
  vector[nGroupsA] dalphaA;
  vector[nGroupsB] dalphaB;
  
  Omega = quad_form_diag(rho, omega[1 : 2]); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + beta[1] * (logPobs[i] - 2.2);
    miu[i, 2] = S1mHat   + beta[2] * (logPobs[i] - 2.2);
    miu[i, 3] = dS1Hat;
  }
  
   for (i in 1 : nAnalytes) {
    for (c in 1 : (nColumns-1)) {
    cmiu[c, i, 1] = clogkwHat[c] + cbeta[c,1] * (logPobs[i] - 2.2);
    cmiu[c, i, 2] = cS1mHat[c]   + cbeta[c,2] * (logPobs[i] - 2.2);
    cmiu[c, i, 3] = cdS1Hat[c];
  }}
  
  for (i in 1 : nAnalytes) { 
   logkwx[i, 1, : ]  =  paramN[i,1]*[1,1,1]';
  for (c in 1 : (nColumns-1)) {
   logkwx[i, c+1, : ]  =  (paramN[i,1]+clogkwN[c,i])*[1,1,1]'; 
  }}
   
  for (d in 1 : nGroupsA) {
    logkwx[idxGroupsA[d,1], 1, 1+idxGroupsA[d,2]] += dlogkwA[d];
    if (idxGroupsA[d,2]==1) {
   logkwx[idxGroupsA[d,1], 1, 3] += dlogkwA[d];
   }
    for (c in 1 : (nColumns-1)) {
   cdlogkwA[c,d] = cdlogkwHat[c,1] + ckappa[c,1]*etacdlogkwA[c,d];
   logkwx[idxGroupsA[d,1], c+1, 1+idxGroupsA[d,2]] += dlogkwA[d]+cdlogkwA[c,d];
   if (idxGroupsA[d,2]==1) {
   logkwx[idxGroupsA[d,1], c+1, 3] += dlogkwA[d]+cdlogkwA[c,d];
   }}}
  
  for (d in 1 : nGroupsB) {
    logkwx[idxGroupsB[d,1], 1, idxGroupsB[d,2]] += dlogkwB[d];
    if (idxGroupsB[d,2]==2) {
    logkwx[idxGroupsB[d,1], 1, 1] += dlogkwB[d];
    }
    for (c in 1 : (nColumns-1)) {
    cdlogkwB[c,d] = cdlogkwHat[c,2]+ ckappa[c,1]*etacdlogkwB[c,d];
    logkwx[idxGroupsB[d,1], c+1, idxGroupsB[d,2]] += dlogkwB[d]+cdlogkwB[c,d];
    if (idxGroupsB[d,2]==2) {
    logkwx[idxGroupsB[d,1], c+1, 1] += dlogkwB[d]+cdlogkwB[c,d];
    }}}
   
  for (i in 1 : nAnalytes) { 
    apHx[i, 1, : ] = [0,0,0]';
    for (c in 1 : (nColumns-1)) {
    apHx[i, c+1, : ] = [0,0,0]'; 
  }}
   
  for (d in 1 : nGroupsA) {
    apHx[idxGroupsA[d,1], 1, 1+idxGroupsA[d,2]] += apH[1];
  if (idxGroupsA[d,2]==1) {
   apHx[idxGroupsA[d,1], 1, 3] += apH[1];
   }
  for (c in 1 : (nColumns-1)) {
   apHx[idxGroupsA[d,1], c+1, 1+idxGroupsA[d,2]] += apH[1]+capH[c,1];
  if (idxGroupsA[d,2]==1) {
   apHx[idxGroupsA[d,1], c+1, 3] += apH[1]+capH[c,1];
   }}}
  
  for (d in 1 : nGroupsB) {
    apHx[idxGroupsB[d,1], 1, idxGroupsB[d,2]] += apH[2];
         if (idxGroupsB[d,2]==2) {
    apHx[idxGroupsB[d,1], 1, 1] += apH[2];
   }
   for (c in 1 : (nColumns-1)) {
    apHx[idxGroupsB[d,1], c+1, idxGroupsB[d,2]] += apH[2]+capH[c,2];
    if (idxGroupsB[d,2]==2) {
    apHx[idxGroupsB[d,1], c+1, 1] += apH[2]+capH[c,2];
   }}}
   
  for (i in 1 : nAnalytes) { 
   S1x[i, 1, 1, : ]  =  paramN[i,2]*[1,1,1]';
   S1x[i, 2, 1, : ]  =  paramN[i,2]*[1,1,1]' + dS1N[i]; 
  for (c in 1 : (nColumns-1)) {
   S1x[i, 1, c+1, : ]  =  (paramN[i,2]+cS1mN[c,i])*[1,1,1]' ;
   S1x[i, 2, c+1, : ]  =  (paramN[i,2]+cS1mN[c,i])*[1,1,1]' + (dS1N[i]+cdS1N[c,i]); 
  }}
   
  for (d in 1 : nGroupsA) {
   S1x[idxGroupsA[d,1], 1, 1, 1+idxGroupsA[d,2]] += dS1mA[d];
   S1x[idxGroupsA[d,1], 2, 1, 1+idxGroupsA[d,2]] += dS1mA[d]+dS1A[d]; 
    
  if (idxGroupsA[d,2]==1) {
   S1x[idxGroupsA[d,1], 1, 1, 3] += dS1mA[d];
   S1x[idxGroupsA[d,1], 2, 1, 3] += dS1mA[d]+dS1A[d];
   }
   
  for (c in 1 : (nColumns-1)) {  
   cdS1mA[c,d]= cdS1mHat[c,1] + ckappa[c,2]*etacdS1mA[c,d];
   cdS1A[c,d] = cddS1Hat[c,1] + ckappa[c,3]*etacdS1A[c,d];
   S1x[idxGroupsA[d,1], 1, c+1, 1+idxGroupsA[d,2]] += dS1mA[d]        +cdS1mA[c,d];
   S1x[idxGroupsA[d,1], 2, c+1, 1+idxGroupsA[d,2]] += dS1mA[d]+dS1A[d]+cdS1mA[c,d]+cdS1A[c,d];
  if (idxGroupsA[d,2]==1) {
   S1x[idxGroupsA[d,1], 1, c+1, 3] += dS1mA[d]        +cdS1mA[c,d];
   S1x[idxGroupsA[d,1], 2, c+1, 3] += dS1mA[d]+dS1A[d]+cdS1mA[c,d]+cdS1A[c,d];
   }}}
  
  for (d in 1 : nGroupsB) {
    S1x[idxGroupsB[d,1], 1, 1, idxGroupsB[d,2]] += dS1mB[d];
    S1x[idxGroupsB[d,1], 2, 1, idxGroupsB[d,2]] += dS1mB[d]+dS1B[d];
    if (idxGroupsB[d,2]==2) {
    S1x[idxGroupsB[d,1], 1, 1, 1] += dS1mB[d];
    S1x[idxGroupsB[d,1], 2, 1, 1] += dS1mB[d]+dS1B[d];
   }
    for (c in 1 : (nColumns-1)) {
    cdS1mB[c,d]= cdS1mHat[c,2] + ckappa[c,2]*etacdS1mB[c,d];
    cdS1B[c,d] = cddS1Hat[c,2] + ckappa[c,3]*etacdS1B[c,d];

    S1x[idxGroupsB[d,1], 1, c+1, idxGroupsB[d,2]] += dS1mB[d]        +cdS1mB[c,d];
    S1x[idxGroupsB[d,1], 2, c+1, idxGroupsB[d,2]] += dS1mB[d]+dS1B[d]+cdS1mB[c,d]+cdS1B[c,d];
  if (idxGroupsB[d,2]==2) {
    S1x[idxGroupsB[d,1], 1, c+1, 1] += dS1mB[d]        +cdS1mB[c,d];
    S1x[idxGroupsB[d,1], 2, c+1, 1] += dS1mB[d]+dS1B[d]+cdS1mB[c,d]+cdS1B[c,d];
   }}} 
  
   S2x[1,1] = 10^(logS2mHat);
   S2x[2,1] = 10^(logS2mHat + dlogS2Hat);
   
   for (c in 1 : (nColumns-1)) {
   S2x[1,c+1] = S2x[1,1];
   S2x[2,c+1] = S2x[2,1];
   }
   
  for (i in 1 : nAnalytes) { 
   dlogkTx[i, 1]  = dlogkT[i];
  for (c in 1 : (nColumns-1)) {
   cdlogkT[c,i] = cdlogkTHat[c]+comegaT[c]*etacdlogkT[c,i];
   dlogkTx[i, c+1]  =  dlogkT[i] + cdlogkT[c,i];
  }}
  
  for (i in 1 : nAnalytes) { 
   pKawx[i, : ]  = [0,0]';
   alphax[i, 1, : ]  = [0,0]'; 
   alphax[i, 2, : ]  = [0,0]'; 
  }
  
  for (d in 1 : nGroupsA) {
    
  alphamA[d] = alphamHat[1] + tau[2]*etaalphamA[d];
  dalphaA[d] = dalphaHat[1] + tau[3]*etadalphaA[d];
   
  pKawx[idxGroupsA[d,1], idxGroupsA[d,2]] = pKawA[d];
  alphax[idxGroupsA[d,1], 1, idxGroupsA[d,2]]= alphamA[d];
  alphax[idxGroupsA[d,1], 2, idxGroupsA[d,2]]= alphamA[d]+dalphaA[d];
   }
  
  for (d in 1 : nGroupsB) {
    alphamB[d] = alphamHat[2] + tau[2]*etaalphamB[d];
    dalphaB[d] = dalphaHat[2] + tau[3]*etadalphaB[d];
    pKawx[idxGroupsB[d,1], idxGroupsB[d,2]] = pKawB[d];
    alphax[idxGroupsB[d,1], 1, idxGroupsB[d,2]]= alphamB[d];
    alphax[idxGroupsB[d,1], 2, idxGroupsB[d,2]]= alphamB[d]+dalphaB[d];
   }
   
  for (i in 1 : nAnalytes) { 
   sigmax[i, 1] = exp(logsigma[i]);
   for (c in 1 : (nColumns-1)) {
   clogsigma[c,i] = clogmsigma[c] + cssigma[c]*etaclogsigma[c,i];
   sigmax[i, c+1] = exp(logsigma[i]+clogsigma[c,i]); 
  }}
}


model {

logkwHat~normal(3.6,0.0793);
S1mHat~normal(4.96,0.0832);
dS1Hat~normal(0.612,0.0512);
dlogkwHat[1]~normal(-0.786,0.0718);
dlogkwHat[2]~normal(-0.97,0.0503);
dS1mHat[1]~normal(0.176,0.121);
dS1mHat[2]~normal(0.11,0.0744);
ddS1Hat[1]~normal(0.279,0.0865);
ddS1Hat[2]~normal(-0.651,0.0564);
logS2mHat~normal(-0.308,0.0143);
dlogS2Hat~normal(0.42,0.0076);
beta[1]~normal(0.831,0.0407);
beta[2]~normal(0.483,0.0449);
alphamHat[1]~normal(2.23,0.155);
alphamHat[2]~normal(-1.36,0.102);
dalphaHat[1]~normal(0.212,0.0981);
dalphaHat[2]~normal(-0.188,0.0744);
dlogkTHat~normal(-0.0898,0.00292);
apH[1]~normal(-0.0275,0.00118);
apH[2]~normal(0.0806,0.0008);
omega[1]~normal(0.925,0.0549);
omega[2]~normal(0.941,0.059);
omega[3]~normal(0.565,0.0343);
omegaT~normal(0.0337,0.00212);
kappa[1]~normal(0.584,0.035);
kappa[2]~normal(0.694,0.0475);
kappa[3]~normal(0.57,0.0378);
msigma~normal(0.386,0.0263);
ssigma~normal(0.813,0.0496);
tau[1]~normal(0.882,0.0494);
tau[2]~normal(0.958,0.0584);
tau[3]~normal(0.793,0.0506);
clogkwHat[1]~normal(0.424,0.0115);
clogkwHat[2]~normal(0.173,0.0132);
clogkwHat[3]~normal(0.105,0.0118);
clogkwHat[4]~normal(0.166,0.0101);
cS1mHat[1]~normal(0.588,0.0198);
cS1mHat[2]~normal(-0.122,0.0218);
cS1mHat[3]~normal(0.362,0.018);
cS1mHat[4]~normal(0.479,0.0167);
cdS1Hat[1]~normal(0.151,0.0134);
cdS1Hat[2]~normal(0.828,0.0278);
cdS1Hat[3]~normal(0.552,0.0535);
cdS1Hat[4]~normal(0.0308,0.0157);
cdlogkwHat[1,1]~normal(0.0412,0.0154);
cdlogkwHat[2,1]~normal(-0.0515,0.0208);
cdlogkwHat[3,1]~normal(-0.0321,0.0167);
cdlogkwHat[4,1]~normal(0.0426,0.0164);
cdlogkwHat[1,2]~normal(0.00532,0.0135);
cdlogkwHat[2,2]~normal(-0.0312,0.0154);
cdlogkwHat[3,2]~normal(-0.0578,0.0127);
cdlogkwHat[4,2]~normal(-0.0205,0.0137);
cdS1mHat[1,1]~normal(-0.024,0.0449);
cdS1mHat[2,1]~normal(-0.202,0.0591);
cdS1mHat[3,1]~normal(-0.501,0.0452);
cdS1mHat[4,1]~normal(-0.218,0.0412);
cdS1mHat[1,2]~normal(0.375,0.0343);
cdS1mHat[2,2]~normal(-0.394,0.0371);
cdS1mHat[3,2]~normal(-0.0862,0.0296);
cdS1mHat[4,2]~normal(0.313,0.033);
cddS1Hat[1,1]~normal(0.0505,0.0267);
cddS1Hat[2,1]~normal(0.134,0.0484);
cddS1Hat[3,1]~normal(0.457,0.0879);
cddS1Hat[4,1]~normal(0.00183,0.0243);
cddS1Hat[1,2]~normal(0.17,0.0155);
cddS1Hat[2,2]~normal(0.573,0.0301);
cddS1Hat[3,2]~normal(0.504,0.0608);
cddS1Hat[4,2]~normal(-0.0301,0.0153);
cbeta[1,1]~normal(0.00719,0.00644);
cbeta[2,1]~normal(-0.00576,0.00746);
cbeta[3,1]~normal(-0.0241,0.0076);
cbeta[4,1]~normal(-0.0102,0.00595);
cbeta[1,2]~normal(-0.0664,0.0107);
cbeta[2,2]~normal(0.113,0.0127);
cbeta[3,2]~normal(-0.0198,0.0103);
cbeta[4,2]~normal(-0.00871,0.00955);
cdlogkTHat[1]~normal(-0.00595,0.00125);
cdlogkTHat[2]~normal(-0.0208,0.0016);
cdlogkTHat[3]~normal(-0.00605,0.00105);
cdlogkTHat[4]~normal(-0.00302,0.00104);
capH[1,1]~normal(-0.0149,0.00159);
capH[2,1]~normal(-0.0219,0.00165);
capH[3,1]~normal(0.00482,0.00146);
capH[4,1]~normal(-0.0226,0.00149);
capH[1,2]~normal(-0.0336,0.00118);
capH[2,2]~normal(-0.044,0.00105);
capH[3,2]~normal(-0.0506,0.000875);
capH[4,2]~normal(-0.013,0.000997);
comega[1,1]~normal(0.107,0.00755);
comega[2,1]~normal(0.129,0.0091);
comega[3,1]~normal(0.114,0.00763);
comega[4,1]~normal(0.0966,0.00656);
comega[1,2]~normal(0.0971,0.0193);
comega[2,2]~normal(0.146,0.0172);
comega[3,2]~normal(0.0901,0.0157);
comega[4,2]~normal(0.0713,0.0169);
comega[1,3]~normal(0.141,0.0108);
comega[2,3]~normal(0.309,0.0203);
comega[3,3]~normal(0.608,0.04);
comega[4,3]~normal(0.172,0.0124);
ckappa[1,1]~normal(0.0688,0.00642);
ckappa[2,1]~normal(0.116,0.00844);
ckappa[3,1]~normal(0.0855,0.00672);
ckappa[4,1]~normal(0.0873,0.00623);
ckappa[1,2]~normal(0.0714,0.0353);
ckappa[2,2]~normal(0.199,0.0302);
ckappa[3,2]~normal(0.0932,0.0314);
ckappa[4,2]~normal(0.0334,0.0238);
ckappa[1,3]~normal(0.0738,0.02);
ckappa[2,3]~normal(0.277,0.0232);
ckappa[3,3]~normal(0.717,0.0417);
ckappa[4,3]~normal(0.0861,0.0164);
comegaT[1]~normal(0.0023,0.00159);
comegaT[2]~normal(0.0115,0.00122);
comegaT[3]~normal(0.00179,0.00122);
comegaT[4]~normal(0.00114,0.000838);
clogmsigma[1]~normal(0.2,0.017);
clogmsigma[2]~normal(-0.017,0.0207);
clogmsigma[3]~normal(-0.26,0.0205);
clogmsigma[4]~normal(-0.1,0.017);
cssigma[1]~normal(0.0927,0.0236);
cssigma[2]~normal(0.169,0.02);
cssigma[3]~normal(0.161,0.019);
cssigma[4]~normal(0.0854,0.0245);
rho~lkj_corr_point_lower_tri(0.864,0.0231);
  
  for (i in 1 : nAnalytes) {
  paramN[i] ~ multi_normal(miu[i,1:2], Omega);
  }
  
  dS1N ~ normal(miu[,3], omega[3]);
  dlogkT ~ normal(dlogkTHat, omegaT);
  dlogkwA ~ normal(dlogkwHat[1], kappa[1]);
  dlogkwB ~ normal(dlogkwHat[2], kappa[1]);
  dS1mA ~ normal(dS1mHat[1], kappa[2]);
  dS1mB ~ normal(dS1mHat[2], kappa[2]);
  dS1A ~ normal(ddS1Hat[1], kappa[3]);
  dS1B ~ normal(ddS1Hat[2], kappa[3]);
  
 for (c in 1 : (nColumns-1)) {   
  clogkwN[c] ~ normal(cmiu[c, ,1], comega[c,1]);
  cS1mN[c] ~ normal(cmiu[c, ,2], comega[c,2]);
  cdS1N[c] ~ normal(cmiu[c, ,3], comega[c,3]);
  etacdlogkT[c] ~ normal(0,1);
  etacdlogkwA[c] ~ normal(0,1);
  etacdlogkwB[c] ~ normal(0,1);
  etacdS1mA[c] ~ normal(0,1);
  etacdS1mB[c] ~ normal(0,1);
  etacdS1A[c] ~ normal(0,1);
  etacdS1B[c] ~ normal(0,1);
 }
 
  pKawA ~ normal(pKaslitA, tau[1]);
  pKawB ~ normal(pKaslitB, tau[1]);
  etaalphamA ~ normal(0,1);
  etaalphamB ~ normal(0,1);
  etadalphaA ~ normal(0,1);
  etadalphaB ~ normal(0,1);

  logsigma  ~ normal(log(msigma),ssigma); 
  
  for (c in 1 : (nColumns-1)) { 
  etaclogsigma[c]  ~ normal(0,1); 
  }

}

generated quantities {
 
array[nObs] real trCond;
array[nObs] real trHatCond;
array[nexpid] vector[nAnalytessim] tr_matrix;
array[nexpid] real maxtr;
array[nexpid] real mintr;
array[nexpid] real mindifftr;

 for(z in 1:nObs){
   
  trHatCond[z] = chromgratrapz(steps[z], 
                                logkwx[analyte[z], column[z],  : ] + dlogkT[analyte[z]]  * hplcparam[z, 11],
                                apHx[analyte[z], column[z], : ],
                                S1x[analyte[z], modifier[z], column[z],  : ],
                                S2x[modifier[z], column[z]],
                                pKawx[analyte[z],  : ],
                                alphax[analyte[z],modifier[z],  : ],
                                R[analyte[z]],
                                hplcparam[z]);

  trCond[z] = student_t_rng(3, trHatCond[z], sigmax[analyte[z],column[z]]);
  tr_matrix[expid[z],analytesim[z]]=trHatCond[z];
  }
  
 
 for(z in 1:nexpid){
  vector[nAnalytessim] sorttr=sort_asc(tr_matrix[z]);
  maxtr[z] = sorttr[nAnalytessim];
  mintr[z] = sorttr[1];
  mindifftr[z] = min(sorttr[2:nAnalytessim]-sorttr[1:nAnalytessim-1]);
  }
  
}


