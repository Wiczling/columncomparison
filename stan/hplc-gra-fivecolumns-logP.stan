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
    cumki2 = td*invki1; 
    
    bo = 0.001 / td;
    
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
      cumki2 += cki_b[1];
      bo = cki_b[2];
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
                   array[,] real dlogkT,
                   array[] vector sigma) {
                     
    real lp = 0;

    for (z in start : end) {
      
    real y_hat = chromgratrapz(steps[z], 
                                logkw[analyte[z], column[z],  : ] + dlogkT[analyte[z],column[z]]  * hplcparam[z, 11],
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
  int nAnalytes;            // number of analytes
  int nColumns;             // number of columns
  int nModifiers;           // number of org. modifiers
  int nObs;                 // number of observations
  array[nObs] int analyte;   // analyte indexes
  array[nObs] int column;   // column indexes
  array[nObs] int modifier;  // modifier indexes
  array[nObs] int<lower=1> steps;  // steps for gradient retention time aproximation
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
}
transformed data {
  int grainsize = 1;
  array[nObs] int ind = rep_array(1, nObs);
  int Nsim = 5;
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
  logkwHat ~ normal(2.2, 2);
  S1mHat ~ normal(4, 1);
  dS1Hat ~ normal(1, 1);
  dlogkwHat ~ normal(-1, 0.125);
  dS1mHat   ~ normal(0, 0.5);
  ddS1Hat   ~ normal(0, 0.25);
  logS2mHat ~ normal(-0.7, 0.125);
  dlogS2Hat ~ normal(1, 0.125);
  beta[{1}] ~ normal(1, 0.125);
  beta[{2}] ~ normal(0.5, 0.5);
  dlogkTHat ~ normal(-0.087, 0.022);
  apH       ~ normal(0, 0.1);
  
  alphamHat[{1}] ~ normal(2, 0.25);
  alphamHat[{2}] ~ normal(-1, 0.25);
  dalphaHat ~ normal(0, 0.125);
  tau[{1}] ~ normal(0, 0.25);
  tau[{2,3}] ~ normal(0, 0.125);
  omega ~ normal(0, 2);
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
  omegaT ~ normal(0, 0.022);
  kappa ~ normal(0, 0.25);
    
  clogkwHat ~ normal(0, 1);
  cS1mHat ~ normal(0, 0.5);
  cdS1Hat ~ normal(0, 0.5);
  to_vector(cdlogkwHat) ~ normal(0, 0.0625);
  to_vector(cdS1mHat) ~ normal(0, 0.25);
  to_vector(cddS1Hat) ~ normal(0, 0.125); 
  to_vector(cbeta) ~ normal(0, 0.25);
  cdlogkTHat ~ normal(0, 0.011);
  to_vector(capH) ~ normal(0, 0.05);
  to_vector(comega) ~ normal(0, 1);
  comegaT ~ normal(0, 0.011);
  to_vector(ckappa) ~ normal(0, 0.125);

  
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
  
  msigma ~ normal(0,1);
  ssigma ~ normal(0,1);
  logsigma  ~ normal(log(msigma),ssigma); 
  
  clogmsigma ~ normal(0,0.125);
  cssigma ~ normal(0,0.125);
  
  for (c in 1 : (nColumns-1)) { 
  etaclogsigma[c]  ~ normal(0,1); 
  }
  
  if (run_estimation == 1) {

    target += reduce_sum(partial_sum, ind, grainsize, trobs, steps,
                         hplcparam, analyte, column, modifier, R, logkwx, apHx,
                         S1x, S2x, pKawx, alphax, dlogkTx, sigmax);
  }
}


generated quantities {
 
vector[Nsim] logPsim =[0,1.5,3,4.5,6]';
 array[Nsim] vector[3] miuPred;
 array[nColumns-1,Nsim] vector[3] cmiuPred;
  
// 1st column
  array[Nsim] vector[2] paramNPred;
  vector[Nsim] dS1NPred;
  vector[Nsim] dlogkTPred;
  vector[Nsim] dlogkwAPred;
  vector[Nsim] dlogkwBPred;
  vector[Nsim] dS1mAPred;
  vector[Nsim] dS1mBPred;
  vector[Nsim] dS1APred;
  vector[Nsim] dS1BPred;
  
// 2nd column
  array[nColumns-1] vector[Nsim] clogkwNPred;
  array[nColumns-1] vector[Nsim] cS1mNPred;
  array[nColumns-1] vector[Nsim] cdS1NPred;
  array[nColumns-1] vector[Nsim] cdlogkTPred;
  array[nColumns-1] vector[Nsim] cdlogkwAPred;
  array[nColumns-1] vector[Nsim] cdlogkwBPred;
  array[nColumns-1] vector[Nsim] cdS1mAPred;
  array[nColumns-1] vector[Nsim] cdS1mBPred;
  array[nColumns-1] vector[Nsim] cdS1APred;
  array[nColumns-1] vector[Nsim] cdS1BPred;

  array[Nsim,nColumns]   vector[3] logkwxPred;
  array[Nsim,nModifiers, nColumns] vector[3] S1xPred;
  array[Nsim,nColumns]  real S2xPred;
  
   for (i in 1 : Nsim) {
    miuPred[i, 1] = logkwHat + beta[1] * (logPsim[i] - 2.2);
    miuPred[i, 2] = S1mHat   + beta[2] * (logPsim[i] - 2.2);
    miuPred[i, 3] = dS1Hat;
    for (c in 1 : (nColumns-1)) {
    cmiuPred[c,i, 1] = clogkwHat[c] + cbeta[c,1] * (logPsim[i] - 2.2);
    cmiuPred[c,i, 2] = cS1mHat[c]   + cbeta[c,2] * (logPsim[i] - 2.2);
    cmiuPred[c,i, 3] = cdS1Hat[c];
  }}
  
  for(i in 1:Nsim){
    paramNPred[i] = multi_normal_rng(miuPred[i,1:2],Omega); 
    dS1NPred[i]   = normal_rng(miuPred[i,3], omega[3]);
     
    dlogkwAPred[i] = normal_rng(dlogkwHat[1],kappa[1]);
    dS1mAPred[i] = normal_rng(dS1mHat[1],kappa[2]);
    dS1APred[i] = normal_rng(ddS1Hat[1],kappa[3]);
    
    dlogkwBPred[i] = normal_rng(dlogkwHat[2],kappa[1]);
    dS1mBPred[i] = normal_rng(dS1mHat[2],kappa[2]);
    dS1BPred[i] = normal_rng(ddS1Hat[2],kappa[3]);
  }
  
  for(i in 1:Nsim){
    for (c in 1 : (nColumns-1)) {
    clogkwNPred[c,i]   = normal_rng(cmiuPred[c,i,1], comega[c,1]);
    cS1mNPred[c,i]     = normal_rng(cmiuPred[c,i,2], comega[c,2]);
    cdS1NPred[c,i]   = normal_rng(cmiuPred[c,i,3], comega[c,3]);
    cdlogkwAPred[c,i] = normal_rng(cdlogkwHat[c,1], ckappa[c,1]);
    cdS1mAPred[c,i] = normal_rng(cdS1mHat[c,1], ckappa[c,2]);
    cdS1APred[c,i] = normal_rng(cddS1Hat[c,1], ckappa[c,3]);
    cdlogkwBPred[c,i] = normal_rng(cdlogkwHat[c,2], ckappa[c,1]);
    cdS1mBPred[c,i] = normal_rng(cdS1mHat[c,2], ckappa[c,2]);
    cdS1BPred[c,i] = normal_rng(cddS1Hat[c,2], ckappa[c,3]);
  }}
  
  for (i in 1 : Nsim) { 
   logkwxPred[i, 1, : ]  =  paramNPred[i,1]*[1,1,1]';
   logkwxPred[i, 1, 3] += dlogkwAPred[i];
   logkwxPred[i, 1, 1] += dlogkwBPred[i];
  for (c in 1 : (nColumns-1)) {
   logkwxPred[i, c+1, : ]  =  (paramNPred[i,1]+clogkwNPred[c,i])*[1,1,1]'; 
   logkwxPred[i, c+1, 3] += dlogkwAPred[i]+cdlogkwAPred[c,i];
   logkwxPred[i, c+1, 1] += dlogkwBPred[i]+cdlogkwBPred[c,i];
   }}
   
  for (i in 1 : Nsim) { 
   S1xPred[i, 1, 1, : ]  =  paramNPred[i,2]*[1,1,1]';
   S1xPred[i, 2, 1, : ]  =  paramNPred[i,2]*[1,1,1]' + dS1NPred[i]; 
   
   S1xPred[i, 1, 1, 3] += dS1mAPred[i];
   S1xPred[i, 2, 1, 3] += dS1mAPred[i]+dS1APred[i]; 
   
   S1xPred[i, 1, 1, 1] += dS1mBPred[i];
   S1xPred[i, 2, 1, 1] += dS1mBPred[i]+dS1BPred[i];
   
  for (c in 1 : (nColumns-1)) {
   S1xPred[i, 1, c+1, : ]  =  (paramNPred[i,2]+cS1mNPred[c,i])*[1,1,1]' ;
   S1xPred[i, 2, c+1, : ]  =  (paramNPred[i,2]+cS1mNPred[c,i])*[1,1,1]' + (dS1NPred[i] + cdS1NPred[c,i]);

   S1xPred[i, 1, c+1, 3] += dS1mAPred[i]        +cdS1mAPred[c,i];
   S1xPred[i, 2, c+1, 3] += dS1mAPred[i]+dS1APred[i]+cdS1mAPred[c,i]+cdS1APred[c,i];

   S1xPred[i, 1, c+1, 1] += dS1mBPred[i]        +cdS1mBPred[c,i];
   S1xPred[i, 2, c+1, 1] += dS1mBPred[i]+dS1BPred[i]+cdS1mBPred[c,i]+cdS1BPred[c,i];
}}
  
   S2xPred[1,1] = 10^(logS2mHat);
   S2xPred[2,1] = 10^(logS2mHat + dlogS2Hat);
   for (c in 1 : (nColumns-1)) {
   S2xPred[1,c+1] = S2xPred[1,1];
   S2xPred[2,c+1] = S2xPred[2,1];
   }

}
