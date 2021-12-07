#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
  
Rcpp::List simulate_Rcpp(

  int max_age,
  int initial_IW, int initial_IM,
  int tmax, int simtlockdown, double ifr, double total_pop,
  NumericVector par, int betafun, // parameters and type of transmission function
  NumericVector siW, NumericVector siM,
  NumericVector tid, NumericVector tit
  
  ){

    Rcpp::NumericVector undetected_infected_byageW(max_age);     // create vector of undetected infected by age of infection
    Rcpp::NumericVector undetected_infected_byageM(max_age);     // create vector of undetected infected by age of infection

    Rcpp::NumericVector all_infected_byageW(max_age);   // create vector of all infected by age of infection
    Rcpp::NumericVector all_infected_byageM(max_age);   // create vector of all infected by age of infection

    Rcpp::NumericVector new_deaths_byage(max_age);     // create vector of newly dead by age of infection

    Rcpp::NumericVector all_deaths_incidence(tmax);    // vector of deaths 

    Rcpp::NumericVector all_detected_incidenceW(tmax);  // vector of detected
    Rcpp::NumericVector all_detected_incidenceM(tmax);  // vector of detected

    Rcpp::NumericVector all_incidenceW(tmax);           // vector of incidence
    Rcpp::NumericVector all_incidenceM(tmax);           // vector of incidence

    Rcpp::NumericVector all_presentingW(tmax);          // vector of people that may present to healthcare
    Rcpp::NumericVector all_presentingM(tmax);          // vector of people that may present to healthcare

    // by default the NumericVector are filled with 0s

    int i, t;
    double new_infectionsW, new_infectionsM, new_deaths, new_detectedW, new_detectedM, new_presentingW, new_presentingM, ni_detectedW, ni_detectedM, ni_deaths;
    double proba_detection;

    // INITIALISE WITH INDIVIDUALS OF ALL AGE OF INFECTION EQUALLY (OK WHEN R=1)
    for(i = 0; i < max_age; i++){

      undetected_infected_byageW[i] = initial_IW; // initialise with initial_IW infected of all ages
      undetected_infected_byageM[i] = initial_IM; // initialise with initial_IW infected of all ages

      all_infected_byageW[i] = initial_IW;
      all_infected_byageM[i] = initial_IM;

    }

    double total_infected = initial_IW * max_age + initial_IM * max_age;
    
    double R01, R02, tR0, kR0, pdetectedmin, pdetectedmax, ktimecoeff, tmidpoint;
    double transmission_advantage; // transmission advantage of mutant
    double* R0vec; R0vec = new double[tmax]; // vector to store R0 as a function of time

    // define parameter values depending on which mode we run the function on:
    if(betafun == 1){ // step function for R0 change
      if(par.length() != 7) exit(EXIT_FAILURE); 
      
      R01 = par[0];
      R02 = par[1];
      pdetectedmin = par[2];
      pdetectedmax = par[3];
      ktimecoeff = par[4];
      tmidpoint = par[5];
      transmission_advantage = par[6]; 

      // fill in R0vec vector
      for(t = 1; t <= tmax; t++){
        if(t < simtlockdown){
          R0vec[t-1] = R01;
        } else {
          R0vec[t-1] = R02;
        }
      }
    }
    if(betafun == 2){ // smooth function
      if(par.length() != 9) exit(EXIT_FAILURE); 

       R01 = par[0];
       R02 = par[1];
       tR0 = par[2];
       kR0 = par[3];
       pdetectedmin = par[4];
       pdetectedmax = par[5];
       ktimecoeff = par[6];
       tmidpoint = par[7];
       transmission_advantage = par[8]; 

       // fill in R0vec vector
      for(t = 1; t <= tmax; t++){
          R0vec[t-1] =  (R01 + (R02 - R01) / (1 + exp(- kR0 * (t - tR0))));
      }
    }
    if(betafun == 3){ // ROvec given as an argument
      if(par.length() != tmax + 5) exit(EXIT_FAILURE); 

      for(t = 1; t <= tmax; t++){
         R0vec[t-1] =  par[t-1];
      }
      pdetectedmin = par[tmax];
      pdetectedmax = par[tmax+1];
      ktimecoeff = par[tmax+2];
      tmidpoint = par[tmax+3];
      transmission_advantage = par[tmax+4]; 
    }
    
    // printf("R01 %f", R01);
    for(t = 1; t <= tmax; t++){

      // 1. compute new infections
      new_infectionsW = 0; new_infectionsM = 0;
      new_deaths = 0;
      new_detectedW = 0; new_detectedM = 0;
      new_presentingW = 0; new_presentingM = 0;

      for(i = 0; i < max_age; i++){ // for each age of infection
        
        // infections W and M:
        new_infectionsW += (total_pop - total_infected) / total_pop * R0vec[t-1] * siW[i] * undetected_infected_byageW[i];
        new_infectionsM += (total_pop - total_infected) / total_pop * R0vec[t-1] * (1 + transmission_advantage) * siM[i] * undetected_infected_byageM[i];

        
        // dead:
        ni_deaths = ifr * tid[i] * (all_infected_byageW[i] + all_infected_byageM[i]); // both detected and undetected infections die at same rate
        new_deaths += ni_deaths;
        
        // detected:
        proba_detection = (pdetectedmin + (pdetectedmax - pdetectedmin) / (1. +  exp(- ktimecoeff * (t - tmidpoint))));
        ni_detectedW = proba_detection * tit[i] * undetected_infected_byageW[i];
        ni_detectedM = proba_detection * tit[i] * undetected_infected_byageM[i];

        new_detectedW += ni_detectedW;
        new_detectedM += ni_detectedM;

        new_presentingW += tit[i] * undetected_infected_byageW[i];
        new_presentingM += tit[i] * undetected_infected_byageM[i];

        undetected_infected_byageW[i] -= ni_detectedW; // remove detected from undetected infected with age of infection i
        undetected_infected_byageM[i] -= ni_detectedM; // remove detected from undetected infected with age of infection i

        // undetected_infected_byage[i] -= ni_deaths; // remove dead from infected with age of infection i -> we neglect this now as deaths are rare and occur after a long time

      } // end of loop on age
      
      // 2. update undetected_infected_byage accordingly
      for(i = max_age - 1; i > 0; i--){
        undetected_infected_byageW[i] = undetected_infected_byageW[i-1]; // increase age of infection
        undetected_infected_byageM[i] = undetected_infected_byageM[i-1]; // increase age of infection
        all_infected_byageW[i] = all_infected_byageW[i-1];
        all_infected_byageM[i] = all_infected_byageM[i-1];
      }
      undetected_infected_byageW[0] = new_infectionsW;       // update new infections at age of infection 0
      undetected_infected_byageM[0] = new_infectionsM;       // update new infections at age of infection 0

      all_infected_byageW[0] = new_infectionsW;   // update new infections at age of infection 0
      all_infected_byageM[0] = new_infectionsM;   // update new infections at age of infection 0

      total_infected += new_infectionsW + new_infectionsM;         // update total infected
      
      // 3. record in vectors
      all_deaths_incidence[t-1] = new_deaths;
      all_presentingW[t-1] = new_presentingW;
      all_presentingM[t-1] = new_presentingM;
      all_detected_incidenceW[t-1] = new_detectedW;
      all_detected_incidenceM[t-1] = new_detectedM;
      all_incidenceW[t-1] = new_infectionsW;
      all_incidenceM[t-1] = new_infectionsM;

    } // end of loop on time

    // return results:
    return(List::create(

    Named("undetected_infected_byageW") = undetected_infected_byageW,
    Named("undetected_infected_byageM") = undetected_infected_byageM,

    Named("all_infected_byageW") = all_infected_byageW,
    Named("all_infected_byageM") = all_infected_byageM,

    Named("all_deaths_incidence") = all_deaths_incidence,

    Named("all_detected_incidenceW") = all_detected_incidenceW,
    Named("all_detected_incidenceM") = all_detected_incidenceM,

    Named("all_presentingW") = all_presentingW,
    Named("all_presentingM") = all_presentingM,

    Named("all_incidenceW") = all_incidenceW,
    Named("all_incidenceM") = all_incidenceM,

    Named("total_infected") = total_infected
    ));
}

