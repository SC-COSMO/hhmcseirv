#' The household model for a multi-compartment SEIR (HH-MC_SEIR) epidemic model.
#'
#' \code{hh_mc_seir_out} implements the household model for a multi-compartment 
#' SEIR (HH-MC_SEIR) epidemic model.
#' @param l_params_all List with all parameters of decision model
#' @return 
#' A data.frame with population for each class over time.
#' @export
hh_mc_seir <- function(l_params_all){
  df_out_hh_mc_seir <- as.data.frame(deSolve::ode(y = l_params_all$state0, 
                                                  times = l_params_all$times, 
                                                  func = hh_mc_seir_dx,
                                                  parms = l_params_all))
  
  # df_out_cosmo <- as.data.frame(deSolve::ode(y = l_params_all$v_states_init, 
  #                                            times = l_params_all$v_times, 
  #                                            func = ifelse(l_params_all$comp, 
  #                                                          yes = cosmo_dXdt_comp, 
  #                                                          no = cosmo_dxdt), 
  #                                            parms = l_params_all))
  l_out_hh_mc_seir <- list(df_out_hh_mc_seir = df_out_hh_mc_seir,
                           l_params_all = l_params_all)
  return(l_out_hh_mc_seir)
}

#### Calculate number of household possibilities ####
hh_mc_seir_out <- function(parameters){
  with(as.list(parameters),{
    v_names_sus    <- c("S") 
    v_names_exp    <- paste("E", letters[seq(1, n_exp_states)], sep = "")
    v_names_inf    <- paste("I", letters[seq(1, n_inf_states)], sep = "")
    v_names_inf_dx <- paste("IDX", letters[seq(1, n_inf_states)], sep = "")
    v_names_rec    <- c("R") 
    v_names_vax    <- c("V") 
    v_names_states <- c("S", 
                        v_names_exp,
                        v_names_inf,
                        "R",
                        "V")
    v_names_states_tot <- c("S", 
                            v_names_exp,
                            v_names_inf,
                            v_names_inf_dx,
                            "R",
                            "V")
    n_states     <- length(v_names_states) # Not including the number of IDX
    n_states_tot <- length(v_names_states_tot) # Including the number of IDX
    
    n_hh_mod <- factorial(n_hhsize + (n_states-1))/(factorial(n_hhsize)*factorial(n_states-1)) #l_params_all
    
    ### Matrix with possible combinations of household members
    df_possibilities <- gen_hh_n(n_hhsize = n_hhsize, 
                                 v_names_states = v_names_states)
    m_possibilities  <- as.matrix(df_possibilities)
    
    #### Naming vectors ####
    v_hh_names <- as.matrix(tidyr::unite(df_possibilities, col = "names", sep = ""))
    # paste("HH",m_possibilities[, 1:n_states], m_possibilities[,2], sep = "")
    ### Names of household members by class names
    v_names_hh  <- paste("HH", v_hh_names, sep = "")
    ### Derivative names of household members by class
    v_names_dhh <- paste("dHH", v_hh_names, sep = "")
    ### All disease states including household states ####
    v_names_states_all <- c(v_names_states_tot, v_names_hh)
    
    #### Initialize state vector ####
    startingI  <- n_inf/n_pop_size
    ### Household state vector
    v_HH0 <- numeric(length = n_hh_mod)
    names(v_HH0) <- v_names_hh
    v_HH0[1] <- 1 - n_hhsize*startingI
    v_HH0[m_possibilities[, "Ia"] == 1][1] <- n_hhsize*startingI
    # v_HH0[3] <- n_hhsize*startingI
    ### Full state vector
    state0 <- c(S = (1-startingI)*n_pop_size,
                rep(0, n_exp_states),
                Ia = startingI*n_pop_size,
                rep(0, (n_inf_states-1)),
                rep(0, n_inf_states),
                R = 0,
                V = 0,
                v_HH0*n_pop_size/n_hhsize,
                Infcomm = 0,
                Infhh   = 0
    )
    names(state0) <- c(v_names_states_all, "Infcomm", "Infhh")
    
    # list2env(as.list(state0), envir = .GlobalEnv)
    
    ### Indexing vectors for household model
    ## Column index for susceptibles
    v_index_hh_sus <- which(v_names_states %in% v_names_sus)
    ## Column index for exposed
    v_index_hh_exp <- which(v_names_states %in% v_names_exp)
    ## Column index for infectious
    v_index_hh_inf <- which(v_names_states %in% v_names_inf)
    ## Column index for recovered
    v_index_hh_rec <- which(v_names_states %in% v_names_rec)
    ## Column index for vaccinated
    v_index_hh_vax <- which(v_names_states %in% v_names_vax)
    
    #### Household matrices ####
    l_transition_matrices <- gen_household_matrices_mc_seir(n_hhsize = n_hhsize,
                                                            n_hh_mod = n_hh_mod,
                                                            v_names_states = v_names_states,
                                                            v_names_exp = v_names_exp,
                                                            v_names_inf = v_names_inf,
                                                            v_names_rec = v_names_rec,
                                                            v_names_vax = v_names_vax,
                                                            v_index_hh_sus = v_index_hh_sus,
                                                            v_index_hh_exp = v_index_hh_exp,
                                                            v_index_hh_inf = v_index_hh_inf,
                                                            v_index_hh_rec = v_index_hh_rec,
                                                            v_index_hh_vax = v_index_hh_vax,
                                                            df_possibilities = df_possibilities,
                                                            r_sigma = r_sigma,
                                                            r_gamma = r_gamma,
                                                            r_omega = r_omega,
                                                            r_vax_omega = r_vax_omega)
    m_comm_trans    <- l_transition_matrices$m_comm_trans
    m_hh_trans      <- l_transition_matrices$m_hh_trans
    m_hh_prog       <- l_transition_matrices$m_hh_prog
    m_hh_recov      <- l_transition_matrices$m_hh_recov
    m_hh_waning     <- l_transition_matrices$m_hh_waning
    m_hh_vax        <- l_transition_matrices$m_hh_vax
    m_hh_waning_vax <- l_transition_matrices$m_hh_waning_vax
    
    l_params <- c(as.list(parameters),
                  list(state0 = state0,
                       v_names_states = v_names_states,
                       v_names_states_tot = v_names_states_tot,
                       v_names_sus    = v_names_sus,
                       v_names_exp    = v_names_exp,
                       v_names_inf    = v_names_inf,
                       v_names_inf_dx = v_names_inf_dx,
                       v_names_rec    = v_names_rec,
                       v_names_vax    = v_names_vax,
                       v_names_hh = v_names_hh,
                       n_states     = n_states, 
                       n_states_tot = n_states_tot, 
                       df_possibilities = df_possibilities,
                       m_possibilities  = m_possibilities,
                       v_index_hh_sus = v_index_hh_sus,
                       v_index_hh_exp = v_index_hh_exp,
                       v_index_hh_inf = v_index_hh_inf,
                       v_index_hh_rec = v_index_hh_rec,
                       m_comm_trans    = m_comm_trans,
                       m_hh_trans      = m_hh_trans  ,
                       m_hh_prog       = m_hh_prog   ,
                       m_hh_recov      = m_hh_recov  ,
                       m_hh_waning     = m_hh_waning,
                       m_hh_vax        = m_hh_vax,
                       m_hh_waning_vax = m_hh_waning_vax
                  )
    )
    l_out <- hh_mc_seir(l_params_all = l_params)
    return(l_out)
  }
  )
  
}

#' Derivatives of the household model for a multi-compartment SEIR (HH-MC_SEIR)
#' epidemic model.
#'
#' \code{hh_mc_seir_dx} computes the derivatives of the household model for 
#' a multi-compartment SEIR (HH-MC_SEIR) epidemic model.
#' 
#' @param time Time at which the derivatives will be computed
#' @param v_pop Vector with the population by class
#' @param l_params_all List with all parameters of SC-COSMO model
#' @return 
#' A list with a vector of derivatives for each class .
#' @export
hh_mc_seir_dx <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # state <- state0
    # list2env(as.list(state), envir = .GlobalEnv)
    S     <- state[v_names_sus]
    v_E   <- state[v_names_exp]
    v_I   <- state[v_names_inf]
    v_IDX <- state[v_names_inf_dx]
    E   <- sum(v_E)
    I   <- sum(v_I)
    IDX <- sum(v_IDX)
    R   <- state[v_names_rec]
    V   <- state[v_names_vax]
    
    ### Total pop
    N = S + E + I + IDX + R + V
    
    ### vector with household members
    # v_HH <- state[(n_states_tot+1):length(state)]
    v_HH <- (state[(n_states_tot+1):(length(state)-2)]/N)*n_hhsize # The "-2" is because we added two compartments: incidence from community and household
    
    # n_err <- sum(v_HH[v_HH<0])
    # v_HH[v_HH<0] = 0
    # v_HH[length(v_HH)] <- -n_err
    
    ### Proportion of infections that are DX
    p_dx <- IDX/(I + IDX)
    
    ### NPI effect
    p_npi <- get_npi(n_time = t, parameters = parameters)
    r_beta_current <- r_beta*n_contacts_comm*p_npi
    ### Vax rate
    r_vax <- get_vax(n_time = t, parameters = parameters)
    
    ### Adjust transmission matrices to account for detection
    v_lambda_hh <- (((I)/N)*(r_beta_current*(1 - p_dx)) + 
                      ((IDX)/N)*(r_beta_current*p_alpha_dx*p_dx))
    m_comm_trans_dx <- m_comm_trans * as.numeric(v_lambda_hh)
    r_tau_avg         <- r_tau*(1 - p_dx) + r_tau*p_alpha_dx*p_dx
    m_hh_trans_dx   <- m_hh_trans*r_tau_avg
    
    ### Adjust vaccination matrix to account for vax rate and effectiveness
    m_hh_vax_rate_eff <- m_hh_vax * (r_vax * eff_vax)
    
    ### Household infection
    household_infection_rate <- gen_household_transmission_mc_seir(r_tau = r_tau_avg, #r_tau, #
                                                                   n_hhsize = n_hhsize,
                                                                   n_contacts_hh = n_contacts_hh,
                                                                   v_HH = v_HH,
                                                                   v_index_hh_sus = v_index_hh_sus,
                                                                   v_index_hh_exp = v_index_hh_exp,
                                                                   v_index_hh_inf = v_index_hh_inf,
                                                                   v_index_hh_rec = v_index_hh_rec,
                                                                   m_possibilities = m_possibilities)
    ## Force of infection from household to community
    n_household_infection_rate <- as.numeric(household_infection_rate)
    print(c(t, n_household_infection_rate))
    
    ### Force of infection
    n_lambda  <- r_beta_current*(I/N) + r_beta_current*p_alpha_dx*(IDX/N)
    
    ### Births
    birth_int <- r_birth*N #+ r_gamma*n_inf_states*p_death_inf*(v_I[n_inf_states] + v_IDX[n_inf_states])
    
    dS <- birth_int +                                       # Inflows (Births)
      r_omega*R +                                           # Inflows 
      r_vax_omega * V -                                     # Coming from V through vaccine waning immunity
      n_lambda*S -                                          # Outflows
      n_household_infection_rate*S -                        # Outflows # previously multiplied by N!!!
      (r_vax*eff_vax) * S -                                 # Outflows
      r_death*S                                             # Outflows
    v_d_E <- c(n_lambda*S + n_household_infection_rate*S,     # Inflows # previously multiplied by N!!!
               (r_sigma*n_exp_states)*v_E[-n_exp_states]) - # Inflows
      (r_sigma*n_exp_states)*v_E -                          # Outflows
      r_death*v_E                                           # Outflows
    v_d_I <- c((r_sigma*n_exp_states)*v_E[n_exp_states],    # Inflows
               (r_gamma*n_inf_states)*v_I[-n_inf_states]) - # Inflows
      (r_gamma*n_inf_states)*v_I -                          # Outflows
      r_dx*v_I -                                            # Outflows
      r_death*v_I
    v_d_IDX <- r_dx*v_I +                                 # Inflows 
      c(0, (r_gamma*n_inf_states)*v_IDX[-n_inf_states]) - # Inflows
      (r_gamma*n_inf_states)*v_IDX -                      # Outflows
      r_death*v_IDX                                       # Outflows
    dR <- (r_gamma*n_inf_states*(1-p_death_inf))*v_I[n_inf_states] + # Inflows
      (r_gamma*n_inf_states*(1-p_death_inf))*v_IDX[n_inf_states] -    # Inflows
      r_omega*R -                                     # Outflows
      r_death*R                                       # Outflows
    dV <- (r_vax*eff_vax) * S - # Incoming S to vaccination
      r_vax_omega * V -         # Leaving V through vaccine waning immunity
      r_death * V               # Leaving V through background mortality
    dD <- N*r_death + 
      r_gamma*n_inf_states*p_death_inf*(v_I[n_inf_states] + v_IDX[n_inf_states])
    
    ## 
    dInfcomm <- n_lambda*S
    dInfhh   <- n_household_infection_rate*N
    
    ### Total death rate
    r_death_tot <- dD/N
    ### Births
    ## Total birth rate
    r_birth_tot <- r_death_tot
    ## Bring offspring to susceptible states
    v_index_keep_sus <- (m_possibilities[, v_index_hh_sus] > 0)
    v_tot_births     <- rep(0, length(v_HH))
    names(v_tot_births) <- v_names_hh
    
    # ### I ADDED THIS 
    # if(sum(v_HH[v_index_keep_sus]) == 0){
    #   v_tot_births[v_names_hh[1]] <- r_birth_tot
    # } else{
    #   v_tot_births[v_index_keep_sus] <- (v_HH[v_index_keep_sus]/sum(v_HH[v_index_keep_sus]))*r_birth_tot
    # }
    # ### I ADDED THIS
    v_tot_births[v_index_keep_sus] <- v_HH[v_index_keep_sus]/sum(v_HH[v_index_keep_sus])*r_birth_tot
    
    # ### I ADDED THIS
    # v_index_keep_alive  <- v_HH > 0 
    # v_tot_deaths        <- rep(0, length(v_HH))
    # names(v_tot_deaths) <- v_names_hh
    # if(sum(v_HH[v_index_keep_alive]) == 0){
    #   v_tot_deaths[v_names_hh[1]] <- r_death_tot
    # } else{
    #   v_tot_deaths[v_index_keep_alive] <- (v_HH[v_index_keep_alive]/sum(v_HH[v_index_keep_alive]))*r_death_tot
    # }
    # ### I ADDED THIS
    
    #### Rates of change in within household epidemics ####
    v_dHH  <- v_tot_births +
      t(m_comm_trans_dx) %*% v_HH + 
      t(m_hh_trans_dx)   %*% v_HH + 
      t(m_hh_prog)       %*% v_HH + 
      t(m_hh_recov)      %*% v_HH +
      t(m_hh_waning)     %*% v_HH +
      t(m_hh_vax_rate_eff)   %*% v_HH +
      t(m_hh_waning_vax) %*% v_HH +
      - r_death_tot       *  v_HH
    # - v_tot_deaths ### I ADDED THIS
    
    ### I ADDED THIS
    v_dHH <- v_dHH*N/n_hhsize
    ### I ADDED THIS
    
    m_test <- cbind(v_tot_births,
                    t(m_comm_trans_dx)   %*% v_HH,
                    t(m_hh_trans_dx)     %*% v_HH,
                    t(m_hh_prog)         %*% v_HH,
                    t(m_hh_recov)        %*% v_HH,
                    t(m_hh_waning)       %*% v_HH,
                    t(m_hh_vax_rate_eff) %*% v_HH,
                    t(m_hh_waning_vax)   %*% v_HH,
                    - r_death_tot       *  v_HH)
    sum(rowSums(m_test))
    
    # return the rate of change
    return(list(c(dS, 
                  v_d_E,
                  v_d_I,
                  v_d_IDX,
                  dR,
                  dV,
                  v_dHH,
                  dInfcomm,
                  dInfhh
    )))
  }) # end with(as.list ...
}

#' Generate household possibilities
#' 
#' \code{gen_hh_n} generates household transmission (i.e., 
#' force of infection).
#' 
#' @param n_hhsize Household size
#' @param v_names_states vector with state names
#' @return 
#' A data.frame with possible combinations of household members within each 
#' epidemic compartment.
gen_hh_n <- function(n_hhsize, v_names_states){
  n_states <- length(v_names_states)
  df_possibilities <- data.frame(id = 1, V1 = 0:n_hhsize)
  for(i in 1:(n_states - 1)){ # i <- 2
    df_possibilities <- left_join(df_possibilities, data.frame(id = 1, temp = 0:n_hhsize), by = "id")
    colnames(df_possibilities)[i+2] <- paste0("V", i+1)
  }
  df_possibilities <-  df_possibilities %>%
    select(-id)
  v_names_cols <- colnames(df_possibilities)
  df_possibilities$new_hhsize <- rowSums(df_possibilities)
  df_possibilities <- df_possibilities %>% 
    filter(new_hhsize == n_hhsize) %>%
    select(-new_hhsize) %>% 
    group_by_all() %>%
    summarise_all(list(mean)) %>%
    ungroup() %>%
    arrange_all(desc)
  colnames(df_possibilities) <- v_names_states
  return(df_possibilities)
}

#' Generate household transmission
#' 
#' \code{gen_household_transmission} generates household transmission (i.e., 
#' force of infection).
#' 
#' @param v_HH State vector of within household epidemics
#' @param v_index_hh_sus Indexing column vector of susceptibles
#' @param v_index_hh_inf Indexing column vector of infectious
#' @param v_index_hh_rec Indexing column vector of recovered
#' @param m_possibilities Matrix with possible combinations of household 
#' members within each epidemic compartment
#' @return 
#' A scalar with the force of infection
gen_household_transmission_mc_seir <- function(r_tau,
                                               n_hhsize,
                                               n_contacts_hh,
                                               v_HH,
                                               v_index_hh_sus,
                                               v_index_hh_exp,
                                               v_index_hh_inf,
                                               v_index_hh_rec,
                                               m_possibilities){
  #### Indexing vectors ####
  ## Indices with active susceptibles and infectious
  v_index_keep_tau <- (m_possibilities[, v_index_hh_sus] > 0 & 
                         rowSums(m_possibilities[, v_index_hh_inf, 
                                                 drop = FALSE]) > 0)
  #### Estimate within HH transmission ####
  ### Vector with combinations that can infect and get infected
  v_hh_size_inf <- m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
    rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop = FALSE])
  
  n_hh_inf_rate <- (r_tau * n_contacts_hh) * 
    (v_hh_size_inf %*% v_HH[v_index_keep_tau])/n_hhsize
  return(n_hh_inf_rate)
}

#' Generate household matrices
#' 
#' \code{gen_household_matrices} generates transition matrices for within 
#' household epidemics.
#' @param n_hhsize Household size
#' @param n_hh_mod Number of household model states
#' @param v_names_states Vector with all names
#' @param v_names_exp Vector with names of exposed states
#' @param v_names_inf Vector with names of infectious states
#' @param v_names_rec Vector with names of recovered states
#' @param v_names_vax Vector with names of vaccination states
#' @param v_index_hh_sus Indexing column vector of susceptibles
#' @param v_index_hh_exp Indexing column vector of exposed
#' @param v_index_hh_inf Indexing column vector of infectious
#' @param v_index_hh_rec Indexing column vector of recovered
#' @param v_index_hh_vax Indexing column vector of vaccination
#' @param m_possibilities Matrix with possible combinations of household 
#' members within each epidemic compartment
#' @param r_sigma  Daily rate of progression of exposed individuals 
#' (Latent period)
#' @param r_gamma Daily rate of recovery of infectious individuals 
#' (Infectiousness period).
#' @param r_omega Waning rate for the household model members within each 
#' epidemic compartment
#' @param r_vax_omega Vaccination waning rate for the household model members
#' @return 
#' A list with three transition matrices: 1) community transmission matrix; 2) 
#' household transmission matrix; and 3) household recovered matrix
gen_household_matrices_mc_seir <- function(n_hhsize,
                                           n_hh_mod,
                                           v_names_states,
                                           v_names_exp,
                                           v_names_inf,
                                           v_names_rec,
                                           v_names_vax,
                                           v_index_hh_sus,
                                           v_index_hh_exp,
                                           v_index_hh_inf,
                                           v_index_hh_rec,
                                           v_index_hh_vax, #
                                           df_possibilities,
                                           r_sigma,
                                           r_gamma,
                                           r_omega,
                                           r_vax_omega #
){
  m_possibilities <- as.matrix(df_possibilities)
  #### Naming vectors ####
  v_hh_names <- as.matrix(unite(df_possibilities, col = "names", sep = ""))
  # paste("HH",m_possibilities[, 1:n_states], m_possibilities[,2], sep = "")
  ### Names of household members by class names
  v_names_hh  <- paste("HH", v_hh_names, sep = "")
  n_hh_mod    <- length(v_names_hh)
  ### Derivative names of household members by class
  v_names_dhh <- paste("dHH", v_hh_names, sep = "")
  
  n_exp_states <- length(v_names_exp)
  n_inf_states <- length(v_names_inf)
  
  v_names_exp_inf <- c(v_names_exp, v_names_inf)
  v_names_inf_rec <- c(v_names_inf, v_names_rec)
  
  #### Indexing vectors ####
  ## Column index for susceptibles and infectious
  v_index_sus_inf_hh <- c(v_index_hh_sus, v_index_hh_inf)
  ## Column index for exposed and infectious
  v_index_exp_inf_hh <- c(v_index_hh_exp, v_index_hh_inf)
  ## Column index for infectious and recovered
  v_index_inf_rec_hh <- c(v_index_hh_inf, v_index_hh_rec)
  ## Column index for susceptibles, exposed and infectious
  v_index_sus_exp_inf_hh <- c(v_index_hh_sus, v_index_hh_exp, v_index_hh_inf)
  ## Column index for susceptibles and vaccinated
  v_index_sus_vax_hh <- c(v_index_hh_sus, v_index_hh_vax)
  
  ## Indices with active susceptibles
  v_index_keep_sus    <- (m_possibilities[, v_index_hh_sus] > 0)
  ## Indices with active exposed
  v_index_keep_exp       <- rowSums((m_possibilities[, v_index_hh_exp, drop=FALSE] > 0)) > 0 
  v_index_keep_exp_first <- (m_possibilities[, v_index_hh_exp[1]] > 0)
  v_index_keep_exp_last  <- (m_possibilities[, v_index_hh_exp[n_exp_states]] > 0)
  m_index_keep_exp       <- as.matrix((m_possibilities[, v_index_hh_exp] > 0), 
                                      ncol = n_exp_states) # matrix of indices to keep exposed by class
  ## Indices with active infectious
  v_index_keep_inf       <- rowSums((m_possibilities[, v_index_hh_inf, drop=FALSE] > 0)) > 0 
  v_index_keep_inf_first <- (m_possibilities[, v_index_hh_inf[1]] > 0)
  v_index_keep_inf_last  <- (m_possibilities[, v_index_hh_inf[n_inf_states]] > 0)
  m_index_keep_inf       <- as.matrix((m_possibilities[, v_index_hh_inf] > 0), 
                                      ncol = n_inf_states) # matrix of indices to keep exposed by class
  
  ## Indices with active exposed and infectious
  m_index_keep_exp_inf <- (m_possibilities[, v_index_exp_inf_hh] > 0) # matrix of indices to keep exposed by classv_index_exp_inf_hh
  ## Indices with active infectious and recovered
  m_index_keep_inf_rec <- (m_possibilities[, v_index_inf_rec_hh] > 0) # matrix of indices to keep exposed by classv_index_inf_rec_hh
  ## Indices with active exposed for HH transmission
  v_index_keep_exp_hh <- ((rowSums(m_possibilities[, v_index_hh_exp, drop=FALSE]) > 0) & 
                            (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  v_index_keep_exp_hh_first <- ((m_possibilities[, v_index_hh_exp[1]] > 0) & 
                                  (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  ## Indices with active infectious for HH transmission
  v_index_keep_inf_hh <- (m_possibilities[, v_index_hh_inf] > 1) #### NOT USED
  ## Indices with active recovered
  v_index_keep_rec    <- (m_possibilities[, v_index_hh_rec] > 0)
  ## Indices with active susceptibles and infectious
  v_index_keep_tau <- ((m_possibilities[, v_index_hh_sus] > 0) & 
                         (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  ## Indices with active vaccinated
  v_index_keep_vax    <- (m_possibilities[, v_index_hh_vax] > 0)
  
  #### Community transmission matrix ####
  m_comm_trans <- matrix(0, 
                         nrow = n_hh_mod, ncol = n_hh_mod, 
                         dimnames = list(v_names_hh, v_names_hh))
  diag(m_comm_trans)[v_index_keep_sus] <- -1* # r_beta*
    m_possibilities[v_index_keep_sus, v_index_hh_sus]
  if(is.null(dim(m_comm_trans[v_index_keep_sus, v_index_keep_exp_first]))){ 
    m_comm_trans[v_index_keep_sus, v_index_keep_exp_first] <- 1* #r_beta*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]  
  } else {
    diag(m_comm_trans[v_index_keep_sus, v_index_keep_exp_first]) <- 1* #r_beta*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]  
  }
  
  
  #### Household matrices ####
  ### Household transmission matrix
  m_hh_trans <- matrix(0, 
                       nrow = n_hh_mod, ncol = n_hh_mod, 
                       dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_trans)[v_index_keep_tau] <- -1* # -r_tau_avg*
    m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
    rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop = FALSE])
  if(is.null(dim(m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first]))){ 
    # If resulting object is a scalar, don't use 'diag'
    m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first] <- 1* # r_tau_avg*
      m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
      rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop = FALSE])  
  } else{ # If resulting object is a matrix, use 'diag'
    diag(m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first]) <- 1* # r_tau_avg*
      m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
      rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop=FALSE])  
  }
  
  
  ### Household progression matrix
  a_hh_prog_out <- array(0, 
                         dim = c(n_hh_mod, n_hh_mod, n_exp_states),
                         dimnames = list(v_names_hh, v_names_hh, v_names_exp))
  m_hh_prog_out_flat <- array(0, 
                              dim = c(n_hh_mod, n_hh_mod),
                              dimnames = list(v_names_hh, v_names_hh))
  for(i in 1:n_exp_states){ # i <- 1
    diag(a_hh_prog_out[, , i])[m_index_keep_exp[, i]] <- diag(a_hh_prog_out[, , i])[m_index_keep_exp[, i]] - 
      (r_sigma*n_exp_states)*m_possibilities[m_index_keep_exp[, i], v_index_hh_exp[i]]
    m_hh_prog_out_flat <- m_hh_prog_out_flat + a_hh_prog_out[, , i]
  }  
  
  
  m_hh_prog_in_flat <- compute_inflows(n_hhsize   = n_hhsize,
                                       n_hh_mod   = n_hh_mod ,
                                       v_names_hh = v_names_hh,
                                       v_names_states  = v_names_states,
                                       v_names_source  = v_names_exp,
                                       v_names_destiny = v_names_inf,
                                       v_names_source_destiny = v_names_exp_inf, # v_names_exp_inf[-length(v_names_exp_inf)], #### CHECK,
                                       v_index_hh_source   = v_index_hh_exp,
                                       v_index_hh_destiny  = v_index_hh_inf,
                                       v_index_keep_source = v_index_keep_exp,
                                       v_index_souce_destiny_hh = v_index_exp_inf_hh,
                                       m_possibilities = m_possibilities,
                                       r_flow = r_sigma,
                                       n_source_states = n_exp_states,
                                       m_hh_out = m_hh_prog_out_flat)
  
  m_hh_prog <- m_hh_prog_out_flat + m_hh_prog_in_flat
  
  ### Household recovered matrix
  a_hh_recov_out <- array(0, 
                          dim = c(n_hh_mod, n_hh_mod, n_inf_states),
                          dimnames = list(v_names_hh, v_names_hh, v_names_inf))
  
  m_hh_recov_out_flat <- array(0, 
                               dim = c(n_hh_mod, n_hh_mod),
                               dimnames = list(v_names_hh, v_names_hh))
  for(i in 1:n_inf_states){ # i <- 1
    diag(a_hh_recov_out[, , i])[m_index_keep_inf[, i]] <- diag(a_hh_recov_out[, , i])[m_index_keep_inf[, i]] - 
      (r_gamma*n_inf_states)*m_possibilities[m_index_keep_inf[, i], v_index_hh_inf[i]]
    
    m_hh_recov_out_flat <- m_hh_recov_out_flat + a_hh_recov_out[, , i]
  }  
  
  m_hh_recov_in_flat <- compute_inflows(n_hhsize = n_hhsize,
                                        n_hh_mod = n_hh_mod ,
                                        v_names_states = v_names_states,
                                        v_names_hh = v_names_hh,
                                        v_names_source = v_names_inf,
                                        v_names_destiny = v_names_rec,
                                        v_names_source_destiny = v_names_inf_rec,
                                        v_index_hh_source = v_index_hh_inf,
                                        v_index_hh_destiny = v_index_hh_rec,
                                        v_index_keep_source = v_index_keep_inf,
                                        v_index_souce_destiny_hh = v_index_inf_rec_hh,
                                        m_possibilities = m_possibilities,
                                        r_flow = r_gamma,
                                        n_source_states = n_inf_states,
                                        m_hh_out = m_hh_recov_out_flat)
  
  m_hh_recov <- m_hh_recov_out_flat + m_hh_recov_in_flat
  
  ### Household waning matrix
  m_hh_waning <- matrix(0, 
                        nrow = n_hh_mod, ncol = n_hh_mod, 
                        dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_waning)[v_index_keep_rec] <- -r_omega*
    m_possibilities[v_index_keep_rec, v_index_hh_rec]
  if(is.null(dim(m_hh_waning[v_index_keep_rec, v_index_keep_sus]))){
    # If resulting object is a scalar, don't use 'diag'
    m_hh_waning[v_index_keep_rec, v_index_keep_sus] <- r_omega *
      m_possibilities[v_index_keep_rec, v_index_hh_rec]
  } else {
    # If resulting object is a matrix, use 'diag'
    diag(m_hh_waning[v_index_keep_rec, v_index_keep_sus]) <- r_omega *
      m_possibilities[v_index_keep_rec, v_index_hh_rec]
  }
  
  ### Household vaccine-induced immunity waning matrix
  m_hh_waning_vax <- matrix(0, 
                            nrow = n_hh_mod, ncol = n_hh_mod, 
                            dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_waning_vax)[v_index_keep_vax] <- -r_vax_omega*
    m_possibilities[v_index_keep_vax, v_index_hh_vax]
  if(is.null(dim(m_hh_waning_vax[v_index_keep_vax, v_index_keep_sus]))){
    # If resulting object is a scalar, don't use 'diag'
    m_hh_waning_vax[v_index_keep_vax, v_index_keep_sus] <- r_vax_omega*
      m_possibilities[v_index_keep_vax, v_index_hh_vax]
  } else {
    # If resulting object is a matrix, use 'diag'
    diag(m_hh_waning_vax[v_index_keep_vax, v_index_keep_sus]) <- r_vax_omega*
      m_possibilities[v_index_keep_vax, v_index_hh_vax]
  }
  
  ### Household vaccination matrix
  m_hh_vax <- matrix(0, 
                     nrow = n_hh_mod, ncol = n_hh_mod, 
                     dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_vax)[v_index_keep_sus] <- -1* # r_vax*r_vax_eff*
    m_possibilities[v_index_keep_sus, v_index_hh_sus]
  if(is.null(dim(m_hh_vax[v_index_keep_sus, v_index_keep_vax]))){ 
    # If resulting object is a scalar, don't use 'diag'
    m_hh_vax[v_index_keep_sus, v_index_keep_vax] <- 1* # r_vax*r_vax_eff*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]
  } else {
    # If resulting object is a matrix, use 'diag'
    diag(m_hh_vax[v_index_keep_sus, v_index_keep_vax]) <- 1* # r_vax*r_vax_eff*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]
  }
  
  # rowSums(m_comm_trans)
  # rowSums(m_hh_trans)
  # rowSums(m_hh_prog)
  # rowSums(m_hh_recov)
  # rowSums(m_hh_waning)
  # rowSums(m_hh_waning_vax)
  # rowSums(m_hh_vax)
  
  return(list(m_comm_trans    = m_comm_trans,
              m_hh_trans      = m_hh_trans,
              m_hh_prog       = m_hh_prog,
              m_hh_recov      = m_hh_recov,
              m_hh_waning     = m_hh_waning,
              m_hh_waning_vax = m_hh_waning_vax,
              m_hh_vax        = m_hh_vax))
}

#' Compute inflows considering competing risks using convolution of binomials

#' \code{compute_inflows} computes inflows 
#' @param n_hhsize
#' @param n_hh_mod,
#' @param v_names_hh,
#' @param v_names_states
#' @param v_names_source,
#' @param v_names_destiny
#' @param v_names_source_destiny,
#' @param v_index_hh_source,
#' @param v_index_hh_destiny,
#' @param v_index_keep_source,
#' @param v_index_souce_destiny_hh,
#' @param m_possibilities,
#' @param r_flow,
#' @param n_source_states,
#' @param m_hh_out
#' @return 
#' A matrix with inflow rates.
#' @export
compute_inflows <- function(n_hhsize,
                            n_hh_mod,
                            v_names_hh,
                            v_names_states,
                            v_names_source,
                            v_names_destiny,
                            v_names_source_destiny,
                            v_index_hh_source,
                            v_index_hh_destiny,
                            v_index_keep_source,
                            v_index_souce_destiny_hh,
                            m_possibilities,
                            r_flow,
                            n_source_states,
                            m_hh_out){
  m_hh_in  <- array(0, 
                    dim = c(n_hh_mod, n_hh_mod),
                    dimnames = list(v_names_hh, v_names_hh))
  
  v_index_keep_source_hh_destiny_first_or <- (rowSums(m_possibilities[, v_index_hh_source[-1], drop=FALSE] > 0) | 
                                                (rowSums(m_possibilities[, v_index_hh_destiny[1], drop=FALSE]) > 0))
  
  v_index_hh_signature <- which(!(v_names_states %in% v_names_source_destiny))
  
  rownames(m_possibilities) <- v_names_hh
  df_possibilities <- as.data.frame(m_possibilities)
  
  m_flow_source  <- m_possibilities[v_index_keep_source, , drop = FALSE] # ADDED
  m_flow_destiny <- m_possibilities[v_index_keep_source_hh_destiny_first_or, , drop = FALSE] # ADDED
  m_flow_destiny <- rbind(m_flow_destiny, 
                          m_flow_source[rownames(m_flow_source)[(!rownames(m_flow_source) %in% rownames(m_flow_destiny))], ]) # these rows represent the states where theres is no progression. It could have been progression but there wasn't for this particular case. i.e., the fraction of people that doesn't move
  
  v_names_HH_flow_source <-  paste("HH", 
                                   as.matrix(tidyr::unite(df_possibilities[v_index_keep_source, ], 
                                                          col = "names", sep = "")), 
                                   sep = "")
  v_names_HH_flow_destiny <-  paste("HH", 
                                    as.matrix(tidyr::unite(df_possibilities[v_index_keep_source_hh_destiny_first_or, ], 
                                                           col = "names", sep = ""))
                                    , sep = "")
  v_names_HH_flow_destiny <- c(v_names_HH_flow_destiny, 
                               v_names_HH_flow_source[!(v_names_HH_flow_source %in% v_names_HH_flow_destiny)])
  
  m_p_flow <- matrix(0, 
                     nrow = length(v_names_HH_flow_source), 
                     ncol = length(v_names_HH_flow_destiny), 
                     dimnames = list(v_names_HH_flow_source, 
                                     v_names_HH_flow_destiny))
  # View(m_p_flow)
  
  p_flow <- 1-exp(-r_flow*n_source_states*1)
  
  for(i in 1:nrow(m_flow_source)){ # i = 41
    # print(rownames(m_flow_source)[i])
    # We want to determine whether the row is binomial or a convolution of binomials
    # and how many columns we are going to need and which columns we need.
    curr_row      <- m_flow_source[i, ]
    name_curr_row <- rownames(m_flow_source)[i]
    # For all the columns that are not first element in the source through the 
    # first element in the destiny, we want to fix their values and total them 
    # up. The household size minus that total is the number of people that can 
    # move.
    signature <- curr_row[v_index_hh_signature] # c(1, (v_index_hh_destiny[1]+1))
    n_movers  <- n_hhsize - sum(signature)
    
    # We want to scan across all elements of the source if ony one of those 
    # numbers is greater than 0, then we are binomial; otherwise, we are 
    # convolution.
    binomial_flag <- FALSE
    if(n_movers == 1){
      binomial_flag <- TRUE
    } else{
      binomial_flag <- sum(curr_row[v_names_source] == n_movers) == 1
    }
    
    # If we are binomial, then the number of states that we need to find is 
    # number of elements in source + 1, which is the diagonal state. 
    if(binomial_flag){
      index_source  <- which(names(curr_row) == names(which(curr_row[v_names_source]>0)))
      index_destiny <- index_source + 1
      # v_index_not_movers <- names(curr_row)[!(names(curr_row) %in% c(names(signature), 
      #                                                                names(curr_row[index_source]), 
      #                                                                names(curr_row[index_destiny])))]
      m_destination_cols <- matrix(0, 
                                   nrow = n_movers + 1,
                                   ncol = length(curr_row), 
                                   dimnames = list(1:(n_movers+1), names(curr_row)))
      t_temp <- t(m_destination_cols)
      t_temp[names(signature), ] <- signature
      m_destination_cols         <- t(t_temp)
      for(mover in 0:n_movers){
        m_destination_cols[mover + 1, index_source]  <- n_movers - mover
        m_destination_cols[mover + 1, index_destiny] <- mover
      }
      v_names_hh_destiny <- paste("HH", 
                                  as.matrix(tidyr::unite(as.data.frame(m_destination_cols),
                                                         col = "names", sep = "")), 
                                  sep = "")
      
      m_p_flow[name_curr_row, v_names_hh_destiny] <- dbinom(0:n_movers, 
                                                            n_movers, 
                                                            prob = p_flow)
    } else{
      # Loop over `m_flow_destiny` columns, gather a list of them that we have to
      # test using the backwards algorithm.
      # Rules:
      #  1. Reject if 
      #     a. first Source column greater than curr_row first source column,
      #     b. destiny signature columns are not equal. to source signature columns
      # For the shorter list, we will use the backwards algorithm
      m_keep_destiny <- cbind(m_flow_destiny, Keep = 0)
      for(j in 1:nrow(m_keep_destiny)){ # j <- 7
        curr_row_destiny <- m_keep_destiny[j, ]
        flag1 <- curr_row_destiny[v_names_source[1]] <= curr_row[v_names_source[1]] 
        flag2 <- sum(curr_row_destiny[names(signature)] == curr_row[names(signature)]) == length(names(signature))
        m_keep_destiny[j, "Keep"] <- flag1 & flag2
      }
      m_keep_destiny_limited <- m_keep_destiny[m_keep_destiny[, "Keep"]==1, ]
      v_names_source_destiny <- c(v_names_source, v_names_destiny[1])       # Documented in SC-COSMO
      for(j in 1:nrow(m_keep_destiny_limited)){ # j <- 2
        v_curr_delta <- vector(mode = "numeric", length = (length(v_names_source_destiny)-1))
        curr_delta <- 0
        for(k in n_source_states:1){# k = 3
          curr_col_source  <- curr_row[v_names_source_destiny[k+1]]
          curr_col_destiny <- m_keep_destiny_limited[j, v_names_source_destiny[k+1]] 
          curr_delta <- curr_col_destiny + curr_delta - curr_col_source
          v_curr_delta[k] <- curr_delta
          if(curr_delta < 0){
            m_keep_destiny_limited[j, "Keep"] <- FALSE
            break
          }
          if(curr_delta > curr_row[v_names_source_destiny[k]]){
            m_keep_destiny_limited[j, "Keep"] <- FALSE
            break
          }
        }
        
        if(m_keep_destiny_limited[j, "Keep"]==TRUE){
          binom_conv <- 1
          for(k in 1:n_source_states){ # k = 2
            # print(curr_row[k])
            # print(v_curr_delta[k])
            binom_conv <- binom_conv*dbinom(x = v_curr_delta[k], 
                                            size = curr_row[v_names_source[k]], 
                                            prob = p_flow)
            # print(binom_conv)
          }
          m_p_flow[name_curr_row, rownames(m_keep_destiny_limited)[j]] <- binom_conv
        }
      }
      
    }
    # t(t(m_keep_destiny[m_keep_destiny[, "Keep"]==1, ]) - c(curr_row, 0))
  }
  
  m_r_flow <- -log(1-m_p_flow)
  if(length(v_names_HH_flow_source)==1){
    m_r_flow[v_names_HH_flow_source, v_names_HH_flow_source] <- 0
  } else{
    diag(m_r_flow[v_names_HH_flow_source, v_names_HH_flow_source]) <- 0
  } 
  
  m_r_prop_flow <- m_r_flow/rowSums(m_r_flow)
  
  if(length(v_names_HH_flow_source)==1){
    m_hh_in[v_names_HH_flow_source, v_names_HH_flow_destiny] <- -1*(m_hh_out[v_names_HH_flow_source, v_names_HH_flow_source] * m_r_prop_flow)  
  } else{
    m_hh_in[v_names_HH_flow_source, v_names_HH_flow_destiny] <- -1*(diag(m_hh_out[v_names_HH_flow_source, v_names_HH_flow_source]) * m_r_prop_flow)  
  }
  
  return(m_hh_in)
}

#' Get non-pharmaceutical intervention (NPI) for a given time
#'
#' \code{get_npi} is used to get the time-varying NPI effect
#'
#' @param n_time A time for which to retrieve the NPI level
#' @param parameters List with all parameters of model
#' @return 
#' A value of the NPI level.
#' @export
get_npi <- function(n_time, parameters){
  p_npi <- parameters$fun_npi(n_time)
  return(p_npi)
}

#' Get vaccination (Vax) rate for a given time
#'
#' \code{get_vax} is used to get the time-varying vaccination rate
#'
#' @param n_time A time for which to retrieve the Vax rate
#' @param parameters List with all parameters of model
#' @return 
#' A value of the Vax rate.
#' @export
get_vax <- function(n_time, parameters){
  r_vax <- parameters$fun_vax(n_time)
  return(r_vax)
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return 
#' A list with all parameters updated.
#' @export
update_param_list <- function(l_params_all, params_updated){
  
  if (typeof(params_updated)!="list"){
    params_updated <- split(unname(params_updated),names(params_updated)) #converte the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}

#' Undetected infections
#'
#' \code{calc_inf_nodx} Calculate number of undetected infections over time.
#'
#' @Param l_out List with output from SC-COSMO and all parameters
#' @return
#' A data.frame with the number of undetected infections as columns over time
#' @export
calc_inf_nodx <- function(l_out){
  df_out <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  df_InfNoDx <- data.frame(time = df_out$time,
                           InfNoDX = rowSums(df_out[, l_params_all$v_names_inf,
                                                    drop=FALSE]),
                           check.names = FALSE)
  return(df_InfNoDx)
}
#' Total infections
#'
#' \code{calc_inf_totals} Calculate total number of infections over time.
#'
#' @Param l_out List with output from SC-COSMO and all parameters
#' @return
#' A data.frame with the total number of infections as columns over time
#' @export
calc_inf_totals <- function(l_out){
  df_out <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  df_Inftot <- data.frame(time = df_out$time,
                          Inftot = rowSums(df_out[, c(l_params_all$v_names_inf,
                                                      l_params_all$v_names_inf_dx),
                                                  drop=FALSE]),
                          check.names = FALSE)
  return(df_Inftot)
}

#' Total exposed
#' 
#' \code{calc_exp_totals} Calculate total number of exposed over time.
#' 
#' @param l_out List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of exposed as columns over time
#' @export
calc_exp_totals <- function(l_out){
  df_out <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  
  df_Exptot <- data.frame(time = df_out$time, 
                          Exptot = rowSums(df_out[, l_params_all$v_names_exp,
                                                  drop=FALSE]), 
                          check.names = FALSE)
  return(df_Exptot)
}



#' Plot SEIRV results
show_MC_SEIRV_model_results <- function(l_out) {
  output <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  df1 <- data.frame(output)
  dfmc <- c()
  dfmc$time <- df1$time
  dfmc$S <- df1$S
  dfmc$E <- rowSums(df1[, l_params_all$v_names_exp, drop = FALSE])
  dfmc$I <- rowSums(df1[, l_params_all$v_names_inf, drop = FALSE])
  dfmc$IDX <- rowSums(df1[, l_params_all$v_names_inf_dx, drop = FALSE])
  dfmc$R <- df1$R
  dfmc$V <- df1$V
  
  dfmc$N <- dfmc$S + dfmc$E + dfmc$I + dfmc$IDX + dfmc$R + dfmc$V
  
  dfmc <- data.frame(dfmc)
  ggplot() + 
    # geom_line(data = dfmc, aes(x = time, y = S), color = "blue"  ) +
    # geom_line(data = dfmc, aes(x = time, y = E), color = "purple") +
    # geom_line(data = dfmc, aes(x = time, y = I), color = "red"   ) +
    # geom_line(data = dfmc, aes(x = time, y = R), color = "green" ) +
    
    # geom_line(data = dfmc, aes(x = time, y = S, color = "S")) +
    geom_line(data = dfmc, aes(x = time, y = E, color = "E")) +
    geom_line(data = dfmc, aes(x = time, y = I, color = "I")) +
    geom_line(data = dfmc, aes(x = time, y = IDX, color = "IDX")) +
    geom_line(data = dfmc, aes(x = time, y = R, color = "R")) +
    geom_line(data = dfmc, aes(x = time, y = V, color = "V")) +
    # scale_color_manual("Compartment",
    #                    values = c("blue", "purple", "red", "green"),
    #                    labels = c("S", "E", "I", "R")) +
    xlab('Time') +
    ylab('Count') +
    theme(legend.position = "bottom")
}

#' Plot SEIR results
show_MC_SEIR_model_results <- function(l_out) {
  output <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  df1 <- data.frame(output)
  dfmc <- c()
  dfmc$time <- df1$time
  dfmc$S <- df1$S
  dfmc$E <- rowSums(df1[, l_params_all$v_names_exp, drop = FALSE])
  dfmc$I <- rowSums(df1[, l_params_all$v_names_inf, drop = FALSE])
  dfmc$R <- df1$R
  dfmc <- data.frame(dfmc)
  ggplot() + 
    # geom_line(data = dfmc, aes(x = time, y = S), color = "blue"  ) +
    # geom_line(data = dfmc, aes(x = time, y = E), color = "purple") +
    # geom_line(data = dfmc, aes(x = time, y = I), color = "red"   ) +
    # geom_line(data = dfmc, aes(x = time, y = R), color = "green" ) +
    geom_line(data = dfmc, aes(x = time, y = S, color = "S")) +
    geom_line(data = dfmc, aes(x = time, y = E, color = "E")) +
    geom_line(data = dfmc, aes(x = time, y = I, color = "I")) +
    geom_line(data = dfmc, aes(x = time, y = R, color = "R")) +
    # scale_color_manual("Compartment",
    #                    values = c("blue", "purple", "red", "green"),
    #                    labels = c("S", "E", "I", "R")) +
    xlab('Time') +
    ylab('Count') +
    theme(legend.position = "bottom")
}

show_MC_EI_model_results <- function(l_out) {
  output <- l_out$df_out_hh_mc_seir
  l_params_all <- l_out$l_params_all
  df1 <- data.frame(output)
  dfmc <- c()
  dfmc$time <- df1$time
  dfmc$S <- df1$S
  dfmc$E   <- rowSums(df1[, l_params_all$v_names_exp, drop = FALSE])
  dfmc$I   <- rowSums(df1[, l_params_all$v_names_inf, drop = FALSE])
  dfmc$IDX <- rowSums(df1[, l_params_all$v_names_inf_dx, drop = FALSE])
  dfmc$R <- df1$R
  dfmc <- data.frame(dfmc)
  ggplot() + 
    # # geom_line(data = dfmc, aes(x = time, y = S), color = "blue") +
    # geom_line(data = dfmc, aes(x = time, y = E), color = "purple") +
    # geom_line(data = dfmc, aes(x = time, y = I), color = "red") +
    # # geom_line(data = dfmc, aes(x = time, y = R), color = "green") +
    # geom_line(data = dfmc, aes(x = time, y = S, color = "S")) +
    geom_line(data = dfmc, aes(x = time, y = E, color = "E")) +
    geom_line(data = dfmc, aes(x = time, y = I, color = "I")) +
    geom_line(data = dfmc, aes(x = time, y = IDX, color = "IDX")) +
    # geom_line(data = dfmc, aes(x = time, y = R, color = "R")) +
    scale_color_manual("Compartment",
                       values = c("purple", "red", "gray"),
                       labels = c("E", "I", "IDX")) +
    xlab('Time') +
    ylab('Count') +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

show_MC_EI_model_results_old <- function(output) {
  # output <- l_out$df_out_hh_mc_seir
  # l_params_all <- l_out$l_params_all
  df1 <- data.frame(output)
  dfmc <- c()
  dfmc$time <- df1$time
  dfmc$S <- df1$S
  dfmc$E <- rowSums(df1[, v_names_exp])
  dfmc$I <- rowSums(df1[, v_names_inf])
  dfmc$R <- df1$R
  dfmc <- data.frame(dfmc)
  ggplot() + 
    # # geom_line(data = dfmc, aes(x = time, y = S), color = "blue") +
    # geom_line(data = dfmc, aes(x = time, y = E), color = "purple") +
    # geom_line(data = dfmc, aes(x = time, y = I), color = "red") +
    # # geom_line(data = dfmc, aes(x = time, y = R), color = "green") +
    # geom_line(data = dfmc, aes(x = time, y = S, color = "S")) +
    geom_line(data = dfmc, aes(x = time, y = E, color = "E")) +
    geom_line(data = dfmc, aes(x = time, y = I, color = "I")) +
    # geom_line(data = dfmc, aes(x = time, y = R, color = "R")) +
    scale_color_manual("Compartment",
                       values = c("purple", "red"),
                       labels = c("E", "I")) +
    xlab('Time') +
    ylab('Count') +
    theme(legend.position = "bottom")
}

#' Daily rate
#'
#' \code{daily_rate} generates the daily rate to achieve a cumulative proportion 
#' over a given duration.
#'
#' @param cum_prop Cumulative proportion
#' @param duration total number of days to achieve proportion
#' @return 
#' A daily rate
#' @export
daily_rate <- function(cum_prop, duration){
  if((cum_prop >= 1) | (cum_prop < 0)){
    stop("cum_prop should be [0,1)")
  }
  return(-log(1 - cum_prop)/duration)
}
