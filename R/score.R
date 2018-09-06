#' Score CREDI response data
#'
#' This function scores response data from administering the CREDI instruments. MAP (default), EAP, or plausible value estimates can be requested.
#' @param data (data.frame) Defaults to NULL. Response data. If NULL, then user is prompted to identify a .csv file with response data. Defaults to NULL.
#' @param reverse_code (Logical) Defaults to TRUE. If TRUE, then reverse coding is automated to appropriately handle the negatively worded items LF9, LF102, LFMH1, LFMH2, LFMH3, LFMH4, LFMH5, LFMH7, LFMH8, & LFMH9. If FALSE, then no reverse coding is applied.
#' @param interactive (Logical) Defaults to TRUE. If TRUE, the user may be prompted with caution messages regarding whether scoring should be continued, where to save the scores, where to save a logfile, etc. If FALSE, continuation is assumed and scores and the user is not prompted to save scores or a logfile.
#' @param EAP (logical) Defaults to FALSE. If TRUE, then EAP estimates are used in place of MAP estimates.
#' @param M (integer) Defaults to 0. Generates M plausible value for each observations.
#' @param MH_iterations (integer) Defaults to 2000. Number of total iterations for the Metropolis-Hastings MCMC sampler if EAP estimates or plausible values are requested. Half of chain is used for burnin and is discarded.
#' @param Chains (integer) Number of MCMC chains used when constructing plausible values. Default is 2 chains.
#' @param Processors (integer) Number of processors for MCMC sampling (usually equal to the number of chains). Default is 2 processors.
#' @param sink_wd (character) Usually not required and defaults to NULL. If specified, a "sink" .txt file is written in the sink_wd directory with screen outputs and messages. Useful for debugging when multiple processors are used for the MCMC sampling.
#' @keywords CREDI
#' @export
#' @examples
#' score()

#
#reverse_code = FALSE
#save_logfile = FALSE
#interactive = FALSE
#data = input_df
#M  = 0
# Chains = 2
# MH_iterations = 2000

score<-function(data = NULL, reverse_code = TRUE, interactive = TRUE, EAP = FALSE, M = 0, MH_iterations = 2000, Chains = 2, Processors = 2, sink_wd = NULL){

  
    # Identify if dialog specifying .csv file should be bypassed.
    bypass = ifelse(is.null(data), FALSE, TRUE)
    if (bypass == TRUE){
      if (!is.data.frame(data)){
        stop("data argument must be type data.frame")
      }
    }

    # Load required pacakges

    require("stats")
    require("svDialogs")
    require("Rcpp")
    require("RcppArmadillo")
    require("numDeriv")
    if (M>0 | EAP == TRUE){
      require("coda")
      require("MHadaptive")
      require("foreach")
      require("doParallel")
      require("doRNG")
      require("mitools")
    }


    # Created log file
    time1 = proc.time()
    log = list(c("------------------------------------"), c("Log for CREDI Scoring Messages"),
            paste("Date:", Sys.time()), c("------------------------------------"))


    # Load in the response data, if not bypassed
    csv_wd = getwd()
    if(bypass == FALSE){
      out_dlgOpen = dlgOpen(title = "Select the .csv file with the CREDI response data",
                            filters = c("csv", "*.csv"))
      csv_file = out_dlgOpen$res
      if (!endsWith(tolower(csv_file), ".csv")){stop("Selected file is not a .csv file.", call. = FALSE)}
      csv_wd = paste(strsplit(csv_file,"/")[[1]][-length(strsplit(csv_file,"/")[[1]])],collapse = "/")
      setwd(csv_wd)
      input_df = read.csv(file = csv_file)
    } else {
      input_df = data
    }

    # Clean the input data
    list_cleaned = clean(input_df = input_df, mest_df = mest_df, reverse_code = reverse_code,
                         interactive = interactive, log = log)
    log = list_cleaned$log
    if(list_cleaned$stop!=0){

      print("*Error: Processing the provided response data resulted in errors. See log for more details.")

      if (interactive == TRUE){
        x<-as.character(readline(prompt = "Would you like to save a log file of warning and error messages? [Y/N]: "))
        x <- toupper(x)
        cut = 0
        while(cut == 0){
          if (x == "Y"){
            cut = 1;
          } else if (x == "N"){
            cut = 1
            stop("Scoring canceled.", call. = FALSE)
          } else {
            x<-as.character(readline(prompt = "Would you like to continue? [Y/N]:"))
            x <- toupper(x)
            cut = 0
          }
        } #end while
        write_log(log = log, folder = csv_wd)
      } #End if interactive

      return(list(log = log))

    } # End if stop != 0
    cleaned_df = list_cleaned$cleaned_df
    items_noresponse = list_cleaned$items_noresponse

    # Crate data matricies
    X = model.matrix(~1 + I( (AGE-18)/10.39 ) + I( ((AGE-18)/10.39)^2 ) + I( ((AGE-18)/10.39)^3 ), data = cleaned_df)
    X_4 = model.matrix(~1 + I( (AGE-18)/10.39 ) + I( ((AGE-18)/10.39)^2 ) + I( ((AGE-18)/10.39)^3 ) + I( ((AGE-18)/10.39)^4 ), data = cleaned_df)
    Y = as.matrix(cleaned_df[,-match(c("ID","AGE",items_noresponse), names(cleaned_df))]); Y[is.na(Y)] = -9L
    MU_LF = X%*%as.matrix(B) #NxK (matrix)
    MU_SF = X%*%as.numeric(beta) #Nx1

    # Obtain necessary parameter matricies
    inds_exclude = match(items_noresponse, mest_df$CREDI_code)
    if (length(inds_exclude)==0){
      LAMBDA = as.matrix(mest_df[,c("MOT","COG","LANG","SEM")])
      TAU = as.vector(mest_df$tau)
      ALPHA = as.vector(mest_df$alpha)
      DELTA = as.vector(mest_df$delta)
    }else{
      LAMBDA = as.matrix(mest_df[-inds_exclude,c("MOT","COG","LANG","SEM")])
      TAU = as.vector(mest_df$tau[-inds_exclude])
      ALPHA = as.vector(mest_df$alpha[-inds_exclude])
      DELTA = as.vector(mest_df$delta[-inds_exclude])
    }


    # Obtain necessary constants
    J = ncol(Y);
    K = 4L
    P = 3L
    N = as.integer(nrow(Y))
    invS = as.matrix(invS)
    SIGMA_SQ= exp(X%*%as.numeric(gamma))


    # initialize the theta values
    THETA0_LF = MU_LF #NxK (matrix)
    THETA0_SF = MU_SF #Nx1 (matrix)

    # Conduct the optimization and plausible values
    if (Processors>1){
      cl<-makePSOCKcluster(Processors)
      registerDoParallel(cl)
    }
    if (M>0){
      pv_array = array(NA, dim = c(N,K,M))
    }
    MAP_LF = 0.*THETA0_LF + NA
    MAP_SF = 0.*THETA0_SF + NA
    if (EAP==TRUE){
      EAP_LF = MAP_LF
      EAP_SF = MAP_SF
    }
    Z_LF = MAP_LF
    Z_SF = MAP_SF
    SE_LF = MAP_LF
    SE_SF = MAP_SF
    writeLines(paste("\nScoring ", N, " observations:"))
    pb<-txtProgressBar(min = 0, max = N, initial = 0, style = 3)
    for (i in 1:N){

      # Obtain the standardized estimates
      center_i = X_4[i,] %*% as.matrix(normcoef_mean)
      scale_i =  X_4[i,] %*% as.matrix(normcoef_sd)

      # Score the long form
      out_LF = optim(par = as.vector(THETA0_LF[i,]),
                     fn = cpp_posterior_density,
                     #gr = cpp_grad_posterior_density,
                     Yi = as.vector(Y[i,]),
                     MUi = as.vector(MU_LF[i,]),
                     invS =invS,
                     TAU = TAU,
                     LAMBDA = LAMBDA,
                     J = J,
                     K = K,
                     method = "BFGS",
                     hessian = TRUE)
      if(out_LF$convergence == 0){
            MAP_LF[i,] = out_LF$par
            fisherInfo = out_LF$hessian
            SE_LF[i,] = sqrt(diag(solve(fisherInfo,diag(K))))
    
            Z_LF[i,1:4] = (MAP_LF[i,]+50-center_i[1,1:4])/scale_i[1,1:4]
    
            # Average of the scores (note that the _SF is misleading b/c it is not on the same scale as short form)
            MAP_SF[i, 1] = mean(out_LF$par)
            SE_SF[i,1] = (1/1)*(1/sqrt(sum(sum(out_LF$hessian))))
            Z_SF[i,1] = weighted.mean(Z_LF[i,], w = (1./SE_LF[i,])^2)
            
            if (EAP == TRUE | M>0){
              
              if (Processors == 1){
                
                sample_i = list(NULL)
                
                for (mcmc in 1:Chains){
    
                  obj_MH = Metro_Hastings(
                    li_func = function(THETAi, Yi, MUi, invS, TAU, LAMBDA, J, K) {
                      return(-1.0 * cpp_posterior_density(THETAi, Yi, MUi, invS, TAU, LAMBDA, J, K))
                    },
                    pars = out_LF$par + runif(length(out_LF$par), -1.5, 1.5),
                    prop_sigma = solve(out_LF$hessian),
                    Yi = as.vector(Y[i, ]),
                    MUi = as.vector(MU_LF[i, ]),
                    invS = invS,
                    TAU = TAU,
                    LAMBDA = LAMBDA,
                    J = J,
                    K = K,
                    iterations = MH_iterations,
                    burn_in = floor(0.5 *MH_iterations),
                    quiet = TRUE
                  )
                  chain_i = as.mcmc(obj_MH$trace)
                  sample_i[[mcmc]] = chain_i
                  
                } # End mcmc in 1:Chains
    
                
                
              } # end if Processors == 1
              
              if (Processors>1){
                  sample_i <- foreach(
                    mcmc = 1:Chains,
                    .packages = c("credicpp","coda","doRNG", "MHadaptive")
                  ) %dorng% {
        
                    if(i==1 & !is.null(sink_wd)){
                      sink(file = paste0(sink_wd,"/chain",mcmc,".txt"))
                    }
        
                    print(paste0("\n --------- ID = ", cleaned_df$ID[i], ", CHAIN = ", mcmc, " ----------"))
                    
                    obj_MH = Metro_Hastings(
                      li_func = function(THETAi, Yi, MUi, invS, TAU, LAMBDA, J, K) {
                        return(-1.0 * cpp_posterior_density(THETAi, Yi, MUi, invS, TAU, LAMBDA, J, K))
                      },
                      pars = out_LF$par,
                      prop_sigma = solve(out_LF$hessian),
                      Yi = as.vector(Y[i, ]),
                      MUi = as.vector(MU_LF[i, ]),
                      invS = invS,
                      TAU = TAU,
                      LAMBDA = LAMBDA,
                      J = J,
                      K = K,
                      iterations = MH_iterations,
                      burn_in = floor(0.5 *MH_iterations),
                      quiet = FALSE
                    )
                    chain_i = obj_MH$trace
                    return(as.mcmc(chain_i))
                    
                    if(!is.null(sink_wd)){sink()}
                  } #End foreach
              } # End if processors
              
              mcmc_i = mcmc.list(sample_i)
              gelman_i= gelman.diag(mcmc_i, autoburnin = FALSE)
              es_i = effectiveSize(mcmc_i)
              converge_i = ifelse(gelman_i$mpsrf<1.1 & min(es_i)>M, TRUE, FALSE)
              if (converge_i == FALSE){
                warning_mcmc = paste0("MCMC chains did not converge when sampling plausible values for ID = ", i, ". Increase the number of iterations of the MH sampler and the number of chains.")
                warning(warning_mcmc)
                log[length(log)+1] = warning_mcmc
                }
              if (M>0){
                matrix_mcmc_i = as.matrix(mcmc_i)
                pv_i = matrix_mcmc_i[round(seq(1,nrow(matrix_mcmc_i), len =  M)), ]
                pv_array[i, , ] = pv_i + 50
              }
              if (EAP == TRUE){
                EAP_LF[i, ]=summary(mcmc_i)$statistics[,1] + 50
                EAP_SF[i, ] = mean(EAP_LF[i,])
              }
              #if(!is.null(sink_wd)){sink()}
              #setwd(wd_tmp)
            } #End if M>0
        

      }#end if out_LF$convergence==0


      setTxtProgressBar(pb, i)
    }# end for i = 1:N

    if (M>0 & Processors>1){stopImplicitCluster()}

    # Clean up the MAP_LF and SE_LF
    MAP_LF = data.frame(round(MAP_LF,3)+50)
    SE_LF = data.frame(round(SE_LF,3)); names(SE_LF) = paste(names(SE_LF),"_SE", sep = "")

    # Clean up the MAP_SF and SE_SF
    MAP_SF = data.frame(OVERALL = round(MAP_SF,3)+50)
    SE_SF = data.frame(OVERALL_SE = round(SE_SF,3))

    #Clean the standardized estimates
    Z_LF = data.frame(round(Z_LF,3))
    names(Z_LF) = paste("z_",names(Z_LF), sep = "")

    Z_SF = data.frame(round(Z_SF, 3))
    names(Z_SF) = "z_OVERALL"

    # Put in the input
    output_df = cbind(data.frame(ID = cleaned_df$ID), Z_LF, Z_SF, MAP_LF, MAP_SF,SE_LF, SE_SF)

    # Write out the data
    if(interactive == TRUE){
        out_dlgDir = dlgSave(default = csv_wd, title = "Save scores as", gui = .GUI)
        out_csv = paste(strsplit(out_dlgDir$res,"/")[[1]],collapse = "/")

        if (!endsWith(out_csv,".csv")){out_csv = paste(out_csv, ".csv", sep = "")}
        write.csv(output_df, file = out_csv, row.names = FALSE)

        log[length(log)+1] = paste("\n Scores written to ", out_csv,".", sep = "")

        txt_wd = paste(strsplit(out_csv,"/")[[1]][-length(strsplit(out_csv,"/")[[1]])],collapse = "/")
        out_txt = paste(txt_wd,"/logfile - CREDI scoring.txt", sep = "")
        write_log(log = log, folder = txt_wd, file = out_txt)

        writeLines("\n")
        writeLines( paste("\n Scores written to ", out_csv,".", sep = "") )


    }

    #Return plausible values
    return_list = list(scores = output_df, log = log)
    if (EAP == TRUE){
      EAP_df = EAP_LF
      names(EAP_df) = c("MOT","COG","LANG","SEM")
      return_list$EAP = EAP_df
    }
    if (M>0){
      imputation_list = lapply(1:M, function(m){df_m = data.frame(pv_array[,,m]); names(df_m) = c("MOT", "COG","LANG","SEM"); df_m = transform(df_m, ID = cleaned_df$ID)})
      imputation_list = imputationList(imputation_list)
      return_list$pvs_list = imputation_list
    }
    return(return_list)
} #End funcion score
