## linear mixed model with repeated measures
metab.mixed.model <-  function(dtst,columns_with_metabolite_data, covariables, dataset_name){

  ## set up results table
  mixed_results = as.data.frame(matrix(data = NA, nrow = length(columns_with_metabolite_data), ncol = 14))
  ## set up reporting var
  s = base::seq(from=1,to=length(columns_with_metabolite_data), by=1)
  for (i in 1:length(columns_with_metabolite_data)){
    
    
    mtb_name <- names(dtst[columns_with_metabolite_data])[i]
    
    rownames(mixed_results)[rownames(mixed_results) == i] <- mtb_name
    if  (is.na(covariables[1])) {
      mtb <- dtst[!is.na(dtst[columns_with_metabolite_data][i]),c(mtb_name,"studyid","Timepoint")]
    } else { mtb <- dtst[!is.na(dtst[columns_with_metabolite_data][i]),c(mtb_name,"studyid","Timepoint", covariables)]
    }
    names(mtb)[1] <- "feature"
    mixed_results[i,1] <- nrow(mtb)
    
    # fit models
    if( is.na( covariables[1] ) ){
      formA = formula(feature ~ (1|studyid))
      formB = formula(feature ~ Timepoint + (1|studyid))
    } else {
      formA = formula(paste0("feature ~ ", paste0(covariables, collapse = " + ") ," + ", "(1|studyid)"))
      formB = formula(paste0("feature ~ ", paste0(covariables, collapse = " + ") ," + ", "Timepoint + (1|studyid)"))  
    }
    
    fitA <- tryCatch(lmerTest::lmer(formA, data = mtb, REML=TRUE), error=function(e) e)
    fitB <- tryCatch(lmerTest::lmer(formB, data = mtb, REML=TRUE), error=function(e) e)
    
    
    # fill mixed_results table
    coef = coef(summary(fitB))
    tp_coef = coef[nrow(coef),]
    mixed_results[i,2:6] <- tp_coef
    mixed_results[i,7] <- toString(fitA@optinfo$conv$lme4$messages)
    mixed_results[i,8] <- toString(fitA@optinfo$conv$lme4$code)
    mixed_results[i,9] <- toString(fitB@optinfo$conv$lme4$messages)
    mixed_results[i,10] <- toString(fitB@optinfo$conv$lme4$code)

  }
  # set column names
  names(mixed_results)[1] <- "n_mixed"
  names(mixed_results)[2] <- "beta_mixed"
  names(mixed_results)[3] <- "se_mixed"
  names(mixed_results)[4] <- "df_mixed"
  names(mixed_results)[5] <- "t_mixed"
  names(mixed_results)[6] <- "P_mixed"
  names(mixed_results)[7] <- "mixed_model_warnings_wo_tp"
  names(mixed_results)[8] <- "mixed_model_conv_wo_tp"
  names(mixed_results)[9] <- "mixed_model_warnings_w_tp"
  names(mixed_results)[10] <- "mixed_model_conv_w_tp"

  
  return(mixed_results)
}






