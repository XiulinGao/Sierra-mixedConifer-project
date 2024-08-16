#### FATES ensemble parameter sets generator ####
# This script adapted approach from Marcos Longo 
# and the Ensemble_FATES_Param.r script Marcos   
# originally wrote. Although, we now also account 
# for different PFT strategies by keeping difference in specific parameters among PFTs
# Author: Xiu Lin Gao 
# Date: 2024-07-07 


####   Reset R before running the script ####
#---~---
# Unload all packages
suppressWarnings({
  plist = names(sessionInfo()$otherPkgs)
  if (length(plist) > 0){
    dummy = sapply(X=paste0("package:",plist),FUN=detach,character.only=TRUE,unload=TRUE)
  }#end if (length(plist) > 0)
  # Remove all variables, reset warnings, close plots and clean up memory
  rm(list=ls())
  options(warn=0)
  invisible(graphics.off())
  invisible(gc())
})#end suppressWarnings


#### Path setting ####
#---~---
home_path      = path.expand("~")
work_path      = file.path(home_path,"Documents","CA-grassland-simuDoc")
parset_path    = file.path(work_path,"ensemble-param-base")
parset_base    = "tunnes-params.csv"
parcorr_base   = NA
ncdf4_in_path  = file.path(home_path,"Documents","CA-grassland-simuDoc","ensemble-param-base")
ncdf4_in_base  = c("p61-tharps-refineN2T14.nc") #for parallel ensemble
ncdf4_out_path = file.path(work_path,"EnsembleParamSet","sierra-test")
ncdf4_out_pref = substr(ncdf4_in_base[1],1,3)

#### Ensemble setting ####
n_pft=4
n_ensemble     = 180 
tasks_per_node = 36  
seed_init      = 6
lhs_eps        = 0.025
lhs_maxIt      = 1000L
verbose = FALSE

#### Package loading and setting checks ####
cat (" + Load packages.\n")
fine.packages = c( ncdf4     = require(ncdf4    ,quietly=TRUE,warn.conflicts=FALSE)
                 , MASS      = require(MASS     ,quietly=TRUE,warn.conflicts=FALSE)
                 , pse       = require(pse      ,quietly=TRUE,warn.conflicts=FALSE)
                 , stringr   = require(stringr  ,quietly=TRUE,warn.conflicts=FALSE )
)#end fine.packages
if (! all(fine.packages)){
  cat (" List of required packages, and the success status loading them:\n")
  print(fine.packages)
  stop(" Some packages are missing and must be installed.")
}#end if (! all(fine.packages))

if (! is.na(seed_init)) dummy = set.seed(seed_init)

parset_file   = file.path(parset_path  ,parset_base  )
ncdf4_in_file = file.path(ncdf4_in_path,ncdf4_in_base[1])
if ( ! file.exists(ncdf4_in_file) ){
  cat ("------------------------------------------------------------------\n")
  cat (" Reference NetCDF file not found!\n"                                 )
  cat (" - NCDF4_IN_PATH = ",ncdf4_in_path,".\n",sep=""                      )
  cat (" - NCDF4_IN_BASE = ",ncdf4_in_base,".\n",sep=""                      )
  cat ("------------------------------------------------------------------\n")
  stop(" This script requires a valid input parameter file in NetCDF format.")
}#end if ( ! file.exists(ncdf4_in_file) )

dummy = dir.create(ncdf4_out_path,recursive=TRUE,showWarnings=FALSE)

n_nodes      = ceiling(n_ensemble/tasks_per_node)
node_digits  = 1L + round(log10(n_nodes))
task_digits  = 1L + round(log10(tasks_per_node))
ens_digits   = 1L + round(log10(n_ensemble))
ens_fmt      = paste0("%",ens_digits,".",ens_digits,"i")
task_fmt     = paste0(  "Node%",node_digits,".", node_digits,"i"
                        ,"_Task%",task_digits,".",task_digits,"i"
)

#### LHS sampling ####
if (file.exists(parset_file)){
  cat (" + Read parameter settings.\n")
  param_config        = read.csv(file=parset_file,header=TRUE,stringsAsFactors=FALSE)
  
  #---~---
  #   Add a few useful columns.
  #---~---
  param_config$pft_var = is.finite(param_config$pft)
  param_config$all_pft = param_config$pft_var & (param_config$pft == 0L)
  param_config$org_var = is.finite(param_config$organ)
  param_config$all_org = param_config$org_var & (param_config$organ == 0L)
  param_config$lhs_spl = param_config$ineq == "FALSE"

  #---~---
  
  #---~---
  #   Name for correlation matrix
  #---~---
  pft_suffix             = ifelse( test = (! param_config$pft_var) | param_config$all_pft 
                                   , yes  = ""
                                   , no   = sprintf("_P%2.2i",param_config$pft)
  )#end ifelse
  organ_suffix           = ifelse( test = (! param_config$org_var) | param_config$all_org
                                   , yes  = ""
                                   , no   = sprintf("_O%2.2i",param_config$organ)
  )#end ifelse
  param_config$corr_name = paste0(param_config$fates_parameter_name,pft_suffix,organ_suffix)
  #---~---
  param_lhs = subset(param_config,lhs_spl)
  
  
  n_param_tol             = nrow(param_config)
  n_param_lhs             = nrow(param_lhs)
  
}else{
  cat ("-----------------------------------------------------------------\n"       )
  cat ("   Parameter configuration file is missing! Check configuration.\n"        )
  cat (" - parset_path = \"",parset_path,"\"\n"                             ,sep="")
  cat (" - parset_base = \"",parset_base,"\"\n"                             ,sep="")
  cat ("-----------------------------------------------------------------\n"       )
  stop(" Parameter instruction file is missing.")
}#end if (file.exists(parset_file))


parcorr_file = ifelse( test = is.na(parcorr_base)
                       , yes  = NA_character_
                       , no   = file.path(parset_path,parcorr_base)
)#end ifelse
if (is.na(parcorr_file)){
  #---~---
  #   Use sample as is, ignoring correlations.
  #---~---
  cat (" + No parameter correlation provided. Assume uncorrelated samples.\n"        )
  param_corr_lhs           = diag(nrow = n_param_lhs)
  dimnames(param_corr_lhs) = list(param_lhs$correlation,param_lhs$correlation)
  #---~---
}else if (file.exists(parcorr_file)){
  #---~---
  #   Read in the correlation information.
  #---~---
  corr_config_lhs = read.csv(file=parcorr_file,header=TRUE,stringsAsFactors=FALSE)
  n_corr_lhs      = nrow(corr_config_lhs)
  #---~---
  
  #---~---
  #   Add a few useful columns.
  #---~---
  corr_config_lhs$pft_var = is.finite(corr_config_lhs$pft)
  corr_config_lhs$all_pft = corr_config_lhs$pft_var & (corr_config_lhs$pft == 0L)
  corr_config_lhs$org_var = is.finite(corr_config_lhs$organ)
  corr_config_lhs$all_org = corr_config_lhs$org_var & (corr_config_lhs$organ == 0L)
  #---~---
  
  
  #---~---
  #   Build the correlation matrix.  Because the same parameter may appear multiple 
  # times in the parameter configuration (due to multiple PFTs or multiple organs),
  # we go through each line of the correlation configuration.
  #---~---
  param_corr_lhs           = matrix(data=0.,nrow=n_param_lhs,ncol=n_param_lhs)
  dimnames(param_corr_lhs) = list(param_lhs$corr_name,param_lhs$corr_name)
  for (m in sequence(n_corr_lhs)){
    #--- Retrieve parameters
    a_parameter = corr_config_lhs$parameter_a[m]
    b_parameter = corr_config_lhs$parameter_b[m]
    ab_pft      = corr_config_lhs$pft        [m]
    ab_organ    = corr_config_lhs$organ      [m]
    ab_corr     = corr_config_lhs$corr       [m]
    ab_pft_var  = corr_config_lhs$pft_var    [m]
    ab_all_pft  = corr_config_lhs$all_pft    [m]
    ab_org_var  = corr_config_lhs$org_var    [m]
    ab_all_org  = corr_config_lhs$all_org    [m]
    #---~---
    
    
    #---~---
    #   Retrieve rows and columns
    #---~---
    rsel = param_lhs$fates_parameter_name %in% a_parameter
    csel = param_lhs$fates_parameter_name %in% b_parameter
    psel = param_lhs$pft       %in% ab_pft
    osel = param_lhs$organ     %in% ab_organ
    #--- Restrict correlation for PFTs
    if (ab_pft_var){
      rsel = rsel & ( psel | (! param_lhs$pft_var ) | ab_all_pft )
      csel = csel & ( psel | (! param_lhs$pft_var ) | ab_all_pft )
    }#end if
    #--- Restrict correlation for organs
    if (ab_org_var){
      rsel = rsel & ( osel | (! param_lhs$org_var ) | ab_all_org )
      csel = csel & ( osel | (! param_lhs$org_var ) | ab_all_org )
    }#end if
    #---~---
    
    
    
    #---~---
    #   Update correlation for rows and columns, taking the correlation matrix symmetry
    # into account.
    #---~---
    param_corr_lhs[rsel,csel] = ab_corr
    param_corr_lhs[csel,rsel] = ab_corr
    #---~---
    
  }#end for (i in sequence(n_corr))
  #--- Update the diagonal, which should always be 1.
  diag(param_corr_lhs) = 1.
  #---~---
}else{
  cat ("-------------------------------------------------------------------\n"       )
  cat ("   Parameter correlation file is missing! In case you do not\n"              )
  cat (" want to provide correlation, set \"parcorr_base = NA_character_\".\n"       )
  cat ("\n"                                                                          )
  cat ("   Current settings:\n"                                                      )
  cat (" - parset_path  = \"",parset_path,"\"\n"                              ,sep="")
  cat (" - parcorr_base = \"",parcorr_base,"\"\n"                             ,sep="")
  cat ("-------------------------------------------------------------------\n"       )
  stop(" Parameter correlation file is missing.")
}#end if
#---~---

cat (" + Build normalised samples for parameters using LHS sample.\n")
if (is.na(parcorr_file)){
  #---~---
  #   Run a default Latin Hypercube sample of uncorrelated variables.
  #---~---
  param_sampleLHS = LHS( factors = param_lhs$corr_name
                       , N       = pmax(n_ensemble,n_param_lhs+2L)
                       , q       = rep(x="qunif",times=n_param_lhs)
                       , method  = "random"
  )#end LHS
  #---~---
}else{
  #---~---
  #   Run a default Latin Hypercube sample whilst accounting for correlation.
  #---~---
  param_sampleLHS = LHS( factors = param_lhs$corr_name
                       , N       = pmax(n_ensemble,n_param_lhs+2L)
                       , q       = rep(x="qunif",times=n_param_lhs)
                       , method  = "HL"
                       , opts    = list( COR   = param_corr_lhs
                                        , eps   = lhs_eps
                                        , maxIt = lhs_maxIt
                      )#end list
  )#end LHS
  #---~---
}#end if (sample_method %in% "default")

cat (" + Build normalised samples for parameters using LHS sampling.\n")
n_sample    = nrow(param_sampleLHS$data)
idx_use     = sort(sample(x=n_sample,size=n_ensemble,replace=FALSE))
param_quant = param_sampleLHS$data[idx_use,,drop=FALSE]
#---~---



#---~---
#   Scale quantiles to parameter units
#---~---
cat (" + Scale quantiles to the parameter range.\n")
add0        = matrix( data     = param_lhs$value_min
                    , nrow     = n_ensemble
                    , ncol     = n_param_lhs
                    , byrow    = TRUE
                    , dimnames = list(NULL,param_lhs$corr_name)
)#end matrix
mult        = matrix( data     = param_lhs$value_max-param_lhs$value_min
                    , nrow     = n_ensemble
                    , ncol     = n_param_lhs
                    , byrow    = TRUE
                    , dimnames = list(NULL,param_lhs$corr_name)
)#end matrix
add0        = as.data.frame(add0)
mult        = as.data.frame(mult)
param_table = add0 + mult * param_quant

#### non-LHS sampling ####
#---~---
#   Generate parameters that are based on reference values from other PFTs and parameters.
#---~---
param_ref = param_config$corr_name[!param_config$lhs_spl]
n_ref     = length(param_ref)
param_pd  = matrix( data = NA_real_
                  , nrow = n_ensemble
                  , ncol = n_param_tol - n_param_lhs
                  , byrow = TRUE
                  , dimnames = list(NULL,param_ref))
param_pd = as.data.frame(param_pd)
param_table = cbind(param_table,param_pd)
rm(param_pd)
opt       = c("lt","gt","eq")
for (r in sequence(nrow(param_table))){
  for(p in sequence(n_ref)){
    targ_col  = param_ref[p]
    ridx      = which(param_config$corr_name==targ_col)
    param_prefix = param_config$fates_parameter_name[ridx]
    tag_min   = param_config$value_min[ridx]
    tag_max   = param_config$value_max[ridx]
    ineq      = param_config$ineq[ridx]
    ref_pft   = param_config$ref_pft[ridx]
    ref_param = param_config$ref_parameter[ridx]
    if(ineq %in% opt){
      ref_col = paste0(ref_param,"_",ref_pft)
      ref_val = param_table[[ref_col]][[r]]
      if(tag_min > ref_val | tag_max < ref_val){
        cat ("-------------------------------------------------------------------\n"       )
        cat ("   Make sure lower and upper bounds are same for\"",targ_col,"\"","and\"",ref_col,"\"\n")
        stop("range error for sampling")
      }
      else{
        if(ineq == "gt"){
          assign_val = runif(1,min=ref_val,max=tag_max)
        }
        if(ineq == "lt"){
          assign_val = runif(1,min=tag_min,max=ref_val)
        }
        if(ineq == "eq"){
          assign_val = ref_val
        }
        param_table[[targ_col]][[r]] = assign_val
      } 
    } #end of  if(fst_cha %in% opt)
    if (ineq == "complex_sla"){
     if(ref_pft=="FALSE"){
       opt_plx = substr(ref_param, 1,2)
       ref_col = substr(ref_param,3,str_length(ref_param))
       ref_val = param_table[[ref_col]][[r]]
       if(opt_plx == "gt"){
         assign_val = runif(1,min=ref_val,max=tag_max)
       }
       if(opt_plx == "lt"){
         assign_val = runif(1,min=tag_min,max=ref_val)
       }
       if(opt_plx == "eq"){
         assign_val = ref_val
       }
       param_table[[targ_col]][[r]] = assign_val
     } #end of if(isFALSE(ref_pft))
     else{
       opt_pft     = substr(ref_pft,1,2)
       ref_pft_col = substr(ref_pft,3,str_length(ref_pft))
       ref_pft_val = param_table[[ref_pft_col]][[r]]
       opt_param   = substr(ref_param,1,2)
       ref_param_col = substr(ref_param,3,str_length(ref_param))
       ref_param_val = param_table[[ref_param_col]][[r]]
       if(is.na(ref_pft_val) | is.na(ref_param_val)){
         stop("make sure sampling of reference value happen ahead of\"",targ_col,"\"")}
       if(opt_pft=="gt" & ref_pft_val < ref_param_val){
         assign_val = runif(1,min=ref_pft_val,max=ref_param_val)}
       else if(opt_pft=="lt" & ref_pft_val > ref_param_val){
         assign_val = runif(1,min=ref_param_val,max=ref_pft_val)}
       else{
         stop("make sure relative difference between salmax and slatop within and among PFTs are correct ")
       }
       param_table[[targ_col]][[r]] = assign_val
     }#end of if(isFALSE(ref_pft))
      
    } #end of if(ineq=="complex")
    
    if(param_prefix == "fates_frag_maxdecomp"){
      ref_col = "fates_frag_maxdecomp_O01"
      ref_val = param_table[[ref_col]][[r]]
      param_table[[targ_col]][[r]] = ref_val * tag_max
    } # end of if(param_prefix == "fates_frag_maxdecomp")
} #end of for(p in sequence(n_ref))
} # end of for (r in sequence(nrow(param_table)))

#### Gnerating parameter files ####
#---~---
#   Loop through all ensemble iterations to make NetCDF files
#---~---
cat (" + Generate NetCDF for each ensemble realisation.\n")
for (e in sequence(n_ensemble)){
  #---~---
  #   Configure ensemble file
  #---~---
  node            = ceiling(e/tasks_per_node)
  task           = 1 + ((e-1)%%tasks_per_node)
  ens_label      = sprintf(ens_fmt,e)
  task_label     = sprintf(task_fmt,node,task)
  ncdf4_out_base = paste0(ncdf4_out_pref,"_",task_label,".nc")
  ncdf4_out_file = file.path(ncdf4_out_path,ncdf4_out_base)
  cat("   - ",ens_label,"/",n_ensemble,": Generate file ",ncdf4_out_base,".\n",sep="")
  #---~---
  
  
  #---~---
  #   Create ensemble file, then open it
  #---~---
  dummy    = file.copy(ncdf4_in_file,ncdf4_out_file,overwrite=TRUE)
  nc_conn  = nc_open(ncdf4_out_file,write=TRUE)
  nc_nvars = nc_conn$nvars
  nc_ndims = nc_conn$ndims
  nc_dlist = rep(NA_character_,times=nc_ndims)
  nc_vlist = rep(NA_character_,times=nc_nvars)
  for (d in sequence(nc_ndims)) nc_dlist[d] = nc_conn$dim[[d]]$name
  for (v in sequence(nc_nvars)) nc_vlist[v] = nc_conn$var[[v]]$name
  #---~---
  
  
  #---~---
  #    We only update values that are in the NetCDF file. In case there is any
  # parameter missing, we warn the user.
  #---~---
  p_change = which(  param_config$fates_parameter_name %in% nc_vlist)
  p_miss   = which(! param_config$fates_parameter_name %in% nc_vlist)
  if (length(p_miss) > 0L){
    cat("     > The following parameters do not exist in the parameter file!:\n")
    for (p in p_miss) cat("       ~ ",param_config$fates_parameter_nameer[p],"\n",sep="")
  }#end if (length(p_miss) > 0L)
  #---~---
  
  #---~---
  #   Loop through parameters and update values
  #---~---
  for (p in p_change){
    #---~---
    #   Shorter names
    #---~---
    p_parameter = param_config$fates_parameter_name [p]
    p_corr_name = param_config$corr_name [p]
    p_pft       = param_config$pft       [p]
    p_organ     = param_config$organ     [p]
    p_pft_var   = param_config$pft_var   [p]
    p_all_pft   = param_config$all_pft   [p]
    p_org_var   = param_config$org_var   [p]
    p_all_org   = param_config$all_org   [p]
    #---~---
    
    #---~---
    #   Retrieve parameter value (original and new value).
    #---~---
    p_out_value = ncvar_get(nc_conn,p_parameter,collapse_degen=FALSE)
    p_new_value = param_table[[p_corr_name]][e]
    #---~---
    
    
    #---~---
    #   Set indices for PFT and organs in case they are needed.
    #---~---
    if (p_pft_var) p_pft_idx = if(p_all_pft){sequence(dim(p_out_value)[1])}else{p_pft  }
    if (p_org_var & !p_pft_var) p_org_idx = if(p_all_org){sequence(dim(p_out_value)[1])}else{p_organ}
    if (p_org_var & p_pft_var)  p_org_idx = if(p_all_org){sequence(dim(p_out_value)[2])}else{p_organ}
    #---~---
    
    
    #---~---
    #   Update parameter
    #---~---
    if (verbose){
      cat("     > Update parameter ",p_parameter," ("
          ,sprintf("%g",signif(p_new_value,4)),").\n",sep="")
    }#end if (verbose)
    #---~---
    if (p_org_var & p_pft_var){
      #---~---
      #   PFT- and organ-specific parameter. Update only the sought PFTs and organs.
      #---~---
      if(n_pft>1){p_out_value[p_pft_idx,p_org_idx] = p_new_value}else{
        p_out_value[p_org_idx] = p_new_value
      }
      
      #---~---
    }else if (p_pft_var & !p_org_var){
      #---~---
      #   PFT-specific parameter. Update only the sought PFTs.
      #---~---
      if(n_pft>1){p_out_value[p_pft_idx] = p_new_value}else{
        p_out_value = p_new_value 
      }
      
      #---~---
    }else if (p_org_var & !p_pft_var){
      #---~----
      # organ-specific parameter. 
      #---~---
      p_out_value[p_org_idx] = p_new_value
      
    }
    
    else{
      #---~---
      #    Global value. We multiply the original value by 0 to preserve the
      # original dimensions.
      #---~---
      p_out_value = p_new_value + 0. * p_out_value
      #---~---
    }#end if (p_global)
    #---~---
    
    
    #---~---
    #   Update value in the NetCDF file
    #---~---
    dummy = ncvar_put(nc_conn,p_parameter,p_out_value)
    #---~---
  }#end for (p in p_change)
  #---~---
  
  
  #---~---
  #   Close file
  #---~---
  dummy = nc_close(nc_conn)
  #---~---
}#end for (e in sequence(n_ensemble))
#---~---


#---~---
#   Write a message confirming success. Please keep this message, it is useful for
# tracking whether or not the script ran fine when this script is called externally.
#---~---
cat("\n")
cat("-----------------------------------------------------\n")
cat(" SUCCESS! All ensemble parameter files were created! \n")
cat("-----------------------------------------------------\n")
cat("\n")
#---~---




