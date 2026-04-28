rm(list=ls())
source("Functions.R")


# -------------------- Soliveres data competition networks ----
## >> Fungi ----

d=read_xls("./Data/matrices.continuous.xls",3)
d_trait=read_xls("./Data/species.list_traits.xls")%>%
  dplyr::filter(., species %in% d$species)%>%
  dplyr::select(., all_of(names(which(sapply(., function(x) !all(is.na(x)))))))


trait_sp=data.frame("Colony"=bcPower(d_trait$colony.growth.24h,.1),"RGR"=bcPower(d_trait$RGR.intra,.1))
rownames(trait_sp)=d_trait$species

d_value=d_AIC=tibble()
for (trait_k in colnames(trait_sp)){
  
  dist_mat = compute_dist_matrix(trait_sp%>%dplyr::select(., all_of(trait_k)),metric = "euclidean")
  
  mat=matrix(1,1,length(d_trait$species))
  colnames(mat)=d_trait$species
  
  trait_value=trait_sp%>%dplyr::select(., all_of(trait_k))%>%dplyr::pull(.)
  uniq_ajd = uniqueness(mat, dist_mat)$Ui
  uniq_val = log(uniq_ajd+.01)
  distinct_ajd= as.numeric(distinctiveness(mat, dist_mat))
  Ajd=igraph::graph_from_adjacency_matrix(as.matrix(d[,-1]))
  
  d_value=rbind(d_value,tibble(Distinctiveness_k=distinct_ajd,
           Uniqueness_k=uniq_val,
           Trait_value=trait_value,
           Trait_k=trait_k,
           Dataset="Fungi"))

  d_AIC=rbind(d_AIC,tibble(AIC=AIC(lm(distinct_ajd~igraph::alpha_centrality(Ajd)),
                                   lm(uniq_val~igraph::alpha_centrality(Ajd)),
                                   lm(trait_value~igraph::alpha_centrality(Ajd)),
                                   lm(distinct_ajd~igraph::eigen_centrality(Ajd)$vector),
                                   lm(uniq_val~igraph::eigen_centrality(Ajd)$vector),
                                   lm(trait_value~igraph::eigen_centrality(Ajd)$vector),
                                   lm(distinct_ajd~apply(as.matrix(d[,-1]),2,mean)),
                                   lm(uniq_val~apply(as.matrix(d[,-1]),2,mean)),
                                   lm(trait_value~apply(as.matrix(d[,-1]),2,mean)))$AIC,
                           Trait_k=trait_k,
                           Functional_index=rep(c("Distinctiveness","Uniqueness","Raw"),3),
                           Network_index = rep(c("Alpha centrality","Eigen centrality", "Mean IS"), each = 3),
                           Dataset="Fungi"
  ))
}




## >> Plants ----

d = read_xls("./Data/matrices.continuous.xls", 4)
d_trait = read_xls("./Data/species.list_traits.xls") %>%
  dplyr::slice(1:nrow(d)) %>%
  dplyr::mutate(., Name_species = sapply(1:nrow(.), function(x) {
    parts = strsplit(.$species[x], "\\.")[[1]]
    return(stringr::str_to_upper(paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))))
  })) %>%
  dplyr::filter(., Name_species %in% d$species) %>%
  dplyr::select(., all_of(names(which(sapply(., function(x) !all(is.na(x)))))))

d_trait$LeafN=as.numeric(d_trait$LeafN)
# Select the specific traits and apply transformations
trait_sp = data.frame(
  "RGR" = d_trait$RGR.intra,
  "Height" = d_trait$height.vegetative,
  "LDMC" = d_trait$LDMC,
  "LeafN" = d_trait$LeafN,
  "SLA" = d_trait$SLA,
  "Seed_mass" = sqrt(d_trait$Seed.mass)
)%>%
  dplyr::mutate(., Name_species=d_trait$Name_species)
rownames(trait_sp) = d_trait$Name_species

PCA_data=as.data.frame(Perform_PCA_variables(trait_sp%>%dplyr::select(., -Name_species))$PCA$scores)%>%
  dplyr::mutate(., Name_species=rownames(Perform_PCA_variables(trait_sp%>%dplyr::select(., -Name_species))$PCA$scores))

trait_sp=trait_sp%>%
  merge(., PCA_data[,c("RC1","RC2","Name_species")],
        by.x="Name_species",by.y="Name_species",all.x=T)
rownames(trait_sp)=trait_sp$Name_species

trait_sp=trait_sp%>%
  dplyr::select(., -Name_species)

for (trait_k in colnames(trait_sp)) {
  
  species_clean = rownames(trait_sp)[!is.na(trait_sp[, trait_k])]
  
  trait_sp_clean = trait_sp[species_clean, trait_k, drop = FALSE]
  
  d_clean = d %>% dplyr::filter(species %in% species_clean)
  
  d_clean=d_clean[,colnames(d_clean)%in%species_clean]
  
  dist_mat = compute_dist_matrix(trait_sp_clean, metric = "euclidean")
  
  mat = matrix(1, 1, length(species_clean))
  colnames(mat) = species_clean
  
  trait_value = trait_sp_clean %>% dplyr::pull(.)
  uniq_ajd = uniqueness(mat, dist_mat)$Ui
  uniq_val = log(uniq_ajd + 0.01)
  distinct_ajd = as.numeric(distinctiveness(mat, dist_mat))
  Ajd = igraph::graph_from_adjacency_matrix(as.matrix(d_clean))
  
  d_value = rbind(d_value, tibble(
    Distinctiveness_k = distinct_ajd,
    Uniqueness_k = uniq_val,
    Trait_value = trait_value,
    Trait_k = trait_k,
    Dataset="Plants"
  ))
  
  d_AIC = rbind(d_AIC, tibble(
    AIC = AIC(
      lm(distinct_ajd ~ igraph::alpha_centrality(Ajd)),
      lm(uniq_val ~ igraph::alpha_centrality(Ajd)),
      lm(trait_value ~ igraph::alpha_centrality(Ajd)),
      lm(distinct_ajd~igraph::eigen_centrality(Ajd)$vector),
      lm(uniq_val~igraph::eigen_centrality(Ajd)$vector),
      lm(trait_value~igraph::eigen_centrality(Ajd)$vector),
      lm(distinct_ajd ~ apply(as.matrix(d_clean), 2, mean)),
      lm(uniq_val ~ apply(as.matrix(d_clean), 2, mean)),
      lm(trait_value ~ apply(as.matrix(d_clean), 2, mean))
    )$AIC,
    Trait_k = trait_k,
    Functional_index = rep(c("Distinctiveness", "Uniqueness", "Raw"), 3),
    Network_index = rep(c("Alpha centrality","Eigen centrality", "Mean IS"), each = 3),
    Dataset="Plants"
  ))
}

d_AIC_heatmap= d_AIC %>%
  dplyr::group_by(.,Trait_k, Network_index, Dataset) %>%
  dplyr::mutate(.,
    Raw_AIC = AIC[Functional_index == "Raw"][1],
    Uniq_AIC = AIC[Functional_index == "Uniqueness"][1],
    AIC_diff_raw = AIC - Raw_AIC,
    AIC_diff_uniq = AIC - Uniq_AIC
  ) %>%
  ungroup(.) %>%
  dplyr::filter(.,Functional_index == "Distinctiveness") %>%
  dplyr::mutate(.,Facet_group = paste(Trait_k, Network_index, Dataset, sep = " | "))%>%
  reshape2::melt(., measure.vars=c("AIC_diff_raw","AIC_diff_uniq"))%>%
  dplyr::mutate(., variable=recode_factor(variable,
                                          "AIC_diff_uniq"="AIC (distinctiveness \n - uniqueness)",
                                          "AIC_diff_raw"="AIC (distinctiveness \n - raw)"
                                          ))

ggplot(d_AIC_heatmap, aes(x = variable, 
                          y = Facet_group, 
                          fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "red",      
                       mid = "white", 
                       high = "blue",
                       midpoint = 0,
                       name = "AIC Difference") +
  labs(title = "",
       x = "", 
       y = "") +
  the_theme2 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "right")


# -------------------- Co-occurence ----


Estimate_association_matrix = function(d, posterior_0 = TRUE, alpha = 0.05, type_model = "GLM_poisson", 
                                       chains = 4, iter = 2000, warmup = 1000, 
                                       control = NULL, prior = NULL,
                                       family_INLA="zeroinflatedpoisson1",
                                       remove_rare_sp=F) {
  
  # 
  # d=t(mat_quadrat)
  # type_model="INLA"
  # family_INLA= "nbinomial0"
  # remove_rare_sp = T
  
  # type_model_options = c(
  #   "GLM_poisson",          
  #   "GLM_TMB",              
  #   "GLM_binom",
  #   "GLM_poisson_Zero",
  #   "GLM_binom_Zero",  
  #   "GLM_truncbinom_Zero",
  #   "STAN_GLM_poisson",
  #   "STAN_GLM_nb",
  #   "INLA_poisson_Zero",
  #   "BRMS_poisson_Zero",
  # )  
  
  if (remove_rare_sp){ #we remove species which have an occurence lower than 4 quadrats
    if (any(colSums(d > 0)<4)){
      which_T=which(colSums(d > 0)<4)
      d=d[,-which_T]
    }
  }
  
  
  # Matrix of interaction. Columns are species
  A_matrix = matrix(0, ncol(d), ncol(d))
  colnames(A_matrix) = rownames(A_matrix) = colnames(d)
  
  S = ncol(d)
  
  for (i in 1:ncol(d)) { # each species
    if (var(d[, i]) > 0) {
      
      # Prepare data for current species
      response = d[, i]
      predictors = d[, -i, drop = FALSE]
      
      if (type_model == "GLM_poisson") {
        mod = arm::bayesglm(response ~ ., 
                            data = as.data.frame(predictors), 
                            family = "poisson",
                            control = list(maxit = 5000, epsilon = 1e-9))
        
        estimcoef = coef(mod)[-1] # get estimation of alpha_{i,j}
        
        if (posterior_0) {
          posterior = coef(arm::sim(mod))[, -1] # get posteriors of beta_ij
          
          for (k in 1:(S - 1)) { # make sure that posteriors don't include 0
            if (sign(quantile(posterior[, k], alpha / 2)) != 
                sign(quantile(posterior[, k], (1 - alpha / 2)))) {
              estimcoef[k] = 0
            }
          }
        }
        alpha_i = rep(0, S)
        alpha_i[-i] = estimcoef
        
      } else if (type_model == "GLM_TMB") {
        data_mod = data.frame(response = response, predictors)
        mod.full = glmmTMB::glmmTMB(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                                    data = data_mod,
                                    family = "poisson",
                                    control = glmmTMB::glmmTMBControl(
                                      optimizer = optim,
                                      optArgs = list(method = "BFGS")
                                    ))
        alpha_i = rep(0, S)
        alpha_i[-i] = glmmTMB::fixef(mod.full)$cond[-1]
        
      } else if (type_model == "GLM_binom") {
        data_mod = data.frame(response = response, predictors)
        mod.full = glmmTMB::glmmTMB(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                                    data = data_mod,
                                    family = "nbinom2",
                                    control = glmmTMB::glmmTMBControl(
                                      optimizer = optim,
                                      optArgs = list(method = "BFGS")
                                    ))
        alpha_i = rep(0, S)
        alpha_i[-i] = glmmTMB::fixef(mod.full)$cond[-1]
        
      } else if (type_model == "GLM_poisson_Zero") {
        data_mod = data.frame(response = response, predictors)
        mod.full = glmmTMB::glmmTMB(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                                    data = data_mod,
                                    family = "poisson", 
                                    ziformula = ~ 1,
                                    control = glmmTMB::glmmTMBControl(
                                      optimizer = optim,
                                      optArgs = list(method = "BFGS")
                                    ))
        alpha_i = rep(0, S)
        alpha_i[-i] = glmmTMB::fixef(mod.full)$cond[-1]
        
      } else if (type_model == "GLM_binom_Zero") {
        data_mod = data.frame(response = response, predictors)
        mod.full = glmmTMB::glmmTMB(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                                    data = data_mod,
                                    family = "nbinom2", 
                                    ziformula = ~ 1,
                                    control = glmmTMB::glmmTMBControl(
                                      optimizer = optim,
                                      optArgs = list(method = "BFGS")
                                    ))
        alpha_i = rep(0, S)
        alpha_i[-i] = glmmTMB::fixef(mod.full)$cond[-1]
        
      } else if (type_model == "GLM_truncbinom_Zero") {
        data_mod = data.frame(response = response, predictors)
        mod.full = glmmTMB::glmmTMB(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                                    data = data_mod,
                                    family = "truncated_nbinom1",
                                    ziformula = ~ 1,
                                    control = glmmTMB::glmmTMBControl(
                                      optimizer = optim,
                                      optArgs = list(method = "BFGS")
                                    ))
        alpha_i = rep(0, S)
        alpha_i[-i] = glmmTMB::fixef(mod.full)$cond[-1]
        
      } else if (type_model == "INLA") {
        require(INLA)
        data_mod = data.frame(response = response, predictors)
        
        mod.full = INLA::inla(as.formula(paste0("response ~", paste0(colnames(data_mod)[-1],collapse ="+"))),
                              data = data_mod,
                              family = family_INLA,
                              control.predictor = list(compute = TRUE),
                              control.compute = list(config = TRUE),
                              control.inla = list(
                                strategy = "adaptive",  
                                control.vb=list(emergency = 50),
                                int.strategy = "ccd",   
                                tolerance = 1e-8,       
                                step.factor = 0.5,      
                                h = 1e-5                
                              ),
                              verbose = FALSE)
        
        # Extract fixed effects
        summary_fixed = mod.full$summary.fixed[-1, ] # remove intercept
        alpha_i = rep(0, S)
        
        
        if (posterior_0) {
          for (k in 1:(S - 1)) { # make sure that posteriors don't include 0
            if (sign(summary_fixed$`0.025quant`[k]) != 
                sign(summary_fixed$`0.975quant`[k])) {
              summary_fixed$mean[k] = 0
            }
          }
        }
        alpha_i[-i] = summary_fixed$mean
        
      } else if (type_model == "BRMS_poisson_Zero") {
        require(brms)
        data_mod = data.frame(response = response, predictors)
        
        # Set default controls if not provided
        if (is.null(control)) {
          control = list(adapt_delta = 0.95, max_treedepth = 12)
        }
        
        mod.full = brms::brm(response ~ .,
                             data = data_mod,
                             family = brms::zero_inflated_poisson(),
                             chains = chains,
                             iter = iter,
                             warmup = warmup,
                             seed = seed,
                             control = control,
                             silent = 2,  # Suppress most output
                             refresh = 0)  # No progress updates
        
        # Extract fixed effects
        posterior_samples = brms::posterior_samples(mod.full, pars = "^b_")
        coef_names = colnames(posterior_samples)
        coef_names = coef_names[!grepl("Intercept", coef_names)]
        
        if (length(coef_names) == (S - 1)) {
          alpha_i = rep(0, S)
          # Get posterior means
          coef_means = apply(posterior_samples[, coef_names], 2, mean)
          alpha_i[-i] = coef_means
          
          # Apply posterior zero check if requested
          if (posterior_0) {
            for (k in 1:(S - 1)) {
              if (sign(quantile(posterior_samples[, k], alpha / 2)) != 
                  sign(quantile(posterior_samples[, k], (1 - alpha / 2)))) {
                alpha_i[-i][k] = 0
              }
            }
          }
        } else {
          warning(paste("BRMS model for species", colnames(d)[i], 
                        "did not return expected number of coefficients"))
          alpha_i = rep(0, S)
        }
        
      } else {
        stop(paste("Unknown model type:", type_model))
      }
      
      A_matrix[i, ] = alpha_i
    }
  }
  
  # # Symmetrize the matrix
  # A_matrix = (A_matrix + t(A_matrix)) / 2
  
  return(A_matrix)
}




df=readxl::read_xlsx("./Data/Final_biodesert.xlsx")
kept_sp=read.table("./Data/Drypop_biodesert_cover.csv",sep=";")%>%
  dplyr::mutate(., Complete_name_sp=paste0(Complete_name,"_",Com))%>%
  dplyr::filter(., Cover>0)%>%
  dplyr::select(., LL,SLA,LDMC,LA,MaxH,MaxLS,Complete_name_sp)%>%
  as.data.frame(.)%>%
  drop_na(.)

all_net=lapply(1:nrow(df),function(plot_id){
  
  all_motif_data=tibble()
  df=readxl::read_xlsx("./Data/Final_biodesert.xlsx")%>%
    add_column(., Complete_name=paste0(gsub(" ","",.$COU),"_",gsub(" ","",.$SITE),"_",.$PLOT))
  
  list_quadrats=list.files("./Data/Cover_quadrats/")
  initial_dataset=readxl::read_xlsx("./Data/Final_biodesert.xlsx")
  
  quadrat_id=paste0(df$Complete_name,".csv")[plot_id]
  
  quadrat_k=read.table(paste0("./Data/Cover_quadrats/",quadrat_id),sep=";")
  Complete_name_quadrat_k=gsub(".csv","",gsub(" ","",quadrat_id))
  
  mat_quadrat=as.matrix(quadrat_k[,-1]) #build the quadrat-abundance matrix
  colnames(mat_quadrat)=paste0("Quadrat_",1:100)
  rownames(mat_quadrat)=quadrat_k[,1]
  mat_quadrat=ceiling(mat_quadrat) #because some observations are integer, we round all values
  
  if (any(is.na(as.numeric(mat_quadrat))==T)){ # in case there is a NA, put a 0 instead
    mat_quadrat[is.na(mat_quadrat)]=0
  }
  if (any(colSums(mat_quadrat)==0)){ # in case there is 0 species abundance in a quadrat, we remove quadrats
    mat_quadrat=mat_quadrat[,-which(colSums(mat_quadrat)==0)]
  }
  if (is.matrix(mat_quadrat)){
    
    if (any(rowSums(mat_quadrat)==0)){ # in case there is 0 species abundance, we remove species
      mat_quadrat=mat_quadrat[-which(rowSums(mat_quadrat)==0),]
    }
    
    if (nrow(mat_quadrat)>4){ #At least five species 
      
      relative_cover_species = rowSums(mat_quadrat)/sum(mat_quadrat)
      nb_species = nrow(mat_quadrat)
      
      #evenness of cover
      evenness_cover = -sum(relative_cover_species*log(relative_cover_species)/log(nb_species))
      
      #inferring association network using GLM
      interaction_network = Estimate_association_matrix(t(mat_quadrat))
    }else{
      interaction_network=NA
    }
  }else{
    interaction_network=NA
  }
  print(plot_id)
  return(list(Interaction_net=interaction_network,Site=df$Complete_name[plot_id]))
})

saveRDS(all_net,"./Data/Motifs/All_interaction_network.rds")


trait_biodesert=read.table("./Data/Drypop_biodesert_cover.csv",sep=";")%>%
  dplyr::mutate(., Complete_name_sp=paste0(Complete_name,"_",Com))%>%
  dplyr::filter(., Cover>0)
dat=trait_biodesert%>%dplyr::select(., Complete_name_sp,Com,Cover)
trait_df=trait_biodesert%>%dplyr::select(., LL,SLA,LDMC,LA,MaxH,MaxLS,Complete_name_sp)%>%
  as.data.frame(.)%>%
  drop_na(.)

rownames(trait_df)=trait_df$Complete_name_sp
dist_mat=compute_dist_matrix(trait_df%>%dplyr::select(., -Complete_name_sp))

colnames(dat)=c("species","site","value")
#test=funrar::funrar_stack(dat%>%dplyr::filter(., species %in% trait_df$Complete_name_sp),
#                          "species",dist_mat)


uniqueness_coms=uniqueness_stack(dat%>%dplyr::filter(., species %in% trait_df$Complete_name_sp),
                                 "species", 
                                 dist_mat)

distinctiveness_coms=distinctiveness_stack(dat%>%dplyr::filter(., species %in% trait_df$Complete_name_sp),
                                 "species",com = "site",abund = "value", 
                                 dist_mat)

all_networks=readRDS("./Data/All_interaction_network.rds")
all_networks[[1]]



trait_biodesert=read.table("./Data/Drypop_biodesert_cover.csv",sep=";")%>%
  dplyr::mutate(., Complete_name_sp=paste0(Complete_name,"_",Com))%>%
  dplyr::filter(., Cover>0)%>%
  dplyr::group_by(Com) %>%
  dplyr::filter(n() > 5) %>%
  ungroup()


# -------------------- Recruit net ----

options(encoding = "latin1")
library(future.apply)
library(igraph)
library(glmmTMB)
library(MuMIn)
library(scales)
library(ciTools)
library(truncnorm)
library(deSolve)
library(lme4)

#
#
# DATA
#
#

cc = read.csv("./Data/Data_S1/CanopyCover.csv")
rn = read.csv("./Data/Data_S1/RecruitNet.csv")

rl = split(rn,rn$Study_site)
cl = split(cc,cc$Study_site)


# COMPUTE X2 test

# compute weighted mean cover
wm.cover = lapply(cl, function(x){
  tot.area = sum(unlist(tapply(x$Sampled_distance_or_area, x$Plot, unique)))
  tapply(x$Cover * x$Sampled_distance_or_area, x$Canopy, function(y) sum(y, na.rm=T)/tot.area)
}
)

# compute total number of recruits
recru.t = lapply(rl, function(x){tapply(x$Frequency, x$Recruit, sum)})


# Select communities with complete data
#1. Remove those that do not provide cover for Canopies
no.cover = sapply(cl,function(x)sum(x$Cover, na.rm=T))
no.cover = names(no.cover)[no.cover==0]
wm.cover = wm.cover[!names(wm.cover) %in% no.cover]
rl = rl[!names(rl) %in% no.cover]
cl = cl[!names(cl) %in% no.cover]


#2. Remove those that do not provide cover for any species that act as nurse
no.cover.any = c()
for(i in 1:length(rl)){
  id.canopy = unique(rl[[i]]$Canopy)
  cover = wm.cover[[i]][id.canopy]
  if(any(is.na(cover))){
    no.cover.any[length(no.cover.any)+1] = names(rl)[i]
  }else{if(any(cover==0)){
    no.cover.any[length(no.cover.any)+1] = names(rl)[i]
  }}
}
wm.cover = wm.cover[!names(wm.cover) %in% no.cover.any]
rl = rl[!names(rl) %in% no.cover.any]
cl = cl[!names(cl) %in% no.cover.any]
recru.t = recru.t[!names(recru.t)%in%c(no.cover.any,no.cover)]

# obtain the life habit.
life.habit = lapply(rl, function(x){
  lh = c(tapply(x$LifeHabit_Canopy, x$Canopy, unique),tapply(x$LifeHabit_Recruit, x$Recruit, unique))
  lh[!duplicated(names(lh))]
})


# obtain the population stability criteria (populations have adults and recruits)
sta.pop = list()
sta.pop.data = data.frame(id = names(cl), precruit= NA, pcover=NA, sps=NA, stable.sps=NA)
for(i in 1:length(cl)){
  canopy.i = cl[[i]]
  canopy.i = tapply(canopy.i$Cover, canopy.i$Canopy, sum, na.rm=T)
  canopy.i = names(canopy.i)[canopy.i>0]
  recrui.i = rl[[i]]
  recrui.i = tapply(recrui.i$Frequency, recrui.i$Recruit, sum, na.rm=T)
  recrui.i = names(recrui.i)[recrui.i>0]
  stable = intersect(canopy.i,recrui.i)
  id.sps = unique(c(recrui.i,canopy.i))
  id.sps = id.sps[!id.sps %in% c("Open","Rock","Fallen_log")]
  stable.info = ifelse(id.sps%in%stable,"S","noS")
  names(stable.info) = id.sps
  sta.pop[[i]] = stable.info
}

# prepare data for X2 test
data.chi = list()
for(i in 1:length(rl)){
  x = rl[[i]]
  wm = wm.cover[[i]]
  x.open = x[x$Canopy=="Open",]
  x = x[x$Canopy!="Open",]
  sps = unique(c(x$Recruit, x$Canopy, names(wm)[names(wm)!="Open"]))
  id.int = paste(rep(sps,each=length(sps)), rep(sps,length(sps)), sep="_div_")
  x$inter_ID = paste(x$Recruit, x$Canopy, sep="_div_")
  freq = tapply(x$Frequency, x$inter_ID, sum)
  freq.open = tapply(x.open$Frequency, x.open$Recruit, sum)
  mis.open = sps [!sps %in% names(freq.open)]
  mis.open.freq = rep(0, length(mis.open))
  names(mis.open.freq) = mis.open
  freq.open =c(freq.open, mis.open.freq)
  data = data.frame(Recruit=sapply(strsplit(id.int,"_div_"), function(x)x[1]),
                     Canopy=sapply(strsplit(id.int,"_div_"), function(x)x[2]), 
                     inter_ID=id.int,
                     Study_site = unique(x$Study_site))
  data$Canopy_Freq = freq[match(id.int,names(freq))]
  data$Canopy_Freq[is.na(data$Canopy_Freq)] = 0
  data$Canopy_cover = wm[match(data$Canopy,names(wm))]
  data$Open_cover = wm["Open"]
  data$Open_Freq = freq.open[match(data$Recruit,names(freq.open))]
  data$recru.t = recru.t[[i]][match(data$Recruit,names(recru.t[[i]]))]
  data = data[data$Canopy_Freq+data$Open_Freq>0 & !is.na(data$Canopy_cover) & data$Canopy_cover>0,]
  data$sta.pop.recruit = sta.pop[[i]][match(data$Recruit,names(sta.pop[[i]]))]
  data$sta.pop.canopy = sta.pop[[i]][match(data$Canopy,names(sta.pop[[i]]))]
  data$life.habit.recruit = life.habit[[i]][match(data$Recruit,names(life.habit[[i]]))]
  data$life.habit.canopy = life.habit[[i]][match(data$Canopy,names(life.habit[[i]]))]
  data.chi[[i]] = data
  print(i)
}
names(data.chi) = names(wm.cover)
data.chi = do.call(rbind,data.chi)
data.chi = split(data.chi, 1:nrow(data.chi))

# conduct X2 test
func.para = function(y){
  chi = chisq.test(x=c(y[,"Canopy_Freq"],y[,"Open_Freq"]),p=c(y[,"Canopy_cover"]/(y[,"Canopy_cover"]+y[,"Open_cover"]),y[,"Open_cover"]/(y[,"Canopy_cover"]+y[,"Open_cover"])),simulate.p.value = T)
  y["stdres"] = chi$stdres[1]
  y["int_p"] = chi$p.value
  y["int_sign"] = ifelse(chi$p.value <= 0.05 & chi$stdres[1] >0, "Positive", "Neutral")
  if(chi$p.value <= 0.05 & chi$stdres[1] <0){y["int_sign"]="Negative"}
  y
}

plan(multisession)  
data.chi = future_lapply(data.chi, func.para)
plan(sequential)
da = do.call(rbind,data.chi)
da = split(da, da$Study_site)

#
#
# RECIPROCAL FACILITATION AND RELATED STATISTICs
#
#

# remove non-woody and transient populations
coms = lapply(da, function(x){
  x = x[!is.na(x$life.habit.recruit),]
  x = x[!is.na(x$life.habit.canopy),]
  x = x[x$Canopy != x$Recruit, ]
  x[x$sta.pop.recruit == "S" & x$sta.pop.canopy == "S" & x$life.habit.recruit == "W" & x$life.habit.canopy == "W",]
})
coms = coms [sapply(coms, nrow)>0]

#compute facilitation  strength (NII)
coms = lapply(coms, function(x){x$w.faci =(x$Canopy_Freq/x$recru.t) / (x$Canopy_cover);x})
coms = lapply(coms, function(x){x$w.open = (x$Open_Freq/x$recru.t) / (x$Open_cover);x})
coms = lapply(coms, function(x) {x$w.faci.rel = x$w.faci - x$w.open; x})


# Detect shortest paths to reciprocity, compute returned beneficts ("reward") and negative interactions with species in shortest path
for(i in 1:length(coms)){
  net = coms[[i]]
  net = net[net[,"int_sign"]=="Positive",1:2]
  net = net[net[,1]!=net[,2],]
  g = graph_from_edgelist(as.matrix(net), directed = TRUE)
  if(nrow(net) >1){
    for(j in 1:nrow(net)){
      path = all_shortest_paths(g, from = net[j,2], to = net[j,1])$res
      if(length(path)>0){
        sps_return = sapply(path,function(x)x$name[2])
        int.canopy.asrecruit = coms[[i]][coms[[i]][,1]==net[j,2] & coms[[i]][,2]%in%sps_return,]
        int.canopy.ascanopy = coms[[i]][coms[[i]][,1]%in%sps_return & coms[[i]][,2]==net[j,2],]
        coms[[i]][coms[[i]][,1]==net[j,1] & coms[[i]][,2]==net[j,2], "short_path"] = length(path[[1]])-2
        coms[[i]][coms[[i]][,1]==net[j,1] & coms[[i]][,2]==net[j,2], "reward_short_path"] = mean(int.canopy.asrecruit[int.canopy.asrecruit$int_sign== "Positive","w.faci.rel"])
        n.neg.loop = c()
        for(l in 1:length(path)){
          n.neg.loop[l] = sum(coms[[i]][,1]%in%path[[l]]$name & coms[[i]][,2]==net[j,2] & coms[[i]]$int_sign=="Negative")	
        }
        coms[[i]][coms[[i]][,1]==net[j,1] & coms[[i]][,2]==net[j,2], "nneg_short_path"] = mean(n.neg.loop)
      }
    }
  }else{coms[[i]][,c("short_path","reward_short_path","nneg_short_path")] = NA}
  if(!any(names(coms[[i]])=="short_path"))coms[[i]][,c("short_path","reward_short_path","nneg_short_path")] = NA
  print(i)
}	

#
#
# RECIPROCAL FACILITATION STATS AT THE POPULATION LEVEL
#
#

# compute degree as nurse, as recruit and reciprocal
res = list()
for (i in 1:length(coms)){
  dai = coms[[i]]
  id.sps = unique(c(dai$Recruit, dai$Canopy))
  resi = data.frame(sps=id.sps, site=names(coms)[i])
  for(j in 1:length(id.sps)){
    recrui = dai[dai$Recruit == id.sps[j] & dai$int_sign== "Positive",]
    nursei = dai[dai$Canopy == id.sps[j] & dai$int_sign== "Positive",]
    resi[j, "dr"] = nrow(recrui)
    resi[j, "dn"] = nrow(nursei)
    resi[j, "drf"] =  sum(!is.na(nursei$short_path), na.rm=T)
    
  }
  res[[i]] = resi
  print(i)	
}	
rest = do.call(rbind, res)

res.pop = data.frame(prop=rep(NA,5),prop.low=NA,prop.high=NA)
rownames(res.pop) = c("facilitated", "facilitator", "both", "none", "reciprocal")
res.pop["facilitated","prop"] = sum(rest$dr>0 & rest$dn==0)/nrow(rest)
res.pop["facilitator","prop"] = sum(rest$dr==0 & rest$dn>0)/nrow(rest)
res.pop["both","prop"] = sum(rest$dr>0 & rest$dn>0)/nrow(rest)
res.pop["none","prop"] = sum(rest$dr==0 & rest$dn==0)/nrow(rest)
res.pop["reciprocal","prop"] = sum(rest$drf>0)/sum(rest$dr>0 & rest$dn>0)

boot = list()
for(i in 1:1000){
  res.b = data.frame(prop=rep(NA,5))
  rownames(res.b) = c("facilitated", "facilitator", "both", "none", "reciprocal")	
  res.boot = res[sample(1:length(res),length(res),replace=T)]
  rest = do.call(rbind, res.boot)
  res.b["facilitated","prop"] = sum(rest$dr>0 & rest$dn==0)/nrow(rest)
  res.b["facilitator","prop"] = sum(rest$dr==0 & rest$dn>0)/nrow(rest)
  res.b["both","prop"] = sum(rest$dr>0 & rest$dn>0)/nrow(rest)
  res.b["none","prop"] = sum(rest$dr==0 & rest$dn==0)/nrow(rest)
  res.b["reciprocal","prop"] = sum(rest$drf>0)/sum(rest$dr>0 & rest$dn>0)
  boot[[i]] = res.b
  print(i)
}
rest = do.call(rbind, res)

boot.p = do.call(cbind,lapply(boot, function(x)x[,1]))
res.pop[,"prop.low"] = apply(boot.p,1,quantile,0.025)
res.pop[,"prop.high"] = apply(boot.p,1,quantile,0.975)


#
#
# RECIPROCAL FACILITATION STATS AT THE INTERACTION LEVEL
#
#


int.lev = do.call(rbind, coms)
direct = sum(int.lev$short_path==0,na.rm=T)/sum(int.lev$int_sign=="Positive")
indirect = sum(int.lev$short_path>0,na.rm=T)/sum(int.lev$int_sign=="Positive")
all = sum(!is.na(int.lev$short_path))/sum(int.lev$int_sign=="Positive")
table(int.lev$short_path)/sum(int.lev$short_path>0,na.rm=T)

res.b = data.frame(direct=NA, indirect=NA,all=NA)
tl = list()
for(i in 1:1000){
  kkboot = do.call(rbind, coms[sample(1:length(coms),length(coms),replace=T)])
  res.b[i,"direct"] = sum(kkboot$short_path==0,na.rm=T)/sum(kkboot$int_sign=="Positive")
  res.b[i,"indirect"] = sum(kkboot$short_path>0,na.rm=T)/sum(kkboot$int_sign=="Positive")
  res.b[i,"all"] = sum(!is.na(kkboot $short_path))/sum(kkboot$int_sign=="Positive")
  tl[[i]] = table(kkboot$short_path)/sum(kkboot$short_path>0,na.rm=T)
  print(i)
}

res.int = rbind(c(direct,indirect,all),apply(res.b,2,quantile,p=c(0.025,0.975)))
rownames(res.int) = c("observed","ci.low","ci.high")
colnames(res.int) = c("Only.direct", "Only.indirect", "All")

tl = lapply(tl,function(x){x = x[match(1:7,names(x))];names(x)=1:7;x[is.na(x)]=0;x})
tm = do.call(rbind,tl)
apply(tm,2,quantile,p=c(0.025,0.975))


#
#
# RECIPROCAL FACILITATION STATS AT THE COMMUNITY LEVEL
#
#


res.com = lapply(coms, function(x){
  direct = sum(x$short_path==0, na.rm=T)/sum(x$int_sign=="Positive")
  indirect = sum(x$short_path>0, na.rm=T)/sum(x$int_sign=="Positive")
  all = sum(!is.na(x$short_path))/sum(x$int_sign=="Positive")
  c(direct,indirect,all)
}
)

res.com = data.frame(do.call(rbind,res.com))
colnames(res.com) = c("direct","indirect","all")
res.com[is.na(res.com)] = 0

sum(res.com$all >0)/nrow(res.com) # proportion of communities with at least one facilitation loop
sum(res.com$all[res.com$all>0])/sum(res.com$all >0)# avg. community level degree of reciprocity 

boot.c = c()
for(i in 1:1000)boot.c[i]=sum(sample(res.com$all,replace=T) >0)/nrow(res.com)
quantile(boot.c,p=c(0.025,0.975))

boot.c1 = c()
for(i in 1:1000){
  res.boot = sample(res.com$all,replace=T)
  boot.c1[i]=sum(res.boot[res.boot>0])/sum(res.boot >0)
}
quantile(boot.c1,p=c(0.025,0.975))

#
#
# RELATIONSHIP BETWEEN FACILITATION STRENGHT AND NUMBER OF INTERMEDIARY SPS (i.e. sorthest path length to reciprocity)
#
#

# Average data per population
dat.avg.sps = lapply(coms, function(x){
  y = x[x$short_path>0 & !is.na(x$short_path),]	
  reward_short_path = tapply(y$reward_short_path,y$Canopy, mean)
  mean.path = tapply(y$short_path,y$Canopy, mean)
  p.neg.loop = tapply(y$nneg_short_path/y$short_path , y$Canopy, mean)
  data.frame(cbind(reward_short_path,mean.path,p.neg.loop))
})

for(i in 1:length(dat.avg.sps))if(nrow(dat.avg.sps[[i]]>0))dat.avg.sps[[i]]$Study_site = names(coms)[i]
dat.avg.sps= dat.avg.sps[sapply(dat.avg.sps,nrow)>0]
for(i in 1:length(dat.avg.sps))dat.avg.sps[[i]]$species = rownames(dat.avg.sps[[i]])
dat.sps = do.call(rbind,dat.avg.sps)

#fit model
m1 = glmmTMB(reward_short_path ~ mean.path  + (1|Study_site) + (1|species), 
              family = Gamma(link = "log"), 
              data = dat.sps)
summary(m1)
r.squaredGLMM(m1)


#plot (Fig. 3A main text)
nd = expand.grid(mean.path=seq(min(dat.sps$mean.path),max(dat.sps$mean.path),length=100))
predictions = predict(m1, newdata = nd, se.fit = TRUE, re.form = NA)
upper = predictions$fit + predictions$se*1.96
lower = predictions$fit - predictions$se*1.96
pred = predictions$fit
with(dat.sps,plot(log(reward_short_path) ~ mean.path, pch=19,cex=0.75,las=1,col=alpha("darkgoldenrod2",0.15), ylab="Log avg. facilitation strength", xlab="Avg. number of intermediate species" ))
polygon(c(seq(min(dat.sps$mean.path),max(dat.sps$mean.path),length=100), rev(seq(min(dat.sps$mean.path),max(dat.sps$mean.path),length=100))), c(upper, rev(lower)), col = alpha("darkgoldenrod2",0.6), border = NA, lwd=2)
lines(pred~seq(min(dat.sps$mean.path),max(dat.sps$mean.path),length=100),col="darkgoldenrod2", lwd=2)
text(expression(paste(R^2,"=",0.32," ; ",italic(P)," < 0.001")), x=4, y=4, col=1,cex=1)
text("A", x=1.2, y=4, col=1,cex=1.5)


#
#
# PROBABILITY NEGATIVE INTERACTIONS WITH SPECIES IN FACILITATION LOOP VS PROBABILITY NEGATIVE INTERACTION IN COMMUNITY
#
#

# empirical probability negative interactions in the community
p.neg.com =sapply(coms, function(x){
  nsps = length(unique(c(x$Recruit, x$Canopy)))
  sum(x$int_sign=="Negative")/((nsps*(nsps-1))- sum(x$int_sign=="Positive"))
})

dat.sps$p.neg.com = p.neg.com[match(dat.sps$Study_site, names(p.neg.com))]

with(dat.sps, wilcox.test(p.neg.loop, p.neg.com, paired=T, alternative = "less"))

#plot (Fig. 3B main text)
hist(asin(sqrt(dat.sps$p.neg.loop)) - asin(sqrt(dat.sps$p.neg.com)), las=1,breaks=50, xlim=c(-asin(sqrt(1)),asin(sqrt(1)))
     ,ylab="Number of populations", xlab="", main="", border="darkgoldenrod4", col="darkgoldenrod2")
xlab=expression(atop(paste(Delta, italic("P"), " negative interactions"), paste("(",italic("P")^"loop"," - ",italic("P")^"community)")))
mtext(side = 1, line = 4, xlab, col = "black")
abline(v=0,col=1,lwd=2, lty=2)
text("B", x=-1.4, y=700, col=1,cex=1.5)

#
#
# RELATIONSHIP RICHNESS AND DEGREE OF INDIRECT RECIPROCITY 
#
#

richness = sapply(coms,function(x)length(unique(c(x$Canopy, x$Recruit))))
degree.ind = sapply(coms,function(kk) sum(kk$short_path>0, na.rm=T)/sum(kk$int_sign=="Positive"))
degree.ind[is.na(degree.ind)] = 0

m1 = glm(richness ~ degree.ind,family = poisson(link = "log"))
summary(m1)
with(summary(m1), 1 - deviance/null.deviance)

#plot (Fig. 4A main text)
pdat = seq(min(degree.ind),max(degree.ind),length=500)
pp = predict(m1, data.frame(degree.ind=pdat))
pp = add_pi(data.frame(richness=pp,degree.ind=pdat), m1, names = c("lwr", "upr"), alpha = 0.05, nsims = 20000)
plot(log(richness) ~ degree.ind, 
     ylab="Log community richness", xlab="Indirect degree of reciprocity", pch="")
lines(pp[,"richness"] ~ pdat, col="darkgoldenrod2", lwd=2)
polygon(c(pdat,rev(pdat)),c(log(pp[,"lwr"]),rev(log(pp[,"upr"]))),col=alpha("darkgoldenrod2",0.2),border=NA)
points(log(richness) ~ degree.ind, pch=20, col="darkgoldenrod2")
text(expression(paste("R"^2,"=",0.75," ; ",italic(P)," < 0.001")), x=0.25, y=6., col=1,cex=1)


#
#
# SIMULATED DYNAMICS IN FACILITATION (+ NEGATIVE) NETWORKS
#
#

# prepare data to generate facilitation and competition matrices 
obs.nets = path.nets = pl.nets = list()
for(i in 1:length(coms)){
  net = coms[[i]]
  net = net[net[,"int_sign"]=="Positive",c("Recruit","Canopy","w.faci.rel","short_path")]
  net = net[net[,1]!=net[,2],]
  
  id.sps = unique(c(coms[[i]][,1],coms[[i]][,2]))
  obs = matrix(0,ncol=length(id.sps),nrow=length(id.sps))
  colnames(obs) = rownames(obs) = id.sps
  pl = paths = obs	
  
  if(nrow(net)>0){
    for(j in 1:nrow(net)){
      obs [net[j,1],net[j,2]] = net[j,3]
    }
    g = graph_from_data_frame(net, directed = TRUE)
    path.net = net.pl = NULL
    for(j in 1:nrow(net)){
      path = all_shortest_paths(g, from = net[j,2], to = net[j,1])$res
      if(length(path)>0){
        sps_loop = unlist(lapply(path,function(x)x$name))
        sps_return = sapply(path,function(x)x$name[2])
        path_length = length(path[[1]])-2
        path.net = rbind(path.net,cbind(sps_loop[sps_loop!=net[j,2]],net[j,2]))
        net.pl = rbind(net.pl,cbind(recruit=net[j,2], canopy=sps_return,short_path=path_length))
      }
    }
    if(!is.null(path.net)){
      path.net = path.net[!duplicated(path.net),]
      net.pl = aggregate(short_path ~ recruit + canopy, data = net.pl, FUN = function(x)mean(as.numeric(x)))
      for(j in 1:nrow(path.net)){
        paths[path.net[j,1],path.net[j,2]] = 1
      }
      for(j in 1:nrow(net.pl)){
        pl[net.pl[j,1],net.pl[j,2]] = net.pl[j,3]
      }
    }
  }
  obs.nets[[i]] = obs
  path.nets[[i]] = paths
  pl.nets[[i]] = pl
  print(i)
}
names(obs.nets) = names(path.nets) = names(pl.nets) = names(coms)

# Observed NII
alphas.pos = unlist(lapply(obs.nets,function(x)x[x>0]))
aerror = qt(0.975,df=length(alphas.pos)-1)*sd(log(alphas.pos))/sqrt(length(alphas.pos))
amean = mean(log(alphas.pos))

# Reproduction rate
rn[rn$Canopy=="Open", "LifeHabit_Canopy"] = "W"
rn = rn[rn$LifeHabit_Canopy =="W" & rn$LifeHabit_Recruit=="W",]
rl = split(rn,rn$Study_site)
rl = rl[names(coms)]
wm.cover = wm.cover[names(coms)]
reprol = list()
for(i in 1:length(rl)){
  x = rl[[i]]
  y = wm.cover[[i]]
  tot = tapply(x$Frequency, x$Recruit, sum)
  open = tapply(x[x$Canopy=="Open","Frequency"], x[x$Canopy=="Open","Recruit"], sum)
  
  recru.open = open[names(open)%in%names(tot)]/tot[names(tot)%in%names(open)]
  reprol[[i]] = recru.open[names(recru.open)%in%names(y)]/y[names(y)%in%names(recru.open)]
}
reprol = unlist(reprol)
error = qt(0.975,df=length(reprol)-1)*sd(log(reprol))/sqrt(length(reprol))
mean = mean(log(reprol))

#
# Lotka y Volterra model 
#

nrun=25
t = seq(0,100,by=1)
prob.competition.basal = 0.95
prob.competition.loop = prob.competition.basal/2.5
slope = 0.075
types = list(c("random.facilitation","random.competition"),c("increased.facilitation","random.competition"),c("random.facilitation","avoiding.loop"),c("increased.facilitation","avoiding.loop"))
data.sim = list()

for(p in 1:length(types)){
  type.fac = types[[p]][1]
  type.com = types[[p]][2]
  dym = list()
  for(i in 1:length(obs.nets)){
    mfac = obs.nets[[i]]
    mpcom = path.nets[[i]]
    mpl = pl.nets[[i]]
    N.sps = nrow(mfac)
    nfaci = sum(mfac>0)
    if(type.com == "random.competition"){
      mpcom[] = prob.competition.basal
    }
    if(type.com == "avoiding.loop"){
      mpcom[mpcom==0] = prob.competition.basal
      mpcom[mpcom==1] = prob.competition.loop
    }
    surv = c()
    for(j in 1:nrun){
      mfac[which(mfac>0)] = exp(runif(nfaci, amean-aerror, amean+aerror))
      if(type.fac == "increased.facilitation"){
        mfac[which(mpl>0)] = exp(amean) + rtruncnorm(sum(mpl>0),0,Inf, slope,0.05)*mpl[which(mpl>0)]
      }
      mcom = runif(N.sps^2) < mpcom
      mcom[mcom] = exp(runif(sum(mcom), amean-aerror, amean+aerror))*-1 		
      m = mfac + mcom
      diag(m) = -1
      r= exp(runif(ncol(m),mean-error,mean+error))
      names(r) = colnames(m)
      par=c(r=r, m=as.matrix(m))
      N = rep(0,ncol(m))
      basal = rep(1,ncol(m))
      names(basal) = colnames(m)
      state = c(N = basal)
      mlv=function(time, N, parameters){
        dN = N*r*(1+as.vector(m %*% N))
        list(c(dN))
      }
      dlv=ode.1D(y=state,times=t,func = mlv,parms = par,nspec = 1,atol = 1e-8, rtol = 1e-8)
      if(nrow(dlv)==100+1){
        ss = sum(dlv[nrow(dlv),-1]>=1)
        surv[j] = ifelse(ss==0,ss/(ncol(dlv)-1),(ss-1)/(ncol(dlv)-2))
      }else{surv[j] = NA}
    }
    dym[[i]] = surv
    print(c(p,i))
  }
  
  data.sim[[p]] = data.frame(
    type.info = paste(types[[p]], collapse="_"),
    N.sps = rep(richness,each=nrun),
    sur = unlist(dym),
    ind = rep(degree.ind,each=nrun),
    id = rep(names(obs.nets),each=nrun)
  )
}
dat = do.call(rbind,data.sim)

# modelling 
m1 = lmer(sur ~ ind * type.info + (1|id), data=dat, REML=F)
m2 = lmer (sur ~ ind + type.info + (1|id), data=dat, REML=F)
AICc(m1)
AICc(m2)
r.squaredGLMM (m1)

# Plot (Fig. 4B main text)
newdat = expand.grid(ind=seq(0,1,0.01),type.info=unique(dat$type.info),id=unique(dat$id))
pred = predict(m1,newdat,re.form=NA)[newdat$id==unique(dat$id)[1]]
boot = bootMer(m1,function(x) predict(x, newdata = newdat, re.form = NA), nsim = 200)
pred.boots = boot$t[1:200,newdat$id==unique(dat$id)[1]]
lower = apply(pred.boots,2,quantile,.025)
upper = apply(pred.boots,2,quantile,.975)
refdata = newdat[newdat$id==unique(dat$id)[1],]
id.type = unique(dat$type.info)
cols= c("chocolate4","gold1","darkolivegreen3","darkseagreen3")
x = refdata$type.info == id.type[1]
with(dat,plot(sur ~ ind, las=1,col=cols[rev(as.numeric(as.factor(dat$type.info)))], pch="",cex=0.5,
              ylab="Proportion of survivor populations", xlab="Degree of indirect reciprocity", ylim=c(0,0.4)))
polygon(c(seq(0,1,0.01), rev(seq(0,1,0.01))), c(upper[x], rev(lower[x])), col = alpha(cols[1],0.2), border = NA, lwd=2)
for(i in 1:length(id.type)){
  x = refdata$type.info == id.type[i]
  lines(pred[x] ~ seq(0,1,0.01), col=cols[i], lwd=2)
  polygon(c(seq(0,1,0.01), rev(seq(0,1,0.01))), c(upper[x], rev(lower[x])), col = alpha(cols[i],0.2), border = NA)
}
legend("topleft", bty = "n", c("Row facilitation network"," facilitation strengh ~ length loop"," competition in loop","facilitation and  competition"), lwd=2, col=cols)
text("B", x=.75, y=.39, col=1,cex=1.5)




# -------------------- Coux ----

all_data=lapply(list.files("./Data/Coux/",full.names = T),function(x){
  return(read.table(x,sep = ",",header = T,row.names = 1))
})
names(all_data)=gsub(".csv","",list.files("./Data/Coux/"))

pol_trait = all_data$pollinator_traits
pol_trait$Carrying_structure=as.factor(as.numeric(as.factor(as.character(pol_trait$Carrying_structure))))
pol_trait$soc_sol=as.factor(as.numeric(as.factor(as.character(pol_trait$soc_sol))))
pol_trait$season=as.factor(as.numeric(as.factor(as.character(pol_trait$season))))
pol_trait$daily=as.factor(as.numeric(as.factor(as.character(pol_trait$daily))))
pol_trait=pol_trait%>%
  dplyr::mutate(., across(c("larv_necpol","larv_plant","larv_carr","larv_paras","larv_dung","larv_dinsct"),as.factor))%>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))}))
  

# assign pollinator trait weights
weight = c(0.5, 0.5,1,0.5,0.5,1, rep(1/6, 6), rep(1,3))

pol_abund = all_data$pollinator_abundances
Compute_functional_importance(trait_data = pol_trait,
                              remove_outlier = T,
                              abundance_data = pol_abund,
                              weight = weight)




## plant traits
plant_trait =  all_data$plant_traits
plant_trait$growth_form=as.factor(as.numeric(as.factor(as.character(plant_trait$growth_form))))
plant_trait$annual_perennial=as.factor(as.numeric(as.factor(as.character(plant_trait$annual_perennial))))
plant_trait$flower_symmetry=as.factor(as.numeric(as.factor(as.character(plant_trait$flower_symmetry))))
plant_trait$inflorescence_symmetry=as.factor(as.numeric(as.factor(as.character(plant_trait$inflorescence_symmetry))))
plant_trait$flower_sex=as.factor(as.numeric(as.factor(as.character(plant_trait$flower_sex))))

# assign plant trait weights
pl.weight = c(rep(1, 8), rep(0.25, 4), rep(1, 3))

p = all_data$plant_abundances_bin
Compute_functional_importance(trait_data = plant_trait,
                              remove_outlier = F,
                              abundance_data = p,
                              weight = pl.weight)




ntw=all_data$interactions %>%
  split(., .$Site) %>%
  lapply(., function(x){
    net = matrix(x$Links, nrow=length(unique(x$Pol_sp)), ncol=length(unique(x$Plant_sp)))
    colnames(net)= unique(x$Plant_sp)
    rownames(net) = unique(x$Pol_sp)
    return(net)
  })


# calculating FD outputs
# 
pol.coords = species.coords(pol_trait, a, weight, na.rm=T) # species scores in the trait space
pol.measures = FD_measures(pol.coords$coords, pol.coords$centr, a2) # keep weighted
# abundance matrix to keep track of abundances in the analysis

pl.coords = species.coords(plant_trait, p, pl.weight)
pl.measures = FD_measures(pl.coords$coords, pl.coords$centr, p)

# with weighted centroid: only originality should change
w.pol.coords = species.coords(pol_trait,a2,weight)
w.pol.measures = FD_measures(w.pol.coords$coords, w.pol.coords$centr, a2)
colnames(w.pol.measures)[which(colnames(w.pol.measures)=="orig")] = "w.orig" 
colnames(w.pol.measures)[which(colnames(w.pol.measures)=="uniq")] = "w.uniq" 


