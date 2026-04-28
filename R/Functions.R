library(tidyverse)
library(readxl)
library(reshape2)
library(funrar)
library(fundiversity)
library(bipartite)
library(magrittr)


the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  # panel.border = element_rect(colour = "black", fill=NA),
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  #legend.box="vertical",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)

Perform_PCA_variables=function(d){
  
  d=d%>%
    dplyr::select(where(~!all(is.na(.))))%>%
    drop_na(.)
  
  fit=psych::principal(d,nfactors = 4,rotate = "varimax")
  
  return(list(plot=ggarrange(plotPCA(fit$scores,fit$loadings,fit$Vaccounted,1,2,c(min(fit$scores[,1]),max(fit$scores[,1])),c(min(fit$scores[,1]),max(fit$scores[,1]))),
                             plotPCA(fit$scores,fit$loadings,fit$Vaccounted,3,4,c(min(fit$scores[,1]),max(fit$scores[,1])),c(min(fit$scores[,1]),max(fit$scores[,1]))),ncol=2),
              PCA=fit))
}


plotPCA= function(fitScores, fitLoadings, fitVaccounted, xIndex, yIndex, xLim, yLim, annotateFactor = 2.6,
                  colorLow = "ivory", colorHigh = "purple", colorMid = "pink", midPoint = 0.09,
                  pointSize = 1, pointAlpha = 1, pointColor = "grey85", xTimeSegment = 2.5,
                  yTimeSegment = 2.5, sizeAnnotation = 5) {
  p=ggplot2::ggplot(data = NULL, ggplot2::aes(x = fitScores[, xIndex], y = fitScores[, yIndex])) +
    ggplot2::theme_classic() +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon") +
    ggplot2::geom_point(size = pointSize, alpha = pointAlpha, color = "grey85") +
    ggplot2::geom_segment(data = NULL, ggplot2::aes(x = 0, y = 0, xend = (fitLoadings[, xIndex] * xTimeSegment),
                                                    yend = (fitLoadings[, yIndex] * yTimeSegment)),
                          arrow = ggplot2::arrow(length = ggplot2::unit(1 / 2, units = "picas")),
                          color = "black") +
    ggplot2::annotate("text", label = rownames(fitLoadings), size = sizeAnnotation,
                      x = (fitLoadings[, xIndex] * annotateFactor), y = (fitLoadings[, yIndex] * annotateFactor)) +
    ggplot2::scale_fill_gradient2(low = colorLow, high = colorHigh, mid = colorMid,  midpoint = midPoint) +
    ggplot2::xlim(xLim) +
    ggplot2::ylim(yLim) +
    ggplot2::xlab(paste0("PC ", xIndex, " (",
                         round(x = fitVaccounted["Proportion Var", xIndex] * 100, digits = 1), "%)")) +
    ggplot2::ylab(paste0("PC ", yIndex, " (",
                         round(x = fitVaccounted["Proportion Var", yIndex] * 100, digits = 1), "%)")) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_text(size = 14))
  return(p)
}

#abundance_data=a2
Compute_functional_importance=function(trait_data, abundance_data, 
                                       weighting_cover = c(F,T),
                                       remove_outlier=T,weight=NULL) {
  
  binary_columns=trait_data %>% #checking trait data to make it work with FD
    select(where(~ n_distinct(., na.rm = TRUE) == 2)) %>%
    names()
  
  # Convert binary columns to 0/1
  trait_data = trait_data %>%
    mutate(across(all_of(binary_columns), 
                  ~ {
                    unique_vals = sort(unique(.[!is.na(.)]))
                    case_when(
                      . == unique_vals[1] ~ 0,
                      . == unique_vals[2] ~ 1,
                      TRUE ~ NA_real_
                    )
                  }))
  
  if (is.null(weight)) weight=rep(1,ncol(trait_data))
  
  #Step 0: check outliers of abundances &  we check data
  
  # outlier check (methods from the R Book, Crawley M.J., 2007, p363)
  leverage=function(x){1/length(x)+(x-mean(x))^2/sum((x-mean(x))^2)}
  all_abundances = abundance_data[abundance_data>0]
  
  if (any(leverage(all_abundances)>2*5/length(all_abundances)) & remove_outlier){
    
    pos_in_vec = which(leverage(all_abundances) > 2*5/length(all_abundances))
    
    positive_positions = which(abundance_data > 0, arr.ind = TRUE)
    
    outlier_matrix_positions = positive_positions[pos_in_vec, , drop = FALSE]
    
    for (i in 1:nrow(outlier_matrix_positions)) {
      row_idx = outlier_matrix_positions[i, 1]
      col_idx = outlier_matrix_positions[i, 2]
      abundance_data[row_idx, col_idx] = 1 #setting to 1, so that there is still presence information
    }    
  }
  
  species_names=rownames(trait_data)
  abundance_data=as.matrix(abundance_data)
  
  if (!all(colnames(abundance_data) %in% species_names)) {
    stop("Missing a trait for a given species") #we stop in case
  }
  
  # Step 1: Calculate uniqueness, Fdiv, Feve, Fdiv & rarity using funrar & fundiversity
  
  dist_matrix=compute_dist_matrix(trait_data) %>% as.matrix()
  dist_matrix_scaled=dist_matrix / max(dist_matrix) #scaling btw 0 & 1
  
  
  d_species=d_all_sites=tibble()
  for (weight_k in weighting_cover) {
    
    if (weight_k){ #true
      rel_abund_matrix=abundance_data / rowSums(abundance_data)
    }else{
      
      abundance_data_01=abundance_data
      abundance_data_01[abundance_data_01>0]=1
      rel_abund_matrix=abundance_data_01 / rowSums(abundance_data_01)
    }
    uniq_species=uniqueness(rel_abund_matrix, dist_matrix_scaled)
    distinc_species=distinctiveness(rel_abund_matrix, dist_matrix_scaled)
    mean_distinct=colMeans(distinc_species, na.rm = TRUE)
    
    d_species=rbind(d_species,
                    tibble(
                      Species=names(mean_distinct),
                      uniqueness = uniq_species$Ui,
                      distinctiveness = mean_distinct,
                      Weighted_abundance=weight_k))
    
    d_all_sites=rbind(d_all_sites,
                      melt(distinc_species,varnames = c("Site","Species"),value.name = "Distinctiveness")%>%
                        dplyr::mutate(., Weighted_abundance=weight_k))
  }
  
  #  Step 2: functional diversity indices
  
  weighted_FD=dbFD(trait_data,abundance_data,w.abun = T,corr = "cailliez",w = weight)
  unweighted_FD=dbFD(trait_data,abundance_data,w.abun = F,corr = "cailliez",w = weight)
  
  d_change_FD=tibble()
  for (i in seq_along(species_names)) { #for each species
    removed_species=species_names[i]
    
    # trait data
    traits_remaining=trait_data[rownames(trait_data) != removed_species, , drop = FALSE]
    
    # abundance data
    abund_remaining=abundance_data[, colnames(abundance_data) != removed_species, drop = FALSE]
    
    # dist_matrix_remaining=compute_dist_matrix(traits_remaining) %>% as.matrix()
    # dist_matrix_remaining_scaled=dist_matrix_remaining / max(dist_matrix_remaining) #scaling btw 0 & 1
    
    Try_FD=tryCatch({
      weighted_FD_remaining=dbFD(traits_remaining,abund_remaining,w.abun = T,correction = "cailliez",weight_vector = weight)
      unweighted_FD_remaining=dbFD(traits_remaining,abund_remaining,w.abun = F,correction = "cailliez",weight_vector = weight)
    },
    error =function(e){
      weighted_FD_remaining = NA
      unweighted_FD_remaining = NA
    })
    
    d_change_FD=rbind(d_change_FD,tibble(
      Species=removed_species,
      Site=names(weighted_FD$FEve),
      Change_wFEve=weighted_FD$FEve-weighted_FD_remaining$FEve,
      Change_wFDis=weighted_FD$FDis-weighted_FD_remaining$FDis,
      Change_wFDiv=weighted_FD$FDiv-weighted_FD_remaining$FDiv,
      Change_FEve=unweighted_FD$FEve-unweighted_FD_remaining$FEve,
      Change_FDis=unweighted_FD$FDis-unweighted_FD_remaining$FDis,
      Change_FDiv=unweighted_FD$FDiv-unweighted_FD_remaining$FDiv))
  }
  
  
  return(list(Change_FD=d_change_FD,
              Species_rarity=d_species,
              Sites=d_all_sites))
}



