### Functions

# reads in strain assignments from mg-stk
read_in_strain_res <- function(intra_phylo_path) {
  lapply(list.files(intra_phylo_path, full.names = TRUE, pattern = "MGS"), function(i) {
    mgs_file <- paste0(i, "/within/IQtree_allsites.strains.txt")
    mgs <- last(strsplit(i, "/")[[1]])
    if(file.exists(mgs_file)) {
      read.delim(mgs_file, header = FALSE) %>%
        setNames(c("Sample", "Strain")) %>% 
        mutate(MGS = mgs)
    } else {
      data.frame()
    }
  }) %>% bind_rows() %>% 
    filter(!is.na(Strain))
}


# Calculates persistence for each individual 
calculate_persistence <- function(strain_res_df_3t) {
  
  lapply(unique(strain_res_df_3t[["MGS"]]), function(ith_mgs) {
    #print(paste0(ith_mgs))
    x <- filter(strain_res_df_3t, MGS == ith_mgs)
    lapply(unique(x[["IndividualID"]]), function(ith_ind) {
      #print(paste0(ith_ind))
      y <- filter(x, IndividualID == ith_ind) %>% 
        mutate(Time_rel = as.numeric(Time_rel)) %>% 
        arrange(Time_rel)
      
      if(nrow(y) > 1) {
        trX = y[["Time_rel"]]
        sCVP = y[["Strain"]]
        
        strainPairs <- sapply(1:(length(sCVP)-1), function(i) c(sCVP[i], sCVP[i+1]), simplify = FALSE)
        timePairs <- sapply(1:(length(trX)-1), function(i) c(trX[i], trX[i+1]), simplify = FALSE)
        
        sameStrains <- sapply(strainPairs, function(i) i[1] == i[2])
        
        sameStrainsFreq <- mean(table(sCVP)/length(sCVP)*100)
        
        Ncts <- sum(sameStrains)
        Nctd <- sum(!sameStrains)
        
        Dcts <- sum(sapply(1:(length(sCVP)-1), function(i) ifelse(sameStrains[i],
                                                                  trX[i+1]-trX[i],
                                                                  (trX[i+1]-trX[i])/2)))
        Dctd <- sum(sapply(1:(length(sCVP)-1), function(i) ifelse(!sameStrains[i],
                                                                  (trX[i+1]-trX[i])/2,
                                                                  0)))
        ind_res <- Ncts/(Nctd+Ncts)*100
        ind_pers <- Dcts/(Dctd+Dcts)*100
        data.frame(MGS = ith_mgs, ind = ith_ind, Ncts = Ncts, Nctd = Nctd, Dcts = Dcts, Dctd = Dctd, 
                   ind_res = ind_res, ind_pers = ind_pers, ind_freq = sameStrainsFreq,
                   IBD = y[["IBD"]][1], IBD2 = y[["IBD2"]][1], Cohort = y[["Cohort"]][1])
      } else {
        data.frame()
      }
      
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# Creates a data frame for correlation plots
get_df_for_correlations <- function(pers_df, ts2) {
  pers_df %>% 
    group_by(MGS) %>% 
    summarise(totalNcts = sum(Ncts),
              totalNctd = sum(Nctd),
              totalN = length(unique(ind)),
              resilience = sum(Ncts)/(sum(Nctd)+sum(Ncts))*100,
              totalDcts = sum(Dcts),
              totalDctd = sum(Dctd),
              persistence = sum(Dcts)/(sum(Dcts)+sum(Dctd))*100) %>% 
    left_join(select(ts2, c("MGS", "meanTE", "enrichedIn", "robustDA")), by = "MGS") %>% 
    filter(!is.na(meanTE))
}

# Generates figure 4B 
plot_ibd_association_correlations <- function(cor_df, min_ind = 2) {
  cor_df %>% 
    filter(totalN >= min_ind) %>% 
    mutate(Type = case_when(robustDA == "non-sig" ~ "no IBD association",
                            robustDA %in% c("variable", "robust") ~ "IBD association")) %>% 
    ggplot(aes(x = Type, y = persistence, fill = Type)) +
    geom_boxplot(show.legend = FALSE) +
    geom_signif(na.rm = TRUE,
                comparisons = list(c("IBD association", "no IBD association")),
                step_increase = 0.1,
                textsize = 3) +
    scale_fill_manual(values = c("IBD association" = "#fba686", "no IBD association" = "white")) +
    theme_bw(base_size = 12) + 
    xlab("") +
    ylab("Intraindividual persistence [%]") + 
    ggtitle("Persistence of species")
}

# Generates figure 4A
plot_persistence_per_ibd_association <- function(cor_df, min_ind = 2) {
  cor_df %>% 
    filter(robustDA %in% c("variable", "robust"), totalN >= min_ind) %>% 
    mutate(enrichedIn = ifelse(enrichedIn == "CTRL", "control", "IBD")) %>% 
    ggplot(aes(x = persistence, y = abs(as.numeric(meanTE)), color = enrichedIn)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_manual("Enriched in", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    theme_bw(base_size = 12) +
    ylab("|IBD-association|") + 
    stat_cor(method = "spearman", cor.coef.name = "rho", show.legend = FALSE, size = 5) +
    theme(legend.position = "inside", 
          legend.position.inside = c(0.125,0.15), 
          legend.background = element_blank(), 
          legend.box.background = element_rect(colour = "black")) +
    xlab("Intraindividual persistence [%]") +
    ggtitle("Persistence of species affected by IBD")
}

# Generates figure S15A
plot_ibd_association_correlations_ibd_healthy <- function(pers_df, ts2, min_ind = 2) {
  cor_df_supp <- pers_df %>%
    group_by(MGS, IBD) %>%
    summarise(totalNcts = sum(Ncts),
              totalNctd = sum(Nctd),
              totalN = length(unique(ind)),
              resilience = sum(Ncts)/(sum(Nctd)+sum(Ncts))*100,
              totalDcts = sum(Dcts),
              totalDctd = sum(Dctd),
              persistence = sum(Dcts)/(sum(Dcts)+sum(Dctd))*100,
              frequency = mean(ind_freq)) %>%
    left_join(select(ts2, c("MGS", "meanTE", "enrichedIn", "robustDA")), by = "MGS") %>%
    filter(!is.na(meanTE))
  
  cor_df_supp %>%
    filter(robustDA %in% c("variable", "robust"), totalN >= min_ind) %>%
    mutate(IBD = ifelse(IBD == "IBD", "IBD", "control"),
           enrichedIn = ifelse(enrichedIn == "IBD", "IBD", "control")) %>%
    ggplot(aes(x = persistence, y = abs(as.numeric(meanTE)), color = enrichedIn)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_manual("Enriched in", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    theme_bw(base_size = 12) +
    ylab("|IBD-association|") +
    stat_cor(method = "spearman", cor.coef.name = "rho", show.legend = FALSE, size = 5, label.y = c(6.25, 5.5), label.x = 50) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.07,0.14),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")) +
    xlab("Persistence [%]") +
    xlim(c(50, 100)) +
    ggtitle("Persistence of species affected by IBD in control and IBD individuals") +
    facet_wrap(~IBD)
}

# Generates figure S15B
plot_ibd_association_correlations_per_cohort <- function(pers_df, ts2, min_ind = 2) {
  
  cor_df_supp_coh <- pers_df %>%
    group_by(MGS, Cohort) %>%
    summarise(totalNcts = sum(Ncts),
              totalNctd = sum(Nctd),
              totalN = length(unique(ind)),
              resilience = sum(Ncts)/(sum(Nctd)+sum(Ncts))*100,
              totalDcts = sum(Dcts),
              totalDctd = sum(Dctd),
              persistence = sum(Dcts)/(sum(Dcts)+sum(Dctd))*100,
              frequency = mean(ind_freq)) %>%
    left_join(select(ts2, c("MGS", "meanTE", "enrichedIn", "robustDA")), by = "MGS") %>%
    filter(!is.na(meanTE))
  
  coh_sizes <- pers_df %>%
    select(ind, Cohort) %>%
    unique() %>%
    group_by(Cohort) %>%
    summarise(n_ind = n())
  
  cor_df_supp_coh %>%
    filter(robustDA %in% c("variable", "robust"), totalN >= min_ind) %>%
    left_join(coh_sizes) %>%
    mutate(Cohort = ifelse(Cohort == "HMP", "HMP + HMP2", Cohort),
           Cohort2 = paste0(Cohort, " (N=", n_ind, ")"),
           enrichedIn = ifelse(enrichedIn == "IBD", "IBD", "control")) %>%
    ggplot(aes(x = persistence, y = abs(as.numeric(meanTE)), color = enrichedIn)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_manual("Enriched in", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    theme_bw(base_size = 12) +
    ylab("|IBD-association|") +
    stat_cor(method = "spearman", cor.coef.name = "rho", show.legend = FALSE, size = 4, label.x = 50, label.y = c(6, 5.25)) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.07,0.07),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")) +
    xlab("Persistence [%]") +
    ggtitle("Persistence of species affected by IBD per cohort") +
    facet_wrap(~Cohort2, nrow = 2) +
    ylim(c(0,6.25)) +
    xlim(c(50,100))
}


# Extracts species that are persistent and enriched in healthy or IBD
get_species_persistent_enriched_in <- function(cor_df, ts2, enriched_in, min_ind = 2, persistence_threshold = 95) {
  cor_df %>% 
    filter(robustDA %in% c("variable", "robust"), totalN >= min_ind) %>% 
    mutate(enrichedIn = ifelse(enrichedIn == "CTRL", "control", "IBD"),
           meanTE = abs(as.numeric(meanTE))) %>% 
    filter(persistence > persistence_threshold, enrichedIn == enriched_in) %>% 
    left_join(select(ts2, c("MGS", "bins")))
}

# Generates figure 4C with definitions of strain changes
plot_strain_changes_definitions <- function() {
  timepoints_plot <- data.frame("Timepoint" = c(1:5),
                                "Individual 1" = c("strain 1","strain 1",NA,NA,NA),
                                "Individual 2" = c("strain 1","strain 2",NA,"strain 1","strain 1"),
                                "Individual 3" = c("strain 1","strain 2","strain 2",NA,"strain 3"),
                                "Individual 4" = c("strain 1","strain 2","strain 3","strain 2","strain 4"),
                                "Individual 5" = c(NA,"strain 1","strain 2","strain 3","strain 1"),
                                "Individual 6" = c("strain 1", NA, "strain 2", NA, NA)) %>% 
    pivot_longer(2:7, names_to = "Type",
                 values_to = "Strain") %>% 
    mutate(Individual = gsub("Individual.", "", Type)) %>% 
    mutate(Individual = factor(Individual, levels = c(6:1))) %>% 
    ggplot(aes(x = Timepoint, y = Individual, shape = Strain)) +
    geom_point(size = 4) +
    theme_bw(base_size = 12) +
    scale_shape_manual(values = c(15, 16, 17, 18), na.translate = FALSE) +
    theme(legend.position = "bottom") +
    scale_x_continuous(position = "top") + 
    guides(shape = guide_legend(nrow = 2)) +
    ggtitle("Evaluated types of strain changes")
  
  types_plot <- data.frame("Individual" = factor(c(1,2,3,4,5,6), levels = c(6:1)),
                           "Vanishing" = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
                           "Retention" = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
                           "Replacement" = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE),
                           "Variability" = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)) %>% 
    pivot_longer(2:5, names_to = "Type", values_to = "Presence") %>% 
    ggplot(aes(x = Type, y = Individual, fill = Presence)) +
    geom_tile(color = "white") +
    theme_bw(base_size = 12) +
    xlab("Type of strain changes") +
    scale_fill_manual(values = c("TRUE" = "#FA8072", "FALSE" = "#FFE0DA")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom") +
    scale_x_discrete(position = "top") 
  
  strain_replacements_defs <- plot_grid(timepoints_plot, types_plot, nrow = 1, rel_widths = c(4, 4), align = "h") 
}



get_first_x_timepoints_samples <- function(metadata, n_timepoints) {
  first5timepoints <- sapply(unique(metadata[["AssmblGrps"]]), function(ith_ind) {
    x <- filter(metadata, AssmblGrps == ith_ind) %>% 
      mutate(Time_rel = as.numeric(Time_rel))
    if(nrow(x) >= n_timepoints) {
      arrange(x, Time_rel)[["X.SmplID"]][1:n_timepoints]
    }
  }, USE.NAMES = FALSE) %>%  
    unlist()
}


# Calculates strain changes in each individual
get_strain_changes_df_5t <- function(strain_res_df_5t, metadata, first5timepoints) {
  lapply(unique(strain_res_df_5t[["MGS"]]), function(ith_mgs) {
    print(paste0(ith_mgs))
    x <- filter(strain_res_df_5t, MGS == ith_mgs)
    lapply(unique(x[["IndividualID"]]), function(ith_person) {
      print(paste0(ith_person))
      y <- filter(x, IndividualID == ith_person)
      
      z <- filter(metadata, `X.SmplID` %in% first5timepoints, IndividualID == ith_person) %>% 
        select(c("X.SmplID", "IndividualID", "Time_rel", "IBD", "IBD2")) %>% 
        mutate(Time_rel = as.numeric(Time_rel)) %>% 
        left_join(., select(y, c("Sample", "Strain")), by = c("X.SmplID" = "Sample")) %>% 
        arrange(Time_rel)
      print(paste0(nrow(z)))
      
      if(all(is.na(z[["Strain"]])) | nrow(z) != 5) {
        data.frame()
      } else {
        first <- z[["Strain"]][1]
        second <- z[["Strain"]][2]
        last <- z[["Strain"]][nrow(z)]
        rest <- z[["Strain"]][2:nrow(z)]
        last_typed <- if(length(which(!is.na(rest))) == 0) {
          NA
        } else {
          last(rest[which(!is.na(rest))])
        }
        
        first_typed_index <- which(!is.na(z[["Strain"]]))[1]
        first_typed <- z[["Strain"]][first_typed_index]
        all_typed <- z[["Strain"]][which(!is.na(z[["Strain"]]))]
        first_typed_positions <- which(z[["Strain"]] == first_typed)
        secondtyped_index <- which(all_typed != first_typed)[1]
        
        pers <- if(!is.na(last_typed)) {
          ifelse(first_typed == last_typed & length(all_typed) > 1, TRUE, FALSE)
        } else {
          FALSE
        }
        
        repl <- if(length(unique(all_typed)) > 1) {
          ifelse(any(first_typed_positions > secondtyped_index), FALSE, TRUE)
        } else {
          FALSE
        }
        
        data.frame(IndividualID = ith_person,
                   MGS = ith_mgs,
                   IBD = z[["IBD"]][1],
                   IBD2 = z[["IBD2"]][1],
                   strain_vanishing = ifelse((!is.na(first) & is.na(last)) | (is.na(first) & !is.na(second) & is.na(last)), TRUE, FALSE),
                   persistent_strain = pers,
                   strain_replacement = repl,
                   strain_introgression = ifelse(length(unique(filter(z, !is.na(Strain))[["Strain"]])) > 2, TRUE, FALSE),
                   n_strains = length(unique(filter(z, !is.na(Strain))[["Strain"]])),
                   n_samples = nrow(z))
      }
      
    }) %>% bind_rows()
  }) %>% bind_rows()
}

# Summarises strain changes to frequencies in controls and IBD
summarise_strain_changes <- function(test_res_new_5t, metadata, by_cohort = FALSE) {
  
  if(by_cohort == TRUE) {
    test_res_new_5t %>% 
      left_join(unique(select(metadata, c("IndividualID", "Cohort")))) %>% 
      filter(Cohort %in% c("LIBD", "HMP")) %>%
      group_by(Cohort, MGS, IBD) %>% 
      summarise(Vanishing = sum(strain_vanishing)/n(),
                Retention = sum(persistent_strain)/n(),
                Replacement = sum(strain_replacement)/n(),
                Variability = sum(strain_introgression)/n(),
                n_ind = n()) %>% 
      mutate(Cohort = as.factor(ifelse(Cohort == "HMP", "HMP + HMP2", Cohort))) 
  } else {
    test_res_new_5t %>% 
      left_join(unique(select(metadata, c("IndividualID", "Cohort")))) %>% 
      filter(Cohort %in% c("LIBD", "HMP")) %>%
      group_by(MGS, IBD) %>%
      summarise(Vanishing = sum(strain_vanishing)/n(),
                Retention = sum(persistent_strain)/n(),
                Replacement = sum(strain_replacement)/n(),
                Variability = sum(strain_introgression)/n(),
                n_ind = n())
  }
}

# Generates figure 4D
plot_strain_changes_boxplots_all <- function(types_df_new_5t_hmp_libd, strain_changes_test_coh_blocked, min_ind = 5) {
  
  lab_df <- data.frame(
    label = round(strain_changes_test_coh_blocked[["pval"]], 4),
    Type = strain_changes_test_coh_blocked[["Statistic"]]
  )
  
  types_df_new_5t_hmp_libd %>% 
    filter(n_ind >= min_ind) %>% 
    pivot_longer(c("Vanishing", "Retention", "Replacement", "Variability"), names_to = "Type", values_to = "fraction") %>% 
    mutate(IBD = ifelse(IBD == "IBD", "IBD", "control")) %>% 
    ggplot(aes(x = IBD, y = fraction, fill = IBD)) +
    geom_boxplot(alpha = 0.5, outliers = FALSE) +
    geom_jitter(size = 1, alpha = 0.3, aes(color = IBD)) +
    facet_wrap(~Type, nrow = 1) +
    geom_bracket(inherit.aes = FALSE, y.position = 1.02, xmin = "control", xmax = "IBD",
                 label = "", tip.length = c(0.01, 0.01)) +
    geom_text(inherit.aes = FALSE, data = lab_df, aes(x = -Inf, y = -Inf, label = label), hjust = -1.9, vjust = -30.75, size = 3) +
    theme_bw(base_size = 12) +
    scale_color_manual("", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    scale_fill_manual("", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    ylab("Fraction of individuals") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12)) +
    xlab("") + 
    ggtitle("Strain changes among controls and IBD")
}

# Generates figure S16
plot_strain_changes_boxplots_per_cohort <- function(types_df_new_5t_coh, min_ind = 5) {
  types_df_new_5t_coh %>%
    filter(n_ind >= min_ind) %>%
    pivot_longer(c("Vanishing", "Retention", "Replacement", "Variability"), names_to = "Type", values_to = "fraction") %>%
    mutate(IBD = ifelse(IBD == "IBD", "IBD", "control"),
           Cohort = as.character(Cohort)) %>%
    ggplot(aes(x = IBD, y = fraction, fill = IBD)) +
    geom_boxplot(alpha = 0.5, outliers = FALSE) +
    geom_jitter(size = 1, alpha = 0.3, aes(color = IBD)) +
    facet_grid(Cohort~Type) +
    geom_signif(na.rm = TRUE, textsize = 3,
                comparisons = list(c("control", "IBD")),
                step_increase = 0.1)+
    theme_bw(base_size = 12) +
    scale_color_manual("", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    scale_fill_manual("", values = c("control" = "#1F78B4", "IBD" = "#E8A400")) +
    ylab("Fraction of individuals") +
    theme(legend.position = "none") +
    xlab("") +
    ylim(c(0, 1.1)) +
    ggtitle("Strain changes among controls and IBD per cohort")
}

# Performs statistical analysis of strain changes between IBD and controls
get_strain_changes_test_results_blocked_by_cohort <- function(types_df, min_ind = 5) {
  figd_test_df <- types_df %>% 
    filter(n_ind >= min_ind) %>% 
    pivot_longer(c("Vanishing", "Retention", "Replacement", "Variability"), names_to = "Type", values_to = "fraction") %>% 
    mutate(IBD = as.factor(ifelse(IBD == "IBD", "IBD", "control")),
           Cohort = as.factor(Cohort))
  
  lapply(c("Replacement", "Retention","Vanishing", "Variability"), function(ith_measure) {
    print(paste0(ith_measure))
    tres <- coin::wilcox_test(fraction ~ IBD | Cohort, data = filter(figd_test_df, Type == ith_measure))
    z <- tres@statistic@teststatistic
    print(paste0("z: ", z))
    data.frame(Statistic = ith_measure,
               pval = coin::pvalue(tres),
               effsize = z/sqrt(nrow(filter(figd_test_df, Type == ith_measure))))
  }) %>% bind_rows()
}

# Combines subfigures of figure 4 into a single figure
combine_fig4 <- function(fig4a, fig4b, fig4c, fig4d) {
  ab <- plot_grid(fig4b, fig4a, labels = c("A", "B"), nrow = 1, rel_widths = c(3, 4))
  plot_grid(plotlist = list(ab, fig4c, fig4d), labels = c("", "C", "D"), ncol = 1, rel_heights = c(3, 3, 3))
}
