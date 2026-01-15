# This function runs lm for case control studies
# It requires:
# * a list of df containing the data. df in the list have to have row = sample names; col = species
# * the meta data matching the samples in the list above. Need a variable cales Cohort and one called Sample_ID
# * meta.summary is a summary of the metadata and it will tell me whether there are covariates in the models or not


calcualte_da_lm<-function( list.df,
                            meta,
                            treat = "IBD",
                           meta.summary){
  require(dplyr)
  require(broom)
  # create list with results
  l.df<-list.df
  l.mod<-list.df
  
  
  
  # loop through all the df in the list
  for(d in 1:length(list.df)) {
    # get names of cohorts and check whether the variables are there
    cohort<-names(list.df)[[d]]
    print(paste0("Running the lm models for the Cohort: ", cohort))
    
    # set the parameter to define the linear model structure
    if(meta.summary[meta.summary$Cohort_precise == cohort,]$age_category == 0 &  
       meta.summary[meta.summary$Cohort_precise == cohort,]$sex == 0){
      model.type<-paste0("spec.ab ~", treat)
    } else  if (meta.summary[meta.summary$Cohort_precise == cohort,]$age_category != 0 &  
                meta.summary[meta.summary$Cohort_precise == cohort,]$sex == 0){
      if(length(unique(subset(meta, Cohort_precise == cohort)$age_category) %>% na.omit) > 1)
      {
        model.type<-paste0("spec.ab ~ age_category + ", treat)
      } else
      { model.type<-paste0("spec.ab ~", treat) }
    } else  if (meta.summary[meta.summary$Cohort_precise == cohort,]$age_category == 0 &  
                meta.summary[meta.summary$Cohort_precise == cohort,]$sex != 0) {
      if(length(unique(subset(meta, Cohort_precise == cohort)$sex) %>% na.omit) > 1)
      {
        model.type<-paste0("spec.ab ~ sex + ", treat) 
      } else {model.type<- paste0("spec.ab ~" ,treat) }
      
    }  else  if (meta.summary[meta.summary$Cohort_precise == cohort,]$age_category != 0 &  
                 meta.summary[meta.summary$Cohort_precise == cohort,]$sex != 0) {
      if(length(unique(subset(meta, Cohort_precise == cohort)$age_category) %>% na.omit) > 1) {age<-"age_category +" } else {age <- ""} 
      if(length(unique(subset(meta, Cohort_precise == cohort)$sex) %>% na.omit) > 1)  {sex<-"sex + " } else {sex <- ""}
      model.type<-paste0("spec.ab ~ ", sex, age, treat)
    }
    
    
    tmp.spec<-list.df[[cohort]]
    l.mod[[d]]<-vector(mode = "list", length = ncol(tmp.spec))
    df<-data.frame( 
                    contrast = vector(),
                    estimate = vector(),
                    SE = vector(),
                    t.ratio = vector(),
                    p.value = vector(),
                    bins = vector())
    
    meta.tmp<-meta[meta$Cohort_precise == cohort, ] # this is hard coded-need to change to make it universal
    # order meta and spec so they are in the same order and they can be merged
    meta.tmp<-meta.tmp[order(meta.tmp$SmplID),] # this is hard coded-need to change to make it universal
    tmp.spec<-tmp.spec[order(rownames(tmp.spec)),]
    names(l.mod[[d]])<-names(tmp.spec)
      
    for(s in 1:ncol(tmp.spec)) {
      print(paste0("This is species: ", s))
      spec<-names(tmp.spec)[s] #  take naem of sepcies  
      tmp.spec<-tmp.spec[order(rownames(tmp.spec)),]
      meta.tmp<-meta.tmp[order(meta.tmp$SmplID),]
      if(all(meta.tmp$SmplID != rownames(tmp.spec))) {
        stop("Metadata and Spec abund table are not in the same order")
      } else {
        spec.ab<-tmp.spec[,s]
        tmp<-cbind.data.frame(spec.ab, meta.tmp) # build df with species and metadata
        # transform the variable in factors
        meta.tmp$age_category<-factor(meta.tmp$age_category, levels = c("child", "adult"))
        meta.tmp$sex<-factor(meta.tmp$sex, levels = c("Female", "Male"))
        x<-try(l.mod[[d]][[s]]<-lm(as.formula(model.type), data = tmp))
        # in case model cannot be fit and crushes
        if(class(x)[1] == "try-error"){
          df[s,]<-c(rep(NA, 6),spec)
        } else {
          mod<-l.mod[[d]][[s]]
          emm<-emmeans::emmeans(x, pairwise~ IBD2, adjust = "none")$contrasts %>% as.data.frame()
          emm$bins<-rep(spec, nrow(emm))
          emm$formula<-rep(model.type, nrow(emm))
          
          df<-rbind.data.frame(df, emm)

          
        } 
        l.df[[d]]<-df
      }
    }
  }
  
  
  # now I need to bind the data frames together
  df.merged<-do.call(rbind, l.df)
  # add cohort to df
  df.merged$Cohort<-gsub("[.].*", "", rownames(df.merged))
  
  setClass("lm.da", slots=list(df.da="data.frame", l.mod.da="list"))  
  glmm.da.results <- new("lm.da", df.da = df.merged, l.mod.da = l.mod)
  return(glmm.da.results)
}



