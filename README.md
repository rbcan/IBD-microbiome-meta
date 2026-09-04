# IBD-microbiome-meta

## List of R scripts used for analyses 

### Misc
Resampling: `Resampling_submit.Rmd` \
Prevalence filtration: `Prevalence_filtration_q0.1_submit.Rmd`,
`Pevalence_functions.r` \
Combine results into overview table:`Combine_MGSresults_tables_q0.1_submit.Rmd` \

### Beta-diversity
Beta-diversity overall: `Beta_IBD_submit.Rmd` \
Longitudinal analyses beta-diversity and inflammation index: `Temporal_variation_Beta_submit.Rmd` \
Make beta-diversity figure and statistics: `MakeFig_Beta_submit.Rmd` \


### Differential abundance
Make differential abundance figure and statistics: `MakeFig_DA_0.1_submit.Rmd` \
Differential abundance meta-analyses IBD subtypes and controls: `DA_meta_IBD_UCvsCD_submit.Rmd`, `DA_function_UCvsCD.R` \
Differential abundance meta-analyses IBD and controls: `DA_meta_IBD_submit.Rmd`, `DA_function.R` \


### Inflamation index
Compare inflammation index and differential abundances: `humanRds_vs_DAq0.1_submit.Rmd` \


### Strains
Make Strain figure 3: `MakeFig_Strains_q0.1_submit.Rmd` \
Analyse strain permanova results and compare with differential abundances: `strain_pmva_R2vsDA_lm_Jun23_DAq0.1_submit.Rmd` \
Strains figure 4: `IBD_strains_functions.R`, `_targets.R` \

### Enterosignatures
Make enterosignature figure and statistics: `MakeFig_ES_submit.Rmd` \
Phyla abundances: `MakeFIg_PhylaAb_hr_submit.Rmd` \


## Data

Generally all required inputs and intermediary files are provided needed to run the scripts. \
DAmeans.ALL.0.1.txt	:	Differential abundance analysis results all comparisons \
DAmeans.CDHC.0.1.txt	:	Differential abundance analysis results CD vs control \
DAmeans.CDUC.0.1.txt	:	Differential abundance analysis results CD vs UC \
DAmeans.IBDHC.0.1.txt	:	Differential abundance analysis results IBD vs control \
DAmeans.UCHC.0.1.txt	:	Differential abundance analysis results UC vs control \
DA.num.overview.0.1.txt	:	Numbers of MGS differentially abundant \
IBDmeta.cur.filt.mc.matched.1TimeP.txt	:	Metadata first timepoint \
IBDmeta.cur.filt.mc.matched.txt	:	Metadata first timepoint \
IBDmeta.cur.filt.mc.overview.txt	:	Metadata summary \
IBD.MGSabundance.txt	:	MGS abundances \
IBD.MGStaxa.txt	:	MGS taxonomy \
IQtree_allsites.treefile	:	MGS tree \
MetagStats.csv	:	Stats from MATAFILER run \
MGScorvsDAibd.q0.1.crosslong.csv	:	Correlation between DA effect and inflammation index  \
MGS.matL0.txt	:	Kingdom-aggegated abundances \
MGS.matL1.txt	:	Phylum-aggegated abundances \
MGS.matL2.txt	:	Class-aggegated abundances \
MGS.matL3.txt	:	Order-aggegated abundances \
MGS.matL4.txt	:	Family-aggegated abundances \
MGS.matL5.txt	:	Genus-aggegated abundances \
MGS.matL6.txt	:	Species-aggegated abundances \
phylum_hr_correlation.txt	:	Correlation between phylum abundances and inflammation index \
prev.filt.MGS.txt	:	prevelence filtered MGS \
resampling.20.txt	:	Samples from resampling procedure IBD, controls \
resampling.CDUC.20.txt	:	Samples from resampling procedure IBD subtypes, controls \
resampling.png	:	Resampling procedure \
DA_analysis/beta.dbrda.data.rds	:	dbRDA results \
DA_analysis/beta.df_all_cross.rds	:	Beta-diversity cross-sectional results dataframe \
DA_analysis/beta.l_spec.long.distances.df.20ibd2.rds	:	Temporal variance beta-diversity longitudinal data IBD, controls \
DA_analysis/beta.l_spec.long.distances.df.20.rds	:	Temporal variance beta-diversity longitudinal data IBD subtypes, controls \
DA_analysis/beta.pair.tempInd.rds	:	Average beta-diversity among 4 consecutive samples per Individual \
DA_analysis/beta.pcoa.data.rds	:	PcoA results \
DA_analysis/beta.pmva.bc.cohort.rds	:	Permanova results cohort effect \
DA_analysis/beta.pmva.bc.ibd.rds	:	Permanova results IBD, controls \
DA_analysis/beta.pmvaPair.bc.ibd2.rds	:	Permanova results IBD subtypes, controls \
DA_analysis/dist_hr_indNest_log2_lm_perCohort.rds	:	Nested (individuls nexted in cohort) linear model on intra-individual microbiome dissimilarity associated to inflammation index  \
DA_analysis/l.meta_cdcontrols.rds	:	Meta-analysis differential abundance CD vs controls \
DA_analysis/l.meta.cdhc.prevfilt.rds	:	Prevalence-filtered meta-analysis differential abundance CD vs controls \
DA_analysis/l.meta.cdhc.prevfilt.sign0.1.20.rds	:	Robust significant (q<0.1) species differentially abundant in all 20 resamplings, CD vs controls \
DA_analysis/l.meta.cdhc.prevfilt.sign0.1.rds	:	Significant (q<0.1) species differentially abundant in 1-20 resamplings, CD vs controls \
DA_analysis/l.meta.cduc.prevfilt.rds	:	Prevalence-filtered meta-analysis differential abundance CD vs UC \
DA_analysis/l.meta.cduc.prevfilt.sign0.1.20.rds	:	Robust significant (q<0.1) species differentially abundant in all 20 resamplings, CD vs UC \
DA_analysis/l.meta.cduc.prevfilt.sign0.1.rds	:	Significant (q<0.1) species differentially abundant in 1-20 resamplings, CD vs UC \
DA_analysis/l.meta.prevfilt.rds	:	Prevalence-filtered meta-analysis differential abundance IBD vs controls \
DA_analysis/l.meta.prevfilt.sign0.1.20.rds	:	Robust significant (q<0.1) species differentially abundant in all 20 resamplings, IBD vs controls \
DA_analysis/l.meta.prevfilt.sign0.1.rds	:	Significant (q<0.1) species differentially abundant in 1-20 resamplings, IBD vs controls \
DA_analysis/l.meta.rds	:	Meta-analysis differential abundance IBD vs controls \
DA_analysis/l.meta_uccd.rds	:	Meta-analysis differential abundance UC vs CD \
DA_analysis/l.meta_uccontrols.rds	:	Meta-analysis differential abundance UC vs. controls \
DA_analysis/l.meta.uchc.prevfilt.rds	:	Prevalence-filtered meta-analysis differential abundance UC vs controls \
DA_analysis/l.meta.uchc.prevfilt.sign0.1.20.rds	:	Robust significant (q<0.1) species differentially abundant in all 20 resamplings, UC vs controls \
DA_analysis/l.meta.uchc.prevfilt.sign0.1.rds	:	Significant (q<0.1) species differentially abundant in 1-20 resamplings, UC vs controls \
DA_analysis/perInd_dist_vs_deltaHR.rds	:	Per-indifidual Bray-Curtis dissimilarity vs change in inflammation index \
ES/h_scaled.tsv	:	Enterosignature abundance weights \
ES/model_fit.tsv	:	Enterosignature model fit scores \
ES/reapply_NMF_IBD.ipynb	:	Enterosignature analysis \
humanReads/dist.HL_hr_cond_gat_intrChrt.rds	:	High and low inflammation index groups per cohort \
humanReads/dist.HL_hr_ind_gat_intraInd.rds	:	High and low inflammation index groups within individuals \
humanReads/humRds_cor_MGS_list_CTRL.rds	:	Correlation between DA effect and inflammation index control only  \
humanReads/humRds_cor_MGS_list_IBD.rds	:	Correlation between DA effect and inflammation index IBD only  \
humanReads/humRds_cor_MGS_list.rds	:	Correlation between DA effect and inflammation index  \
humanReads/MGScor_blInd_DA_gat.q0.1.rds	:	Significant correlation between DA effect and inflammation index blocked in individuals \
humanReads/MGScorvsDAibd.crosslong.CTRLonly.q0.1.rds	:	Significant correlation between DA effect and inflammation index control only \ 
humanReads/MGScorvsDAibd.crosslong.IBDonly.q0.1.rds	:	Significant correlation between DA effect and inflammation index IBD only  \
humanReads/MGScorvsDAibd.crosslong.q0.1.rds	:	Significant correlation between DA effect and inflammation index  \
rds/cohort_summary.rds	:	Cohort summary \
rds/meta_hr.rds	:	Metadata with human read percent \
rds/MGS.toptax_DAFigq0.1.rds	:	Top taxa shown in Figure 2 \
rds/overview.MGS.DA_q0.1_strains_cor_genes_PosNegSeparate.rds	:	Overview Data of various analyses combined \
strains/ibd2.strains.pmva.all.June23.rds	:	Permanova results un cophenetic distances of strains IBD, controls \
strains/ibd.strains.pmva.all.June23.rds	:	Permanova results un cophenetic distances of strains IBD subtypes, controls \
strains/StrainClust.overlapIBD_IBD2.June23.rds	:	Overlap of robust strain clusters IBD/controls vs CD/UC/controls \
strains/StrainClust.overlapIBD_IBD2.variable.June23.rds	:	Overlap of variable strain clusters IBD/controls vs CD/UC/controls \
strains/strainR2vsDA.lm.June23_q0.1.rds	:	Linear model comparing strain clusterind and species differential abundances \
strains/strainR2vsTE.inclNonSig.ibd_ibd2.June23.rds	:	Results from linear model comparing strain clusterind and species differential abundances IBD subtypes, controls \
strains/strainR2vsTE.inclNonSig_q0.1.June23.rds	:	Results from linear model comparing strain clusterind and species differential abundances IBD, controls \
strains/strainR2vsTE.incTips_q0.1.June23.rds	:	Results from linear model comparing strain clusterind and species differential abundances IBD, controls (only significant) \
strains/tree.df.prev.June23.rds	:	Numbr of tips per MGS-tree \
