# DEG Gene List [CCBR] (f7fc094e-d625-4053-8f61-b44df3178260): v45
ccbr1223_mouse_DEGGeneList <- function(ccbr1223_mouse_DEG) {

## This function filters DEG table

## --------- ##
## Libraries ##
## --------- ##

library(tidyverse)
library(dplyr)
library(tidyselect)
library(tibble)
library(ggplot2)
library(plotrix)
    
## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Input parameters
deg_table <- ccbr1223_mouse_DEG 

# Basic parameters
gene_names_column <- "Gene"
significance_column <- "pval"
significance_cutoff <- 0.001
change_column <- "logFC"
change_cutoff <- 1
filtering_mode <- "in any contrast"

# Advanced parameters
include_estimates <- c("FC","logFC","tstat","pval","adjpval")
round_estimates <- TRUE
rounding_decimal_for_percent_cells = 0

# Filter parameters
contrast_filter = "none"
contrasts = c()
groups = c()
groups_filter = "none"

# Label parameters
label_font_size = 6
label_distance = 1

# Visualization parameters
y_axis_expansion = 0.08
fill_colors =c("steelblue1","whitesmoke")
pie_chart_in_3d = TRUE
bar_width = 0.4
draw_bar_border = TRUE
force_barchart = TRUE

## -------------------------------- ##
## Parameter Misspecifation Errors  ##
## -------------------------------- ##

## -------------------------------- ##
## Functions                        ##
## -------------------------------- ##

## --------------- ##
## Main Code Block ##
## --------------- ##

## If include_estimates param is empty due to template upgrade,
## then fill it with default values.
if (length(include_estimates) == 0){
    include_estimates <- c("FC","logFC","tstat","pval","adjpval")
}
## select DEG stat columns
estimates = paste0("_",include_estimates)
signif = paste0("_", significance_column)
change = paste0("_", change_column)
deg_table <- deg_table %>% dplyr::select(gene_names_column, ends_with(c(estimates, signif, change)))

contrasts_name = deg_table %>% dplyr::select(ends_with(signif)) %>% colnames() 
contrasts_name = unlist(strsplit(contrasts_name, signif))
if ( contrast_filter == "keep") {
    contrasts_name = intersect(contrasts_name, contrasts)    
} else if ( contrast_filter == "remove") {
    contrasts_name = setdiff(contrasts_name, contrasts)    
}
contrasts_name = paste0(contrasts_name, "_")

groups_name = deg_table %>% dplyr::select(ends_with(c("_mean","_sd"))) %>% colnames() 
groups_name = unique(gsub("_mean|_sd", "", groups_name))
if ( groups_filter == "keep") {
    groups_name = intersect(groups_name, groups)    
} else if ( contrast_filter == "remove") {
    groups_name = setdiff(groups_name, groups)    
}
groups_name = paste0(groups_name,"_")

deg_table <- deg_table %>% dplyr::select(gene_names_column, starts_with(c(groups_name,contrasts_name)))

## select filter variables
datsignif <- deg_table %>% 
    dplyr::select(gene_names_column, ends_with(signif)) %>%
    tibble::column_to_rownames(gene_names_column)
datchange <- deg_table %>%
    dplyr::select(gene_names_column, ends_with(change)) %>% tibble::column_to_rownames(gene_names_column)
genes <- deg_table[,gene_names_column]

## filter genes
significant <- datsignif <= significance_cutoff
changed <- abs(datchange) >= change_cutoff
if (filtering_mode == "in any contrast") {
    selgenes <- apply(significant & changed, 1, any)
    select_genes <- genes[selgenes]
    } else {
    selgenes <- apply(significant & changed, 1, all)
    select_genes <- genes[selgenes]
    
}
# stop if 0 genes selected with the selection criteria
if (length(select_genes) == 0) {
    stop("ERROR: Selection criteria select no genes - change stringency of the Significance cutoff and/or Change cutoff parameters")
}
cat(sprintf("Total number of genes selected with %s ≤ %g and |%s| < %g is %g", significance_column, significance_cutoff, change_column, change_cutoff, sum(selgenes)))

##.output dataset
out <- deg_table %>% dplyr::filter(get(gene_names_column) %in% select_genes)
if (round_estimates) {
    out <- out %>% mutate_if(is.numeric, ~signif(., 3))
}

## do plot
significant <- apply( datsignif, 2, function(x) x <= significance_cutoff )  
changed <- apply( datchange, 2, function(x) abs(x) >= change_cutoff )
dd <- significant & changed
if (draw_bar_border){
    bar_border = 'black'
} else {
    bar_border = NA
}

## If fill_colors is blank due to template upgrade, then
## give it default values.
if(length(fill_colors) == 0){
    fill_colors <- c("steelblue1","whitesmoke")
}

if (filtering_mode == "in any contrast") {
    say_contrast = paste(colnames(dd), collapse=" | ")
    say_contrast = gsub("_pval|_adjpval","", say_contrast)

Var2df <- reshape2::melt(apply(dd, 2, table))
if ("L1" %in% names(Var2df)) {
  Var2df <- Var2df %>%
    rename(Var2 = L1)
}

    tab <-  Var2df %>% 
        dplyr::mutate(Significant=ifelse(Var1, "TRUE", "FALSE")) %>% 
        dplyr::mutate(Significant=factor(Significant, levels=c("TRUE","FALSE")), Count=value, Count_format=format(round(value, 1), nsmall=0, big.mark=",")) %>% 
        dplyr::mutate(Var2=gsub("_pval|_adjpval", "", Var2)) %>%
        group_by(Var2) %>% 
        dplyr::mutate(Percent=round(Count / sum(Count)*100,rounding_decimal_for_percent_cells)) %>% 
        dplyr::mutate(Label=sprintf("%s (%g%%)", Count_format, Percent))

    pp <- ggplot(tab, aes(x="", y=Count, labels=Significant, fill=Significant)) +
        geom_col(width=bar_width, position="dodge", col=bar_border) + facet_wrap(~Var2) +
        scale_fill_manual(values=fill_colors) + theme_bw(base_size=20) +
        xlab("Contrast") + ylab("Number of Genes") +
        geom_text(aes(label=Label), color=c("black"), size=label_font_size, position=position_dodge(width=bar_width), vjust=-label_distance) +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        ggtitle(sprintf("%s<%g & |%s|>%g %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode)) + 
        theme(legend.key.size = unit(3,"line"), legend.position='top') +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        theme(strip.background = element_blank(), strip.text = element_text(size=26)) +
        xlab("") + 
        scale_y_continuous(name="", expand=c(y_axis_expansion, 0))
print(pp) 

} else {
    
    say_contrast = paste(colnames(dd), collapse=" & ")
    say_contrast = gsub("_pval|_adjpval","", say_contrast)
    dd <- apply(dd, 1, function(x) all(x==TRUE))
            
    if (force_barchart) {
        dd = data.frame(dd)   
        colnames(dd) = say_contrast                     
        tab <- reshape2::melt(apply(dd, 2, table))  %>% 
        dplyr::mutate(Significant=ifelse(Var1, "TRUE", "FALSE")) %>% 
        dplyr::mutate(Significant=factor(Significant, levels=c("TRUE","FALSE")), Count=value, Count_format=format(round(value, 1), nsmall=0, big.mark=",")) %>% 
        group_by(Var2) %>% 
        dplyr::mutate(Percent=round(Count / sum(Count)*100,rounding_decimal_for_percent_cells)) %>% 
        dplyr::mutate(Label=sprintf("%s (%g%%)", Count_format, Percent))

    pp <- ggplot(tab, aes(x="", y=Count, labels=Significant, fill=Significant)) +
        geom_col(width=bar_width, position="dodge", col=bar_border) + facet_wrap(~Var2) +
        scale_fill_manual(values=fill_colors) + theme_bw(base_size=20) +
        xlab("Contrast") + ylab("Number of Genes") +
        geom_text(aes(label=Label), color=c("black"), size=label_font_size, position=position_dodge(width=bar_width), vjust=-label_distance) +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        ggtitle(sprintf("%s<%g & |%s|>%g %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode)) + 
        theme(legend.key.size = unit(3,"line"), legend.position='top') +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        theme(strip.background = element_blank(), strip.text = element_text(size=26)) +
        xlab("") + 
        scale_y_continuous(name="", expand=c(y_axis_expansion, 0))
print(pp) 
    } else {
    
    N = c( sum(dd), length(dd)-sum(dd))
    Nk = format(round(as.numeric(N), 1), nsmall=0, big.mark=",")
    P = round(N/sum(N)*100,rounding_decimal_for_percent_cells)
    if (label_font_size > 0) {
        labs = c(sprintf("Significant\n%s (%g%%)", Nk[1], P[1]) , sprintf("Non-Significant\n%s (%g%%)", Nk[2], P[2]))
    } else { 
        labs = NULL
    }
    if (pie_chart_in_3d) {
        pie3D(N, radius=0.8, height=0.06, col = fill_colors, theta=0.9, start=0, explode=0, labels=labs, labelcex=label_font_size, shade=0.7, sector.order=1:2, border=FALSE)
    title(main=sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast), cex.main=4, line=-2)
    } else {
        labs = gsub("\n", ": ", labs)
        pie3D(N, radius=0.8, height=0.06, col = fill_colors, theta=0.9, start=45, explode=0, labels=labs, labelcex=label_font_size, shade=0.7, sector.order=1:2, border=NULL)
    title(main=sprintf("%s<%g & |%s|>%g %s: %s", significance_column, significance_cutoff, change_column, change_cutoff, filtering_mode, say_contrast), cex.main=4, line=-2)
    }
}

} 

return(out)
    
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in the code for any table.

#######################
## End of Template   ##
#######################

######### Node Execution Steps ##########
print("template_ccbr1223_mouse_DEGGeneList.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ccbr1223_mouse_DEG
var_ccbr1223_mouse_DEG<-readRDS(paste0(rds_output,"/var_ccbr1223_mouse_DEG.rds"))

if (!('Seurat' %in% class(var_ccbr1223_mouse_DEG))) { if (!(class(var_ccbr1223_mouse_DEG) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ccbr1223_mouse_DEG <- as.data.frame(var_ccbr1223_mouse_DEG)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ccbr1223_mouse_DEGGeneList<-ccbr1223_mouse_DEGGeneList(var_ccbr1223_mouse_DEG)
saveRDS(var_ccbr1223_mouse_DEGGeneList, paste0(rds_output,"/var_ccbr1223_mouse_DEGGeneList.rds"))
