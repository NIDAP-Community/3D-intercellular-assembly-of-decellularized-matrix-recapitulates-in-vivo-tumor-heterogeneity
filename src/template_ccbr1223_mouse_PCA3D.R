# PCA 3D [CCBR] (0309860a-90f0-490d-8e61-69b483923da1): v38
ccbr1223_mouse_PCA3D <- function(ccbr1223_mouse_BatchCorrectedCounts, ccbr1223_mouse_metadata) {

## This function generates an interactive 3D PCA plot

## --------- ##
## Libraries ##
## --------- ##

library(edgeR)
library(colorspace)
library(plotly)
library(htmlwidgets)
    
## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Basic parmaters
Counts_Matrix <- ccbr1223_mouse_BatchCorrectedCounts
Sample_Metadata_Table <- ccbr1223_mouse_metadata
FeatureID_Name_Column <- "Gene"
Sample_Names_Column <- "Sample"
Samples_to_Include <- c("CA_MC38_1","CA_MC38_2","CA_MC38_3","CA_MC38_4","CA_MC38_5","ECM_MC38_1","ECM_MC38_3","ECM_MC38_4","ECM_MC38_5","ECM_MC38_6","CA_CT26_1","CA_CT26_2","CA_CT26_3","CA_CT26_4","CA_CT26_5","CA_CT26_6","ECM_CT26_1","ECM_CT26_2","ECM_CT26_3","ECM_CT26_4","ECM_CT26_5","ECM_CT26_6")
Group_Column <- "Treatment_CellLine"
Plot_Labels_Column <-"Label"
Samples_to_Rename_Manually <- c()
# Advanced parameters
Outlier_Samples_to_Remove <- c()
Use_CPM <- FALSE  
# Visualization parameters
Points_Size <-8
Point_Color <- c()
Plot_Title <- "PCA 3D"    

##--------------- ##
## Error Messages ##
## -------------- ##

## --------- ##
## Functions ##
## --------- ##

# generate random colors
getourrandomcolors <- function(k) {
  seed = 10
  n <- 2e3
  ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
  ourColorSpace <- as(ourColorSpace, "LAB")
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  set.seed(seed)  
  km <- kmeans(currentColorSpace, k, iter.max = 20)
  return(unname(hex(LAB(km$centers))))
}

## --------------- ##
## Main Code Block ##
## --------------- ##

###########################################
#This code block does input data validation

# get the sample data

# remove outliers
Samples_to_Include <-
  Samples_to_Include[!Samples_to_Include %in% Outlier_Samples_to_Remove]
# include samples
Samples_to_Include <-
  Samples_to_Include[Samples_to_Include != FeatureID_Name_Column]
  cat(sprintf("Total number of samples included in the PCA plot: %g", length(Samples_to_Include)))
Counts_Matrix[, append(FeatureID_Name_Column, Samples_to_Include)] -> df
Sample_Metadata_Table <-
  Sample_Metadata_Table[Sample_Metadata_Table[, Sample_Names_Column]  %in% Samples_to_Include, ]
Sampinfo <- Sample_Metadata_Table
colnames(df)[colnames(df) == FeatureID_Name_Column] <- "Gene"
df -> edf.orig

###### PCA plot ##############

# evaluate and display PCA plot
rownames(Sampinfo) <- Sampinfo[, Sample_Names_Column]
Sampinfo <-
  Sampinfo[match(colnames(edf.orig[, -1]), Sampinfo[, Sample_Names_Column]),]
Sampinfo = Sampinfo[complete.cases(Sampinfo[, Sample_Names_Column]), ]
cat(paste0(
  "\nTotal number of genes included in the PCA plot: ",
  nrow(edf.orig)
))
edf <- edf.orig[, match(Sampinfo[, Sample_Names_Column], colnames(edf.orig))]
idx = !rowMeans(edf) == 0
if (Use_CPM) {
  edf = edgeR::cpm(edf[idx, ])
}
edf.orig = edf.orig[idx, ]
tedf <- t(edf)
colnames(tedf) <- edf.orig[, 1]
tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
tedf <- tedf[, apply(tedf, 2, var) != 0]
pca <- prcomp(tedf, scale. = T)
pca.df <- as.data.frame(pca$x)
pca.df$group <- Sampinfo[, Group_Column]
pca.df$sample <- Sampinfo[, Plot_Labels_Column]

# pca stats
stats <-
  data.frame(id = paste0("PC", 1:length(pca$sdev)),
             eigenvalue = (pca$sdev) ^ 2) %>% mutate(variance = eigenvalue * 100 / sum(eigenvalue)) %>% mutate(cumvariance =
                                                                                                                 cumsum(variance)) %>% mutate(
                                                                                                                   variance_label = sprintf("%s (%.1f%% variance)", id, variance),
                                                                                                                   cumvariance_label = sprintf("%s (%.1f%% cumulative variance)", id, cumvariance)
                                                                                                                 )
axis_title = sub(" variance", "", stats$variance_label)
# if rename samplea
if (!is.null(Samples_to_Rename_Manually)) {
  if (Samples_to_Rename_Manually != c("")) {
    for (x in Samples_to_Rename_Manually) {
      old <- strsplit(x, ": ?")[[1]][1]
      new <- strsplit(x, ": ?")[[1]][2]
      pca.df$sample <-
        ifelse(pca.df$sample == old, new, pca.df$sample)
    }
  }
}
# set the colors to be used in the plot
if (length(unique(Sampinfo[, Group_Column])) > length(Point_Color)) {
  ## Original color-picking code.
  k = length(unique(Sampinfo[, Group_Column])) - length(Point_Color)
  more_cols <- getourrandomcolors(k)
  Point_Color <- c(Point_Color , more_cols)
} else {
  Point_Color <- Point_Color[1:length(unique(Sampinfo[, Group_Column]))]
}
names(Point_Color) <- unique(Sampinfo[, Group_Column])
# set label size (may add to user-defined parameters in the future)
Label_Font_Size <- 24

# plot PCA
cat("\nRunning PCA...")
fig <- plot_ly(
    pca.df,
    x = ~ PC1,
    y = ~ PC2,
    z = ~ PC3,
    color = ~ group,
    colors = Point_Color,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = Points_Size),
    hoverinfo = 'text',
    text = ~ sample,
    size = Label_Font_Size
  )
legend = TRUE
if (legend == FALSE) {
    fig <- hide_legend(fig)
}
fig <-
    fig %>% layout(
      title = list(text = Plot_Title, size = 5),
      scene = list(
        xaxis = list(title = axis_title[1]),
        yaxis = list(title = axis_title[2]),
        zaxis = list(title = axis_title[3])
    ),
      legend = list(
        itemsizing = 'constant',
        size = 12,
        y = 0.5
    )
)
  
###### OUTPUT ##############
# save in dataset container
# auto removed: output <- new.output()
# auto removed: output_fs <- output$fileSystem()
# html widget
html_file <- sprintf("Plot_%s.html", gsub(" ", "_", Plot_Title))
htmlwidgets::saveWidget(fig, html_file)
# auto removed: output_fs$upload(html_file, html_file)
cat(sprintf("\nFiles saved in the output dataset:\n• %s", html_file))
# components table
component_file <-
    sprintf("Components_%s.csv", gsub(" ", "_", Plot_Title))
    write.csv(stats, row.names = FALSE, file  (component_file, 'w'))
cat(sprintf("\n• %s", component_file))
# coordinates table
coordinate_file <-
    sprintf("Coordinates_%s.csv", gsub(" ", "_", Plot_Title))
write.csv(pca.df,
        row.names = FALSE,
        file(coordinate_file, 'w'))
cat(sprintf("\n• %s", coordinate_file))
  
  ## keep for later solution to quality of the downloaded static image
  #fig = config(fig, toImageButtonOptions = list(format= 'svg', # one of png, svg, jpeg, webp
  #filename= sub(".html","", html_file), height= 300, width= 500,           scale= 1))  
  
# display visualization
print(fig)

}

## ---------------------------- ##
## Global Imports and Functions ##
## ---------------------------- ##

## Functions defined here will be available to call in
## the code for any table.

## --------------- ##
## End of Template ##
## --------------- ##

######### Node Execution Steps ##########
print("template_ccbr1223_mouse_PCA3D.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ccbr1223_mouse_BatchCorrectedCounts
var_ccbr1223_mouse_BatchCorrectedCounts<-readRDS(paste0(rds_output,"/var_ccbr1223_mouse_BatchCorrectedCounts.rds"))

if (!('Seurat' %in% class(var_ccbr1223_mouse_BatchCorrectedCounts))) { if (!(class(var_ccbr1223_mouse_BatchCorrectedCounts) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ccbr1223_mouse_BatchCorrectedCounts <- as.data.frame(var_ccbr1223_mouse_BatchCorrectedCounts)}}
#############################


# Processing input variable: var_ccbr1223_mouse_metadata
var_ccbr1223_mouse_metadata<-readRDS(paste0(rds_output,"/var_ccbr1223_mouse_metadata.rds"))

if (!('Seurat' %in% class(var_ccbr1223_mouse_metadata))) { if (!(class(var_ccbr1223_mouse_metadata) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ccbr1223_mouse_metadata <- as.data.frame(var_ccbr1223_mouse_metadata)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ccbr1223_mouse_PCA3D<-ccbr1223_mouse_PCA3D(var_ccbr1223_mouse_BatchCorrectedCounts,var_ccbr1223_mouse_metadata)
saveRDS(var_ccbr1223_mouse_PCA3D, paste0(rds_output,"/var_ccbr1223_mouse_PCA3D.rds"))
