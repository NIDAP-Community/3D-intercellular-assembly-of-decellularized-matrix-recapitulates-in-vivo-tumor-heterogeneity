# Expression Heatmap [CCBR] [scRNA-seq] [Bulk] (89a32987-10d9-4233-91c3-e9adf3dcc517): v578
ccbr1223_human_ExpressionHeatmap <- function(ccbr1223_human_BatchCorrectedCounts,ccbr1223_human_metadata) {
    ## This function uses pheatmap to draw a heatmap, scaling first by rows
    ## (with samples in columns and genes in rows)

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(colorspace)
    library(dendsort)
    library(ComplexHeatmap)
    library(dendextend)
    library(tibble)
    library(stringr)
    library(RColorBrewer)
    library(dplyr)
    library(grid)
    library(gtable)
    library(gridExtra)
    library(gridGraphics)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    counts_matrix <- ccbr1223_human_BatchCorrectedCounts
    sample_metadata <- ccbr1223_human_metadata
    gene_column_name <- "Gene"
    sample_name_column <- "Sample"
    sample_label_column = "Label"
    samples_to_include = c("ECM_HT29_1","ECM_HT29_2","ECM_HT29_3","ECM_HT29_4","ECM_HT29_5","ECM_HT29_6","CA_HT29_1","CA_HT29_2","CA_HT29_3","CA_HT29_4","CA_HT29_5","CA_HT29_6")
    
    #Gene Parameters
    include_all_genes <- FALSE
    filter_top_genes_by_variance = TRUE
    top_genes_by_variance_to_include <- 500
    specific_genes_to_include_in_heatmap = "None" 
    cluster_genes <- TRUE
    gene_distance_metric <- "correlation"
    gene_clustering_method <- "average"
    display_gene_dendrograms = TRUE
    display_gene_names <- FALSE
    center_and_rescale_expression <- TRUE

    #Sample Parameters
    cluster_samples <- TRUE
    arrange_sample_columns <- FALSE
    order_by_gene_expression <- FALSE
    gene_to_order_columns <- " "
    gene_expression_order <- "low_to_high"
    smpl_distance_metric <- "correlation"
    smpl_clustering_method <- "complete"
display_smpl_dendrograms <- TRUE
    reorder_dendrogram <- FALSE
    reorder_dendrogram_order <- c()
    display_sample_names <- TRUE
    manually_rename_samples <- FALSE
    samples_to_rename <- c()

    #Annotation
    group_columns <- c("Treatment_CellLine","Cell_Line","Treatment","Extraction_Batch")
    assign_group_colors <- FALSE
    assign_color_to_sample_groups <- c()
    group_colors <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")

    #Visual Parameters:
    heatmap_color_scheme <- "Default"
    autoscale_heatmap_color <- TRUE
    set_min_heatmap_color <- -2
    set_max_heatmap_color <- 2
    aspect_ratio <- "Auto"
    legend_font_size <- 10 
    gene_name_font_size <- 4
    sample_name_font_size <- 8
    display_numbers <- FALSE

    #Advanced Parameters
    return_z_scores <- FALSE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(include_all_genes == TRUE && filter_top_genes_by_variance == TRUE){
        stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
    }

    if((cluster_samples == TRUE && arrange_sample_columns == TRUE) | (arrange_sample_columns == TRUE && order_by_gene_expression == TRUE) | 
    (order_by_gene_expression == TRUE && cluster_samples == TRUE)) {
     stop("ERROR: Only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' may be set as TRUE at one time. Leaving all three of these FALSE will result in the sample column order matching the 'Samples to Include'.")   
    }

    ## --------- ##
    ## Functions ##
    ## --------- ##

    getourrandomcolors<-function(k){
        seed=10
        n <- 2e3
        ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
        ourColorSpace <- as(ourColorSpace, "LAB")
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(currentColorSpace, k, iter.max=20)
        return( unname(hex(LAB(km$centers))))
    }

    ## Begin pal() color palette function∂:
    pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
        if (n < 1L) {
            return(character(0L))
        }
        h <- rep(h, length.out = 2L)
        c <- c[1L]
        l <- rep(l, length.out = 2L)
        power <- rep(power, length.out = 2L)
        rval <- seq(1, -1, length = n)
        rval <- hex(
            polarLUV(
                L = l[2L] - diff(l) * abs(rval)^power[2L], 
                C = c * abs(rval)^power[1L],
                H = ifelse(rval > 0, h[1L], h[2L])
            ),
            fixup=fixup, ...
        )
        if (!missing(alpha)) {
            alpha <- pmax(pmin(alpha, 1), 0)
            alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                width = 2L, upper.case = TRUE)
            rval <- paste(rval, alpha, sep = "")
        }
        return(rval)
    } 
    # End pal() color palette function:

    ## Begin doheatmap() function:
    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
        #require(pheatmap)
        #require(dendsort)
        col.pal <- np[[col]]
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # Define metrics for clustering
        drows1 <- gene_distance_metric
        dcols1 <- smpl_distance_metric
        minx = min(dat)
        maxx = max(dat)
        if (autoscale_heatmap_color) {
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=100)
            legbreaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)
        # Run cluster method using 
        hcrow = hclust(dist(dat), method=gene_clustering_method)
        hc = hclust(dist(t(dat)), method=smpl_clustering_method)

        if (FALSE) {
            sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
        } else {
            sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        }
        if (clus) {
            colclus <- sort_hclust(hc)
        } else {
            colclus = FALSE
        }
        if (clus2) {
            rowclus <- sort_hclust(hcrow)
        } else {
            rowclus = FALSE
        }
        if (display_smpl_dendrograms) {
            smpl_treeheight <- 25
        } else {
            smpl_treeheight <- 0
        }
        if (display_gene_dendrograms) {
            gene_treeheight <- 25
        } else {
            gene_treeheight <- 0
        }
        hm.parameters <- list(
            dat, 
            color=col.pal,
            legend_breaks=legbreaks,
            legend=TRUE,
            scale="none",
            treeheight_col=smpl_treeheight,
            treeheight_row=gene_treeheight,
            kmeans_k=NA,
            breaks=breaks,
            display_numbers=dispnum,
            number_color = "black",
            fontsize_number = 8,
            height=80,
            cellwidth = NA, 
            cellheight = NA, 
            fontsize= legend_font_size,   
            fontsize_row=gene_name_font_size,
            fontsize_col=sample_name_font_size,
            show_rownames=rn, 
            show_colnames=cn,
            cluster_rows=rowclus, 
            cluster_cols=clus,
            clustering_distance_rows=drows1, 
            clustering_distance_cols=dcols1,
            annotation_col = annotation_col,
            annotation_colors = annot_col,
            labels_col = labels_col
        )
        mat = t(dat)
        callback = function(hc, mat) {
            dend = rev(dendsort(as.dendrogram(hc)))
            if(reorder_dendrogram == TRUE) {
                dend %>% dendextend::rotate(reorder_dendrogram_order) -> dend
            } else {
                dend %>% dendextend::rotate(c(1:nobs(dend))) 
            }
            as.hclust(dend)
        }

        ## Make Heatmap
        phm <- do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
        
    }
    # End doheatmap() function.

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Build different color spectra options for heatmap:
    np0 = pal(100) 
    np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1) # Blue to Red
    np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)) # Red to Vanilla
    np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) # Violet to Pink
    np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #Red to yellow to blue
    np5 = colorRampPalette(c("steelblue","white", "red"))(100)  # Steelblue to White to Red

    ## Gather list of color spectra and give them names for the GUI to show.
    np = list(np0, np1, np2, np3, np4, np5)
    names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")

    ## Parse input counts matrix. Subset by samples.
    df1 <- counts_matrix
    # Swap out Gene Name column name, if it's not 'Gene'.
    if(gene_column_name != "Gene"){
        # Drop original Gene column
        df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
        # Rename column to Gene
        colnames(df1)[which(colnames(df1) == gene_column_name)] <- 'Gene'
    }
    # Get sample columns
    samples_to_include <- samples_to_include[samples_to_include != gene_column_name]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

    # Build new counts matrix containing only sample subset chosen by user.
    df1 <- df1[,append("Gene", samples_to_include)]
    df.orig = df1
    df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
    df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
    #df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
    row.names(df.mat) <- df$Gene
    rownames(df.mat) <- str_wrap(rownames(df.mat),30) #for really long geneset names
    df.mat <- as.data.frame(df.mat)

    ## Subset counts matrix by genes.
    # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
    if (include_all_genes == FALSE) {
        # Add user-submitted gene list (optional).
        genes_to_include_parsed = c()
        genes_to_include_parsed = strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
        #genes_to_include_parsed = gsub("_"," ",genes_to_include_parsed)
        df.mat[genes_to_include_parsed,] -> df.final.extra.genes
        if(filter_top_genes_by_variance == TRUE) {
            # Want to filter all genes by variance.
            df.final = as.matrix(df.mat)
            var <- matrixStats::rowVars(df.final)
            df <- as.data.frame(df.final)
            rownames(df) <- rownames(df.final)
            df.final <- df
            df.final$var <- var
            df.final %>% rownames_to_column("Gene") -> df.final 
            df.final %>% dplyr::arrange(desc(var)) -> df.final
            df.final.extra.genes = dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
            df.final = df.final[1:top_genes_by_variance_to_include,]
            df.final = df.final[complete.cases(df.final),]
            # Rbind user gene list to variance-filtered gene list and deduplicate.
            df.final <- rbind(df.final, df.final.extra.genes)
            df.final <- df.final[!duplicated(df.final),] 
            rownames(df.final) <- df.final$Gene
            df.final$Gene <- NULL
            df.final$var <- NULL
        } else {
            # Want to use ONLY user-provided gene list.
            df.final <- df.final.extra.genes
            df.final <- df.final[!duplicated(df.final),]
            # Order genes in heatmap by user-submitted order of gene names.
            df.final <- df.final[genes_to_include_parsed,]
            #df.final$Gene <- NULL
        }
    } else {
        df.final <- df.mat
        df.final$Gene <- NULL
    }
    
        ## Optionally apply centering and rescaling (default TRUE).
    if (center_and_rescale_expression == TRUE) {
            tmean.scale = t(scale(t(df.final)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
    } else {
            tmean.scale = df.final
    }

    if(order_by_gene_expression == TRUE){
        gene_to_order_columns <- gsub(" ","",gene_to_order_columns)
        if(gene_expression_order == "low_to_high"){
        tmean.scale <- tmean.scale[,order(tmean.scale[gene_to_order_columns,])] #order from low to high 
        } else{
        tmean.scale <- tmean.scale[,order(-tmean.scale[gene_to_order_columns,])] #order from high to low  
        }
    }

    df.final <- as.data.frame(tmean.scale)

    ## Parse input sample metadata and add annotation tracks to top of heatmap.
    annot <- sample_metadata
    # Filter to only samples user requests.
    annot=annot %>% dplyr::filter(.data[[sample_name_column]] %in% samples_to_include) 
    
    # Arrange sample options.
    if(arrange_sample_columns) {
      annot = annot[match(samples_to_include,annot[[sample_name_column]]),]
      for(x in group_columns){
                annot[,x]= factor(annot[,x],levels=unique(annot[,x]))  
        }       
      annot = annot %>% dplyr::arrange_(.dots=group_columns,.by_group = TRUE)
      df.final <- df.final[,match(annot[[sample_name_column]],colnames(df.final))] 
    }
    
    # Build subsetted sample metadata table to use for figure.
    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    group_colors <- colorlist[group_colors]

    annot %>% dplyr::select(group_columns) -> annotation_col    
    annotation_col = as.data.frame(unclass(annotation_col))
    annotation_col[] <- lapply(annotation_col,factor)
    x <- length(unlist(lapply(annotation_col,levels)))
    if(x>length(group_colors)){
        k=x-length(group_colors)
        more_cols<- getourrandomcolors(k) 
        group_colors <- c(group_colors, more_cols)
    }
    rownames(annotation_col) <- annot[[sample_label_column]]
    annot_col = list()
    b=1
    i=1
    while (i <= length(group_columns)){
        nam <- group_columns[i]
        grp <- as.factor(annotation_col[,i])
        c <- b+length(levels(grp))-1
        col = group_colors[b:c]
        names(col) <- levels(grp)
        assign(nam,col)
        annot_col = append(annot_col,mget(nam))
        b = c+1
        i=i+1
    }

    if(assign_group_colors == TRUE){
            colassign <- assign_color_to_sample_groups
            groupname <- c()
            groupcol <- c() 
            for (i in 1:length(colassign)) {
                groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
                groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
            }
            annot_col[[1]][groupname] <- groupcol
    }

    ## Setting labels_col for pheatmap column labels.
    
    ## Set order of columns based on smaple name input
    #colnames(df.final)%>%print
    #df.final=df.final[,c(samples_to_include)]

    if (manually_rename_samples == TRUE) {
        # Use user-provided names to rename samples.
        replacements = samples_to_rename
        old <- c()
        new <- c()
        labels_col <- colnames(df.final)
        for (i in 1:length(replacements)) {
            old <- strsplit(replacements[i], ": ?")[[1]][1]
            new <- strsplit(replacements[i], ": ?")[[1]][2]
            old=gsub("^[[:space:]]+|[[:space:]]+$","",old)
            new=gsub("^[[:space:]]+|[[:space:]]+$","",new)
            labels_col[labels_col==old]=new           
        }
    } else {
        
        #annot=annot[match(colnames(df.final),annot[[sample_name_column]]),]
        #print(annot[[sample_label_column]])
        #print(colnames(df.final))
        #colnames(df.final)=annot[[sample_label_column]]

        old=annot[[sample_name_column]]
        new=annot[[sample_label_column]]
            names(old)=new
            df.final=rename(df.final,any_of(old))
        labels_col <- colnames(df.final)
    }

    ## Print number of genes to log.
    print(paste0("The total number of genes in heatmap: ", nrow(df.final)))

    ## Make the final heatmap.
    p <- doheatmap(dat=df.final, clus=cluster_samples, clus2=cluster_genes, ht=50, rn=display_gene_names, cn=display_sample_names, col=heatmap_color_scheme, dispnum=display_numbers)
    p@matrix_color_mapping@name <- " "
    p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
    p@column_title_param$gp$fontsize <- 10
    print(p)

    ## If user sets toggle to TRUE, return Z-scores.
    ## Else return input counts matrix by default (toggle FALSE).
    ## Returned matrix includes only genes & samples used in heatmap.
    if(return_z_scores){
        df.new <- data.frame(tmean.scale) # Convert to Z-scores.
        df.new %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    } else {
        df.final %>% rownames_to_column("Gene") -> df.new
        return(df.new)
    }
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
print("template_ccbr1223_human_ExpressionHeatmap.R #########################################################################")
# Setup Compute Environment
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
############################

# Processing input variable: var_ccbr1223_human_BatchCorrectedCounts
var_ccbr1223_human_BatchCorrectedCounts<-readRDS(paste0(rds_output,"/var_ccbr1223_human_BatchCorrectedCounts.rds"))

if (!('Seurat' %in% class(var_ccbr1223_human_BatchCorrectedCounts))) { if (!(class(var_ccbr1223_human_BatchCorrectedCounts) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ccbr1223_human_BatchCorrectedCounts <- as.data.frame(var_ccbr1223_human_BatchCorrectedCounts)}}
#############################


# Processing input variable: var_ccbr1223_human_metadata
var_ccbr1223_human_metadata<-readRDS(paste0(rds_output,"/var_ccbr1223_human_metadata.rds"))

if (!('Seurat' %in% class(var_ccbr1223_human_metadata))) { if (!(class(var_ccbr1223_human_metadata) %in% c('data.frame', 'RFilePath', 'character', 'list'))) {var_ccbr1223_human_metadata <- as.data.frame(var_ccbr1223_human_metadata)}}
#############################

# Saving Output Dataset
invisible(graphics.off())
var_ccbr1223_human_ExpressionHeatmap<-ccbr1223_human_ExpressionHeatmap(var_ccbr1223_human_BatchCorrectedCounts,var_ccbr1223_human_metadata)
saveRDS(var_ccbr1223_human_ExpressionHeatmap, paste0(rds_output,"/var_ccbr1223_human_ExpressionHeatmap.rds"))
