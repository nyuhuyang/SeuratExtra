library(dplyr)
library(magrittr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gplots)
library(Matrix)
library(qlcMatrix)
library(cowplot)
library(cowplot)

# Set a default value if an object is nullf
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
        if (!is.null(x = lhs)) {
                return(lhs)
        } else {
                return(rhs)
        }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
        return(rhs)
    } else {
        return(lhs)
    }
}

# Generic GetAssayData functions will generate Error in t.default(x = GetAssayData(object = object, slot = slot)[features,  : 
# argument is not a matrix
#GetAssayData <- function(object, slot = "data", assay = NULL){
#    assay <- assay %||% DefaultAssay(object = object)
#    assay.data <- GetAssay(object = object, assay = assay)
#    mtrix <- slot(object = assay.data, name = slot)
#    if(class(mtrix) =="dgCMatrix") mtrix = as.matrix(mtrix)
#    return(mtrix)
#}

#' .AddMetaColor: prepare meta.data data frame to store color code
#' @param mat factor at column 1, rownames is cell.names
#' @param colors vector of hex colors
#' @param df_colors data frame to add to meta data directly
#' @example 
# singlerDF = .AddMetaColor(mat = singler$singler[[2]]$SingleR.single$labels,
#                        colors = singler_colors[3:26])
.AddMetaColor <- function(mat, colors){
    if(class(mat) != "data.frame") mat = as.data.frame(mat)
    #mat$cell.names = rownames(mat)
    if(class(mat[,1]) == "factor") mat[,1] = droplevels(mat[,1])
    #mat$index <- as.numeric(as.factor(mat[,1]))
    if(length(colors) < length(unique(mat[,1]))) {
        stop(paste("Not enough colors! Provide at least", 
            length(unique(mat[,1])),"different colors"))
        } else 
    if(length(colors) > length(unique(mat[,1]))) {
        colors = colors[1:length(unique(mat[,1]))]
        }
    if(is.null(names(colors))) names(colors) = levels(as.factor(mat[,1]))
    mat$colors = plyr::mapvalues(mat[,1], from = names(colors),to = colors)
    colnames(mat)[2] = paste0(colnames(mat)[1],".", colnames(mat)[2])
   
     # remove non-color if there is any
    mat[!(mat$cell_types.colors %in% colors),2] = NA
    return(mat)
}


#' AddMetaColor: convert one MetaData label into color scheme and store into MetaData
#' @object seurat object
#' @label colname in metadata
#' @colors vector of hex colors, same format as the result of ExtractMetaColor(object)
# MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = Singler.colors)
AddMetaColor<- function(object, label = NULL, colors = NULL){
    
    if(is.null(label)) label <- FindIdentLabel(object)
    if(is.null(colors)) colors <- SingleR:::singler.colors
    mat = object[[label]]
    colnames(mat) = get("label")
    newMetaData = .AddMetaColor(mat = mat, colors = colors)
    object <- AddMetaData(object = object, metadata = newMetaData)
    
    return(object)
}

#' find Alias gene names
#' @param df marker data frame with Alias notation
#' @param gene gene name
#' @example Alias(df = df_markers, gene = "Cd19")
Alias <- function(df, gene = "HLA-DRB1"){
        
        df = as.data.frame(df)
        df_remove_alias = df[,-grep("Alias",colnames(df))]
        species <- CheckSpecies(gene)
        gene1 = toupper(gene)
        ind <- which(df_remove_alias == gene1, arr.ind = TRUE)
        if(nrow(ind) == 0) return(NULL) # no record
        # unique gene name has alias
        ind <- which(df == gene1, arr.ind = TRUE)
        if(nrow(ind) == 1) {
                cell_type <- colnames(df)[ind[2]]
                gene.alias <-  df[ind[1],paste0(cell_type,".Alias")]
        } else {
                # duplicate gene name and has no alias
                if(nrow(ind) > 1 & (ind[2,2]-ind[1,2]==1)) return(NULL)
        
                # duplicate gene name and has alias
                if(nrow(ind) > 1 & (ind[2,2]-ind[1,2]!=1)) {
                        cell_type <- colnames(df)[ind[1,2]]
                        gene.alias = df[ind[1,1],paste0(cell_type,".Alias")]
                }
        }
        if(is.null(gene.alias)) return(NULL)
        if(CheckSpecies(gene)=="Mouse") {
            gene.alias = paste0(" (",Hmisc::capitalize(tolower(gene.alias)),")")
        }
        if(CheckSpecies(gene)=="Human") {
            gene.alias = paste0(" (",toupper(tolower(gene.alias)),")")
        }
        return(gene.alias)
        }


#' Support FeaturePlot.2
#' @param breaks parameter pass to cut function. greater than or equal to 2. 
#'               when breaks = 0, label as with/without expression
BlendPlot <- function (data.use, features, data.plot, pt.size, pch.use,alpha,
                       cols.use, dim.codes, min.cutoff, max.cutoff, breaks=0,
                       no.axes, no.legend,
                       dark.theme)
{
    num.cols <- length(x = cols.use)
    cols.not.provided <- colors(distinct = TRUE)
    cols.not.provided <- cols.not.provided[!(grepl(pattern = paste(cols.use,
    collapse = "|"), x = cols.not.provided, ignore.case = TRUE))]
    if (num.cols > 4) {
        cols.use <- cols.use[c(1:4)]
    } else if ((num.cols == 2) || (num.cols == 3)) {
        blend <- BlendColors(cols.use[c(num.cols - 1, num.cols)])
        cols.use <- c(cols.use, blend)
        if (num.cols == 2) {
            cols.use <- c(sample(x = cols.not.provided, size = 1),
            cols.use)
        }
    } else if ((num.cols == 1)) {
        if (cols.use %in% rownames(x = brewer.pal.info)) {
            palette <- brewer.pal(n = 3, name = cols.use)
            cols.use <- c(palette, BlendColors(palette[c(2,
            3)]))
        }
        else {
            cols.high <- sample(x = cols.not.provided, size = 2,
            replace = FALSE)
            cols.use <- c(cols.use, cols.high, BlendColors(cols.high))
        }
    } else if (num.cols <= 0) {
        cols.use <- c("yellow", "red", "blue", BlendColors("red",
        "blue"))
    }
    names(x = cols.use) <- c("low", "high1", "high2", "highboth")
    length.check <- vapply(X = list(features, min.cutoff,
    max.cutoff), FUN = function(x) {
        return(length(x = x) != 2)
    }, FUN.VALUE = logical(length = 1))
    if (any(length.check)) {
        stop("An overlayed FeaturePlot only works with two features and requires two minimum and maximum cutoffs")
    }
    data.gene <- stats::na.omit(object = data.frame(data.use[features,]))
    min.cutoff <- c(SetQuantile(cutoff = min.cutoff[1], data = data.gene[features[1],]),
                    SetQuantile(cutoff = min.cutoff[2], data = data.gene[features[2],]))
    max.cutoff <- c(SetQuantile(cutoff = max.cutoff[1], data = data.gene[features[1],]),
                    SetQuantile(cutoff = max.cutoff[2], data = data.gene[features[2],]))
    max.cutoff[max.cutoff == 0] = 0.00001
    cell.names <- colnames(x = data.gene)
    data.gene <- matrix(data = vapply(X = data.gene, 
                                      FUN = function(x) ifelse(test = x < min.cutoff, yes = min.cutoff, no = x), 
                                      FUN.VALUE = c(1,1)), nrow = 2)
    data.gene <- matrix(data = vapply(X = as.data.frame(x = data.gene),
                                      FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff, no = x),
                                      FUN.VALUE = c(1, 1)), nrow = 2)
    data.gene <- as.data.frame(x = data.gene)
    rownames(x = data.gene) <- features
    colnames(x = data.gene) <- cell.names
    
    if (all(data.gene == 0)) {
            stop("All cells have zero expression!")
    }
    if (breaks < 2) {
        
        cuts = as.matrix(data.gene)
        cuts[1, data.gene[1,] > 0] = 1
        cuts[1, data.gene[1,] <= 0] = 0
        cuts[2, data.gene[2,] > 0] = 1
        cuts[2, data.gene[2,] <= 0] = 0
        cuts = as.data.frame(cuts)
        
        } else if ( breaks >= 2){
            
        cuts <- apply(X = data.gene, MARGIN = 1, FUN = cut, breaks = breaks, labels = FALSE)
        cuts.dim <- dim(x = cuts) 
        if (cuts.dim[1] > cuts.dim[2]) {
            cuts <- t(x = cuts)
        }
        }

    # check if minimal is not 0
    #if(any(apply(cuts,1,min) !=0)) {
    #    cuts[apply(cuts,1,min) !=1,] = 1
    #}
    data.cut = apply(X = cuts, MARGIN = 2, FUN = function(x) {
            return(if ((x[1] == 0) && (x[2] > 0)) {
                    "high2"
            } else if ((x[1] > 0) && (x[2] == 0)) {
                    "high1"
            } else if ((x[1] > 0) && (x[2] > 0)) {
                    "highboth"
            } else {
                    "low"
            })
    })
    data.cut <- as.factor(x = data.cut)
    data.cut <- factor(data.cut, levels = c("low","high1","high2","highboth"))
    data.plot$colors <- data.cut

    if(length(pt.size) == 1) {
        data.plot$pt.size <- pt.size
    } else if(length(pt.size) > 1 & length(pt.size) < 5){
        pt.size = c(pt.size, rep(pt.size[length(pt.size)], 4-length(pt.size)))
        data.plot$pt.size <- plyr::mapvalues(x = base::as.character(data.plot$colors),
                                            from = c("low","high1","high2","highboth"),
                                            to = pt.size) %>% as.numeric()
    }

    if(length(alpha) == 1) {
        data.plot$alpha <- alpha
    } else if(length(alpha) > 1 & length(alpha) < 5){
        alpha = c(alpha, rep(alpha[length(alpha)], 4-length(alpha)))
        data.plot$alpha <- plyr::mapvalues(x = base::as.character(data.plot$colors),
                                             from = c("low","high1","high2","highboth"),
                                             to = alpha) %>% as.numeric()
    }
    data.plot = arrange(data.plot, colors)
    
    if(breaks < 2){
        legend.names <- c(high1 = features[1],
                          high2 = features[2], highboth = "Both")
    } else if(breaks < 2){
        legend.names <- c(high1 = paste("High", features[1]),
                          high2 = paste("High", features[2]), highboth = "High both")
    }
    title <- paste0(features, collapse = " + ")
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y,color = colors))
    p <- p + geom_point(size = data.plot$pt.size,
                        alpha = data.plot$alpha, shape = pch.use)
    p <- p + scale_color_manual(values = cols.use, limits = c("low","high1",
    "high2", "highboth"), labels = legend.names, guide = guide_legend(title = NULL,
    override.aes = list(size = 2)))
    if (no.axes) {
        p <- p + labs(title = title, x = "", y = "") + theme(axis.line = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank())
    } else {
        p <- p + labs(title = title, x = dim.codes[1], y = dim.codes[2])
    }
    if (no.legend) {
        p <- p + theme(legend.position = "none")
    }
    p <- p + theme_bw() + theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
    if (dark.theme) {
        p <- p + DarkTheme()
    }
    return(p)
}


BlendColors <- function (..., as.rgb = FALSE) 
{
    colors <- base::as.character(x = c(...))
    if (length(x = colors) < 2) {
        stop("Please provide two or more colors to blend")
    }
    alpha.value <- 255
    if (sum(sapply(X = colors, FUN = grepl, pattern = "^#")) != 
        0) {
        hex <- colors[which(x = grepl(pattern = "^#", x = colors))]
        hex.length <- sapply(X = hex, FUN = nchar)
        if (9 %in% hex.length) {
            hex.alpha <- hex[which(x = hex.length == 9)]
            hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, 
                               stop = 9)
            dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
            dec.vals <- dec.vals/255
            alpha.value <- dec.vals[1]
            for (val in dec.vals[-1]) {
                alpha.value <- alpha.value + (val * (1 - alpha.value))
            }
            alpha.value <- alpha.value * 255
        }
    }
    rgb.vals <- sapply(X = colors, FUN = grDevices::col2rgb)
    if (nrow(x = rgb.vals) != 3) {
        rgb.vals <- t(x = rgb.vals)
    }
    blend <- apply(X = rgb.vals, MARGIN = 1, FUN = mean)
    if (as.rgb) {
        result <- matrix(data = blend, nrow = 3, dimnames = list(c("red", 
                                                                   "green", "blue"), "blend"))
    }
    else {
        result <- grDevices::rgb(matrix(data = blend, ncol = 3), 
                                 alpha = alpha.value, maxColorValue = 255)
    }
    return(result)
}


#' a supprting function for FeaturePlot.1
#' Change ggplot color scale to increase contrast gradient
#' #https://github.com/satijalab/seurat/issues/235
#' @param p ggplot object
#' @param alpha.use Define transparency of points
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
ChangeColorScale <- function(p1, alpha.use = 1,
                             gradient.use = c("yellow", "red"),
                             scaled.expression.threshold = 0.001) {
    # Order data by scaled gene expresion level
    # Compute maximum value in gene expression
    if (length(p1$data$feature_values)>0){                   # SingleFeaturePlot.1 
        p1$data = p1$data[order(p1$data$feature_values),] 
        max.scaled.exp <- max(p1$data$feature_values) 
    }
    
    # Define lower limit of scaled gene expression level
    if (scaled.expression.threshold == 0) {
        scaled.expression.threshold <- min(p1$data$feature_values)+0.001
    }
    
    
    # Fill points using the scaled gene expression levels
    p1$layers[[1]]$mapping$fill <- p1$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p1$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p1 = p1 + guides(colour = FALSE)
    p1 = p1 + scale_colour_gradientn(colours = gradient.use, guide = F,
                                     limits = c(scaled.expression.threshold,
                                                max.scaled.exp),
                                     na.value = "gray97") +
        scale_fill_gradientn(colours = gradient.use,
                             name = expression(atop(expression)),
                             limits = c(scaled.expression.threshold,
                                        max.scaled.exp),
                             na.value = "gray97") +
        scale_alpha_continuous(range = alpha.use, guide = F)
    
    # Return plot
    return(p1)
}


#' check sepcies
#' @param subject could be Seurat object, or gene names
#' @example CheckSpecies("CD8A")
#' @example CheckSpecies(object)
CheckSpecies <- function(subject){
        if(class(subject)=="Seurat"){
                subject = rownames(subject)[1]
        }
        # Human gene
        if(subject == toupper(subject))
                return("Human")
        # Mouse gene
        if(subject == Hmisc::capitalize(tolower(subject)))
                return("Mouse") 

}


# Make color scheme vector
# names must be ordered vector or factor
color_scheme <- function(color,names){
    df_names <- as.data.frame(table(names))
    df_names_Var <- df_names$Freq
    color_code <- base::as.character(unlist(mapply(rep, color, df_names_Var)))
    names(color_code) = names
    return(color_code)
}


# Combine and print multiple PNG,
#' @param ... PNG path
#' @param ncol grid.arrange argments
#' @param do.print TRUE/FALSE print in device
#' @param do.save TRUE/FALSE save in output/date folder
#' @param bottom_text grid.arrange argments
#' @example CombPngs(GSEA.plots.path.list, ncol = 3)
CombPngs <- function(...,ncol = 2, do.print = FALSE, do.save = TRUE, bottom_text = NULL){
    Img <- lapply(..., function(x) grid::rasterGrob(as.raster(png::readPNG(x)),interpolate = TRUE))
    
    path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    if(do.print) {
        do.call(gridExtra::grid.arrange ,c(Img, ncol = ncol, bottom_text = NULL))
    } else if(do.save) {
        jpeg(paste0(path,deparse(substitute(...)),"_CombPngs.jpeg"), units=units, width=10, height=7,res=600)
        do.call(gridExtra::grid.arrange ,c(Img, ncol = ncol, bottom_text = NULL))
        dev.off()
    }
}


#' Produce a contingency table with certain gene expression at different ident
#'
#' @param object seurat object
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, ... 
#' Any argument that can be retreived using FetchData
#' @export Cells_ident nX2 dataframe with orig.ident at column1, Freq at column2
#' @examples
#' CountsbyIdent(SSCs,"Gfra1")
CountsbyIdent <- function(object,subset.name,...){
    "Generate nX2 dataframe with orig.ident at column1, Freq at column2"
    cells.use <- WhichCells(object=object,subset.name=subset.name,...)
    Cells.use <- data.frame("cells"=cells.use,"ident"=sub('_.*', '', cells.use))
    Cells_ident <- as.data.frame(table(Cells.use$ident))
    colnames(Cells_ident)[2] <- subset.name
    return(Cells_ident)
}


#' Convert data frame to list
#'
#' This function will convert a data frame to a list, even if they are unequal length
#'
#' @param df
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
#' genelist <- df2list(brainTxDbSets_df)
df2list <- function(df){
    if(is.matrix(df)) df <- as.data.frame(df)
    list <- lapply(df, as.vector) # as.vector! not as.character
    list <- lapply(list, function(x) x[!is.na(x)])
    list <- lapply(list, function(x) x[!(x == "")])
    list <- lapply(list, function(x) gsub("\\s","", x)) #remove space
    names(list) <- names(df)
    return(list)
}


#' DoHeatmap.1, automatically group top DE genes from FindMakers output
#' @param dge_markers FindAllMarkers results
#' @param features specify the genes
#' @param cells specify the cells
#' @param group.colors hex color character vector matching to group.by
#' @param colors hex color character vector for heatmap
#example   DoHeatmap.1(object,top,Top_n = 15, 
#                   group.order = major_cells,ident.use = "all cell types",
#                   group.label.rot = T,cex.row = 5,remove.key =T)
#
DoHeatmap.1 <- function(object, dge_markers = NULL,features = NULL, cells = NULL,
                        no.legend =F,unique.name= T, Top_n = 10, group.by = "ident",
                        group.bar = TRUE, group.colors = NULL, colors=NULL, disp.min = -2.5, disp.max = NULL, slot = "scale.data",
                        assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45,
                        raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                        combine = TRUE,title = "",title.size = 14,do.print = FALSE,
                        position = "right",save.path = NULL,file.name = NULL,
                        cex.row=12,legend.size = NULL,units="in", width=10, height=7,res=600,...){
    if(is.null(file.name)){
        v <- UniqueName(object = object, fileName = deparse(substitute(object)), unique.name = unique.name)
        v = paste0(v,"_",FindIdentLabel(object))
        if(!no.legend) v = paste0(v, "_Legend")
        file.name = paste0("Heatmap_top",Top_n,"_",v,".jpeg")
    }
    
    if(class(title) != "character") stop("Title is incorrect")
    if (!is.null(x = dge_markers)) {
        
        colnames(dge_markers)[grep("cluster*.",colnames(dge_markers))[1]] = "cluster"
        top <-  dge_markers %>%
            group_by(cluster) %>%
            top_n(Top_n, avg_logFC)
        features = c(base::as.character(top$gene),features)
    }
    if((length(group.colors) != length(unique(Idents(object)))) &
       !is.null(group.colors) ) stop("length of colors do not match!")
    heatmap <- DoHeatmap(object = object, features = features, cells = cells, group.by = group.by,
                         group.bar = group.bar, group.colors, disp.min = disp.min, disp.max = disp.max,
                         slot = slot, assay = assay, label = label, size = size,
                         hjust = hjust, angle = angle, raster = raster, draw.lines = draw.lines,
                         lines.width = lines.width, group.bar.height = group.bar.height,
                         combine = combine)+
        scale_y_discrete(position = position)
    if(pal_gsea) heatmap = heatmap + scale_fill_gradientn(colors = ggsci::pal_gsea()(12))
    if(!is.null(colors)) heatmap = heatmap + scale_fill_gradientn(colors = colors)
    if(!is.null(title)) {
        heatmap = heatmap+ ggtitle(title)+
            theme(plot.title = element_text(size=title.size, hjust = 0.5,face="plain"))
    }
    heatmap = heatmap + theme(axis.text.y = element_text(size = cex.row))
    if(!is.null(legend.size)) {
        heatmap = heatmap + theme(legend.text = element_text(size = legend.size),
                                  legend.title = element_text(size = legend.size*2),
                                  legend.key.size = unit(legend.size/10,"cm"))
    }
    if(no.legend) heatmap = heatmap + NoLegend()
    if(do.print){
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path, "/", file.name),units=units, width=width, height=height,res=res)
        print(heatmap)
        dev.off()
    } else return(heatmap)
}

#' DoHeatmap.2, add more than one annotation bar to the heatmap
#' @param dge_markers FindAllMarkers results
#' @param features specify the genes
#' @param cells specify the cells
#' @param group.by modify the group.by, when length(group.by) >1,
#' instead of adding additional heatmap, produce additonal annotation bar(s) on top of heatmap, 
#' @param group1.colors hex color character vector matching to group.by[1]
#' @param group2.colors hex color character vector matching to group.by[2]
#' @param colors hex color character vector for heatmap, default is Seurat::PurpleAndYellow(), or ggsci::pal_gsea()(12)
#' @param position gene name position
#' @param nrow pass to wrap_plots
#' @param ncol pass to wrap_plots
#' @param design pass to wrap_plots
DoHeatmap.2 <- function(object, dge_markers = NULL,features = NULL, cells = NULL,
                        no.legend =F,unique.name= T, 
                        group.by = c("X4clusters","orig.ident"),
                        group.bar = TRUE, group1.colors = Singler.colors, 
                        group2.colors=Singler.colors, colors = Seurat::PurpleAndYellow(),
                        disp.min = -2.5, disp.max = NULL, slot = "scale.data",
                        assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45,
                        raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
                        combine = TRUE,title = "",title.size = 14,do.print = FALSE,
                        position = "right",save.path = NULL,file.name = NULL,
                        cex.row=12,legend.size = 12,
                        nrow = 5, ncol = 6, design = c(patchwork::area(1, 1, 5, 5),
                                                       patchwork::area(3, 6, 3, 6)),
                        units="in", width=10, height=7,res=600,...){
        if(is.null(file.name)){
                v <- UniqueName(object = object, fileName = deparse(substitute(object)), 
                                unique.name = T)
                v = paste0(v,"_",FindIdentLabel(object))
                if(!no.legend) v = paste0(v, "_Legend")
                file.name = paste0("Heatmap2_top",v,".jpeg")
        }
        
        if(class(title) != "character") stop("Title is incorrect")
        if (!is.null(x = dge_markers)) {
                
                colnames(dge_markers)[grep("cluster*.",colnames(dge_markers))[1]] = "cluster"
                top <-  dge_markers %>%
                        group_by(cluster) %>%
                        top_n(Top_n, avg_logFC)
                features = c(base::as.character(top$gene),features)
        }
        
        for (i in 1:length(x = group.by)) {
                len = length(unique(object@meta.data[,group.by[i]]))
                assign("group.colors.len",length(get(paste0("group",i,".colors"))))
                if(len <= group.colors.len) {
                        assign(paste0("group",i,".colors"),
                               get(paste0("group",i,".colors"))[1:len])
                } else stop(paste("Need",len, "colors for",group.by[i]))
        }
                
        cells <- cells %||% colnames(x = object)
        if (is.numeric(x = cells)) {
                cells <- colnames(x = object)[cells]
        }
        assay <- assay %||% DefaultAssay(object = object)
        DefaultAssay(object = object) <- assay
        features <- features %||% VariableFeatures(object = object)
        features <- rev(x = unique(x = features))
        disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                         yes = 2.5, no = 6)
        possible.features <- rownames(x = GetAssayData(object = object, 
                                                       slot = slot))
        if (any(!features %in% possible.features)) {
                bad.features <- features[!features %in% possible.features]
                features <- features[features %in% possible.features]
                if (length(x = features) == 0) {
                        stop("No requested features found in the ", slot, 
                             " slot for the ", assay, " assay.")
                }
                warning("The following features were omitted as they were not found in the ", 
                        slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                                         collapse = ", "))
        }
        data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                                   slot = slot)[features, cells, drop = FALSE])))
        object <- suppressMessages(expr = StashIdent(object = object, 
                                                     save.name = "ident"))
        if(!all(group.by %in% colnames(object@meta.data))) {
                stop("group.by=",paste(group.by[!(group.by %in% colnames(object@meta.data))], "doesn't exsit in meta.data!"))
        }
        group.by <- group.by %||% "ident"
        groups.use <- object[[group.by]][cells, , drop = FALSE]
        
        # main heatmap
        data.group <- data
        group.use <- groups.use[, 1, drop = TRUE]
        if (!is.factor(x = group.use)) {
                group.use <- factor(x = group.use)
        }
        names(x = group.use) <- cells
        if (draw.lines) {
                lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                                0.0025)
                placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) * 
                                                           lines.width), FUN = function(x) {
                                                                   return(Seurat:::RandomName(length = 20))
                                                           })
                placeholder.groups <- rep(x = levels(x = group.use), 
                                          times = lines.width)
                group.levels <- levels(x = group.use)
                names(x = placeholder.groups) <- placeholder.cells
                group.use <- as.vector(x = group.use)
                names(x = group.use) <- cells
                group.use <- factor(x = c(group.use, placeholder.groups), 
                                    levels = group.levels)
                na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                                        ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                                     colnames(x = data.group)))
                data.group <- rbind(data.group, na.data.group)
        }
        lgroup <- length(levels(group.use))
        plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                         colors = colors,
                                         disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                         cell.order = names(x = sort(x = group.use)), group.by = group.use)
        if(!no.legend) {
                plot = plot + guides(color = guide_legend(override.aes = list(shape = 15,alpha=1,size = legend.size*0.8)))+
                        theme(legend.text = element_text(size = legend.size),
                              legend.title = element_text(size = legend.size*1.2))
                }
        plot <- plot + theme(line = element_blank())
        if(!is.null(title)) {
                plot = plot+ ggtitle(title)+
                        theme(plot.title = element_text(size=title.size, hjust = 0.5,face="plain"))
        }
        plot = plot + scale_y_discrete(position = position)
        plot = plot + theme(axis.text.y = element_text(size = cex.row,colour = "black"))

        if(no.legend) plot = plot + NoLegend()
        if (group.bar) {
                default.colors <- c(scales::hue_pal()(length(x = levels(x = group.use))))
                if (!is.null(x = names(x = group1.colors))) {
                        cols <- unname(obj = group1.colors[levels(x = group.use)])
                } else {
                        cols <- group1.colors[1:length(x = levels(x = group.use))] %||% 
                                default.colors
                }
                if (any(is.na(x = cols))) {
                        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
                        cols <- Seurat:::Col2Hex(cols)
                        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                                        start = 1, stop = 7)))))
                        through <- length(x = default.colors)
                        while (length(x = col.dups) > 0) {
                                pal.max <- length(x = col.dups) + through
                                cols.extra <- hue_pal()(pal.max)[(through + 
                                                                          1):pal.max]
                                cols[col.dups] <- cols.extra
                                col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                                                start = 1, stop = 7)))))
                        }
                }
                group.use2 <- sort(x = group.use)
                if (draw.lines) {
                        na.group <- Seurat:::RandomName(length = 20)
                        levels(x = group.use2) <- c(levels(x = group.use2), 
                                                    na.group)
                        group.use2[placeholder.cells] <- na.group
                        cols <- c(cols, "#FFFFFF")
                }
                names(x = cols) <- levels(x = group.use2)
                
                # extract coordicates of ggplot2
                #' @export x.min,x.max,y.min,y.max, coordinate values anti-clockwise from left bottom
                Extract_coord <- function(g){
                        pbuild <- ggplot_build(plot = g)
                        x.min <- min(pbuild$layout$panel_params[[1]]$x.range)
                        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                        y.min <- min(pbuild$layout$panel_params[[1]]$y.range)
                        y.max <- max(pbuild$layout$panel_params[[1]]$y.range)
                        
                        return(list("x.min" = x.min, "x.max" = x.max, "y.min" = y.min,"y.max" = y.max))
                }
                coord <- Extract_coord(plot)
                y.range = coord$y.max - coord$y.min
                y.pos = coord$y.max + y.range * 0.015
                coord$y.max = y.pos + group.bar.height * y.range
                coord$x.min <- coord$x.min + 0.1
                coord$x.max <- coord$x.max - 0.1
                
                plot <- plot + 
                        annotation_raster(raster = t(x = cols[group.use2]), 
                                          xmin = coord$x.min, xmax = coord$x.max, 
                                          ymin = y.pos, ymax = coord$y.max) + 
                        coord_cartesian(ylim = c(0, coord$y.max), clip = "off") + 
                        scale_color_manual(values = cols)
                
                for (i in 2:ncol(x = groups.use)) {
                        group.use3 = groups.use[names(group.use2),i]
                        if (!is.factor(x = group.use3)) group.use3 %<>% factor()
                        group.use3 %<>% droplevels()
                        if(length(levels(group.use3)) > length(group2.colors)) stop("Not enough color!")
                        group3_colors = group2.colors[1:length(unique(group.use3))]
                        if (draw.lines) group3_colors = c(group3_colors[1:(length(group3_colors)-1)],
                                                          "#FFFFFF")
                        names(x = group3_colors) <- levels(x = group.use3)
                        coord <- Extract_coord(plot)
                        
                        bar.ymin = coord$y.max+ y.range * 0.01
                        bar.ymax = bar.ymin + group.bar.height * y.range
                        
                        plot <- plot + 
                                annotation_raster(raster = t(x = group3_colors[group.use3]), 
                                                  xmin = coord$x.min, xmax = coord$x.max, 
                                                  ymin = bar.ymin, 
                                                  ymax = bar.ymax)+
                                coord_cartesian(ylim = c(0, bar.ymax), clip = "off")
                        if(!no.legend){
                                df <- as.data.frame(table(group.use3))
                                df$color = group3_colors[df$group.use3]
                                colnames(df)[1] = group.by[i]
                                legend <- ggplot(df,aes_string(x = "Freq", fill = group.by[i]))+ 
                                        geom_bar() + scale_fill_manual(values=df$color)+ 
                                        theme(legend.text = element_text(size = legend.size),
                                              legend.title = element_text(size = legend.size*1.2))
                                legend <- cowplot::get_legend(legend)
                                legend <- ggpubr::as_ggplot(legend)

                                plot <- patchwork::wrap_plots(plot,legend, ncol = ncol, nrow = nrow,design = design)
                        }
                }

                if (label) {
                        y.max = coord$y.max
                        x.max = coord$x.max
                        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% 
                                attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(), 
                                     which = "pos")
                        x <- data.frame(group = sort(x = group.use), 
                                        x = x.divs)
                        label.x.pos <- tapply(X = x$x, INDEX = x$group, 
                                              FUN = function(y) {
                                                      if (isTRUE(x = draw.lines)) {
                                                              mean(x = y[-length(x = y)])
                                                      }
                                                      else {
                                                              mean(x = y)
                                                      }
                                              })
                        label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                                  label.x.pos)
                        plot <- plot + geom_text(stat = "identity", 
                                                 data = label.x.pos, aes_string(label = "group", 
                                                                                x = "label.x.pos"), y = y.max + y.max * 
                                                         0.03 * 0.5, angle = angle, hjust = hjust, 
                                                 size = size)
                        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * 
                                                                                         size), clip = "off"))
                }
        }
        if(do.print){
                if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
                if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
                jpeg(paste0(save.path,"/", file.name),units=units, width=width, height=height,res=res)
                print(plot)
                dev.off()
        } else return(plot)
}


#' DoHeatmap.matrix, generates Seurat like heatmap using matrix
#' @param data.use expression data.frame
#' @param group.by factor/character vector of column group labels
#' @param group1.colors hex color character vector matching to group.by[1]
#' @param group2.colors hex color character vector matching to group.by[2]
#' @param colors hex color character vector for heatmap, default is Seurat::PurpleAndYellow(), or ggsci::pal_gsea()(12)
#' @param save.path folder to save
DoHeatmap.matrix <- function (data.use, features = NULL, cells = NULL, 
                              no.legend =F, group.by = "ident", group.bar = TRUE, 
                              group.colors = NULL,disp.min = -2.5, disp.max = NULL, slot = "scale.data", 
                              assay = NULL, label = TRUE, colors = Seurat::PurpleAndYellow(), size = 5.5, hjust = 0, angle = 45, 
                              raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, 
                              combine = TRUE,title = "",title.size = 14,do.print = FALSE,
                              pal_gsea = TRUE,position = "right",save.path = NULL, file.name = NULL,
                              unique.name= T, cex.row=12,legend.size = NULL,units="in", 
                              width=10, height=7,res=600,...) 
{
    if(is.null(file.name)){
        v <- UniqueName(object = object, fileName = deparse(substitute(object)), unique.name = unique.name)
        v = paste0(v,"_",FindIdentLabel(object))
        if(!no.legend) v = paste0(v, "_Legend")
        file.name = paste0("Heatmap_top",Top_n,"_",v,".jpeg")
    }
    if(class(title) != "character") stop("Title is incorrect")
    cells <- cells %||% colnames(x = data.use)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = data.use)[cells]
    }
    features <- rev(x = unique(x = features))
    disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                     yes = 2.5, no = 6)
    possible.features <- rownames(x = data.use)
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if (length(x = features) == 0) {
            stop("No requested features found in the ", slot, 
                 " slot for the ", assay, " assay.")
        }
        warning("The following features were omitted as they were not found in the ", 
                slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                                 collapse = ", "))
    }
    data <- as.data.frame(x = as.matrix(x = t(x = data.use[features, cells, drop = FALSE])))
    groups.use <- data.frame("ident" = group.by, row.names = cells)
    plots <- vector(mode = "list", length = ncol(x = groups.use))
    for (i in 1:ncol(x = groups.use)) {
        data.group <- data
        group.use <- groups.use[, i, drop = TRUE]
        if (!is.factor(x = group.use)) {
            group.use <- factor(x = group.use)
        }
        names(x = group.use) <- cells
        if (draw.lines) {
            lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                        0.0025)
            placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) * 
                                                   lines.width), FUN = function(x) {
                                                       return(Seurat:::RandomName(length = 20))
                                                   })
            placeholder.groups <- rep(x = levels(x = group.use), 
                                      times = lines.width)
            group.levels <- levels(x = group.use)
            names(x = placeholder.groups) <- placeholder.cells
            group.use <- as.vector(x = group.use)
            names(x = group.use) <- cells
            group.use <- factor(x = c(group.use, placeholder.groups), 
                                levels = group.levels)
            na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                                    ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                                 colnames(x = data.group)))
            data.group <- rbind(data.group, na.data.group)
        }
        lgroup <- length(levels(group.use))
        plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                cell.order = names(x = sort(x = group.use)), group.by = group.use)
        if (group.bar) {
            default.colors <- c(scales::hue_pal()(length(x = levels(x = group.use))))
            cols <- group.colors[1:length(x = levels(x = group.use))] %||% 
                default.colors
            if (any(is.na(x = cols))) {
                cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
                cols <- Seurat:::Col2Hex(cols)
                col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                                start = 1, stop = 7)))))
                through <- length(x = default.colors)
                while (length(x = col.dups) > 0) {
                    pal.max <- length(x = col.dups) + through
                    cols.extra <- scales::hue_pal()(pal.max)[(through + 
                                                          1):pal.max]
                    cols[col.dups] <- cols.extra
                    col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                                    start = 1, stop = 7)))))
                }
            }
            group.use2 <- sort(x = group.use)
            if (draw.lines) {
                na.group <- Seurat:::RandomName(length = 20)
                levels(x = group.use2) <- c(levels(x = group.use2), 
                                            na.group)
                group.use2[placeholder.cells] <- na.group
                cols <- c(cols, "#FFFFFF")
            }
            pbuild <- ggplot_build(plot = plot)
            names(x = cols) <- levels(x = group.use2)
            y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
            y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 
                y.range * 0.015
            y.max <- y.pos + group.bar.height * y.range
            plot <- plot + annotation_raster(raster = t(x = cols[group.use2]), 
                                             xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                coord_cartesian(ylim = c(0, y.max), clip = "off") + 
                scale_color_manual(values = cols)
            if (label) {
                x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% 
                    pbuild$layout$panel_params[[1]]$x$break_positions()
                x <- data.frame(group = sort(x = group.use), 
                                x = x.divs)
                label.x.pos <- tapply(X = x$x, INDEX = x$group, 
                                      FUN = median) * x.max
                label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                          label.x.pos)
                plot <- plot + geom_text(stat = "identity", 
                                         data = label.x.pos, aes_string(label = "group", 
                                                                        x = "label.x.pos"), y = y.max + y.max * 
                                             0.03 * 0.5, angle = angle, hjust = hjust, 
                                         size = size)
                plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                         y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * 
                                                                             size), clip = "off"))
            }
        }
        plot <- plot + theme(line = element_blank())
        plots[[i]] <- plot
    }
    if (combine) {
        plots <- patchwork::wrap_plots(plots)
    }
    plots = plots + scale_y_discrete(position = position)
    if(!is.null(colors)) plots = plots + scale_fill_gradientn(colors = colors)
    if(!is.null(title)) {
        plots = plots+ ggtitle(title)+ 
            theme(plot.title = element_text(size=title.size, hjust = 0.5,face="plain"))
    }
    plots = plots + theme(axis.text.y = element_text(size = cex.row))
    if(!is.null(legend.size)) {
        plots = plots + theme(legend.text = element_text(size = legend.size),
                                  legend.title = element_text(size = legend.size*2),
                                  legend.key.size = unit(legend.size/10,"cm"))
    }
    if(no.legend) plots = plots + NoLegend()
    if(do.print){
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path,"/",file.name), units=units, width=width, height=height,res=res)
        print(plots)
        dev.off()
    } else return(plots)
}

#' wrapper for DotPlot from Seurat
#' Allow log expression

#' @param log.data = log2, 
#' @param cluster.features = FALSE, cluster features
#' @param exp.min = NA, Set minimal expression level, 
#' @param exp.max = NA, Set maximal expression level 
#' @param n.breaks = NULL, number of breaks for expression color bar
#' @examples
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features = cd_genes)
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' DotPlot(object = pbmc_small, features = cd_genes, split.by = 'groups')
#'
DotPlot.1 <- function(
        object,
        assay = NULL,
        features,
        log.data = log2,
        cols = c("lightgrey", "blue"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        idents = NULL,
        group.by = NULL,
        split.by = NULL,
        cluster.idents = FALSE,
        cluster.features = FALSE,
        scale = TRUE,
        scale.by = 'radius',
        scale.min = NA,
        scale.max = NA,
        exp.min = NA,
        exp.max = NA,
        n.breaks = NULL
) {
        assay <- assay %||% DefaultAssay(object = object)
        DefaultAssay(object = object) <- assay
        split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
        scale.func <- switch(
                EXPR = scale.by,
                'size' = scale_size,
                'radius' = scale_radius,
                stop("'scale.by' must be either 'size' or 'radius'")
        )
        feature.groups <- NULL
        if (is.list(features) | any(!is.na(names(features)))) {
                feature.groups <- unlist(x = sapply(
                        X = 1:length(features),
                        FUN = function(x) {
                                return(rep(x = names(x = features)[x], each = length(features[[x]])))
                        }
                ))
                if (any(is.na(x = feature.groups))) {
                        warning(
                                "Some feature groups are unnamed.",
                                call. = FALSE,
                                immediate. = TRUE
                        )
                }
                features <- unlist(x = features)
                names(x = feature.groups) <- features
        }
        cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
        data.features <- FetchData(object = object, vars = features, cells = cells)
        data.features$id <- if (is.null(x = group.by)) {
                Idents(object = object)[cells, drop = TRUE]
        } else {
                object[[group.by, drop = TRUE]][cells, drop = TRUE]
        }
        if (!is.factor(x = data.features$id)) {
                data.features$id <- factor(x = data.features$id)
        }
        id.levels <- levels(x = data.features$id)
        data.features$id <- as.vector(x = data.features$id)
        if (!is.null(x = split.by)) {
                splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
                if (split.colors) {
                        if (length(x = unique(x = splits)) > length(x = cols)) {
                                stop("Not enough colors for the number of groups")
                        }
                        cols <- cols[1:length(x = unique(x = splits))]
                        names(x = cols) <- unique(x = splits)
                }
                data.features$id <- paste(data.features$id, splits, sep = '_')
                unique.splits <- unique(x = splits)
                id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
        }

        data.plot <- lapply(
                X = unique(x = data.features$id),
                FUN = function(ident) {
                        data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
                        avg.exp <- apply(
                                X = data.use,
                                MARGIN = 2,
                                FUN = function(x) {
                                        return(mean(x = expm1(x = x)))
                                }
                        )
                        PercentAbove <- function(x, threshold) {
                                return(length(x = x[x > threshold]) / length(x = x))
                        }
                        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
                        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
                }
        )
        names(x = data.plot) <- unique(x = data.features$id)
        if (cluster.idents) {
                mat <- do.call(
                        what = rbind,
                        args = lapply(X = data.plot, FUN = unlist)
                )
                mat <- scale(x = mat)
                id.levels <- id.levels[hclust(d = dist(x = mat))$order]
        }
        if (cluster.features) {
                mat <- do.call(
                        what = rbind,
                        args = lapply(X = data.plot, FUN = unlist)
                )
                mat <- scale(x = mat)
                features <- features[hclust(d = dist(x = t(mat[,1:(ncol(mat)/2)])))$order]
        }
        data.plot <- lapply(
                X = names(x = data.plot),
                FUN = function(x) {
                        data.use <- as.data.frame(x = data.plot[[x]])
                        data.use$features.plot <- rownames(x = data.use)
                        data.use$id <- x
                        return(data.use)
                }
        )
        data.plot <- do.call(what = 'rbind', args = data.plot)
        #if(is.function(log.data)) data.plot$avg.exp = log.data(data.plot$avg.exp+1)
        if (!is.null(x = id.levels)) {
                data.plot$id <- factor(x = data.plot$id, levels = id.levels)
        }
        if (length(x = levels(x = data.plot$id)) == 1) {
                scale <- FALSE
                warning(
                        "Only one identity present, the expression values will be not scaled",
                        call. = FALSE,
                        immediate. = TRUE
                )
        }
        avg.exp.scaled <- sapply(
                X = unique(x = data.plot$features.plot),
                FUN = function(x) {
                        data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
                        if (scale) {
                                data.use <- scale(x = data.use)
                                data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                        } else if(is.function(log.data)){
                                data.use <- log.data(x = data.use+1)
                        } else stop("No scale and log.data is not function")
                        return(data.use)
                }
        )
        avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
        if (split.colors) {
                avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
        }
        data.plot$avg.exp.scaled <- avg.exp.scaled
        data.plot$features.plot <- factor(
                x = data.plot$features.plot,
                levels = features
        )
        data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
        data.plot$pct.exp <- data.plot$pct.exp * 100
        if (split.colors) {
                splits.use <- vapply(
                        X = as.character(x = data.plot$id),
                        FUN = gsub,
                        FUN.VALUE = character(length = 1L),
                        pattern =  paste0(
                                '^((',
                                paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
                                ')_)'
                        ),
                        replacement = '',
                        USE.NAMES = FALSE
                )
                data.plot$colors <- mapply(
                        FUN = function(color, value) {
                                return(colorRampPalette(colors = c('grey', color))(20)[value])
                        },
                        color = cols[splits.use],
                        value = avg.exp.scaled
                )
        }
        color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
        color.by <- ifelse(test = is.function(log.data), yes = 'avg.exp.scaled', no = color.by)
        
        if (!is.na(x = scale.min)) {
                data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
        }
        if (!is.na(x = scale.max)) {
                data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
        }
        if (!is.na(x = exp.min)) {
                data.plot[data.plot[,color.by] < exp.min, color.by] <- exp.min
        }
        if (!is.na(x = exp.max)) {
                data.plot[data.plot[,color.by] > exp.max, color.by] <- exp.max
        }
        if (!is.null(x = feature.groups)) {
                data.plot$feature.groups <- factor(
                        x = feature.groups[data.plot$features.plot],
                        levels = unique(x = feature.groups)
                )
        }
        plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
                geom_point(mapping = aes_string(size = 'pct.exp', fill = color.by),
                           color = "black", pch=21) +
                scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                guides(size = guide_legend(title = 'Percent Expressed')) +
                labs(
                        x = 'Features',
                        y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
                ) +
                cowplot::theme_cowplot()
        if (!is.null(x = feature.groups)) {
                plot <- plot + facet_grid(
                        facets = ~feature.groups,
                        scales = "free_x",
                        space = "free_x",
                        switch = "y"
                ) + theme(
                        panel.spacing = unit(x = 1, units = "lines"),
                        strip.background = element_blank()
                )
        }
        if (split.colors) {
                plot <- plot + scale_color_identity()
        } else if (length(x = cols) == 1) {
                plot <- plot + scale_color_distiller(palette = cols)
        } else if (length(x = cols) == 2){
                plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
        } else {
                plot <- plot + scale_fill_gradientn(
                        colours=cols,
                        n.breaks = n.breaks,
                        name = ifelse(test = is.function(log.data),
                                      yes = expression(atop("Mean non-zero\n expression",log[2](UMI+1))),
                                      no = expression(atop("Mean non-zero\n expression",(UMI)))), 
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "fill"
                        )
        }
        if (!split.colors) {
                plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
        }
        return(plot)
}


#' Extend ElbowPlot function to select best resolution/cluster numbers
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
#' @param save.path path to save jpeg file when do.print = T. Default path is under output folder. 
ElbowPlot.1 <- function (object, graph.name = "integrated_snn",reductions = "umap",
                         check.by = c("Resolutions","Cluster.numbers"),
                         title = "",
                         unique.name= T, units="in", width=10, height=7,res=600,
                         do.print = F, do.return = T, save.path = NULL,file.name =NULL,...) 
{
    if (is.null(file.name)) {
        file.name <- UniqueName(object,fileName = NULL,unique.name = unique.name)
        file.name = paste0(file.name,"_",check.by[1])
        file.name = paste0("ElbowPlot_",VarName,".jpeg")
        }
    clustering.name <- grep(graph.name, colnames(object@meta.data), value = T)
    if (!any(clustering.name %in% colnames(object@meta.data))) {
        stop("Provided graph.name not present in Seurat object")
    }
    res = gsub(paste0(graph.name,"_res."),"",clustering.name) %>% as.numeric() %>% sort
    cell.embeddings = object[[reductions]]@cell.embeddings
    k.num <- c()
    total.dist <- c()
    for(i in seq_along(res)){
        cluster_k = unique(object[[paste0(graph.name,"_res.",res[i])]])[,1]
        k.num[i] = length(cluster_k)
        # calculating the total distance from center to each data points
        Cluster_center <- function(df){
            cnt = colMeans(df)
            Dist = apply(df,1,function(x,cnt) {(sqrt((x[1] - cnt[1])^2+(x[2]-cnt[2])^2))},cnt)
            return(sum(Dist))
        }
        total.dist[i] <- sapply(cluster_k, function(k){
                            cells = object[[paste0(graph.name,"_res.",res[i])]][,1] %in% k
                            Cluster_center(cell.embeddings[cells,])
                            }) %>% sum
        Progress(i, length(res))
    }

    if(check.by[1] == "Resolutions") data.use = data.frame(Resolutions = res, total.dist = total.dist)
    if(check.by[1] == "Cluster.numbers") data.use = data.frame(Cluster.numbers = k.num, total.dist = total.dist)
    plot <- ggscatter(data = data.use,
                      x = check.by[1], y = "total.dist",
                      ylab = "Total distance",
                      title = title) + TitleCenter()
    
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path,"/",file.name), units=units, width=width, height=height,res=600)
        print(plot)
        dev.off()
    }
    if(do.return & Sys.info()[['sysname']] != "Linux") return(plot)
}

#' ExtractMetaColor: extract color code from meta.data
#' @param object seurat object, meta.data slot must have "color"
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' @export cell.colors color code
#' @example 
# ExtractMetaColor(object = MCL, group.by = "singler1sub")
ExtractMetaColor <- function(object, group.by = NULL){
    meta.data = object@meta.data
    group.by <- group.by %||% FindIdentLabel(object)
    color_index <- paste0(group.by ,".colors")
    if(!any(color_index %in% colnames(meta.data))) {
        return(gg_color_hue(length(unique(as.character(object@meta.data[,group.by])))))
        } 
    else {
            label = sub("\\.colors","",color_index)
            meta.data = meta.data[,c(label,color_index)]
            #meta.data$index <- as.numeric(as.factor(meta.data[,1]))
            df_colors = meta.data[!duplicated(meta.data[,group.by]),]
            df_colors = df_colors[order(df_colors[,group.by]),]
        }
    if(all(is.na(df_colors[,color_index]))) df_colors[,color_index] = gg_color_hue(nrow(df_colors))
    cell.colors = base::as.character(df_colors[,color_index])
    names(cell.colors) = df_colors[,group.by]
    cell.colors = cell.colors[order(names(cell.colors))]
    return(cell.colors)
}


#' use Findmark results to generate eulerr data frame for venn diagram
#' @param df Seruat::FindMarkers results.
#' @param shape shape of venn diagram
#' @param key Names of venn diagram catergroy, which must be the sub components of df$cluster.
#'  defaul NULL means all components of df$cluster.
#' @param cut_off choose within c("p_val","p_val_adj","avg_logFC")
#' @param cut_off_value corrsponding cut off value.
#' @param eulerr::euler table
#' @param do.legend TRUE/FALSE
#' @param do.return TRUE/FALSE return plot
#' @param return.raw TRUE/FALSE return pos.share_genes list
#' @param do.print print figures
#' @param save.path path to save
#' @export g plot object
#' @export pos_genes positive shared gene list
#' @example eulerr(T_cells_markers,shape =  "ellipse",cut_off = "avg_logFC", cut_off_value = 0.01)
eulerr <- function(df, key = NULL, cut_off = "avg_logFC",cut_off_value = 0.05, 
                   do.lenged = TRUE,do.return = TRUE, return.raw = FALSE,
                   do.print = FALSE,save.path = NULL,file.name =NULL,
                   shape = c("circle", "ellipse"),units="in", width=width, height=height,res=res,...){
        df$cluster <- as.vector(df$cluster)
        df$gene <- as.vector(df$gene)
        if(!is.null(key)) df <- df[(df$cluster %in% key),]
        df_list <- split(df,df$cluster)
        
        if(cut_off == "avg_logFC"){
            pos_genes <- sapply(df_list, function(df) df[(df$avg_logFC > -cut_off_value),"gene"])
        }  
        if(any(cut_off %in% c("p_val","p_val_adj"))){
                pos_genes1 <- lapply(df_list, function(df) df[(df$avg_logFC > 0),"gene"])
                shared_genes <- lapply(df_list, function(df) df[(abs(df[,cut_off]) > cut_off_value),"gene"])
                pos_genes <- mapply(function(x,y) unique(c(x,y)), pos_genes1, shared_genes)
        }
        euler_df <- eulerr::euler(pos_genes,shape = shape,...)
        
        g <- plot(euler_df, quantities = TRUE, lty = 1:6,
                  legend = do.lenged, main = paste(cut_off," : ",cut_off_value))
        if(do.print) {
            if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
            if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
            if(is.null(file.name)) file.name = paste0("Venn_",cut_off,"_",cut_off_value,".jpeg")
            jpeg(paste0(save.path,"/",file.name), units=units, width=10, height=7,res=600)
            print(g)
            dev.off()
        }
        if(return.raw) return(pos_genes)
        if(do.return & Sys.info()[['sysname']] != "Linux") return(g)
        
}


#' modify FeaturePlot 
#' @param cut_off_value corrsponding cut off value.
#' @param no.AxesLabel NoAxesLabel
#' @param ChangeColorScale ChangeColorScale
#' @export save.path folder to save 
FeaturePlot.1 <- function (object, features, dims = c(1, 2), cells = NULL, 
                           cols = ifelse(test = c(blend,blend), yes = c("#ff0000", "#00ff00"), no = c("lightgrey", "red")),
                           pt.size = NULL, order = FALSE, min.cutoff = NA, 
                           max.cutoff = NA, reduction = NULL, split.by = NULL, shape.by = NULL, 
                           slot = "data", threshold = 0.001, blend = FALSE, blend.threshold = 0.5, label = FALSE, alpha =1,
                           label.size = 4,text.size = 15, repel = FALSE, ncol = NULL, combine = TRUE, border = FALSE,
                           coord.fixed = FALSE, by.col = TRUE,do.print = F,do.return =T,
                           title = NULL,unique.name = F,legend.title = NULL,no.AxesLabel = F,
                           no.legend = F,units="in", width=10, height=7,res = 600,
                           save.path = NULL, file.name=NULL) {
    if(is.null(file.name)){
        file.name <- UniqueName(object,fileName = deparse(substitute(object)),unique.name = unique.name)
        file.name = paste0(file.name,"_",FindIdentLabel(object),
                         ifelse(!is.null(split.by), yes = paste0("_",split.by), no =""))
    }
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                      axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",size = 14, margin = margin(r = 7)))
    if (is.null(reduction)) {
        default_order <- c("umap", "tsne", "pca")
        reducs <- which(default_order %in% names(object@reductions))
        reduction <- default_order[reducs[1]]
    }
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    if (blend && length(x = features) != 2) {
        stop("Blending feature plots only works with two features")
    }
    if (blend) {
        if (length(x = cols) > 2) {
            warning("Blending feature plots only works with two colors; using first two colors", 
                    call. = FALSE, immediate. = TRUE)
        } else if (length(x = cols) < 2) {
            warning("Blended feature plots require two colors, using default colors", 
                    call. = FALSE, immediate. = TRUE)
            cols <- c("#ff0000", "#00ff00")
        }
    }
    if (blend && length(x = cols) != 2) {
        stop("Blending feature plots only works with two colors")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", features), 
                      cells = cells, slot = slot)
    if(length(features)>1 ) {
        Minimal <- apply(data[,features],2,min)
        data[,features] = sweep(data[,features], 2, Minimal,"-")
        range <- apply(data[,features],2,max) - apply(data[,features],2,min)
        data[,features] = sweep(data[,features], 2, range,"/")*2.5
    }

    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", 
             paste(features, collapse = ", "), " in slot ", slot, 
             call. = FALSE)
    } else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[,feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[,feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                                max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
                                                                              ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                       FUN = function(index) {
                                           data.feature <- as.vector(x = data[, index])
                                           min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 
                                                                                                   3], data.feature)
                                           max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 
                                                                                                   3], data.feature)
                                           data.feature[data.feature < min.use] <- min.use
                                           data.feature[data.feature > max.use] <- max.use
                                           if (brewer.gran == 2) {
                                               return(data.feature)
                                           }
                                           data.cut <- if (all(data.feature == 0)) {
                                               0
                                           } else {
                                               as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                                breaks = brewer.gran)))
                                           }
                                           return(data.cut)
                                       })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split.by <- if (is.null(x = split.by)) {
        Seurat:::RandomName()
    } else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells], 
               object[[split.by, drop = TRUE]][cells])
    }
    if (!is.factor(x = data$split.by)) {
        data$split.by <- factor(x = data$split.by)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    data <- reshape2::melt(data, id.vars = c(dims,"ident","split.by"),
                           variable.name = "features", 
                           value.name = "feature_values")
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
    if (blend) {
        ncol <- 4
        color.matrix <- BlendMatrix(two.colors = cols, col.threshold = blend.threshold)
        colors <- list(color.matrix[, 1], color.matrix[1, ], 
                       as.vector(x = color.matrix))
    }
    
    if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, 
                                                              feature])) + 1
        cols.use <- colors[[1]][sort(x = unique(x = cols.use))]
    } else {
        cols.use <- NULL
    }
    
    plot <- Seurat:::SingleDimPlot(data = data, 
                            dims = dims, col.by = "feature_values", 
                            order = order, pt.size = pt.size, cols = cols.use,
                            shape.by = shape.by, label = FALSE) + scale_x_continuous(limits = xlims) + 
        scale_y_continuous(limits = ylims) +  cowplot::theme_cowplot()
    if (label) {
        plot <- LabelClusters(plot = plot, id = "ident", 
                              repel = repel, size = label.size)
    }
    if (!is.null(x = split.by) & is.null(ncol)) {
        plot <- plot + facet_grid(facets = features ~ split.by)
    }
    if (!is.null(x = split.by) & !is.null(ncol) & length(features) > 1) {
        plot <- plot + facet_wrap(facets = features ~ split.by, ncol = ncol)
    }
    if (!is.null(x = split.by) & !is.null(ncol) & length(features) == 1) {
        plot <- plot + facet_wrap(facets = . ~ split.by, ncol = ncol)
    }
    if (is.null(x = split.by) & length(features) > 1) {
        plot <- plot + facet_wrap(facets = features ~.,ncol = ncol)
    }
    plot <- plot + theme(strip.background = element_blank(),
              strip.text =  element_text(face="plain",size=text.size),
              axis.title = element_text(size=text.size),
              axis.text = element_text(size=text.size))
    if(border) plot = plot + theme(panel.border = element_rect(size = 0.5,colour = "grey"),
                                   axis.line = element_blank())
    if(no.AxesLabel) plot =  plot + NoAxesLabel()
    if (!blend) {
        plot <- suppressMessages(ChangeColorScale(plot, gradient.use = cols, alpha.use = alpha,
                             scaled.expression.threshold = threshold))
    }
    if (coord.fixed) {
        plot <- plot + coord_fixed()
    }
    if(no.legend) {
        plot = plot + NoLegend()
        L = ""
    } else L = "_Legend"
    if(!is.null(title)) {
        plot = plot + 
            ggtitle(title)+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    }
    if(!is.null(legend.title)) plot = plot + labs(fill = legend.title)
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path,"/",file.name), units=units, width=width, height=height,res=res)
        print(plot)
        dev.off()
    }
    if(do.return & Sys.info()[['sysname']] != "Linux") return(plot)
}


# FeaturePlot in Seurat 3 doesn't support two gene overlay.
# modify FeaturePlot in Seurat 2 and support Seurat 3 object
FeaturePlot.2 <- function (object, features, min.cutoff = NA, max.cutoff = NA,
                           dims = c(1,2), cells = NULL, pt.size = NULL,slot = "data",
                           reduction = "tsne", split.by = NULL, shape.by = NULL, alpha = 1,
                           cols.use = c("yellow",  "red"), pch.use = 16, overlay = TRUE, do.hover = FALSE,
                           data.hover = "ident", do.identify = FALSE, breaks =2,
                           use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE,
                           dark.theme = FALSE, do.return = TRUE, vector.friendly = FALSE,
                           unique.name = F, do.print = FALSE, file.name = NULL,
                           units="in", width=10, height=7,res=600,...)
{
    if(is.null(file.name)){
        VarName <- UniqueName(object,fileName = deparse(substitute(object)),unique.name = unique.name)
        VarName = paste0(VarName,"_",FindIdentLabel(object),
                         ifelse(!is.null(split.by), yes = paste0("_",split.by), no =""))
        file.name = paste0("FeaturePlot_",VarName,"_",paste(features,collapse = "-"),
        "_",reduction,"_",L,".jpeg")
    }

    cells <- cells %||% colnames(x = object)
    if (is.null(x = nCol)) {
        nCol <- 2
        if (length(x = features) == 1) {
            nCol <- 1
        }
        if (length(x = features) > 6) {
            nCol <- 3
        }
        if (length(x = features) > 9) {
            nCol <- 4
        }
    }
    num.row <- floor(x = length(x = features)/nCol - 1e-05) +1
    if (overlay | do.hover) {
        num.row <- 1
        nCol <- 1
    }
    par(mfrow = c(num.row, nCol))
    dims <- paste0(Key(object = object[[reduction]]), dims)
    data.plot <- FetchData(object = object, vars = dims, cells = cells, slot = slot)
    x1 <- dims[1]
    x2 <- dims[2]
    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]

    data.use <- t(x = FetchData(object = object, vars = features, cells = cells, slot = slot))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data.use[feature,
        ]), no = cutoff)
    }, cutoff = min.cutoff, feature = features)
    
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data.use[feature,
        ]), no = cutoff)
    }, cutoff = max.cutoff, feature = features)
    check_lengths = unique(x = vapply(X = list(features,
    min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check_lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    if (overlay) {
        pList <- list(BlendPlot(data.use = data.use, features = features,
        data.plot = data.plot, pt.size = pt.size, pch.use = pch.use,alpha = alpha,
        cols.use = cols.use, dim.codes = dims, min.cutoff = min.cutoff,
        max.cutoff = max.cutoff, no.axes = no.axes, no.legend = no.legend,
        dark.theme = dark.theme, breaks = breaks))
    }
    else {
        pList <- mapply(FUN = FeaturePlot, feature = features,
        min.cutoff = min.cutoff, max.cutoff = max.cutoff, alpha = alpha,
        MoreArgs = list(data.use = data.use, data.plot = data.plot,
        pt.size = pt.size, pch.use = pch.use, cols.use = cols.use,
        dim.codes = dims, no.axes = no.axes, no.legend = no.legend,
        dark.theme = dark.theme, vector.friendly = vector.friendly),
        SIMPLIFY = FALSE)
    }
    pList[[1]] = pList[[1]] + theme(text = element_text(size=12),
                                    plot.title = element_text(size = 16,
                                                              hjust = 0.5))
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path, "/", file.name), units=units, width=width, height=height,res=res)
        print(pList[[1]])
        dev.off()
    }
    if (do.return & Sys.info()[['sysname']] != "Linux") {
        return(pList[[1]])
    }
}

#' @param stats Seurat findAllmarker output
#' @param pathways pathway list
#' @param title sample names in title
FgseaBarplot <- function(stats=res, pathways=hallmark, nperm=1000,cluster = 1,
                         title="", pathway.name = "Hallmark", hjust=0.5, 
                         width=10, height = 7, no.legend = FALSE,
                         cut.off = c("pval","padj"),cut.off.value = 0.25,
                         do.print = TRUE, do.return = FALSE, save.path = NULL, file.name = NULL){
    
    res = stats[order(stats["avg_logFC"]),]
    geneRank = res[res$cluster == cluster,c("gene","avg_logFC")] %>% tibble::deframe
    
    fgseaRes <- fgseaMultilevel(pathways=pathways, stats=geneRank, nperm=nperm)
    print(dim(fgseaRes))
    
    if(cut.off == "pval"){
        (topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=25), pathway])
        (topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=25), pathway])
    }
    if(cut.off == "padj"){
        (topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=25), pathway])
        (topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=25), pathway])
    }
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    out_write=fgseaRes[match(topPathways,pathway)]
    colnames(out_write)[1] <- 'Pathways'
    out_write$To_Plot_value <- -log10(out_write$pval)
    out_write$sign <- ifelse(out_write$NES >0,1,-1)
    out_write$To_Plot_value <- out_write$To_Plot_value*out_write$sign
    out_write$sign <- ifelse(out_write$NES >0,"Upregulated", "Downregulated")
    
    legend_title = ifelse(cut.off == "padj",yes = "Adjusted p value", no = "P value")
    
    p<-ggbarplot(out_write,
                 x = "Pathways",
                 y = "NES",
                 #fill = "sign",           # change fill color by mpg_level
                 color = "white",            # Set bar border colors to white
                 palette = "jco",            # jco journal color palett. see ?ggpar
                 sort.val = "asc",          # Sort the value in descending order
                 sort.by.groups = FALSE,     # Don't sort inside each group
                 ylab = 'Normalized Enrichment Score',
                 legend.title = paste(legend_title,"<",cut.off.value),
                 rotate = TRUE,
                 title = paste(pathway.name,"pathways in",sample,cluster),
                 ggtheme = theme_minimal(base_size = 15))+
        guides(fill = guide_legend(reverse = TRUE))+
        theme(text = element_text(size=12),
              plot.title = element_text(size = 16,hjust = hjust))
    if(cut.off == "pval") p = p + geom_col(aes(fill = pval < cut.off.value))
    if(cut.off == "padj") p = p + geom_col(aes(fill = padj < cut.off.value))
    if(no.legend) p = p + NoLegend()
    if(do.print){
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        if(is.null(file.name)) file.name =  paste0("Barplot_",sample,"_",cluster,"-",pathway.name,".jpeg")
        jpeg(paste0(save.path, "/", file.name), units=units, width=width, height=height,res=600)
        print(p)
        dev.off()
    }
    if(do.return & Sys.info()[['sysname']] != "Linux") return(p)

}

#' FgseaDotPlot generate Dot plot using findmarker results based on FGSEA
#' @param stats findmarker results
#' @param pathways pathway list
#' @param Rowv determines if and how the row dendrogram should be reordered. 
#' By default, NULL or FALSE, then no dendrogram is computed and no reordering is done.
#' If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
#' @param Colv 	determines if and how the column dendrogram should be reordered.
#' Has the options as the Rowv argument above.
#' @param title add to title names
#' @param padj padj cut off
#' @param pval pval cut off
#' @param order.yaxis.by c(1,"pval") means order y axis by pval in cluster 1
#' @param order.xaxis specify order of x axis
#' @param do.return return fgsea data frame
#' @param return.raw return fgsea raw data
#' @export save.path folder to save
#' @param ... ggballoonplot param
#' @example FgseaDotPlot(stats=res, pathways=hallmark,title = "each B_MCL clusters")
FgseaDotPlot <- function(stats=results, pathways=NULL,
                         size = "-log10(pval)", Rowv = NULL,Colv = NULL,
                         font.ytickslab = 15,
                         fill = "NES", title="", order.yaxis.by = c(1,"pval"),
                         order.xaxis = NULL,decreasing = T,
                         pathway.name = "Hallmark",padj = 0.25, pval=0.05,
                         do.return = F,return.raw = F,font.main = 18,
                         verbose=T,save.path = NULL, file.name = NULL,
                         units = "in",width=10, height=7,hjust=0.5,...){
    
    clusters = unique(as.character(stats$cluster))
    message("Calculate fgsea for each cluster.")
    fgseaRes <- list()
    for(i in seq_along(clusters)){
        geneRank = stats[stats$cluster == clusters[i],]
        geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
        fgseaRes[[i]] <- fgseaMultilevel(pathways=pathways, stats=geneRank)
        fgseaRes[[i]] = as.data.frame(fgseaRes[[i]])
        fgseaRes[[i]] = fgseaRes[[i]][,c("pathway","pval","padj","NES")]
        if(clusters[i] == order.yaxis.by[1]) {
            order.yaxis = fgseaRes[[i]][order(fgseaRes[[i]][,order.yaxis.by[2]],
                                              decreasing = decreasing), "pathway"]
        }
        if(!is.null(pval)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$pval < pval,]
        if(!is.null(padj)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$padj < padj,]
        if(nrow(fgseaRes[[i]]) > 0 ) {
            fgseaRes[[i]]$cluster = clusters[i]
        } else fgseaRes[[i]] =NULL
        Progress(i, length(clusters))
    }
    df_fgseaRes <- data.table::rbindlist(fgseaRes) %>% as.data.frame()
    if(nrow(df_fgseaRes) == 0) stop("No significant pathway! Try higher p-value!")
    df_fgseaRes = df_fgseaRes[!is.na(df_fgseaRes[, "pathway"]),]
    df_fgseaRes[," -log10(pval)"] = -log10(df_fgseaRes$pval)
    df_fgseaRes[," -log10(padj)"] = -log10(df_fgseaRes$padj)
    if(verbose) print(round(dim(df_fgseaRes)/length(clusters)))
    
    if(isTRUE(Rowv) | isTRUE(Colv)) {
        mtx_fgseaRes <- df_fgseaRes[,c("pathway","NES","cluster")]
        mtx_fgseaRes %<>% tidyr::spread(cluster,NES)
        rownames(mtx_fgseaRes) = mtx_fgseaRes[,"pathway"]
        mtx_fgseaRes = mtx_fgseaRes[,-grep("pathway",colnames(mtx_fgseaRes))]
        mtx_fgseaRes %<>% as.matrix()
        mtx_fgseaRes[is.na(mtx_fgseaRes)] = 0
    }
    if(isTRUE(Rowv)) {
        hcr <- hclust(as.dist(1-cor(t(mtx_fgseaRes), method="spearman")),
                      method="ward.D2")
        ddr <- as.dendrogram(hcr)
        rowInd <- order.dendrogram(ddr)
        order.yaxis = rownames(mtx_fgseaRes)[rowInd]
    } else {
        order.yaxis = order.yaxis[order.yaxis %in% df_fgseaRes[,"pathway"]]
    }
    if(isTRUE(Colv)) {
        hcc <- hclust(as.dist(1-cor(mtx_fgseaRes, method="spearman")),
                      method="ward.D2")
        ddc <- as.dendrogram(hcc)
        colInd <- order.dendrogram(ddc)
        order.xaxis = colnames(mtx_fgseaRes)[colInd]
    }
    df_fgseaRes[,"pathway"] %<>% as.factor
    df_fgseaRes[,"pathway"] %<>% factor(levels = order.yaxis)
    if(!is.null(order.xaxis)) {
        df_fgseaRes[,"cluster"] %<>% as.factor()
        df_fgseaRes[,"cluster"] %<>% factor(levels = order.xaxis)
        df_fgseaRes %<>% with(df_fgseaRes[order(pathway,cluster),])
    }
    # generate color pal_gsea scale based on NES range.
    rescale_colors <- function(colors = pal_gsea(), Range = range(df_fgseaRes$NES)){
        if(Range[1]>0) return(pal_gsea()(12)[7:12])
        if(Range[2]<0) return(pal_gsea()(12)[1:6])
        if(Range[1]<0 & Range[2]>0) {
            remove <- (Range[2] +Range[1]) / (Range[2] -Range[1])
            if(remove>0) return(pal_gsea()(12)[max(1,12*remove):12])
            if(remove<0) return(pal_gsea()(12)[1:(12+min(-1,12*remove)+1)])
        }
    }
    #font.ytickslab= min(font.ytickslab,round(height*300/dim(df_fgseaRes)[1]))
    plots <- ggballoonplot(df_fgseaRes, x = "cluster", y = "pathway",
                           size = size, fill = fill,
                           size.range = c(1, 5),
                           font.ytickslab= font.ytickslab,
                           title = title,
                           legend.title = ifelse(fill =="NES",
                                                 "Normalized\nenrichment\nscore",
                                                 NULL),
                           xlab = "", ylab = "",...) +
        scale_fill_gradientn(colors = rescale_colors())+ #RPMG::SHOWPAL(ggsci::pal_gsea()(12))
        theme(plot.title = element_text(hjust = hjust,size = font.main))
    if(size == "padj") plot = plot + scale_size(breaks=c(0,0.05,0.10,0.15,0.2,0.25),
                                                labels=rev(c(0,0.05,0.10,0.15,0.2,0.25)))
    if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    if(is.null(file.name)) file.name =  paste0("Dotplot_",title,"_",pathway.name,
                                               "_",padj,"_",pval,".jpeg")
    jpeg(paste0(save.path, "/", file.name),units=units, width=width, height=height,res=600)
    print(plots)
    dev.off()
    if(do.return & return.raw) {
        return(fgseaRes)
    } else if(do.return) return(df_fgseaRes)
}


#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @param p.adjust.methods c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' #' @export gde.all data frame
#' @example FindAllMarkers.UMI(object)
FindAllMarkers.UMI <- function (object, assay = NULL, features = NULL, logfc.threshold = 0.25, 
                              test.use = "MAST", slot = "data", min.pct = 0.1, min.diff.pct = -Inf, 
                              p.adjust.methods = "bonferroni",
                              node = NULL, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                              random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
                              min.cells.group = 3, pseudocount.use = 1, return.thresh = 0.01, 
                              ...) 
    {
        MapVals <- function(vec, from, to) {
            vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
            vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
            return(unname(obj = vec2))
        }
        if ((test.use == "roc") && (return.thresh == 0.01)) {
            return.thresh <- 0.7
        }
        if (is.null(x = node)) {
            idents.all <- sort(x = unique(x = Idents(object = object)))
        } else {
            tree <- Tool(object = object, slot = "BuildClusterTree")
            if (is.null(x = tree)) {
                stop("Please run 'BuildClusterTree' before finding markers on nodes")
            }
            descendants <- DFT(tree = tree, node = node, include.children = TRUE)
            all.children <- sort(x = tree$edge[, 2][!tree$edge[, 
                                                               2] %in% tree$edge[, 1]])
            descendants <- MapVals(vec = descendants, from = all.children, 
                                   to = tree$tip.label)
            drop.children <- setdiff(x = tree$tip.label, y = descendants)
            keep.children <- setdiff(x = tree$tip.label, y = drop.children)
            orig.nodes <- c(node, as.numeric(x = setdiff(x = descendants, 
                                                         y = keep.children)))
            tree <- drop.tip(phy = tree, tip = drop.children)
            new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
            idents.all <- (tree$Nnode + 2):max(tree$edge)
        }
        genes.de <- list()
        messages <- list()
        for (i in 1:length(x = idents.all)) {
            if (verbose) {
                message("Calculating cluster ", idents.all[i])
            }
            genes.de[[i]] <- tryCatch(expr = {
                FindMarkers.UMI(object = object, assay = assay, ident.1 = if (is.null(x = node)) {
                    idents.all[i]
                }
                else {
                    tree
                }, ident.2 = if (is.null(x = node)) {
                    NULL
                }
                else {
                    idents.all[i]
                }, features = features, logfc.threshold = logfc.threshold, 
                test.use = test.use, slot = slot, min.pct = min.pct, 
                p.adjust.methods = p.adjust.methods,
                min.diff.pct = min.diff.pct, verbose = verbose, 
                only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
                random.seed = random.seed, latent.vars = latent.vars, 
                min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                pseudocount.use = pseudocount.use, ...)
            }, error = function(cond) {
                return(cond$message)
            })
            if (class(x = genes.de[[i]]) == "character") {
                messages[[i]] <- genes.de[[i]]
                genes.de[[i]] <- NULL
            }
        }
        gde.all <- data.frame()
        for (i in 1:length(x = idents.all)) {
            if (is.null(x = unlist(x = genes.de[i]))) {
                next
            }
            gde <- genes.de[[i]]
            if (nrow(x = gde) > 0) {
                if (test.use == "roc") {
                    gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                                         myAUC < (1 - return.thresh)))
                }
                else if (is.null(x = node) || test.use %in% c("bimod", 
                                                              "t")) {
                    gde <- gde[order(gde$p_val, -gde[, 2]), ]
                    gde <- subset(x = gde, subset = p_val < return.thresh)
                }
                if (nrow(x = gde) > 0) {
                    gde$cluster <- idents.all[i]
                    gde$gene <- rownames(x = gde)
                }
                if (nrow(x = gde) > 0) {
                    gde.all <- rbind(gde.all, gde)
                }
            }
        }
        if ((only.pos) && nrow(x = gde.all) > 0) {
            return(subset(x = gde.all, subset = gde.all[, 2] > 0))
        }
        rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
        if (nrow(x = gde.all) == 0) {
            warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
        }
        if (length(x = messages) > 0) {
            warning("The following tests were not performed: ", 
                    call. = FALSE, immediate. = TRUE)
            for (i in 1:length(x = messages)) {
                if (!is.null(x = messages[[i]])) {
                    warning("When testing ", idents.all[i], " versus all:\n\t", 
                            messages[[i]], call. = FALSE, immediate. = TRUE)
                }
            }
        }
        if (!is.null(x = node)) {
            gde.all$cluster <- MapVals(vec = gde.all$cluster, from = new.nodes, 
                                       to = orig.nodes)
        }
        return(gde.all)
    }


#' FindIdentLabel: Find identical label between ident and metadata
#' @object seurat object
#' @label colname in metadata
FindIdentLabel <- function(object){
    ident.label <- Idents(object)
    labels <- apply(object@meta.data,2,
                     function(x) all(ident.label == x)) %>% .[.] %>% .[!is.na(.)]
    label <- names(labels[labels])
    label =  label[!(label %in% c("seurat_clusters","ident"))][1]
    return(label)
}


#' Seurat 3
#' Calculate average UMI and attach to FindMarkers results 
#' Modified Seurat::FindMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @param p.adjust.methods p.adjust.methods
#' @export gde.all data frame
#' @example FindMarkers.UMI(object,ident.1 = "1")
FindMarkers.UMI <- function (object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, 
                              subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL, 
                              features = NULL, logfc.threshold = 0.25, test.use = "wilcox", 
                             p.adjust.methods = "fdr",
                              min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, 
                              max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                              min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1, 
                              ...){
    if (!is.null(x = group.by)) {
        if (!is.null(x = subset.ident)) {
            object <- subset(x = object, idents = subset.ident)
        }
        Idents(object = object) <- group.by
    }
    if (!is.null(x = assay) && !is.null(x = reduction)) {
        stop("Please only specify either assay or reduction.")
    }
    data.slot <- ifelse(
        test = test.use %in% c("negbinom", "poisson", "DESeq2"),
        yes = 'counts',
        no = slot
    )
    if (is.null(x = reduction)) {
        assay <- assay %||% DefaultAssay(object = object)
        data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)
    } else {
        if (data.slot == "counts") {
            stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
        }
        data.use <- t(x = Embeddings(object = object, reduction = reduction))
    }
    if (is.null(x = ident.1)) {
        stop("Please provide ident.1")
    } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
        if (is.null(x = ident.2)) {
            stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
        }
        tree <- if (is(object = ident.1, class2 = 'phylo')) {
            ident.1
        } else {
            Tool(object = object, slot = 'BuildClusterTree')
        }
        if (is.null(x = tree)) {
            stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
        }
        ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
        ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
    }
    if (length(x = as.vector(x = ident.1)) > 1 &&
        any(as.character(x = ident.1) %in% colnames(x = data.use))) {
        bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
        if (length(x = bad.cells) > 0) {
            stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
        }
    } else {
        ident.1 <- WhichCells(object = object, idents = ident.1)
    }
    # if NULL for ident.2, use all other cells
    if (length(x = as.vector(x = ident.2)) > 1 &&
        any(as.character(x = ident.2) %in% colnames(x = data.use))) {
        bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
        if (length(x = bad.cells) > 0) {
            stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
        }
    } else {
        if (is.null(x = ident.2)) {
            ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
        } else {
            ident.2 <- WhichCells(object = object, idents = ident.2)
        }
    }
    if (!is.null(x = latent.vars)) {
        latent.vars <- FetchData(
            object = object,
            vars = latent.vars,
            cells = c(ident.1, ident.2)
        )
    }
    counts <- switch(
        EXPR = data.slot,
        'scale.data' = GetAssayData(object = object[[assay]], slot = "counts"),
        numeric()
    )
    de.results <- FindMarkers.Assay(
        object = data.use,
        slot = data.slot,
        counts = counts,
        cells.1 = ident.1,
        cells.2 = ident.2,
        features = features,
        reduction = reduction,
        logfc.threshold = logfc.threshold,
        test.use = test.use,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        verbose = verbose,
        only.pos = only.pos,
        max.cells.per.ident = max.cells.per.ident,
        random.seed = random.seed,
        latent.vars = latent.vars,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        pseudocount.use = pseudocount.use,
        ...)
    if(p.adjust.methods != "bonferroni") {
        de.results$p_val_adj = p.adjust(
            p = de.results$p_val,
            method = p.adjust.methods,
            n = nrow(x = object)
        )
    }
    #de.results$avg_logFC = log2(exp(1)) * de.results$avg_logFC
    if(slot == "data"){
        avg_UMI.1 <- Matrix::rowMeans(expm1(x = data.use[, ident.1]))
        avg_UMI.2 <- Matrix::rowMeans(expm1(x = data.use[, ident.2]))
    }
    else if(slot == "counts"){
        avg_UMI.1 <- Matrix::rowMeans(x = data.use[, ident.1])
        avg_UMI.2 <- Matrix::rowMeans(x = data.use[, ident.2])
    }
    avg_UMI <-data.frame(avg_UMI.1, avg_UMI.2)
    de.results <- cbind(de.results,avg_UMI[match(rownames(de.results),rownames(avg_UMI)),])
    
    return(de.results)
}



#' find marker across by conditions
#' Modified FindMarkers.UMI function, compare the same ident across conditions
#' @param ident.1 dent.1 list
#' @param ident.2 dent.2 list
#' @param p.adjust.methods c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @export save.path folder to save
#' @example FindPairMarkers(object, ident.1 = 1:8, ident.2 = c(5:8,1:4))
FindPairMarkers <- function(object, ident.1, ident.2 = NULL, features = NULL,
                            group.by = NULL, return.thresh = 0.05,
                            logfc.threshold = 0.05, test.use = "MAST", min.pct = 0.1,
                            min.diff.pct = -Inf, only.pos = FALSE,
                            p.adjust.methods = "fdr",
                            max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
                            min.cells.feature = 3,min.cells.group=3, slot = "data",
                            assay = NULL,save.path = NULL,save.files = FALSE,...){
    #prepare save folder name
    if(save.files){
        if(class(ident.1) == "numeric" & class(ident.2) == "numeric") {
            ident1="";ident2=""
        } else {
            ident1 <- unique(gsub('\\_.*', '', ident.1))
            ident2 <- unique(gsub('\\_.*', '', ident.2))
        }
    }
    if(length(ident.1) != length(ident.1)) stop("pair lenth don't match")
    assay = DefaultAssay(object) %||% "RNA"
    gde <- list()
    for(i in 1:length(ident.1)) {
        ident.1vs2 <- paste(ident.1[i],"/", ident.2[i])
        print(ident.1vs2)
        ident_1 <- if(class(ident.1)=="list") ident.1[[i]] else ident.1[i]
        ident_2 <- if(class(ident.2)=="list") ident.2[[i]] else ident.2[i]
        gde[[i]] <- FindMarkers.UMI(object = object, 
                                    ident.1 = ident_1,
                                    ident.2 = ident_2, 
                                    assay = assay, group.by = group.by,
                                    features = features, slot = slot,
                                    p.adjust.methods = p.adjust.methods,
                                    logfc.threshold = logfc.threshold, test.use = test.use, 
                                    min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                    only.pos = only.pos, min.cells.feature = min.cells.feature, 
                                    min.cells.group = min.cells.group, latent.vars = latent.vars, 
                                    max.cells.per.ident = max.cells.per.ident,...)
        gde[[i]] <- gde[[i]][order(-gde[[i]]$avg_logFC,gde[[i]]$p_val),]
        gde[[i]] <- subset(x = gde[[i]], subset = p_val < return.thresh)
        gde[[i]]$cluster1.vs.cluster2 <- ident.1vs2
        gde[[i]]$gene <- rownames(x = gde[[i]])
        if(save.files){
            if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
            if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
            write.csv( gde[[i]], paste0(save.path,ident2,"_vs_",ident1,".csv"))
        }
    }
    return(bind_rows(gde))
}


# FilterGenes
#' filter gene names according to seurat object, produce uniformed gene format
#' pass non-existing genes or ill-formaed genes names to downstream analysis
#' will generate error
#' @param object Seurat object version 2
#' @param marker.genes gene names, marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
#' @param unique return unique name or not
#' @export marker.genes filtered, well-format gene names
#' @example FilterGenes(object, c("cdh5,PECAM1,flt1,Vwf,Plvap,Kdr"))
FilterGenes <- function(object, marker.genes, unique= TRUE, verbose = T){
        if(missing(object)) 
                stop("A seurat object must be provided first")
        if(class(object) != "Seurat")
                stop("A seurat object must be provided first")
        if(missing(marker.genes)) 
                stop("A list of marker genes must be provided")

        marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        marker.genes <- gsub(" ","",marker.genes)
        species = CheckSpecies(object)
        if(species == "Human") marker.genes <- toupper(marker.genes)
        if(species == "Mouse") marker.genes <- Hmisc::capitalize(tolower(marker.genes)) 

        if(verbose) print(paste("Before filtration:",length(marker.genes)))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(object))
        if(unique) marker.genes <- unique(marker.genes)
        if(verbose) print(paste("After filtration:",length(marker.genes)))
        return(as.character(marker.genes))
}


FPKM <- function(counts, lengths) {
    rownames(counts) = tolower(rownames(counts))
    names(lengths) = tolower(names(lengths))
    A = intersect(rownames(counts), names(lengths))
    counts = counts[A, ]
    lengths = lengths[A]
    rate = counts/lengths
    sweep(rate,2,FUN="/",STATS=colSums(counts)) * 1e+06
}

#=====Clean memory======================
GC <- function()
{
    while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1,
                                                         4]) {
    }
}


#' Extract ColorHexa from Seurat TSNE plot
#' @param object aligned seurat object with ident
#' @param return.vector TRUE/return full color vector, FALSE/return color levels
#' @param cells.use only return ColorHexa of selected cells
#' @param ... other TSNEPlot inputs 
#' @export colors: color vector named by cell ID
gg_colors <- function(object = object, return.vector=FALSE, cells.use = NULL,
                      print.plot = T, no.legend = TRUE, do.label = TRUE,colors.use =NULL,
                      do.return = TRUE, label.size = 6, gg_title="", ...){
    
    g1 <- Seurat::TSNEPlot(object = object, no.legend = no.legend,
                           do.label = do.label,do.return = do.return,
                           label.size = label.size, colors.use= colors.use,...)
    if(print.plot) print(g1)
    g <- ggplot2::ggplot_build(g1)
    #        print(unique(g$data[[1]]["colour"]))
    colors_df <- g$data[[1]][,c("colour","group")]
    
    #select color by cells ID
    cells <- Seurat::WhichCells(object)
    colors_df$cell.names <- cells
    if(!is.null(cells.use)) {
        colors_df <- colors_df[(colors_df$cell.names %in% cells.use),]
    }
    colors_df = colors_df[order(colors_df$group),]
    if(return.vector) {
        print(head(colors));print(length(colors))
        return(colors)
    } else {
        colors <- unique(colors_df$colour)
        print(colors)
        return(as.character(colors))
    }
}



# select ggplot color
gg_color_hue <- function(n=5) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
#gg_color_hue(4)
#library(scales)
#show_col(hue_pal()(6))

# Basic function to convert human to mouse gene names
#' @example genes <- Human2Mouse(humGenes)
Human2Mouse <- function(x){
    
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                     values = x , mart = human, attributesL = c("mgi_symbol"), 
                     martL = mouse, uniqueRows=T)
    
    humanx <- unique(genesV2[, 2])
    
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}

#' Convert list to data frame
#'
#' This function will convert a list to a data frame, even if they are unequal length
#'
#' @param fname list
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
list2df <- function(list){
        cbind.fill <- function(...){
                nm <- list(...) 
                nm <- lapply(nm, as.matrix)
                n <- max(sapply(nm, nrow)) 
                do.call(cbind, lapply(nm, function (x) 
                        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
        }

        df <- do.call(cbind.fill, list)
        colnames(df) = names(list)
        return(df)
}


# make corlor bar for DoHeatmap
#' @param object Seurat object
#' @param dge_markers FindAllMarkers results
#' @param Top_n top_n(Top_n, avg_logFC)
#' @param features extra genes to add beyound FindAllMarkers results 
#' @param remove.legend TRUE/FALSE
#' @param color color scheme
#' @export g vertical ggplot
#' @example MakeCorlorBar(EC, top, Top_n=40)
MakeCorlorBar <- function(object, dge_markers = NULL, Top_n = NULL, features = NULL, color =NULL,unique.name = F,
                          no.legend =F, legend.size = 10,do.print = TRUE,save.path =NULL,
                          do.return=FALSE,width=10, height=7){
    if(unique.name) {
        v <- paste(unique(object$orig.ident),collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
    if(!no.legend) v = paste0(v, "_Legend")
    
    if(!is.null(features)){
        marker_bar <- data.frame("gene" = features,
                                 "cluster" = "marker.genes")
        rownames(marker_bar) = marker_bar$gene 
    } else marker_bar <- NULL
    
    if(!is.null(dge_markers)){
        colnames(dge_markers)[grep("cluster",colnames(dge_markers))] ="cluster"
        dge_markers %<>% group_by(cluster)
        if(!is.null(Top_n)) dge_markers %<>% top_n(Top_n, avg_logFC)
        dge_markers = rbind(dge_markers[,c("gene","cluster")],
                                     marker_bar)
    } else dge_markers <- marker_bar

    gene.use = rownames(object@assays$RNA@data)
    dge_markers = dge_markers[!duplicated(dge_markers$gene),]
    dge_markers = dge_markers[dge_markers$gene %in% gene.use,]
    
    dge_markers$x <- 1:nrow(dge_markers)
    
    wide <- table(dge_markers$cluster) %>% as.data.frame
    colnames(wide) = c("cluster","w")
    
    dge_markers %<>%   dplyr::left_join(wide, by = "cluster") %>%
            group_by(cluster) %>%
            dplyr::summarise(median = median(x, na.rm = TRUE)) %>%
            dplyr::inner_join(wide, by = "cluster")
    
    g <- ggplot(data = dge_markers, aes(xmin = median - w / 2, xmax = median + w / 2,
                                      ymin = 0, ymax = 1, fill = cluster)) +
        geom_rect(colour = "white")+
        theme_bw() + 
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_blank(), 
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
    coord_flip() + scale_x_reverse()
        
    if(!no.legend) g = g + theme(legend.title = element_text(size=legend.size*1.2),
                                 legend.text = element_text(size=legend.size),
                                 legend.key.size = unit(legend.size/5,"line"))
    if(no.legend) g = g + theme(legend.position="none")
    if(!is.null(color)) {
        colors_fill = color_scheme(color = color,
                                   names = unique(dge_markers$cluster))
        g = g + scale_fill_manual(values = colors_fill)
    }
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path,"Doheatmap_vcolorbar_",v,".jpeg"),
             units=units, width=width, height=height,res=600)
        print(g)
        dev.off()
    } else if(do.return & Sys.info()[['sysname']] != "Linux") {
        return(g)
    }
}


# Support function for TSNEPlot.1, modified from Seurate:::LabelClusters function
LabelRepel <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, 
                       repel = TRUE,color = NULL,alpha= NULL, ...){
    xynames <- unlist(x = Seurat:::GetXYAesthetics(plot = plot), use.names = TRUE)
    if (!id %in% colnames(x = plot$data)) {
        stop("Cannot find variable ", id, " in plotting data")
    }
    if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
        warning("Cannot find splitting variable ", id, " in plotting data")
        split.by <- NULL
    }
    data <- plot$data[, c(xynames, id, split.by)]
    g <- ggplot_build(plot)
    data_color <- g$data[[1]]
    data = cbind(data, color = data_color$colour[match(data[,1], data_color$x)])
    possible.clusters <- base::as.character(x = na.omit(object = unique(x = data[,id])))
    groups <- clusters %||% base::as.character(x = na.omit(object = unique(x = data[,id])))
    if (any(!groups %in% possible.clusters)) {
        stop("The following clusters were not found: ", 
             paste(groups[!groups %in% 
                              possible.clusters], collapse = ","))
    }
    labels.loc <- lapply(X = groups, FUN = function(group) {
        data.use <- data[data[, id] == group, , drop = FALSE]
        data.medians <- if (!is.null(x = split.by)) {
            do.call(what = "rbind", 
                    args = lapply(X = unique(x = data.use[, split.by]), 
                                  FUN = function(split) {
                                      medians <- apply(X = data.use[data.use[, split.by] == 
                                                                        split, xynames, drop = FALSE], 
                                                       MARGIN = 2, 
                                                       FUN = median, na.rm = TRUE)
                                      medians <- as.data.frame(x = t(x = medians))
                                      medians[, split.by] <- split
                                      return(medians)
                                      }))
        } else {
            as.data.frame(x = t(x = apply(X = data.use[, xynames, 
                                                       drop = FALSE], MARGIN = 2, 
                                          FUN = median, na.rm = TRUE)))
        }
        data.medians[, id] <- group
        return(data.medians)
    })
    labels.loc <- do.call(what = "rbind", args = labels.loc)
    labels <- labels %||% groups
    if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
        stop("Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (", 
             nrow(x = labels.loc), ").")
    }
    names(x = labels) <- groups
    for (group in groups) {
        labels.loc[labels.loc[, id] == group, id] <- labels[group]
    }
    labels.loc = labels.loc[order(labels.loc[,id]),]
    labels.loc = cbind.data.frame(labels.loc, 
                                 color = data$color[match(labels.loc[,id], data[,id])])
    geom.use <- ifelse(test = repel, yes = ggrepel::geom_label_repel, 
                       no = geom_label)
    plot = plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames["x"], 
                                                                   y = xynames["y"], 
                                                                   label = id),
                           colour = labels.loc$color,
                           alpha = alpha,...)
    return(plot)
    
}

#Make rownamess with duplicated values unique in a dataframe
#'@param object Seurat object or expresion matrix
#'@param features genes to be added in scale.datas
MakeUniqueGenes <- function(object, features, verbose = T){
    if(class(object) == "Seurat") assay = DefaultAssay(object)
    featuresNum <- make.unique(features, sep = ".")
    dupFeatures <- duplicated(features)
    dupFeaturesIndex = which(dupFeatures)
    
    if(length(dupFeaturesIndex) == 0 ) {
        message("All genes are unique")
        return(object)
    }
    
    if(class(object) == "Seurat") {
        scale.data = object[[assay]]@scale.data
    } else if(class(object) %in% c("data.frame","matrix")){
        scale.data = object
    } else stop("expression profile must be provided")

    for (i in 1:length(dupFeaturesIndex)){
        if(verbose) svMisc::progress(i/length(dupFeaturesIndex)*100)
        scale.data <- rbind(scale.data, 
                            "dupGene" = scale.data[features[dupFeaturesIndex[i]], ])
        rownames(scale.data)[nrow(scale.data)] = featuresNum[dupFeaturesIndex[i]]
        if(i == length(dupFeaturesIndex)) cat("Done!\n")
    }
    if(class(object) == "Seurat") {
        object[[assay]]@scale.data = scale.data
        return(object)
    } else if(class(object) %in% c("data.frame","matrix")){
        return(scale.data)
    }
}


NoAxesLabel <- function (..., keep.text = FALSE, keep.ticks = FALSE) 
{
    blank <- element_blank()
    if (!keep.text) {
        no.axes.theme <- theme(axis.text.x = blank,
                               axis.text.y = blank, axis.title.x = blank, axis.title.y = blank,
                               validate = TRUE, ...)
    }
    if (!keep.ticks) {
        no.axes.theme <- no.axes.theme + theme(axis.ticks.x = blank, 
                                               axis.ticks.y = blank, validate = TRUE, ...)
    }
    return(no.axes.theme)
}

Progress <- function(i,len){
    svMisc::progress(i/len*100)
}


# modify hamonry commit ee0877a  for Seurat v3 
RunHarmony.1 <- function (object, group.by, dims.use, group.by.secondary = NULL, 
                          theta = 1, theta2 = 1, sigma = 0.1, alpha = 0.1, nclust = 100, 
                          tau = 0, block.size = 0.05, max.iter.harmony = 10, max.iter.cluster = 200, 
                          epsilon.cluster = 1e-05, epsilon.harmony = 1e-04, burn.in.time = 10, 
                          plot_convergence = FALSE)
{
    if (!"Seurat" %in% class(object)) {
        stop("This Function is meant to be run on a Seurat object!")
    }
    if (!"pca" %in% names(object@reductions)) {
        stop("PCA must be computed before running Harmony.")
    }
    if (missing(dims.use)) {
        dims.use <- 1:ncol(object@reductions$pca@cell.embeddings)
    }
    else if (!all(dims.use %in% 1:ncol(object@reductions$pca@cell.embeddings))) {
        stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
    }
    if (length(dims.use) == 1) {
        stop("only specified one dimension in dims.use")
    }
    if (!group.by %in% colnames(object@meta.data)) {
        stop(sprintf("ERROR: Primary integration variable [%s] is not in meta.data"))
    }
    missing.vars <- setdiff(group.by, colnames(object@meta.data))
    if (length(missing.vars) > 0) {
        msg <- gettextf(ngettext(length(missing.vars),
                                 "trying to integrate over missing variable: %s",
                                 "trying to integrate over missing variables: %s",
                                 domain = "R-base"),
                        paste(missing.vars, collapse = ", "))
        stop(msg)
    }
    message("Starting harmony")
    message(sprintf("Using top %d PCs", length(dims.use)))
    if (!is.null(group.by.secondary)) {
        batches_secondary <- object@meta.data[[group.by.secondary]]
    }
    else {
        batches_secondary <- NULL
    }
    harmonyEmbed <- HarmonyMatrix(object@reductions$pca@cell.embeddings,
                                  object@meta.data[[group.by]], batches_secondary, theta, 
                                  theta2, sigma, alpha, nclust, tau, block.size, max.iter.harmony, 
                                  max.iter.cluster, epsilon.cluster, epsilon.harmony, 
                                  burn.in.time, plot_convergence)
    rownames(harmonyEmbed) <- row.names(object@meta.data)
    colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
    
    #data.use <- PrepDR(object = object,features = NULL,verbose = F)
    #feature.loadings <- (as.matrix(x = data.use) %*% as.matrix(x = harmonyEmbed))
    
    object[["harmony"]] <- CreateDimReducObject(
        embeddings = harmonyEmbed,
        #loadings = feature.loadings,
        key = "harmony_",
        assay = DefaultAssay(object = object)
    )
    return(object)
}


paramSweep_v4 <- function (seu, PCs = 1:10, sct = FALSE) 
{
    require(Seurat)
    require(fields)
    require(parallel)
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
    pK.test <- round(pK * min.cells)
    pK <- pK[which(pK.test >= 1)]
    orig.commands <- seu@commands
    assay <- DefaultAssay(seu)
    if (nrow(seu@meta.data) > 10000) {
        real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                     10000, replace = FALSE)]
        data <- seu[[assay]]@counts[, real.cells]
        n.real.cells <- ncol(data)
    }
    if (nrow(seu@meta.data) <= 10000) {
        real.cells <- rownames(seu@meta.data)
        data <- seu[[assay]]@counts
        n.real.cells <- ncol(data)
    }
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    output2 <- list()
    for (n in 1:length(pN)){
        output2[[n]] <- parallel_paramSweep_v3(n, n.real.cells, real.cells, pK, pN, data, orig.commands, 
                                              PCs, sct)
        message(paste0("Complete ",n,":",length(pN)," ----------------------"))
    }
    stopCluster(cl)
    sweep.res.list <- list()
    list.ind <- 0
    for (i in 1:length(output2)) {
        for (j in 1:length(output2[[i]])) {
            list.ind <- list.ind + 1
            sweep.res.list[[list.ind]] <- output2[[i]][[j]]
        }
    }
    name.vec <- NULL
    for (j in 1:length(pN)) {
        name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                      sep = "_"))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)
}


# generate expression txt file for GSEA analysis
#' @param object Seurat object
#' @param k an integer for the number of folds. createFolds argment
#' @param do.return TRUE/FALSE
#' @param continuous.label NULL/continuous label #http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_Phenotype_Labels
#' @example PrepareGSEA(object, k = 50, continuous.label = major_cells)
PrepareGSEA <- function(object, k = 1, do.return = FALSE, continuous.label = NULL){
    
    set.seed(201)
    object@scale.data = NULL
    ident = FindIdentLabel(object)[1]
    
    if(k > 1){
        split_object <- SplitSeurat(object, split.by = ident)
        #Split object@ident to k group =====
        for(i in 1:length(split_object)){
            meta_index <- caret::createFolds(split_object[[i]]@meta.data[,ident],
                                         k = k, list = TRUE, 
                                         returnTrain = FALSE)
            for(n in 1:length(meta_index)){
                split_object[[i]]@meta.data[meta_index[[n]],"GSEA"] = 
                    paste(split_object[[i]]@meta.data[meta_index[[n]],ident],
                          n, sep = "_")
                }
        }
        new_object <-  Reduce(MergeSeurat,split_object) %>%
            SetAllIdent(id = "GSEA")
        if(!is.null(continuous.label))
            new_object@ident <- factor(x = new_object@ident,
                                    levels = paste0(rep(continuous.label,each = k),
                                                    "_",rep(1:k)))
        } else new_object <- object
    
        if(k == 1 & !is.null(continuous.label)) {
            new_object@ident <- factor(x = new_object@ident,
                                       levels = continuous.label)
        }

                
    print("#====Calculate Average Expression======")
    GSEA_expr <- AverageExpression(new_object)
    GSEA_name <- data.frame("NAME" = rownames(GSEA_expr),
                            "DESCRIPTION" = rep(NA,nrow(GSEA_expr)),
                            stringsAsFactors = F)
    GSEA_expr <- cbind.data.frame(GSEA_name,GSEA_expr)
    
    samples <- gsub("_([0-9]+).*$", "", colnames(GSEA_expr)[-c(1,2)])
    
    if(is.null(continuous.label)) {
        cls_list <- list(c(length(samples),length(unique(samples)), 1),
                     paste(c("#",unique(samples)), collapse = " "),
                     paste(samples, collapse = " "))
    } else 
        if(all(continuous.label %in% unique(object@ident))){
            numeric <- match(samples,continuous.label)
            cls_list <- list("#numeric",
                             paste(c("#",unique(samples)), collapse = "."),
                             paste(numeric, collapse = " "))
    } else stop("Incorrect continuous.label!")
    
    if(rownames(GSEA_expr)[1] == 
        Hmisc::capitalize(tolower(rownames(GSEA_expr)[1]))){
        print("#====Replace gene names ======")
        rownames.GSEA_expr = rownames(GSEA_expr)
        human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), #filters = "mgi_symbol",
                         values = rownames(GSEA_expr) , mart = mouse,
                         attributesL = c("hgnc_symbol"), 
                         martL = human, uniqueRows=T)
        rm = duplicated(genesV2[,1])
        genesV2 = genesV2[!rm,]
        colnames(genesV2) = c("gene","NAME")
        colnames(GSEA_expr)[1] = "gene"
        GSEA_expr <- merge(genesV2,GSEA_expr,by = "gene")
        GSEA_expr = GSEA_expr[,-1]
    }
    file.name = paste(unique(object@ident), collapse = "_")
    if(!is.null(continuous.label)) file.name <- paste(continuous.label,collapse = "_")
    if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    write.table(GSEA_expr, file = paste0(save.path, file.name,
                                         "_",k,".txt"),
                sep = "\t", quote = FALSE,row.names = FALSE)
    
    fn = paste0(save.path, file.name,"_",k,".cls")
    if (file.exists(fn)) file.remove(fn)
    lapply(cls_list, cat, "\n", file=fn,append=TRUE)
    
    if(do.return) return(GSEA_expr)
    
}


# modify the RNA velocity function to avoid error
# Error in seq.default(rx[1], rx[2], length.out = grid.n) : 'from' must be a finite number
# provide Seurat object to remove cells
prepare.velocity.on.embedding.cor <- function (object, reduction = "umap", slot = "RunVelocity",
                                               n = 100,scale = "log",  n.cores = parallel::detectCores(logical=F), 
                                               cc = NULL) 
{
    emb = Embeddings(object = object, reduction = reduction)
    vel = Tool(object = object, slot = slot)
    randomize <- FALSE
    em <- as.matrix(vel$current)
    ccells <- intersect(rownames(emb), colnames(em))
    em <- em[, ccells]
    emb <- emb[ccells, ]
    nd <- as.matrix(vel$deltaE[, ccells])
    cgenes <- intersect(rownames(em), rownames(nd))
    nd <- nd[cgenes, ]
    em <- em[cgenes, ]
    if (randomize) {
        nd <- t(apply(nd, 1, function(x) (rbinom(length(x), 
                                                 1, 0.5) * 2 - 1) * abs(sample(x))))
    }
    colDeltaCorSqrt <- function(e, d, nthreads = 1L) {
        .Call('_velocyto_R_colDeltaCorSqrt', PACKAGE = 'velocyto.R', e, d, nthreads)
    }
    if (is.null(cc)) {
        cat("delta projections ... ")
        if (scale == "log") {
            cat("log ")
            cc <- colDeltaCorLog10(em, (log10(abs(nd) + 1) * 
                                            sign(nd)), nthreads = n.cores)
        }
        else if (scale == "sqrt") {
            cat("sqrt ")
            cc <- colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)), 
                                  nthreads = n.cores)
        }
        else if (scale == "rank") {
            cat("rank ")
            cc <- colDeltaCor((apply(em, 2, rank)), (apply(nd, 
                                                           2, rank)), nthreads = n.cores)
        }
        else {
            cat("linear ")
            cc <- colDeltaCor(em, nd, nthreads = n.cores)
        }
        colnames(cc) <- rownames(cc) <- colnames(em)
        diag(cc) <- 0
    }
    if(anyNA(cc)){
        mtx = which(is.na(cc), arr.ind = TRUE)
        object <- subset(object, cells = unique(rownames(mtx)), invert = T)
    }
    return(object)
}



#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#' @param barcodes.fileName = 'barcodes.tsv',
#' @param gene.fileName = 'genes.tsv',
#' @param features.fileName = 'features.tsv.gz',
#' @param matrix.fileName = 'matrix.mtx',
#' 
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For output from CellRanger < 3.0
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#'
#' # For output from CellRanger >= 3.0 with multiple data types
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#' data <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
#' seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
#' }
#'
Read10X.1 <- function(
        data.dir = NULL,
        gene.column = 2,
        unique.features = TRUE,
        strip.suffix = FALSE,
        barcodes.fileName = 'barcodes.tsv',
        gene.fileName = 'genes.tsv',
        features.fileName = 'features.tsv.gz',
        matrix.fileName = 'matrix.mtx') {
        full.data <- list()
        for (i in seq_along(along.with = data.dir)) {
                run <- data.dir[i]
                if (!dir.exists(paths = run)) {
                        stop("Directory provided does not exist")
                }
                barcode.loc <- file.path(run, barcodes.fileName)
                gene.loc <- file.path(run, gene.fileName)
                features.loc <- file.path(run, features.fileName)
                matrix.loc <- file.path(run, matrix.fileName)
                # Flag to indicate if this data is from CellRanger >= 3.0
                pre_ver_3 <- file.exists(gene.loc)
                if (!pre_ver_3) {
                        addgz <- function(s) {
                                return(paste0(s, ".gz"))
                        }
                        barcode.loc <- addgz(s = barcode.loc)
                        matrix.loc <- addgz(s = matrix.loc)
                }
                if (!file.exists(barcode.loc)) {
                        stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
                }
                if (!pre_ver_3 && !file.exists(features.loc) ) {
                        stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
                }
                if (!file.exists(matrix.loc)) {
                        stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
                }
                data <- readMM(file = matrix.loc)
                cell.names <- readLines(barcode.loc)
                if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
                        cell.names <- as.vector(x = as.character(x = sapply(
                                X = cell.names,
                                FUN = ExtractField,
                                field = 1,
                                delim = "-"
                        )))
                }
                if (is.null(x = names(x = data.dir))) {
                        if (i < 2) {
                                colnames(x = data) <- cell.names
                        } else {
                                colnames(x = data) <- paste0(i, "_", cell.names)
                        }
                } else {
                        colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
                }
                feature.names <- read.delim(
                        file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
                        header = FALSE,
                        stringsAsFactors = FALSE
                )
                if (any(is.na(x = feature.names[, gene.column]))) {
                        warning(
                                'Some features names are NA. Replacing NA names with ID from the opposite column requested',
                                call. = FALSE,
                                immediate. = TRUE
                        )
                        na.features <- which(x = is.na(x = feature.names[, gene.column]))
                        replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
                        feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
                }
                if (unique.features) {
                        fcols = ncol(x = feature.names)
                        if (fcols < gene.column) {
                                stop(paste0("gene.column was set to ", gene.column,
                                            " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                                            " Try setting the gene.column argument to a value <= to ", fcols, "."))
                        }
                        rownames(x = data) <- make.unique(names = feature.names[, gene.column])
                }
                # In cell ranger 3.0, a third column specifying the type of data was added
                # and we will return each type of data as a separate matrix
                if (ncol(x = feature.names) > 2) {
                        data_types <- factor(x = feature.names$V3)
                        lvls <- levels(x = data_types)
                        if (length(x = lvls) > 1 && length(x = full.data) == 0) {
                                message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
                        }
                        expr_name <- "Gene Expression"
                        if (expr_name %in% lvls) { # Return Gene Expression first
                                lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
                        }
                        data <- lapply(
                                X = lvls,
                                FUN = function(l) {
                                        return(data[data_types == l, , drop = FALSE])
                                }
                        )
                        names(x = data) <- lvls
                } else{
                        data <- list(data)
                }
                full.data[[length(x = full.data) + 1]] <- data
        }
        # Combine all the data from different directories into one big matrix, note this
        # assumes that all data directories essentially have the same features files
        list_of_data <- list()
        for (j in 1:length(x = full.data[[1]])) {
                list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
                # Fix for Issue #913
                list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
        }
        names(x = list_of_data) <- names(x = full.data[[1]])
        # If multiple features, will return a list, otherwise
        # a matrix.
        if (length(x = list_of_data) == 1) {
                return(list_of_data[[1]])
        } else {
                return(list_of_data)
        }
}


# Run GSEA and generate reports
#'@example ReportGSEA(file = "c5.all.827-D14.827-D16.827-D28.827-D0.Gsea.1552707417126",pos=T,ncol=3)
ReportGSEA <- function(file, pos=T,ncol = 3){
        (gsea_path <- paste("~/gsea_home/output",tolower(format(Sys.Date(), "%b%d")),
                            file,sep ='/'))
        #(pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*pos.*xls"))
        p<-1; if(pos) p <-2
        (pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*xls"))
        GSEA_output <- readr::read_delim(paste0(gsea_path,"/",pos.xls.path[p]),"\t", 
                                         escape_double = FALSE, trim_ws = TRUE)
        print(GSEA_output %>% .[-c(2,3,12)] %>% head(50) %>% kable() %>% kable_styling())
        
        (GSEA.plots <- sapply(GSEA_output$NAME[1:9], function(name) {
                paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                        sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                        .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                        paste(gsea_path, ., sep = "/")) 
        CombPngs(GSEA.plots, ncol = ncol)
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        GSEA_path <- paste0(path, "GSEA_xls")
        if(!dir.exists(GSEA_path)) dir.create(GSEA_path, recursive = T)
        file.copy(paste0(gsea_path,"/",pos.xls.path),GSEA_path)
        file.rename(paste0(path,list.files(path,pattern=paste0("GSEA.plots_CombPngs"))), 
                    paste0(path,file,"-",p,".jpeg"))
        
}


# draw rectangle by four corner coordinates
#' @param x_left,x_right,y_bottom,y_top, coordinate values anti-clockwise from left bottom
#' @param ... argments pass to geom_segment
#' @example rectangle(-4.6, -4.1, -4.7, -5.5,colour = "blue")
rectangle <- function(x_left = -Inf, x_right = Inf, y_bottom = -Inf, y_top = Inf, ...){
    list(geom_segment(aes(x = x_left, xend = x_right, y = y_top, yend = y_top),...),
    geom_segment(aes(x = x_right, xend = x_right, y = y_top, yend = y_bottom),...),
    geom_segment(aes(x = x_left, xend = x_right, y = y_bottom, yend = y_bottom),...),
    geom_segment(aes(x = x_left, xend = x_left, y = y_top, yend = y_bottom),...))
}


SetQuantile <- function (cutoff, data) {
    if (grepl(pattern = "^q[0-9]{1,2}$", x = base::as.character(x = cutoff),
    perl = TRUE)) {
        this.quantile <- as.numeric(x = sub(pattern = "q", replacement = "",
        x = base::as.character(x = cutoff)))/100
        data <- unlist(x = data)
        data <- data[data > 0]
        cutoff <- quantile(x = data, probs = this.quantile)
    }
    return(as.numeric(x = cutoff))
}


# add alpha
SingleDimPlot.1 <- function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL, alpha = 1,
                             shape.by = NULL, order = NULL, label = FALSE, repel = FALSE, 
                             label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                             sizes.highlight = 1, na.value = "gray97") 
{
    pt.size <- pt.size %||% AutoPointSize(data = data)
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    if (!is.data.frame(x = data)) {
        data <- as.data.frame(x = data)
    }
    if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
        stop("Cannot find dimensions to plot in data")
    }
    else if (is.numeric(x = dims)) {
        dims <- colnames(x = data)[dims]
    }
    if (!is.null(x = cells.highlight)) {
        highlight.info <- SetHighlight(cells.highlight = cells.highlight, 
                                       cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||% 
                                           pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||% 
                                           "#C3C3C3", pt.size = pt.size)
        order <- highlight.info$plot.order
        data$highlight <- highlight.info$highlight
        col.by <- "highlight"
        pt.size <- highlight.info$size
        cols <- highlight.info$color
    }
    if (!is.null(x = order) && !is.null(x = col.by)) {
        if (typeof(x = order) == "logical") {
            if (order) {
                data <- data[order(data[, col.by]), ]
            }
        }
        else {
            order <- rev(x = c(order, setdiff(x = unique(x = data[, 
                                                                  col.by]), y = order)))
            data[, col.by] <- factor(x = data[, col.by], levels = order)
            new.order <- order(x = data[, col.by])
            data <- data[new.order, ]
            if (length(x = pt.size) == length(x = new.order)) {
                pt.size <- pt.size[new.order]
            }
        }
    }
    if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
        warning("Cannot find ", col.by, " in plotting data, not coloring plot")
        col.by <- NULL
    }
    else {
        col.index <- match(x = col.by, table = colnames(x = data))
        if (grepl(pattern = "^\\d", x = col.by)) {
            col.by <- paste0("x", col.by)
        }
        else if (grepl(pattern = "-", x = col.by)) {
            col.by <- gsub(pattern = "-", replacement = ".", 
                           x = col.by)
        }
        colnames(x = data)[col.index] <- col.by
    }
    if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
        warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
    }
    plot <- ggplot(data = data) + geom_point(mapping = aes_string(x = dims[1], 
                                                                  y = dims[2], color = paste0("`", col.by, "`"), 
                                                                  shape = shape.by), alpha = alpha,
                                             size = pt.size) + guides(color = guide_legend(override.aes = list(size = 3))) + 
        labs(color = NULL)
    if (label && !is.null(x = col.by)) {
        plot <- LabelClusters(plot = plot, id = col.by, repel = repel, 
                              size = label.size)
    }
    if (!is.null(x = cols)) {
        plot <- plot + if (length(x = cols) == 1 && (is.numeric(x = cols) || 
                                                     cols %in% rownames(x = brewer.pal.info))) {
            scale_color_brewer(palette = cols, na.value = na.value)
        }
        else {
            scale_color_manual(values = cols, na.value = na.value)
        }
    }
    plot <- plot + cowplot::theme_cowplot()
    return(plot)
}


#' re-order seurat idents factors
sortIdent <- function(object,numeric=F){
    
    levels  <- Idents(object) %>% unique %>% base::as.character()
    if(numeric) levels %<>% as.numeric 
    #if(!is.null(FindIdentLabel(object))) {
    #    object@meta.data[,FindIdentLabel(object)] %<>% factor(levels = sort(levels))
    #}
    Idents(object) %<>% factor(levels = sort(levels))
    
    return(object)
}

#' Modified UMAPPlot
#' @param label.repel 
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
PCAPlot.1 <- function(object,dims = c(1, 2),cells = NULL,cols = NULL, pt.size = NULL,
                       reduction = "pca",group.by = NULL,split.by = NULL,shape.by = NULL,
                       order = NULL,label = FALSE,label.repel= FALSE, label.size = 4,
                       repel = TRUE,alpha = 0.85, text.size = 15,title.size = 15,
                       legend.size = 15,
                       cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
                       na.value = 'gray97',combine = TRUE,ncol = NULL,title = NULL,legend.title = NULL,
                       no.legend = F,do.print = F,do.return = T,unique.name=F,
                       units= "in",width=10, height=7,hjust = 0.5,border = FALSE,
                       save.path = NULL,file.name = NULL, ...) {
    if(is.null(file.name)){
        VarName <- UniqueName(object,fileName = deparse(substitute(object)),unique.name = unique.name)
        VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object))
        if(!no.legend) VarName = paste0(VarName, "_Legend")
        VarName = paste0(VarName, ifelse(!is.null(split.by),
                                         yes = paste0("_",split.by),
                                         no =""))
        L = ifelse(no.legend, "", "_Legend")
        if(!label) L = paste0(L,"_noLabel")
        file.name = paste0("PCAPlot_",VarName,L,".jpeg")
    }
    cols = cols %||% ExtractMetaColor(object, group.by = group.by)
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$pca@cell.embeddings[cells,dims]
    data <- as.data.frame(x = data)
    dims <- paste0(object@reductions$pca@key, dims)
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- Seurat:::SingleDimPlot(data = data[, c(dims, x, split.by, shape.by)], 
                                       dims = dims, col.by = x, cols = cols,
                                       pt.size = pt.size, shape.by = shape.by, order = order,
                                       label = FALSE, cells.highlight = cells.highlight,
                                       cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                                       na.value = na.value,...)
        if (label & label.repel == F ) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
                                  size = label.size, split.by = split.by)
        }
        if (label & label.repel) {
            plot <- LabelRepel(plot = plot, id = x, repel = repel,
                               size = label.size, split.by = split.by, color= cols,
                               alpha = alpha)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + theme(strip.background = element_blank(),
                                 strip.text = element_text(face="plain",size=text.size))+ 
                facet_wrap(facets = vars(!!sym(x = split.by)),
                           ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                               length(x = unique(x = data[, split.by]))
                           } else ncol )
            
        }
        if(border == TRUE) plot = plot + NoAxesLabel() + 
                theme(panel.border = element_rect(colour = "black"))
        return(plot)
    })
    if (combine) {
        plots <- patchwork::wrap_plots(plots = plots, ncol = if (!is.null(x = split.by) &&
                                                        length(x = group.by) > 1) {
            1
        }
        else {
            ncol
        }, ...)
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
            theme(plot.title = element_text(hjust = hjust,size=title.size,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
    } else {
        plots = plots + theme(legend.text = element_text(size=legend.size))+
            guides(colour = guide_legend(override.aes = list(size=legend.size/5))) 
    }
    if (!is.null(legend.title)) plots = plots + labs(fill = legend.title)
    if(border) plots = plots + NoAxesLabel() + NoGrid()+
        theme(panel.border = element_rect(colour = "black"))

    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path, "/", file.name), 
             units=units, width=width, height=height,res=600)
        print(plots)
        dev.off()
    }
    #if(do.return & Sys.info()[['sysname']] != "Linux") return(plots)
    if(do.return) return(plots)
}


#' prepare exp and tsne file
#' @param split.by split objecty by. Colname of meta.data. Set to FALSE if want to use all 
#' @example samples =  c("All_samples","nt", "hgg","lgg"),
#' Rshiny_path <- "Rshiny/Malignant_Transformation/"
#' PrepareShiny(object, samples, Rshiny_path, verbose = T)
#' 
PrepareShiny <- function(object, samples, Rshiny_path, split.by = "orig.ident",reduction = "tsne",
                         verbose = F, assay = NULL){
        if(missing(object) | class(object) != "Seurat") stop("samples is not provided")
        if(missing(samples)) stop("samples is not provided")
        if(missing(Rshiny_path)) stop("Rshiny_path is not provided")
        assay = DefaultAssay(object) %||% assay
        Idents(object) <-  split.by
        avaible_samples <- samples %in% c("All_samples",object@meta.data[,split.by])
        if (!all(avaible_samples))
                stop(paste(paste(samples[!avaible_samples],collapse = " "),
                           "are not exist in the data."))
        object <- subset(object, idents = samples[-which("All_samples" %in% samples)])
        
        # prepare max_exp
        max_exp = rowMax(object[[assay]]@data) %>% as.vector()
        max_exp = max_exp/log(2)
        names(max_exp) = rownames(object)
        
        exp <- list()
        tsne <- list()
        for (i in seq_along(samples)){
                sample <- samples[i]
                if(sample == "All_samples") {
                        single_object <- object
                } else single_object <- subset(object, idents = sample)
                #============== exp csv===============
                data <- GetAssayData(single_object)
                data <- as(data, "sparseMatrix")
                data = data/log(2)
                #bad <- rowMax(data) == 0
                #data = data[!bad,]
                
                if(verbose) {
                        print(sample)
                        print(format(object.size(data),units="MB"))
                }
                exp[[i]] = data
                #============== tsne csv===============
                tsne[[i]] = Embeddings(single_object, reduction = reduction)
                
                svMisc::progress(i/length(samples)*100)
        }
        names(exp) = samples
        names(tsne) = samples
        
        all_genes = sort(rownames(object[[assay]]@data))
        shiny_data_path <- paste0(Rshiny_path, "data/")
        if(!dir.exists(shiny_data_path)) dir.create(shiny_data_path, recursive = T)
        save(exp,tsne,max_exp,all_genes, file = paste0(shiny_data_path,basename(Rshiny_path),".Rda"))
}



TitleCenter <- function(size = 15, hjust = 0.5){
    theme(text = element_text(size=size),
          plot.title = element_text(hjust = hjust))
}


#' Modified TSNEPlot
#' @param label.repel
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
#' @export save.path folder to save
TSNEPlot.1 <- function(object,dims = c(1, 2),cells = NULL, cols = NULL, pt.size = NULL,
                       reduction = "tsne",group.by = NULL,split.by = NULL,shape.by = NULL,
                       order = NULL,label = FALSE,label.repel= FALSE, label.size = 4,
                       repel = TRUE,alpha = 0.85, text.size = 15,title.size = 15,
                       legend.size = 15,
                       cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
                       na.value = 'gray97',combine = TRUE,ncol = NULL,title = NULL,legend.title = NULL,
                       no.legend = F,do.print = F,do.return = T,unique.name=F,
                       units= "in",width=10, height=7,hjust = 0.5,border = FALSE,
                       save.path = NULL, file.name = NULL,...) {
    
    if(is.null(file.name)){
        VarName <- UniqueName(object,fileName = deparse(substitute(object)),unique.name = unique.name)
        VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object))
        if(!no.legend) VarName = paste0(VarName, "_Legend")
        VarName = paste0(VarName, ifelse(!is.null(split.by),
                                         yes = paste0("_",split.by),
                                         no =""))
        L = ifelse(no.legend, "", "_Legend")
        if(!label) L = paste0(L,"_noLabel")
        file.name = paste0("TSNEPlot_",VarName,L,".jpeg")
    }
    
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    cols = cols %||% ExtractMetaColor(object, group.by = group.by)
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$tsne@cell.embeddings[cells,dims]
    data <- as.data.frame(x = data)
    dims <- paste0(object@reductions$tsne@key, dims)
    object[["ident"]] <- Idents(object = object)
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- Seurat:::SingleDimPlot(data = data[, c(dims, x, split.by, shape.by)],
                                       dims = dims, col.by = x, cols = cols,
                                       pt.size = pt.size, shape.by = shape.by, order = order,
                                       label = FALSE, cells.highlight = cells.highlight,
                                       cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                                       na.value = na.value,...)
        if (label & label.repel == F ) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
                                  size = label.size, split.by = split.by)
        }
        if (label & label.repel) {
            plot <- LabelRepel(plot = plot, id = x, repel = repel,
                               size = label.size, split.by = split.by, color= cols,
                               alpha = alpha)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + theme(strip.background = element_blank(),
                                 strip.text = element_text(face="plain",size=text.size))+
                facet_wrap(facets = vars(!!sym(x = split.by)),
                           ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                               length(x = unique(x = data[, split.by]))
                           } else ncol )
            
        }
        if(border == TRUE) plot = plot + NoAxesLabel() +
                theme(panel.border = element_rect(colour = "black"))
        
        return(plot)
    })
    if (combine) {
        plots <- patchwork::wrap_plots(plots = plots, ncol = if (!is.null(x = split.by) &&
                                                                 length(x = group.by) > 1) {
            1
        }
        else {
            ncol
        }, ...)
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
            theme(plot.title = element_text(hjust = hjust,size=title.size,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
        
    } else {
        plots = plots + theme(legend.text = element_text(size=legend.size))+
            guides(colour = guide_legend(override.aes = list(size=legend.size/5)))
    }
    
    if (!is.null(legend.title)) plots = plots + labs(fill = legend.title)
    if(border) plots = plots + NoAxesLabel() + NoGrid()+
        theme(panel.border = element_rect(colour = "black"))
    
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path, "/", file.name), units=units, width=width, height=height,res=600)
        print(plots)
        dev.off()
    }
    #if(do.return & Sys.info()[['sysname']] != "Linux") return(plots)
    if(do.return) return(plots)
}


# VolcanoPlots to demonstrate Differential expressed genes
# https://zhuanlan.zhihu.com/p/82785739?utm_source=ZHShareTargetIDMore&utm_medium=social&utm_oi=642996063045423104
VolcanoPlots <- function(data, cut_off = c("p_val_adj","p_val"), cut_off_value = 0.05, cut_off_logFC = 0.25,top = 15,
                         sort.by = "p_val_adj",
                         cols = c("#ba2832","#d2dae2","#2a71b2"),
                         cols.order = c('Upregulated','Stable','Downregulated'),
                         alpha=0.8, size=2,
                         legend.size = 12, legend.position = "bottom", ...) {
    data[,paste0("log10_",cut_off[1])] = -log10(data[,cut_off[1]])
    data$change = ifelse(data[,cut_off[1]] < cut_off_value &
                             abs(data$avg_logFC) >= cut_off_logFC, 
                         ifelse(data$avg_logFC > cut_off_logFC ,'Upregulated','Downregulated'),
                         'Stable')
    cols.order = switch (legend.position,
                         "bottom" = rev(cols.order),
                         "right" = cols.order
    )
    cols = switch (legend.position,
                 "bottom" = rev(cols),
                 "right" = cols.order
    )
    if(!is.null(cols.order)) data$change %<>% as.factor() %>% factor(levels = cols.order)

    colnames(data)[grep("cluster",colnames(data))]="cluster"
    # 将需要标记的基因放置在单独的数组
    Up <- data[data$change %in% "Upregulated",]
    Down <- data[data$change %in% "Downregulated",]
    if(sort.by == "p_val_adj") {
        Up_gene_index <- rownames(Up)[Up[,sort.by] <= tail(head(sort(Up[,sort.by],decreasing = F),top),1)]
        Down_gene_index <- rownames(Down)[Down[,sort.by] <= tail(head(sort(Down[,sort.by],decreasing = F),top),1)]
    }
    if(sort.by == "avg_logFC") {
        Up_gene_index <- rownames(Up)[Up[,sort.by] >= tail(head(sort(Up[,sort.by],decreasing = T),top),1)]
        Down_gene_index <- rownames(Down)[Down[,sort.by] <= tail(head(sort(Down[,sort.by],decreasing = F),top),1)]
    }
    p<-ggplot(
        #设置数据
        data, 
        mapping = aes_string(x = "avg_logFC", 
                             y = paste0("log10_",cut_off[1]),
                             fill = "change"))+
       geom_point(mapping = aes_string(color = "change"), alpha=alpha, size=size,...)
        # 辅助线
    p = p + geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8)
    p = p + geom_hline(yintercept = -log10(cut_off_value),lty=4,col="black",lwd=0.8)
        
        # 坐标轴
    p = p + theme_bw()+
        labs(x="log2(fold change)",
             y= paste("-log10 (",ifelse(cut_off[1] == "p_val_adj", "adjusted p-value","p-value"),")"))+
        
        # 图例
        theme(plot.title = element_text(hjust = 0.5), 
              axis.title=element_text(size=12),
              legend.position=legend.position, 
              legend.title = element_blank(),
              legend.text = element_text(size = legend.size),
        )
    if(!is.null(cols.order)) names(cols) = cols.order
    
    p = p + ggrepel::geom_text_repel(data = data[c(Down_gene_index, Up_gene_index),], 
                                     aes(label = gene),
                                     size = size,box.padding = unit(0.5, "lines"),
                                     point.padding = unit(0.8, "lines"), 
                                     segment.color = "black", 
                                     show.legend = FALSE)
    #p = p + scale_colour_manual(values=cols[unique(data$change)])
    p = p + scale_fill_manual(values=cols[unique(data$change)])
    return(p)
}


#' Modified UMAPPlot
#' @param label.repel
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
UMAPPlot.1 <- function(object,dims = c(1, 2),cells = NULL,cols = NULL, pt.size = NULL,
                       reduction = "tsne",group.by = NULL,split.by = NULL,shape.by = NULL,
                       order = NULL,label = FALSE,label.repel= FALSE, label.size = 4,
                       repel = TRUE,alpha = 1, text.size = 15,title.size = 15,
                       legend.size = 15,
                       cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
                       na.value = 'gray97',combine = TRUE,ncol = NULL,title = NULL,legend.title = NULL,
                       no.legend = F,do.print = F,do.return = T,unique.name=F,
                       units= "in",width=10, height=7,hjust = 0.5,border = FALSE,
                       save.path = NULL,file.name = NULL,...) {
    if(is.null(file.name)){
        VarName <- UniqueName(object,fileName = deparse(substitute(object)),unique.name = unique.name)
        VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object))
        if(!no.legend) VarName = paste0(VarName, "_Legend")
        VarName = paste0(VarName, ifelse(!is.null(split.by),
                                         yes = paste0("_",split.by),
                                         no =""))
        L = ifelse(no.legend, "", "_Legend")
        if(!label) L = paste0(L,"_noLabel")
        file.name = paste0("UMAPPlot_",VarName,L,".jpeg")
    }
    
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    cols = cols %||% ExtractMetaColor(object, group.by = group.by)
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$umap@cell.embeddings[cells,dims]
    data <- as.data.frame(x = data)
    dims <- paste0(object@reductions$umap@key, dims)
    object[["ident"]] <- Idents(object = object)
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- Seurat:::SingleDimPlot(data = data[, c(dims, x, split.by, shape.by)],
                                       dims = dims, col.by = x, cols = cols,
                                       pt.size = pt.size, shape.by = shape.by, order = order,
                                       label = FALSE, cells.highlight = cells.highlight,
                                       cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                                       na.value = na.value,...)
        if (label & label.repel == F ) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
                                  size = label.size, split.by = split.by)
        }
        if (label & label.repel) {
            plot <- LabelRepel(plot = plot, id = x, repel = repel,
                               size = label.size, split.by = split.by, color= cols,
                               alpha = alpha)
        }
        if(alpha != 1 ) plot <- plot + scale_fill_manual(values = alpha(.3))
        if (!is.null(x = split.by)) {
            plot <- plot + theme(strip.background = element_blank(),
                                 strip.text = element_text(face="plain",size=text.size))+
                facet_wrap(facets = vars(!!sym(x = split.by)),
                           ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                               length(x = unique(x = data[, split.by]))
                           } else ncol )
            
        }
        if(border == TRUE) plot = plot + NoAxesLabel() +
                theme(panel.border = element_rect(colour = "black"))
        return(plot)
    })
    if (combine) {
        plots <- patchwork::wrap_plots(plots = plots, ncol = if (!is.null(x = split.by) &&
                                                                 length(x = group.by) > 1) {
            1
        }
        else {
            ncol
        }, ...)
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
            theme(plot.title = element_text(hjust = hjust,size=title.size,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
    } else {
        plots = plots + theme(legend.text = element_text(size=legend.size))+
            guides(colour = guide_legend(override.aes = list(size=legend.size/5)))
    }
    if (!is.null(legend.title)) plots = plots + labs(fill = legend.title)
    if(border) plots = plots + NoAxesLabel() + NoGrid()+
        theme(panel.border = element_rect(colour = "black"))
    
    if(do.print) {
        if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        jpeg(paste0(save.path, "/", file.name), units=units, width=width, height=height,res=600)
        print(plots)
        dev.off()
    }
    #if(do.return & Sys.info()[['sysname']] != "Linux") return(plots)
    if(do.return) return(plots)
}


Singler.colors <- c("#7FC97F","#BEAED4","#FDC086","#386CB0","#F0027F",
                    "#BF5B17","#666666","#1B9E77","#7570B3","#66A61E",
                    "#E6AB02","#A6761D","#A6CEE3","#B2DF8A","#FB9A99",
                    "#E31A1C","#FF7F00","#6A3D9A","#8DA0CB",
                    "#4DAF4A","#984EA3","#c6c386","#999999","#66C2A5",
                    "#FC8D62","#A6D854","#FFD92F","#BEBADA",
                    "#FB8072","#80B1D3","#FDB462","#BC80BD","#B3B3B3",
                    "#33A02C","#B3DE69","#4038b0","#ee7576","#e94749","#E78AC3","#ff0000",
                    "#A65628","#d80172","#F781BF","#D95F02","#E7298A",
                    "#1F78B4","#FDBF6F","#CAB2D6","#B15928","#FBB4AE",
                    "#B3CDE3",
                    '#0173b2','#de8f05','#029e73','#d55e00','#cc78bc','#ca9161','#fbafe4','#949494','#ece133','#56b4e9', # seaborn.color_palette colorblind
                    "#00AFBB", "#E7B800", "#FC4E07",
                    "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
sns.RdBu_r_199 = c('#063264', '#073467', '#08366a', '#0a3b70', '#0c3d73', '#0d3f76', '#0e4179',
               '#10457e', '#114781', '#124984', '#144e8a', '#15508d', '#175290', '#185493',
               '#1a5899', '#1b5a9c', '#1c5c9f', '#1e61a5', '#1f63a8', '#2065ab', '#2267ac',
               '#246aae', '#266caf', '#276eb0', '#2a71b2', '#2b73b3', '#2c75b4', '#2e77b5',
               '#307ab6', '#327cb7', '#337eb8', '#3480b9', '#3783bb', '#3885bc', '#3a87bd',
               '#3c8abe', '#3e8cbf', '#3f8ec0', '#408fc1', '#4393c3', '#4695c4', '#4997c5', 
               '#4f9bc7', '#529dc8', '#569fc9', '#59a1ca', '#5fa5cd', '#62a7ce', '#65a9cf',
               '#6bacd1', '#6eaed2', '#71b0d3', '#75b2d4', '#7bb6d6', '#7eb8d7', '#81bad8',
               '#84bcd9', '#8ac0db', '#8dc2dc', '#90c4dd', '#96c7df', '#98c8e0', '#9bc9e0',
               '#9dcbe1', '#a2cde3', '#a5cee3', '#a7d0e4', '#acd2e5', '#aed3e6', '#b1d5e7',
               '#b3d6e8', '#b8d8e9', '#bbdaea', '#bddbea', '#c2ddec', '#c5dfec', '#c7e0ed',
               '#cae1ee', '#cfe4ef', '#d1e5f0', '#d2e6f0', '#d4e6f1', '#d7e8f1', '#d8e9f1',
               '#dae9f2', '#ddebf2', '#deebf2', '#e0ecf3', '#e1edf3', '#e4eef4', '#e6eff4',
               '#e7f0f4', '#eaf1f5', '#ecf2f5', '#edf2f5', '#eff3f5', '#f2f5f6', '#f3f5f6',
               '#f5f6f7', '#f7f6f6', '#f7f5f4', '#f8f4f2', '#f8f3f0', '#f8f1ed', '#f9f0eb',
               '#f9efe9', '#f9eee7', '#f9ebe3', '#faeae1', '#fae9df', '#fae7dc', '#fbe6da',
               '#fbe5d8', '#fbe4d6', '#fce2d2', '#fce0d0', '#fcdfcf', '#fdddcb', '#fddcc9',
               '#fddbc7', '#fdd9c4', '#fcd5bf', '#fcd3bc', '#fbd0b9', '#fbccb4', '#facab1',
               '#fac8af', '#f9c6ac', '#f9c2a7', '#f8bfa4', '#f8bda1', '#f8bb9e', '#f7b799',
               '#f7b596', '#f6b394', '#f6af8e', '#f5ac8b', '#f5aa89', '#f5a886', '#f3a481',
               '#f2a17f', '#f19e7d', '#ef9979', '#ee9677', '#ec9374', '#eb9172', '#e98b6e',
               '#e8896c', '#e6866a', '#e48066', '#e37e64', '#e27b62', '#e17860', '#de735c',
               '#dd7059', '#dc6e57', '#db6b55', '#d86551', '#d7634f', '#d6604d', '#d35a4a',
               '#d25849', '#d05548', '#cf5246', '#cc4c44', '#cb4942', '#c94741', '#c6413e',
               '#c53e3d', '#c43b3c', '#c2383a', '#bf3338', '#be3036', '#bd2d35', '#ba2832',
               '#b82531', '#b72230', '#b61f2e', '#b3192c', '#b1182b', '#ae172a', '#ab162a',
               '#a51429', '#a21328', '#9f1228', '#991027', '#960f27', '#930e26', '#900d26',
               '#8a0b25', '#870a24', '#840924', '#7f0823', '#7c0722', '#790622', '#760521',
               '#700320', '#6d0220', '#6a011f' )
# maturation score colors
my.cols.RYG <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b",
"#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"))(11)
sns.RdBu_r_15 = c('#175290', '#2a71b2', '#3f8ec0', '#6bacd1', '#9bc9e0', '#c2ddec', '#e0ecf3',
                  '#f7f6f6', '#fbe5d8', '#fbccb4', '#f5aa89', '#e48066', '#d05548', '#ba2832', '#930e26')
#' produce unique string baes on Seurat object's meta data
UniqueName <- function(object, fileName = NULL, unique.name = T){
    if(is.logical(unique.name)){
        if(unique.name) {
            idents <- base::as.character(unique(object$orig.ident))
            fileName <- paste(c(fileName,idents[1:min(5,length(idents))]),collapse = "_")
        } 
    } else if(unique.name %in% colnames(object@meta.data)) {
        fileName <- paste(c(fileName, base::as.character(unique(object@meta.data[,unique.name]))),collapse = "_")
    }
    return(fileName)
}

# unlist and use.names without additional index in the names
Unlist <- function(List){
    for(i in seq_along(List)){
        names(List[[i]]) = rep(names(List)[i],length(List[[i]]))
    }
    flatten_List <- base::unlist(List,use.names=TRUE)
    names(flatten_List) %<>% gsub("\\..*","",.)
    return(flatten_List)
}


