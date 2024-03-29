# TRUE Global Environment
#' modified function and pass ref.list arg to CreateSinglerObject
#' @param N sub sample size
CreateBigSingleRObject.1 <- function (counts, annot = NULL, project.name, xy = NULL, clusters = NULL, 
                                    N = 10000, min.genes = 200, technology = "10X", species = "Human", 
                                    citation = "", ref.list = list(), normalize.gene.length = F, 
                                    variable.genes = "de", fine.tune = T, reduce.file.size = T, 
                                    do.signatures = F, do.main.types = T, temp.dir = getwd(), 
                                    numCores = SingleR.numCores) 
{
    n = ncol(counts)
    s = seq(1, n, by = N)
    dir.create(paste0(temp.dir, "/singler.temp/"), showWarnings = FALSE)
    for (i in s) {
        print(i)
        A = seq(i, min(i + N - 1, n))
        singler = CreateSinglerObject(counts[, A], annot = annot[A], 
                                      project.name = project.name, min.genes = min.genes, 
                                      ref.list =ref.list, normalize.gene.length = normalize.gene.length,
                                      technology = technology, species = species, citation = citation,
                                      variable.genes =variable.genes, fine.tune = fine.tune, 
                                      do.signatures = do.signatures, do.main.types = do.main.types,
                                      clusters = NULL, numCores = numCores)
        save(singler, file = paste0(temp.dir, "/singler.temp/", 
                                    project.name, ".", i, ".RData"))
    }
    singler.objects.file <- list.files(paste0(temp.dir, "/singler.temp/"), 
                                       pattern = "RData", full.names = T)
    singler.objects = list()
    for (i in 1:length(singler.objects.file)) {
        load(singler.objects.file[[i]])
        singler.objects[[i]] = singler
    }
    singler = SingleR.Combine(singler.objects, order = NULL, 
                              clusters = clusters, xy = xy)
    singler
}


# TRUE Global Environment
#' modified SCINA function to split big data into smaller size
#' @param N split into smaller sample subsets. Suggest > 10000 to avoid null signature

BigSCINA <- function (exp, signatures, N = 10000, max_iter = 100, convergence_n = 10, 
                      convergence_rate = 0.99, sensitivity_cutoff = 1, rm_overlap = 1, 
                      allow_unknown = 1, log_file = "SCINA.log") 
{
    n = ncol(exp)
    N1 = N + round((n %% N)/(n %/% N))+1 # increase the N by add back rest of remainder plus 1
    s = seq(1, n, by = N1)
    scina_list <- list()
    for (i in seq_along(s)) {
        A = seq(s[i], min(s[i] + N1 - 1, n))
        print(paste(head(A,1),"<->",tail(A,1)))
        print(length(A))
        system.time(scina_list[[i]] <- SCINA.1(exp = exp[, A], signatures, max_iter = max_iter, convergence_n = convergence_n, 
                          convergence_rate = convergence_rate, sensitivity_cutoff = sensitivity_cutoff,
                          rm_overlap = rm_overlap, allow_unknown = allow_unknown, 
                          log_file = log_file))
        Progress(i, length(s))
    }
    scina <- list(cell_labels = NULL, 
                  probabilities = matrix(1:length(signatures),
                                         dimnames = list(names(signatures),
                                                         "gene")))
    message("merging scina results")
    for(m in seq_along(scina_list)){
        scina$cell_labels = c(scina$cell_labels, scina_list[[m]]$cell_labels)
        scina$probabilities = merge(scina$probabilities, 
                                    scina_list[[m]]$probabilities,
                                    by = "row.names", all = TRUE)
        rownames(scina$probabilities) = scina$probabilities$Row.names
        scina$probabilities = scina$probabilities[,-grep("Row.names|gene",colnames(scina$probabilities))]
        Progress(m, length(s))
    }
    return(scina)
}

# modify SCINA to remove empty signatures. 
# SCINA already add this part in commit 88fda75 on Feb 5, 202
# But I can't install that version in linux
SCINA.1 <- function (exp, signatures, max_iter = 100, convergence_n = 10, 
                     convergence_rate = 0.99, sensitivity_cutoff = 1, rm_overlap = 1, 
                     allow_unknown = 1, log_file = "SCINA.log") 
{
        cat("Start running SCINA.", file = log_file, append = F)
        cat("\n", file = log_file, append = T)
        all_sig = unique(unlist(signatures))
        invert_sigs = grep("^low_", all_sig, value = T)
        if (!identical(invert_sigs, character(0))) {
                cat("Converting expression matrix for low_genes.", file = log_file, 
                    append = T)
                cat("\n", file = log_file, append = T)
                invert_sigs_2add = unlist(lapply(invert_sigs, function(x) strsplit(x, 
                                                                                   "_")[[1]][2]))
                invert_sigs = invert_sigs[invert_sigs_2add %in% row.names(exp)]
                invert_sigs_2add = invert_sigs_2add[invert_sigs_2add %in% 
                                                            row.names(exp)]
                sub_exp = -exp[invert_sigs_2add, , drop = F]
                row.names(sub_exp) = invert_sigs
                exp = rbind(exp, sub_exp)
                rm(sub_exp, all_sig, invert_sigs, invert_sigs_2add)
        }
        quality = SCINA:::check.inputs(exp, signatures, max_iter, convergence_n, 
                                       convergence_rate, sensitivity_cutoff, rm_overlap, log_file)
        if (quality$qual == 0) {
                cat("EXITING due to invalid parameters.", file = log_file, 
                    append = T)
                cat("\n", file = log_file, append = T)
                stop("SCINA stopped.")
        }
        signatures = quality$sig
        max_iter = quality$para[1]
        convergence_n = quality$para[2]
        convergence_rate = quality$para[3]
        sensitivity_cutoff = quality$para[4]
        exp = as.matrix(exp)
        exp = exp[unlist(signatures), , drop = F]
        labels = matrix(0, ncol = convergence_n, nrow = dim(exp)[2])
        unsatisfied = 1
        theta = list()
        for (i in 1:length(signatures)) {
                theta[[i]] = list()
                theta[[i]]$mean = t(apply(exp[signatures[[i]], , drop = F], 
                                          1, function(x) quantile(x, c(0.7, 0.3))))
                tmp = apply(exp[signatures[[i]], , drop = F], 1, var)
                theta[[i]]$sigma1 = diag(tmp, ncol = length(tmp))
                theta[[i]]$sigma2 = theta[[i]]$sigma1
        }
        empty = c()
        for (i in 1:length(theta)) {
                if(dim(theta[[i]]$sigma1)[1] == 0) empty = c(empty,i)
        }
        theta[empty] = NULL
        signatures[empty] = NULL
        
        if (allow_unknown == 1) {
                tao = rep(1/(length(signatures) + 1), length(signatures))
        } else {
                tao = rep(1/(length(signatures)), length(signatures))
        }
        sigma_min = min(sapply(theta, function(x) min(c(diag(x$sigma1), 
                                                        diag(x$sigma2)))))/100
        remove_times = 0
        while (unsatisfied == 1) {
                prob_mat = matrix(tao, ncol = dim(exp)[2], nrow = length(tao))
                row.names(prob_mat) = names(signatures)
                iter = 0
                labels_i = 1
                remove_times = remove_times + 1
                while (iter < max_iter) {
                        iter = iter + 1
                        for (i in 1:length(signatures)) {
                                #if(i %in% empty) next
                                theta[[i]]$inverse_sigma1 = theta[[i]]$inverse_sigma2 = chol2inv(chol(theta[[i]]$sigma1))
                        }
                        for (r in 1:dim(prob_mat)[1]) {
                                #if(r %in% empty) next
                                prob_mat[r, ] = tao[r] * SCINA:::density_ratio(exp[signatures[[r]], 
                                                                                   , drop = F], theta[[r]]$mean[, 1], theta[[r]]$mean[, 
                                                                                                                                      2], theta[[r]]$inverse_sigma1, theta[[r]]$inverse_sigma2)
                        }
                        prob_mat = t(t(prob_mat)/(1 - sum(tao) + colSums(prob_mat)))
                        tao = rowMeans(prob_mat)
                        for (i in 1:length(signatures)) {
                                #if(i %in% empty) next
                                
                                theta[[i]]$mean[, 1] = (exp[signatures[[i]], 
                                                            ] %*% prob_mat[i, ])/sum(prob_mat[i, ])
                                theta[[i]]$mean[, 2] = (exp[signatures[[i]], 
                                                            ] %*% (1 - prob_mat[i, ]))/sum(1 - prob_mat[i, ])
                                keep = theta[[i]]$mean[, 1] <= theta[[i]]$mean[, 2]
                                if (any(keep)) {
                                        theta[[i]]$mean[keep, 1] = rowMeans(exp[signatures[[i]][keep], , drop = F])
                                        theta[[i]]$mean[keep, 2] = theta[[i]]$mean[keep, 1]
                                }
                                tmp1 = t((exp[signatures[[i]], , drop = F] - theta[[i]]$mean[, 1])^2)
                                tmp2 = t((exp[signatures[[i]], , drop = F] - theta[[i]]$mean[, 2])^2)
                                diag(theta[[i]]$sigma1) = diag(theta[[i]]$sigma2) = colSums(tmp1 * prob_mat[i, ] + tmp2 * (1 - prob_mat[i, ]))/dim(prob_mat)[2]
                                diag(theta[[i]]$sigma1)[diag(theta[[i]]$sigma1) < 
                                                                sigma_min] = sigma_min
                                diag(theta[[i]]$sigma2)[diag(theta[[i]]$sigma2) < 
                                                                sigma_min] = sigma_min
                        }
                        labels[, labels_i] = apply(rbind(1 - colSums(prob_mat), 
                                                         prob_mat), 2, which.max) - 1
                        if (mean(apply(labels, 1, function(x) length(unique(x)) == 
                                       1)) >= convergence_rate) {
                                cat("Job finished successfully.", file = log_file, 
                                    append = T)
                                cat("\n", file = log_file, append = T)
                                break
                        }
                        labels_i = labels_i + 1
                        if (labels_i == convergence_n + 1) {
                                labels_i = 1
                        }
                        if (iter == max_iter) {
                                cat("Maximum iterations, breaking out.", file = log_file, 
                                    append = T)
                                cat("\n", file = log_file, append = T)
                        }
                }
                colnames(prob_mat) = colnames(exp)
                row.names(prob_mat) = names(signatures)
                row.names(labels) = colnames(exp)
                dummytest = sapply(1:length(signatures), function(i) mean(theta[[i]]$mean[, 1] - theta[[i]]$mean[, 2] == 0))
                if (all(dummytest <= sensitivity_cutoff)) {
                        unsatisfied = 0
                }
                else {
                        rev = which(dummytest > sensitivity_cutoff)
                        cat(paste("Remove dummy signatures:", rev, sep = " "), 
                            file = log_file, append = T, sep = "\n")
                        signatures = signatures[-rev]
                        tmp = 1 - sum(tao)
                        tao = tao[-rev]
                        tao = tao/(tmp + sum(tao))
                        theta = theta[-rev]
                }
        }
        return(list(cell_labels = c("unknown", names(signatures))[1 +labels[, labels_i]], probabilities = prob_mat))
}

CreatGeneSetsFromSingleR <- function(object = object, cell.type = "B_cell",main.type = FALSE,
                                     species = "Human"){
        "Find Marker genes list for one cell type"
        #=== Create Score matrix======
        if(main.type) {
                de_gene_mode = "de.genes.main"; types_mode = "main_types"
        } else {
                de_gene_mode = "de.genes"; types_mode = "types"
        }
        genes = unique(unlist(object[[de_gene_mode]][[cell.type]]))
        types = unique(object[[types_mode]])
        S_matrix = matrix(data = genes, ncol = 1,
                          dimnames = list(NULL, "gene"))
        #=== Calculate ranking score=====
        for (object_type in types){
                Gen <- object[[de_gene_mode]][[cell.type]][[object_type]]
                Gen_mat = matrix(data = c(Gen,length(Gen):1), ncol = 2,
                                 dimnames = list(NULL, c("gene",object_type)))
                S_matrix =  merge(S_matrix,Gen_mat, by="gene",all = T)
        }
        
        S_matrix = S_matrix[!is.na(S_matrix$gene),]
        # Store the gene name
        if(tolower(species) == "human") Genes = toupper(S_matrix$gene)
        if(tolower(species) == "mouse") Genes = Hmisc::capitalize(tolower(S_matrix$gene))
        S_matrix$gene = 0
        S_matrix = as.matrix(S_matrix)
        S_matrix[is.na(S_matrix)] <- 0
        S_matrix = apply(S_matrix,2, as.numeric)
        S_matrix[,"gene"] = rowSums(S_matrix)
        rownames(S_matrix) = Genes
        markers <-rownames(S_matrix)[order(S_matrix[,"gene"],decreasing = T)]
        
        return(markers)
}


CreatGenePanelFromSingleR <- function(object = object, main.type = TRUE, species = "human"){
        "Find Marker genes list for All cell types"
        if(main.type) {
                types = unique(object$main_types)
        } else { types = unique(object$types)
        }
        marker_list = lapply(types, function(x) {
                CreatGeneSetsFromSingleR(object = object, cell.type = x, main.type = main.type,
                                         species = species)
        })
        names(marker_list) <- types
        marker_df <- list2df(marker_list)
        # arrange columns
        marker_df= marker_df[,order(colnames(marker_df))]
        return(marker_df)
}


CreateSinglerReference <- function(name, expr, types, main_types,
de.num = 200,de.main.num = 300){
    ref = list(name=name, data = expr, types=types, main_types=main_types)
    
    # if using the de method, we can predefine the variable genes
    ref$de.genes = CreateVariableGeneSet(expr,types,de.num)
    ref$de.genes.main = CreateVariableGeneSet(expr,main_types,de.main.num)
    
    # if using the sd method, we need to define an sd threshold
    sd = rowSds(expr)
    sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
    ref$sd.thres = sd.thres
    
    return(ref)
    
}

CreateVariableGeneSet <- function (ref_data, types, n) 
{
        mat = medianMatrix(ref_data, types)
        genes = lapply(1:ncol(mat), function(j) {
                lapply(1:ncol(mat), function(i) {
                        s = sort(mat[, j] - mat[, i], decreasing = T)
                        s = s[s > 0]
                        names(s)[1:min(n, length(s))]
                })
        })
        names(genes) = colnames(mat)
        for (i in 1:length(genes)) {
                names(genes[[i]]) = colnames(mat)
        }
        genes
}


FineTune <- function(x, main.type = FALSE){
        # for both main types and sub-types
        x = gsub(" ","_",x)
        x = gsub("B-cells","B_cells",x)
        x = gsub("cells","cell",x)
        x = gsub("cell","cells",x)
        x = gsub("cytes","cyte",x)
        x = gsub("cyte","cytes",x)
        x = gsub("blasts","blast",x)
        x = gsub("blast","blasts",x)
        x = gsub("Macrophages","Macrophage",x)
        x = gsub("Macrophage","Macrophages",x)
        x = gsub("Neutrophils","Neutrophil",x)
        x = gsub("Neutrophil","Neutrophils",x)
        x = gsub("Smooth_muscle_cells","Smooth_muscle",x)
        x = gsub("Smooth_muscle","Smooth_muscle_cells",x)
        if(main.type){
                x = gsub("CD4\\+_T-cells","T_cells",x)
                x = gsub("CD8\\+_T-cells","T_cells",x)
                x = gsub("HSC_-G-CSF","HSC",x)
                x = gsub("HSC_CD34\\+","HSC",x)
                x = gsub("Memory_B_cells","B_cells",x)
                x = gsub("naive B-cells","B_cells",x)
                x = gsub("Pre-B_cells_CD34-","B_cells",x)
                x = gsub("Pro-B_cells_CD34\\+","B_cells",x)
                x = gsub("Pro-Myelocytes","Myelocytes",x) 
                x = gsub("Tregs","T_cells",x)
                
        } else {
                x = gsub("CD14\\+_Monocytes","Monocytes:CD14\\+",x)
                x = gsub("CD4\\+_T-cells","T_cells:CD4\\+",x)
                x = gsub("CD4\\+_T_cells","T_cells:CD4\\+",x)
                x = gsub("CD8\\+_T-cells","T_cells:CD8\\+",x)
                x = gsub("CD8\\+_T_cells","T_cells:CD8\\+",x)
                x = gsub("CD4\\+_Tcm","T_cells:CD4\\+",x)
                x = gsub("CD4\\+_Tem","T_cells:CD4\\+",x)
                x = gsub("CD8\\+_Tcm","T_cells:CD8\\+",x)
                x = gsub("CD8\\+_Tem","T_cells:CD8\\+",x)
                x = gsub("Class-switched_memory_B_cells","B_cells:Class-switched_memory",x)
                x = gsub("FCGR3A\\+_Monocytes","Monocytes:FCGR3A\\+",x)
                x = gsub("Macrophages_M1","Macrophages:M1",x)
                x = gsub("Macrophages_M2","Macrophages:M2",x)
                x = gsub("Memory_B_cells","B_cells:Memory",x)
                x = gsub("naive_B_cells","B_cells:Naive_B_cells",x)
                x = gsub("Plasma_cells","B_cells:Plasma_cells",x)
                x = gsub("Pre-B_cells_CD34-","B_cells:Pre-B_cells_CD34-",x)
                x = gsub("Pro-B_cells_CD34\\+","B_cells:Pro-B_cells_CD34\\+",x)
                x = gsub("Tregs","T_cells:Tregs",x)
        }
        return(x)
}


#=====Clean memory======================
GC <- function()
{
    while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1,
                                                         4]) {
    }
}


#' search cell type markers database. 
#' Provide gene names and return the most likely cell types
#' @param df data.frame cell type markers database
#' @param marker gene name
#' @export result_df the most likely cell types and rank
SearchMarker <- function(df, marker){
        result = apply(df, 2, function(x) which(grepl(marker, x))[1])
        result_df = data.frame(sort(unlist(result)))
        colnames(result_df) = marker
        return(result_df)
}

#' search cell type markers database. 
#' Provide gene names and return the most likely cell types
#' @param df data.frame cell type markers database
#' @param markers gene names
#' @export result_df the most likely cell types and rank
#' @example 
# markers_All <- read.csv("output/20190311/markers.All.csv")
# neutril <- FilterGenes(object, markers_All[markers_All$cluster == 4,"gene"][1:20])
# neutril_df <- SearchAllMarkers(Hpca_Blueprint_encode_main,neutril)
# neutril_df["Neutrophils",]
SearchAllMarkers <- function(df, markers){
        results <- lapply(markers,function(x) SearchMarker(df,x))
        names(results) <- markers
        
        # remove empty elements
        idx <- lapply(results, nrow) %>% unlist
        non_empty <- which(idx>0)
    
        # add row.names as the new column
        results_list <- lapply(results[non_empty], function(x) {
                x$cell_type = rownames(x)
                x
        })
        merged_results <- Reduce(function(x, y) merge(x, y, all=TRUE,by="cell_type"),
                                 results_list)
        rownames(merged_results) = merged_results$cell_type
        merged_results = merged_results[,-1]
        merged_results = removeNA(merged_results)
        return(merged_results)
}

#' Plot a heatmap of the scores for all the single cells
#'
#' @param SingleR the output from the SingleR function
#' @param cells.use single cells to present, if NULL all single cells presented
#' @param types.use cell types to present, if NULL all cell types presented
#' @param clusters a clustering to present as annotation in the heatmap
#' @param top.n number of cell types to presents. Default is 40. This can have an effect on the clustering which is performed only on the cell types presented.
#' @param normalize if TRUE scores are normalized to a 0-1 scale.
#' @param order.by.clusters if TRUE columns are ordered by the input clusters, and are not clustered again
#' @param cells_order an input order for the column
#' @param silent if TRUE do not draw the plot  
SingleR.DrawHeatmap = function(SingleR,cells.use = NULL, types.use = NULL,
                               clusters=NULL,top.n=40,normalize=T,
                               order.by.clusters=F,cells_order=NULL,silent=F,
                               fontsize_row=9,...) {
    scores = SingleR$scores
    if (!is.null(cells.use)) {
        scores = scores[cells.use,]
    }
    if (!is.null(types.use)) {
        scores = scores[,types.use]
    }
    
    m = apply(t(scale(t(scores))),2,max)
    
    thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]
    
    data = as.matrix(scores)
    
    if (normalize==T) {
        mmax = rowMaxs(data)
        mmin = rowMins(data)
        data = (data-mmin)/(mmax-mmin)
        data = data^3
    }
    data = data[,m>(thres-1e-6)]
    
    
    data = t(data)
    
    if (!is.null(clusters)) {
        clusters = as.data.frame(clusters)
        colnames(clusters) = 'Clusters'
        rownames(clusters) = colnames(data)
        
    }
    additional_params = list(...)
    if (is.null(additional_params$annotation_colors)) {
        annotation_colors = NA
    } else {
        annotation_colors = additional_params$annotation_colors
    }
    clustering_method = 'ward.D2'
    if (order.by.clusters==T) {
        data = data[,order(clusters$Clusters)]
        clusters = clusters[order(clusters$Clusters),,drop=F]
        pheatmap(data,border_color=NA,show_colnames=FALSE,
                 clustering_method=clustering_method,fontsize_row=fontsize_row,
                 annotation_col = clusters,cluster_cols = F,silent=silent, 
                 annotation_colors=annotation_colors)
    } else if (!is.null(cells_order)) {
        data = data[,cells_order]
        clusters = clusters[cells_order,,drop=F]
        pheatmap(data,border_color=NA,show_colnames=FALSE,
                 clustering_method=clustering_method,fontsize_row=fontsize_row,
                 annotation_col = clusters,cluster_cols = F,silent=silent, 
                 annotation_colors=annotation_colors)
    } else {
        if (!is.null(clusters)) {
            pheatmap(data,border_color=NA,show_colnames=FALSE,
                     clustering_method=clustering_method,fontsize_row=fontsize_row,
                     annotation_col = clusters,silent=silent, 
                     annotation_colors=annotation_colors)
        } else {
            pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
                     clustering_method=clustering_method,fontsize_row=fontsize_row,
                     silent=silent, annotation_colors=annotation_colors)
            
        }
    }
}


SingleR.PlotTsne.1 <- function (SingleR, xy, labels = SingleR$labels, score.thres = 0,
                clusters = NULL, do.letters = TRUE, dot.size = 1, do.labels = FALSE,
                do.legend = TRUE, label.size = 3, title = "", colors = singler.colors,
                font.size = NULL, alpha = 0.5, label.repel = TRUE, text.repel = FALSE, force = 1)
{
    if (do.labels == TRUE)
    do.letters = FALSE
    df = data.frame(row.names = SingleR$cell.names)
    df$x = xy[, 1]
    df$y = xy[, 2]
    if (SingleR$method == "cluster") {
        df$ident = clusters.map.values(clusters, labels)
    }
    else {
        df$ident = labels
    }
    if (score.thres > 0) {
        max.score = apply(SingleR$scores, 1, max)
        df$ident[max.score < score.thres] = "X"
    }
    df$ident = factor(df$ident)
    SYMBOLS = c(LETTERS, tolower(letters), c(0:9))
    df$initIdent = SYMBOLS[as.numeric(df$ident)]
    num.levels = length(levels(df$ident))
    p = ggplot(df, aes(x = x, y = y,color = ident))
    p = p + geom_point(aes(color = ident), size = dot.size,
    alpha = alpha, stroke = 0)
    if (do.letters == TRUE) {
        symbols = SYMBOLS[1:num.levels]
        names(symbols) = lev = levels(df$ident)
        p = p + geom_point(aes(shape = ident), size = 2 * dot.size/5,
        color = "black")
        p = p + scale_shape_manual(values = symbols)
    }
    if (do.labels == TRUE) {
        centers <- df %>% dplyr::group_by(ident) %>% dplyr::summarize(x = median(x),
        y = median(y))
        p = p + geom_point(data = centers, aes(x = x, y = y),
        size = 0, alpha = 0)
        if (label.repel == TRUE) {
            p = p + ggrepel::geom_label_repel(data = centers,
            aes(label = ident),
            size = label.size,
            force = force)
        }
        else if (text.repel == TRUE){
            p = p + ggrepel::geom_text_repel(data = centers,
            aes(label = ident),
            force = force,
            size = label.size,
            color = "black")
        }
        if (label.repel == FALSE & text.repel == FALSE) {
            p = p + geom_text(data = centers,
            aes(label = ident), size = label.size, color = "black")
        }
        p = p + guides(colour = guide_legend(override.aes = list(alpha = 1)))
        x.range = layer_scales(p)$x$range$range
        add_to_x = sum(abs(x.range)) * 0.03
        p = p + xlim(x.range[1] - add_to_x, x.range[2] + add_to_x)
    }
    else {
        if (is.null(font.size)) {
            font.size = 250 * (1/num.levels)
            font.size = max(font.size, 5)
            font.size = min(font.size, 10)
        }
        if (num.levels > 35 & num.levels < 60) {
            p = p + theme(legend.position = "bottom", legend.direction = "vertical",
            legend.text = element_text(size = 6), legend.title = element_blank()) +
            guides(col = guide_legend(ncol = 5, override.aes = list(size = 2,
            alpha = 1)))
        }
        else if (num.levels > 60) {
            p = p + theme(legend.position = "bottom", legend.direction = "vertical",
            legend.text = element_text(size = 6), legend.title = element_blank()) +
            guides(col = guide_legend(ncol = 9, override.aes = list(size = 2,
            alpha = 1)))
        }
        else {
            p = p + theme(legend.text = element_text(size = font.size),
            legend.title = element_blank()) + guides(color = guide_legend(ncol = 1,
            override.aes = list(size = 3, alpha = 1)))
        }
    }
    lev = levels(df$ident)
    cols = colors[1:length(lev)]
    names(cols) = lev
    cols[names(cols) == "X"] = "black"
    p = p + scale_color_manual(values = cols)
    p = p + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle(title)
    if (do.legend == FALSE) {
        p = p + theme(legend.position = "none")
    }
    p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
}


#' Plot a feature on the t-SNE plot, modified for Seurat 3
#'
#' @param SingleR the output from the SingleR function
#' @param seurat the seurat object version 3
#' @param plot.feature the feature to plot - 'MaxScore' for coloring according to the top score per cell, 'nGene' for number of non-zero genes, 'nUMI' for number of UMIs, or a gene name.
#' @param dot.size size of the dot in the plot.
#'
#' @return ggplot2 object
SingleR.PlotFeature = function(SingleR, seurat, plot.feature='nCount_RNA', 
                               dot.size=1, reduction = "tsne", title=NULL) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (packageVersion('Seurat')>=3) {
        xy = seurat@reductions[[reduction]]@cell.embeddings
    } else {
        xy = seurat@dr[reduction]@cell.embeddings
    }
    df = as.data.frame(xy[SingleR$cell.names,])
    if (length(plot.feature)==nrow(df)) {
        df$Feature=plot.feature
        tit = 'Feature'
    }   else if (plot.feature=='MaxScore') {
        df$Feature = apply(SingleR$scores,1,max) 
        tit = 'Max Score'
    }   else if (plot.feature %in% colnames(object@meta.data)) {
        df$Feature = seurat@meta.data[[plot.feature]]
        tit = plot.feature
    }   else {
        DefaultAssay(object)
        df$Feature = seurat@assays[[DefaultAssay(object)]]@data[plot.feature,]
        tit = plot.feature
    }
    if (is.null(title)) {
        title = tit
    }
    dims <- paste0(Key(object = object[[reduction]]), 1:2)
    
    ggplot(df, aes_string(x=dims[1], y=dims[2])) + 
        geom_point(aes_string(color="Feature"), size=dot.size)+
        scale_colour_gradient(low='gray',high='blue')+
        ggtitle(title) + theme_classic()
}


removeNA <- function(df){
        AllNA = apply(df,2,function(x) length(which(!is.na(x))))
        df = df[,AllNA !=0]
        return(df)
}


SplitSingler <- function(singler = singler, split.by = "conditions"){
    "
    split singler object by certein criteria
    
    Inputs
    -------------------
    singler: singler object with seurat
    split.by: the criteria to split
    
    Outputs
    --------------------
    Singler.subsets: list of subseted singler object by certein conditions,
    plus levels of the conditions
    "
    if(is.null(singler$seurat))
    stop("A seurat object must be provided first, add singler$seurat = your_seurat_object")
    cell.subsets <- SplitCells(object = singler$seurat, split.by = split.by)
    
    Singler.subsets <- list()
    for(i in 1:(length(cell.subsets)-1)){
        cell.index <- which(singler$seurat@cell.names %in% cell.subsets[[i]])
        Singler.subsets[[i]] <- SingleR.Subset.1(singler=singler,
        subsetdata =cell.index)
    }
    Singler.subsets[[i+1]] <- cell.subsets[[i+1]] # record conditions in the last return
    return(Singler.subsets)
}

SplitSingleR.PlotTsne <- function(singler = singler, split.by = "conditions",
select.plots = NULL, return.plots = FALSE,
do.label= TRUE,do.letters = FALSE,main =FALSE,
show.2nd =FALSE,label.size = 4, dot.size = 3,
legend.size = NULL,... ){
    "
    split singler by certein criteria, and generate TSNE plot
    
    Inputs
    -------------------
    singler: singler object with seurat
    split.by: the criteria to split
    select.plots：change order to output, such as c(2,1)
    return.data: TRUE/FASLE, return splited ojbect or not.
    show.subtype: TRUE/FASLE, show sub cell type or not.
    
    Outputs
    --------------------
    return.plots: if return.data = TRUE
    "
    object.subsets <- SplitSingler(singler = singler, split.by = "conditions")
    levels <- object.subsets[[length(object.subsets)]]
    
    out <- list()
    p <- list()
    if(is.null(select.plots)) select.plots <- 1:length(levels)
    sp = select.plots
    if(main) main_or_sub = "SingleR.single.main" else main_or_sub = "SingleR.single"
    if(show.2nd) st =2 else st =1
    for(i in 1:length(select.plots)){
        out[[i]] <- SingleR.PlotTsne.1(SingleR= object.subsets[[sp[i]]]$singler[[st]][[main_or_sub]],
        xy = object.subsets[[sp[i]]]$meta.data$xy,
        do.label=do.label,do.letters = do.letters,
        labels = object.subsets[[sp[i]]]$singler[[st]][[main_or_sub]]$labels,
        label.size = label.size, dot.size = dot.size,...)
        out[[i]] = out[[i]] +
        ggtitle(levels[sp[i]])+
        theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5))
        if(!is.null(legend.size))
        out[[i]] = out[[i]] + theme(legend.text = element_text(size=legend.size))
    }
    out <- out[lapply(out,length)>0] # remove NULL element
    if(return.plots) return(out) else print(do.call(plot_grid, out))
}



#' test a matirx's max/median/min of row/column'
testMMM <- function(x,MARGIN = 2) {
    par(mfrow=c(3,1))
    Max <- apply(x,MARGIN,max)
    Median <- apply(x,MARGIN,median)
    Min <- apply(x,MARGIN,min)
    xlim = c(min(x),max(x))
    hist(Max, xlim = xlim)
    hist(Median, xlim = xlim)
    hist(Min, xlim = xlim)
    par(mfrow=c(1,1))
}

#cell.use <- TrueCells(MCL@meta.data[,c("singler1sub","singler1main","singler2sub",
#            "singler2main","kang")], celltype = c("B_cells","MCL"))
TrueCells <- function(df, celltype = c("B_cells","MCL")){
    celltype = paste(paste0(celltype,".*"),collapse = "|")
    cell_id <- apply(df,2,function(cell) cell[grepl(celltype,cell)])
    cell_id <- lapply(cell_id, as.matrix)
    cell_df <- Reduce(function(x, y) merge(x, y, by = "row.names",all.x= TRUE), cell_id)
    return(cell_df[,4])
}
