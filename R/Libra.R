library(DESeq2)
library(Matrix)
library(Seurat)
library(blme)
library(dplyr)
library(edgeR)
library(forcats)
library(glmmTMB)
library(limma)
library(lme4)
library(lmerTest)
library(magrittr)
library(matrixStats)
library(methods)
library(parallel)
library(pbmcapply)
library(purrr)
library(stats)
library(tester)
library(tibble)
library(Libra)
#' Add loop for comparison when more than 2 groups differenial analysis.
#' @examples 
#' data("hagai_toy")
#' hagai_toy@meta.data[hagai_toy$label == "unst" & hagai_toy$replicate == "mouse1","label"] <- "new_label"
#' DE = run_de(hagai_toy)
#' Error in `group_by()`:
#' ! Must group by variables found in `.data`.
#' ✖ Column `cell_type` is not found.
#' 
#' 
#' 
#' Run differential expression
#'
#' Perform differential expression on single-cell data. Libra implements a total
#' of 22 unique differential expression methods that can all be accessed from
#' one function. These methods encompass traditional single-cell methods as well
#' as methods accounting for biological replicate including pseudobulk and
#' mixed model methods. The code for this package has been largely inspired
#' by the Seurat and Muscat packages. Please see the documentation of these
#' packages for further information.
#'
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocle3}, or
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}.
#' @param ident.1 Identity class to define markers for; if NULL, use all groups for comparison
#' @param ident.2 A second identity class for comparison; if NULL, use all other cells for comparison
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.   
#' @param de_family the differential expression family to use. Available options
#' are:
#' \itemize{
#' \item{"singlecell"}: Uses traditionally methods implemented by Seurat to
#' test for differentially expressed genes. These methods do not take biological
#' replicate into account. For single cell methods there are \code{six} options
#' for \code{de_method} that can be used, while no input for
#' \code{de_test} is required:
#' \itemize{
#' \item{"wilcox"}: Wilcoxon Rank-Sum test. The default.
#' \item{"bimod"}: Likelihood ratio test
#' \item{"t"}: Student's t-test
#' \item{"negbinom"}: Negative binomial linear model
#' \item{"LR"}: Logistic regression
#' \item{"MAST"}: MAST (requires installation of the \code{MAST} package).
#' }
#'
#' \item{"pseudobulk"}: These methods first convert the single-cell expression
#' matrix to a so-called 'pseudobulk' matrix by summing counts for each gene
#' within biological replicates, and then performing differential expression
#' using bulk RNA-seq methods. For pseudobulk methods there are \code{six}
#' different methods that can be accessed by combinations of \code{de_method}
#' and \code{de_type}. First specify \code{de_method} as one of the following:
#' \itemize{
#' \item{"edgeR"}: The edgeR method according to
#' Robinson et al, Bioinformatics, 2010. For this method please specify
#' \code{de_type} as either \code{"LRT"} or \code{"QLF"} as the null hypothesis
#' testing approach. See http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
#' for further information. The default.
#' \item{"DESeq2"}: The DESeq2 method according to
#' Love et al, Genome Biology, 2014. For this method please specify
#' \code{de_type} as either \code{"LRT"} or \code{"Wald"} as the null hypothesis
#' testing approach. See https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' for further information.
#' \item{"limma"}: The limma method according to
#' Ritchie et al, Nucleic Acids Research, 2015. For this method please specify
#' \code{de_type} as either \code{"voom"} or \code{"trend"} as the precise
#' normalization and null hypothesis testing approach within limma. See
#' https://bioconductor.org/packages/release/bioc/html/limma.html for
#' further information.
#' }
#' \item{"mixedmodel"}: Mixed model methods also take biological replicate into
#' account by modelling it as a random effect. Please note that these methods
#' are generally extremely computationally intensive. For mixed model methods
#' there are \code{ten} different methods that can be accessed by
#' combinations of \code{de_method} and \code{de_type}.
#' First specify \code{de_method} as one of the following:
#' \itemize{
#' \item{"negbinom"}: Negative binomial generalized linear mixed model.
#' The default.
#' \item{"linear"}: Linear mixed model.
#' \item{"poisson"}: Poisson generalized linear mixed model.
#' \item{"negbinom_offset"}: Negative binomial generalized linear mixed model
#' with an offset term to account for sequencing depth differences between cells.
#' \item{"poisson_offset"}: Poisson generalized linear mixed model with an
#' offset term to account for sequencing depth differences between cells.
#' }
#' For each of these options the user has the option to use either a Wald or
#' Likelihood ratio testing method by setting \code{de_type} to \code{"Wald"}
#' or \code{"LRT"}. Default is LRT.
#' }
#' @param de_method the specific differential expression testing method to use.
#' Please see the documentation under \code{de_family} for precise usage options,
#' or see the documentation at https://github.com/neurorestore/Libra. This
#' option will default to \code{wilcox} for \code{singlecell} methods, to
#' \code{edgeR} for \code{pseudobulk} methods, and \code{negbinom} for
#' \code{mixedmodel} methods.
#' @param de_test the specific mixed model test to use. Please see the
#' documentation under \code{de_family} for precise usage options,
#' or see the documentation at https://github.com/neurorestore/Libra. This
#' option defaults to \code{NULL} for \code{singlecell} methods, to \code{LRT}
#' for \code{pseudobulk} and \code{mixedmodel} methods.
#' @param n_threads number of threads to use for parallelization in mixed models.
#'
#' @return a data frame containing differential expression results with the
#' following columns:
#' \itemize{
#' \item{"cell type"}: The cell type DE tests were run on. By default Libra
#' will run DE on all cell types present in the original meta data.
#' \item{"gene"}: The gene being tested.
#' \item{"avg_logFC"}: The average log fold change between conditions. The
#' direction of the logFC can be controlled using factor levels of \code{label_col}
#' whereby a positive logFC reflects higher expression in the first level of
#' the factor, compared to the second.
#' \item{"p_val"}: The p-value resulting from the null hypothesis test.
#' \item{"p_val_adj"}: The adjusted p-value according to the Benjamini
#' Hochberg method (FDR).
#' \item{"de_family"}: The differential expression method family.
#' \item{"de_method"}: The precise differential expression method.
#' \item{"de_type"}: The differential expression method statistical testing type.
#' }
#'
#' @importFrom magrittr  %<>%
#' @importFrom forcats fct_recode
#' @importFrom dplyr group_by mutate select ungroup arrange
#' @examples data("hagai_toy")
#' hagai_toy@meta.data[hagai_toy$label == "unst" & hagai_toy$replicate == "mouse1","label"] <- "new_condition"
#' DE <- run_de.1(hagai_toy,min_reps = 1)
#' DE <- run_de.1(hagai_toy,min_reps = 1, ident.2 = "unst")
#' DE <- run_de.1(hagai_toy,min_reps = 1, ident.1 = "lps4")
#' DE <- run_de.1(hagai_toy,min_reps = 1, ident.2 = "new_condition")
run_de.1 = function(input,
                  meta = NULL,
                  replicate_col = 'replicate',
                  cell_type_col = 'cell_type',
                  label_col = 'label',
                  de_family = 'pseudobulk',
                  de_method = 'edgeR',
                  de_type = 'LRT',
                  ident.1 = NULL,
                  ident.2 = NULL,
                  min_cells = 3,
                  min_reps = 1,
                  min_features = 0,

                  n_threads = 2) {
    
    # first, make sure inputs are correct
    inputs = Libra:::check_inputs(
        input = input,
        meta = meta,
        replicate_col = replicate_col,
        cell_type_col = cell_type_col,
        label_col = label_col)
    input = inputs$expr
    meta = inputs$meta
    
    # run differential expression
    DE = switch(de_family,
                pseudobulk = pseudobulk_de.1(
                    input = input,
                    meta = meta,
                    replicate_col = replicate_col,
                    cell_type_col = cell_type_col,
                    label_col = label_col,
                    ident.1 = ident.1, 
                    ident.2 = ident.2,
                    min_cells = min_cells,
                    min_reps = min_reps,
                    min_features = min_features,
                    de_family = 'pseudobulk',
                    de_method = de_method,
                    de_type = de_type
                ),
                mixedmodel = Libra:::mixedmodel_de(
                    input = input,
                    meta = meta,
                    replicate_col = replicate_col,
                    cell_type_col = cell_type_col,
                    label_col = label_col,
                    min_features = min_features,
                    de_family = 'mixedmodel',
                    de_method = de_method,
                    de_type = de_type,
                    n_threads = n_threads
                ),
                singlecell = Libra:::singlecell_de(
                    input = input,
                    meta = meta,
                    cell_type_col = cell_type_col,
                    label_col = label_col,
                    min_features = min_features,
                    de_method = "MAST"
                )
    )
    
    # clean up the output
    suppressWarnings(
        colnames(DE) %<>%
            fct_recode('p_val' = 'p.value',  ## DESeq2
                       'p_val' = 'pvalue',  ## DESeq2
                       'p_val' = 'p.value',  ## t/wilcox
                       'p_val' = 'P.Value',  ## limma
                       'p_val' = 'PValue'  , ## edgeR
                       'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                       'p_val_adj' = 'adj.P.Val',      ## limma
                       'p_val_adj' = 'FDR',            ## edgeER
                       'avg_logFC' = 'log2FoldChange', ## DESEeq2
                       'avg_logFC' = 'logFC', ## limma/edgeR
                       'avg_logFC' = 'avg_log2FC' # Seurat V4
            )
    ) %>%
        as.character()
    
    DE %<>%
            # calculate adjusted p values
            group_by(cell_type) %>%
            mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
            # make sure gene is a character not a factor
            mutate(gene = as.character(gene)) %>%
        # invert logFC to match Seurat level coding
        #mutate(avg_logFC = avg_logFC * -1) %>%
        dplyr::select(gene,
                      avg_logFC,
                      p_val,
                      p_val_adj,
                      log2UMI.1,
                      log2UMI.2,
                      cluster,
                      cell_type,
                      de_method
        ) %>%
        ungroup() %>%
        group_by(cell_type,cluster) %>%
        arrange(desc(avg_logFC), .by_group = TRUE) %>%
        ungroup()
        
    return(DE)
}

#' Run pseudobulk differential expression methods
#' 
#' Run pseudobulk differential expression methods on single-cell data
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocle3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate 
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.
#' @param de_method the specific differential expression testing method to use.
#'   Defaults to edgeR.
#' @param de_type the specific parameter of the differential expression testing
#'   method. Defaults to LRT for edgeR, LRT for DESeq2, and trend for limma.
#' @param ident.1 Identity class to define markers for; if NULL, use all groups for comparison
#' @param ident.2 A second identity class for comparison; if NULL, use all other cells for comparison
#' @return a data frame containing differential expression results.
#'  
#' @importFrom magrittr %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate n_distinct
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest
#'   glmFit glmLRT topTags cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom purrr map 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' 
#' 
pseudobulk_de.1 = function(input, 
                         meta = NULL, 
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         ident.1 = NULL,
                         ident.2 = NULL,
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0,
                         de_family = 'pseudobulk',
                         de_method = 'edgeR',
                         de_type = 'LRT') {
    # check args
    if (de_method == 'limma') {
        if (de_type != 'voom') {
            # change default type to use
            de_type = 'trend'  
        }
    }
    
    # get the pseudobulks list
    print("Calculating the pseudobulks")
    pseudobulks = to_pseudobulk.1(
        input = input,
        meta = meta,
        replicate_col = replicate_col,
        cell_type_col = cell_type_col,
        label_col = label_col,
        min_cells = min_cells,
        min_reps = min_reps,
        min_features = min_features,
        external = F
    )
    
    results = pbapply::pblapply(pseudobulks, function(x) {
            # create targets matrix
            targets = data.frame(group_sample = colnames(x)) %>%
                    mutate(group = gsub(".*\\:", "", group_sample))
            if(!is.null(ident.1) & !any(ident.1 %in% targets$group)) return(NULL)
            if(!is.null(ident.2) & !any(ident.2 %in% targets$group)) return(NULL)
            
            ## optionally, carry over factor levels from entire dataset
            if (is.factor(meta$label)) {
                    targets$group %<>% factor(levels = levels(meta$label))
                    targets$group %<>% droplevels()
            }
            if (n_distinct(targets$group) >= 2) {
                    find_all_pseudobulk_de(x, targets,de_family = de_family,
                                       de_method = de_method,
                                       de_type = de_type,
                                       ident.1 = ident.1, ident.2 = ident.2)
            } else return(NULL)
            
            
    })
    results %<>% bind_rows(.id = 'cell_type')
    if(nrow(results) == 0){
            print("no gene is significantly differentially expressed.")
            return(NULL)
    } else return(results[,c("gene","logFC","logCPM","LR","PValue","FDR","log2UMI.1","log2UMI.2","cluster","cell_type","de_method")])
}



#' Run pseudobulk differential expression methods when group number == 2
#' @param x element of pseudobulks results
#' @param targets matrix
#' @param de_method the specific differential expression testing method to use.
#'   Defaults to edgeR.
#' @param de_type the specific parameter of the differential expression testing
#'   method. Defaults to LRT for edgeR, LRT for DESeq2, and trend for limma.
#' @param ident.1 Identity class to define markers for; if NULL, use all groups for comparison
#' @param ident.2 A second identity class for comparison; if NULL, use all other cells for comparison
#' @return a data frame containing differential expression results.
#'  
#' @importFrom magrittr %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate n_distinct
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest
#'   glmFit glmLRT topTags cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom purrr map 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @examples find_pseudobulk_de(x, targets,de_family = 'pseudobulk',de_method = 'edgeR',de_type = 'LRT')
find_pseudobulk_de <- function(x, targets, de_family = 'pseudobulk',
                               de_method = 'edgeR',
                               de_type = 'LRT',
                               ident.1 = NULL,
                               ident.2 = NULL){
        if (n_distinct(targets$group) != 2) return(NULL)
        # create design
        design = model.matrix(~ group, data = targets)
        
        DE <- switch(de_method,
                    edgeR = {
                            tryCatch({
                                    y = edgeR::DGEList(counts = x, group = targets$group) %>%
                                            edgeR::calcNormFactors(method = 'TMM') %>%
                                            edgeR::estimateDisp(design)
                                    test = switch(de_type,
                                                  QLF = {
                                                          fit = edgeR::glmQLFit(y, design)
                                                          test = edgeR::glmQLFTest(fit, coef = -1)
                                                  },
                                                  LRT = {
                                                          fit = edgeR::glmFit(y, design = design)
                                                          test = edgeR::glmLRT(fit)
                                                  })
                                    res = topTags(test, n = Inf) %>%
                                            as.data.frame() %>%
                                            rownames_to_column('gene') %>%
                                            # flag metrics in results
                                            mutate(cluster = paste(ident.1,"vs",ident.2),
                                                   de_method = paste0('pseudobulk ',de_method,"-",de_type))
                            }, error = function(e) {
                                    message(e)
                                    data.frame()
                            })
                    },
                    DESeq2 = {
                            tryCatch({
                                    dds = DESeqDataSetFromMatrix(countData = x,
                                                                 colData = targets,
                                                                 design = ~ group)
                                    dds = switch(de_type,
                                                 Wald = {
                                                         dds = try(DESeq2::DESeq(dds,
                                                                                 test = 'Wald',
                                                                                 fitType = 'parametric',
                                                                                 sfType = 'poscounts',
                                                                                 betaPrior = F))
                                                 },
                                                 LRT = {
                                                         dds = try(DESeq2::DESeq(dds,
                                                                                 test = 'LRT',
                                                                                 reduced = ~ 1,
                                                                                 fitType = 'parametric',
                                                                                 sfType = 'poscounts',
                                                                                 betaPrior = F))
                                                 }
                                    )
                                    res = results(dds)
                                    # write
                                    res = as.data.frame(res) %>%
                                            mutate(gene = rownames(x)) %>%
                                            # flag metrics in results
                                            mutate(cluster = paste(ident.1,"vs",ident.2),
                                                   de_method = paste0('pseudobulk ',de_method,"-",de_type))
                            }, error = function(e) {
                                    message(e)
                                    data.frame()
                            })
                    },
                    limma = {
                            tryCatch({
                                    x = switch(de_type,
                                               trend = {
                                                       trend_bool = T
                                                       dge = DGEList(as.matrix(x), group = targets$group)
                                                       dge = calcNormFactors(dge)
                                                       x = new("EList")
                                                       x$E = cpm(dge, log = TRUE, prior.count = 3)
                                                       x
                                               },
                                               voom = {
                                                       counts = all(as.matrix(x) %% 1 == 0)
                                                       if (counts) {
                                                               trend_bool = F
                                                               x = voom(as.matrix(x), design)
                                                               x
                                                       }
                                               })
                                    # get fit
                                    fit = lmFit(x, design) %>%
                                            eBayes(trend = trend_bool, robust = trend_bool)
                                    # format the results
                                    res = fit %>%
                                            # extract all coefs except intercept
                                            topTable(number = Inf, coef = -1) %>%
                                            tibble::rownames_to_column('gene') %>%
                                            # flag metrics in results
                                            mutate(cluster = paste(ident.1,"vs",ident.2),
                                                   de_method = paste0('pseudobulk ',de_method,"-",de_type))
                            }, error = function(e) {
                                    message(e)
                                    data.frame()
                            })
                    }
        )
        rowMeans_mean <- function(mtx){
                if(is.vector(mtx)) return(mtx)
                if(is.data.frame(mtx)) return(rowMeans(mtx))
        }
        DE %<>% mutate(log2UMI.1 = log2(rowMeans_mean(x[DE$gene,targets$group == ident.1])+1),
                        log2UMI.2 = log2(rowMeans_mean(x[DE$gene,targets$group == ident.2])+1))
        return(DE)
}




#' Run pseudobulk differential expression methods when group number > 2
#' @param x element of pseudobulks results
#' @param targets matrix
#' @param de_method the specific differential expression testing method to use.
#'   Defaults to edgeR.
#' @param de_type the specific parameter of the differential expression testing
#'   method. Defaults to LRT for edgeR, LRT for DESeq2, and trend for limma.
#' @param ident.1 Identity class to define markers for; if NULL, use all groups for comparison
#' @param ident.2 A second identity class for comparison; if NULL, use all other cells for comparison
#' @return a data frame containing differential expression results.
#'  
#' @importFrom magrittr %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate n_distinct
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest
#'   glmFit glmLRT topTags cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom purrr map 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @examples find_all_pseudobulk_de(x, targets,de_family = 'pseudobulk',de_method = 'edgeR',de_type = 'LRT')

find_all_pseudobulk_de <- function(x, targets, de_family = 'pseudobulk',
                               de_method = 'edgeR',
                               de_type = 'LRT', ident.1 = NULL, ident.2 = NULL){
        if (n_distinct(targets$group) < 2) return(NULL)
        
        groups <- switch(class(targets$group),
                         factor = levels(targets$group),
                         unique(targets$group))
        
        if(!is.null(ident.1)) {
                ident.1.grps <- groups[groups %in% ident.1]
                } else ident.1.grps <- groups
        DE <- pbapply::pblapply(ident.1.grps, function(grp) {
                targets0 <- targets
                if(is.null(ident.2)) {
                        rest_groups <- paste(groups[!groups %in% grp],collapse = "|")
                        targets0$group %<>% gsub(rest_groups, "rest", .)
                        ident.2 <- "rest"
                } else {
                        if(grp == ident.2) return(NULL)
                        targets0 %<>% filter(group %in% c(grp,ident.2))
                }
                targets0$group %<>% factor(levels = c(ident.2,grp))
                find_pseudobulk_de(x[,targets0$group_sample],targets0,de_family = ,de_family,de_method = de_method,
                                   de_type =  de_type,ident.1 = grp, ident.2 = ident.2)
                
        })
        names(DE) <- ident.1.grps
        
        if(is.null(ident.2)) {
                DE %<>% bind_rows(.id = "cluster")
        } else DE %<>% bind_rows
        return(DE)
        
}

#' Create a pseudobulk matrix
#' 
#' Convert a single-cell expression matrix (i.e., genes by cells)
#' to a pseudobulk matrix by summarizing counts within biological replicates
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocole3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate 
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.
#' @return a list of pseudobulk matrices, for each cell type.
#'  
#' @importFrom magrittr %<>% extract
#' @importFrom dplyr %>% rename_ count group_by filter pull n_distinct distinct
#'   summarise
#' @importFrom purrr map map_int
#' @importFrom Matrix rowSums colSums
#' @importFrom stats setNames
#' @export
#' 
to_pseudobulk.1 = function(input, 
                         meta = NULL, 
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0,
                         external = T) {
    if (external) {
        # first, make sure inputs are correct
        inputs = Libra:::check_inputs(
            input, 
            meta = meta,
            replicate_col = replicate_col,
            cell_type_col = cell_type_col,
            label_col = label_col)
        expr = inputs$expr
        meta = inputs$meta
    } else {
        expr = input
    }
    
    # convert to characters
    meta %<>% mutate(replicate = as.character(replicate),
                     cell_type = as.character(cell_type),
                     label = as.character(label),
                     cell_type_label = paste0(as.character(cell_type),"_",
                                              as.character(label)))
    
    # keep only cell types ~ group with enough cells
    
    keep_meta <- meta %>%
            dplyr::count(cell_type_label) %>%
            filter(n >= min_cells) %>%
            pull(cell_type_label) %>%
            unique()
    meta %<>% filter(cell_type_label %in%  keep_meta)
    # keep only cell types with enough cells
    keep = meta %>%
            dplyr::count(cell_type, label) %>%
            group_by(cell_type) %>%
            filter(all(n >= min_cells)) %>%
            pull(cell_type) %>%
            unique()
    
    # process data into gene x replicate x cell_type matrices
    pseudobulks = keep %>%
            purrr::map( ~ {
            print(.)
            cell_type = .
            meta0 = meta %>% filter(cell_type == !!cell_type)
            expr0 = expr %>% magrittr::extract(, meta0$cell_barcode)
            # catch cell types without replicates or conditions
            if (n_distinct(meta0$label) < 2) return(NA)
            replicate_counts = distinct(meta0, label, replicate) %>%
                group_by(label) %>%
                summarise(replicates = n_distinct(replicate)) %>%
                pull(replicates)
            if (any(replicate_counts < min_reps)) return(NA)
            
            # process data into gene X replicate X cell_type matrice
            mm = model.matrix(~ 0 + replicate:label, data = meta0)
            mat_mm = expr0 %*% mm
            keep_genes = Matrix::rowSums(mat_mm > 0) >= min_features
            mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
            mat_mm %<>% as.data.frame()
            colnames(mat_mm) %<>% gsub("^replicate", "", .) %>%
                                  gsub(":label", ":", .)
            # drop empty columns
            keep_samples = colSums(mat_mm) > 0
            mat_mm %<>% magrittr::extract(, keep_samples)
            return(mat_mm)
        }) %>%
        setNames(keep)
    
    # drop NAs
    pseudobulks %<>% magrittr::extract(!is.na(.))
    
    # also filter out cell types with no retained genes
    min_dim = purrr::map(pseudobulks, as.data.frame) %>% purrr::map(nrow)
    pseudobulks %<>% magrittr::extract(min_dim > 1)
    
    # also filter out types without replicates
    min_repl = purrr::map_int(pseudobulks, ~ {
        # make sure we have a data frame a not a vector
        tmp = as.data.frame(.)
        targets = data.frame(group_sample = colnames(tmp)) %>%
            mutate(group = gsub(".*\\:", "", group_sample))
        if (n_distinct(targets$group) == 1)
            return(as.integer(0))
        min(table(targets$group))
    })
    pseudobulks %<>% magrittr::extract(min_repl >= min_reps)
    return(pseudobulks)
}


