# returns tests that require count data
DEmethods_counts <- function() {
        c("negbinom", "poisson", "DESeq2")
}

# returns tests that do not support feature pre-filtering
DEmethods_noprefilter <- function() {
        c("DESeq2")
}

# returns tests that support latent variables (latent.vars)
DEmethods_latent <- function() {
        c('negbinom', 'poisson', 'MAST', "LR")
}

# returns tests that require CheckDots
DEmethods_checkdots <- function() {
        c('wilcox', 'MAST', 'DESeq2')
}

# returns tests that do not use Bonferroni correction on the DE results
DEmethods_nocorrect <- function() {
        c('roc')
}

# returns tests that require count data
DEmethods_counts <- function() {
        c("negbinom", "poisson", "DESeq2")
}

#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @param p.adjust.methods c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' #' @export gde.all data frame
#' @examples
#' data("pbmc_small")
#' # Find markers for all clusters
#' all.markers <- FindAllMarkers_UMI(object = pbmc_small)
#' head(x = all.markers)
#' 
FindAllMarkers_UMI <- function(
        object,
        assay = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = 'wilcox',
        p.adjust.methods = "bonferroni",
        slot = 'data',
        min.pct = 0.1,
        min.diff.pct = -Inf,
        node = NULL,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        mean.fxn = NULL,
        fc.name = NULL,
        base = 2,
        return.thresh = 1e-2,
        densify = FALSE,
        ...
) {
        MapVals <- function(vec, from, to) {
                vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
                vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
                return(unname(obj = vec2))
        }
        if ((test.use == "roc") && (return.thresh == 1e-2)) {
                return.thresh <- 0.7
        }
        if (is.null(x = node)) {
                idents.all <- sort(x = unique(x = Idents(object = object)))
        } else {
                if (!PackageCheck('ape', error = FALSE)) {
                        stop(cluster.ape, call. = FALSE)
                }
                tree <- Tool(object = object, slot = 'BuildClusterTree')
                if (is.null(x = tree)) {
                        stop("Please run 'BuildClusterTree' before finding markers on nodes")
                }
                descendants <- DFT(tree = tree, node = node, include.children = TRUE)
                all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
                descendants <- MapVals(
                        vec = descendants,
                        from = all.children,
                        to = tree$tip.label
                )
                drop.children <- setdiff(x = tree$tip.label, y = descendants)
                keep.children <- setdiff(x = tree$tip.label, y = drop.children)
                orig.nodes <- c(
                        node,
                        as.numeric(x = setdiff(x = descendants, y = keep.children))
                )
                tree <- ape::drop.tip(phy = tree, tip = drop.children)
                new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
                idents.all <- (tree$Nnode + 2):max(tree$edge)
        }
        genes.de <- list()
        messages <- list()
        for (i in 1:length(x = idents.all)) {
                if (verbose) {
                        message("Calculating cluster ", idents.all[i])
                }
                genes.de[[i]] <- tryCatch(
                        expr = {
                                FindMarkers_UMI(
                                        object = object,
                                        assay = assay,
                                        ident.1 = if (is.null(x = node)) {
                                                idents.all[i]
                                        } else {
                                                tree
                                        },
                                        ident.2 = if (is.null(x = node)) {
                                                NULL
                                        } else {
                                                idents.all[i]
                                        },
                                        features = features,
                                        logfc.threshold = logfc.threshold,
                                        test.use = test.use,
                                        p.adjust.methods = p.adjust.methods,
                                        slot = slot,
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
                                        mean.fxn = mean.fxn,
                                        fc.name = fc.name,
                                        base = base,
                                        densify = densify,
                                        ...
                                )
                        },
                        error = function(cond) {
                                return(cond$message)
                        }
                )
                if (is.character(x = genes.de[[i]])) {
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
                                gde <- subset(
                                        x = gde,
                                        subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
                                )
                        } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
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
                warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
                for (i in 1:length(x = messages)) {
                        if (!is.null(x = messages[[i]])) {
                                warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
                        }
                }
        }
        if (!is.null(x = node)) {
                gde.all$cluster <- MapVals(
                        vec = gde.all$cluster,
                        from = new.nodes,
                        to = orig.nodes
                )
        }
        return(gde.all)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics, with UMI average
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FindMarkers_UMI <- function(object, ...) {
        UseMethod(generic = 'FindMarkers_UMI', object = object)
}

#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param counts Count matrix if using scale.data for DE tests. This is used for
#' computing pct.1 and pct.2 and for filtering features based on fraction
#' expressing
#' @param features Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  \item{"negbinom"} : Identifies differentially expressed genes between two
#'   groups of cells using a negative binomial generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"poisson"} : Identifies differentially expressed genes between two
#'   groups of cells using a poisson generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"LR"} : Uses a logistic regression framework to determine differentially
#'  expressed genes. Constructs a logistic regression model predicting group
#'  membership based on each feature individually and compares this to a null
#'  model with a likelihood ratio test.
#'  \item{"MAST"} : Identifies differentially expressed genes between two groups
#'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
#'  package to run the DE testing.
#'  \item{"DESeq2"} : Identifies differentially expressed genes between two groups
#'  of cells based on a model using DESeq2 which uses a negative binomial
#'  distribution (Love et al, Genome Biology, 2014).This test does not support
#'  pre-filtering of genes based on average difference (or percent detection rate)
#'  between cell groups. However, genes may be pre-filtered based on their
#'  minimum detection rate (min.pct) across both cell groups. To use this method,
#'  please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param verbose Print a progress bar once expression testing begins
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'LR', 'negbinom', 'poisson', or 'MAST'
#' @param min.cells.feature Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param fc.results data.frame from FoldChange
#'
#' @importFrom Matrix rowMeans
#' @importFrom stats p.adjust
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers default
#'
FindMarkers_UMI.default <- function(
        object,
        slot = "data",
        counts = numeric(),
        cells.1 = NULL,
        cells.2 = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = 'wilcox',
        p.adjust.methods = "bonferroni",
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        fc.results = NULL,
        densify = FALSE,
        ...
) {
        Seurat:::ValidateCellGroups(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                min.cells.group = min.cells.group
        )
        features <- features %||% rownames(x = object)
        # reset parameters so no feature filtering is performed
        if (test.use %in% DEmethods_noprefilter()) {
                features <- rownames(x = object)
                min.diff.pct <- -Inf
                logfc.threshold <- 0
        }
        data <- switch(
                EXPR = slot,
                'scale.data' = counts,
                object
        )
        # feature selection (based on percentages)
        alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
        names(x = alpha.min) <- rownames(x = fc.results)
        features <- names(x = which(x = alpha.min >= min.pct))
        if (length(x = features) == 0) {
                warning("No features pass min.pct threshold; returning empty data.frame")
                return(fc.results[features, ])
        }
        alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
        features <- names(
                x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
        )
        if (length(x = features) == 0) {
                warning("No features pass min.diff.pct threshold; returning empty data.frame")
                return(fc.results[features, ])
        }
        # feature selection (based on logFC)
        if (slot != "scale.data") {
                total.diff <- fc.results[, 1] #first column is logFC
                names(total.diff) <- rownames(fc.results)
                features.diff <- if (only.pos) {
                        names(x = which(x = total.diff >= logfc.threshold))
                } else {
                        names(x = which(x = abs(x = total.diff) >= logfc.threshold))
                }
                features <- intersect(x = features, y = features.diff)
                if (length(x = features) == 0) {
                        warning("No features pass logfc.threshold threshold; returning empty data.frame")
                        return(fc.results[features, ])
                }
        }
        # subsample cell groups if they are too large
        if (max.cells.per.ident < Inf) {
                set.seed(seed = random.seed)
                if (length(x = cells.1) > max.cells.per.ident) {
                        cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
                }
                if (length(x = cells.2) > max.cells.per.ident) {
                        cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
                }
                if (!is.null(x = latent.vars)) {
                        latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
                }
        }
        de.results <- Seurat:::PerformDE(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                test.use = test.use,
                verbose = verbose,
                min.cells.feature = min.cells.feature,
                latent.vars = latent.vars,
                densify = densify,
                ...
        )
        print("PerformDE done !!")
        de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
        if (only.pos) {
                de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
        }
        if (test.use %in% DEmethods_nocorrect()) {
                de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
        } else {
                de.results <- de.results[order(de.results[,"p_val"], -de.results[, 1]), ]
                de.results[,"p_val_adj"] = p.adjust(
                        p = de.results[,"p_val"],
                        method = p.adjust.methods,
                        n = nrow(x = object)
                )
        }
        if(slot == "data"){
                avg_UMI.1 <- Matrix::rowMeans(expm1(x = object[features, cells.1]))
                avg_UMI.2 <- Matrix::rowMeans(expm1(x = object[features, cells.2]))
        }
        else if(slot == "counts"){
                avg_UMI.1 <- Matrix::rowMeans(x = object[features, cells.1])
                avg_UMI.2 <- Matrix::rowMeans(x = object[features, cells.2])
        }
        avg_UMI <-data.frame(avg_UMI.1, avg_UMI.2)
        de.results <- cbind(de.results,avg_UMI[match(rownames(de.results),rownames(avg_UMI)),])
        
        return(de.results)
}

#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers Assay
#'
FindMarkers_UMI.Assay <- function(
        object,
        slot = "data",
        cells.1 = NULL,
        cells.2 = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = 'wilcox',
        p.adjust.methods = "bonferroni",
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        mean.fxn = NULL,
        fc.name = NULL,
        base = 2,
        densify = FALSE,
        recorrect_umi = TRUE,
        ...
) {

        data.slot <- ifelse(
                test = test.use %in% DEmethods_counts(),
                yes = 'counts',
                no = slot
        )
        data.use <-  GetAssayData(object = object, slot = data.slot)
        counts <- switch(
                EXPR = data.slot,
                'scale.data' = GetAssayData(object = object, slot = "counts"),
                numeric()
        )
        fc.results <- FoldChange(
                object = object,
                slot = data.slot,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                pseudocount.use = pseudocount.use,
                mean.fxn = mean.fxn,
                fc.name = fc.name,
                base = base
        )
        de.results <- FindMarkers_UMI(
                object = data.use,
                slot = data.slot,
                counts = counts,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                logfc.threshold = logfc.threshold,
                p.adjust.methods = p.adjust.methods,
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
                fc.results = fc.results,
                densify = densify,
                ...
        )
        return(de.results)
}


#' @param recorrect_umi Recalculate corrected UMI counts using minimum of the median UMIs when performing DE using multiple SCT objects; default is TRUE
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers SCTAssay
#'
FindMarkers_UMI.SCTAssay <- function(
        object,
        slot = "data",
        cells.1 = NULL,
        cells.2 = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = 'wilcox',
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        mean.fxn = NULL,
        fc.name = NULL,
        base = 2,
        densify = FALSE,
        recorrect_umi = TRUE,
        ...
) {
        data.slot <- ifelse(
                test = test.use %in% DEmethods_counts(),
                yes = 'counts',
                no = slot
        )
        if (recorrect_umi && length(x = levels(x = object)) > 1) {
                cell_attributes <- SCTResults(object = object, slot = "cell.attributes")
                observed_median_umis <- lapply(
                        X = cell_attributes,
                        FUN = function(x) median(x[, "umi"])
                )
                model.list <- slot(object = object, "SCTModel.list")
                median_umi.status <- lapply(X = model.list,
                                            FUN = function(x) { return(tryCatch(
                                                    expr = slot(object = x, name = 'median_umi'),
                                                    error = function(...) {return(NULL)})
                                            )})
                if (any(is.null(unlist(median_umi.status)))){
                        stop("SCT assay does not contain median UMI information.",
                             "Run `PrepSCTFindMarkers()` before running `FindMarkers()` or invoke `FindMarkers(recorrect_umi=FALSE)`.")
                }
                model_median_umis <- SCTResults(object = object, slot = "median_umi")
                min_median_umi <- min(unlist(x = observed_median_umis))
                if (any(unlist(model_median_umis) != min_median_umi)){
                        stop("Object contains multiple models with unequal library sizes. Run `PrepSCTFindMarkers()` before running `FindMarkers()`.")
                }
        }
        
        data.use <-  GetAssayData(object = object, slot = data.slot)
        counts <- switch(
                EXPR = data.slot,
                'scale.data' = GetAssayData(object = object, slot = "counts"),
                numeric()
        )
        fc.results <- FoldChange(
                object = object,
                slot = data.slot,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                pseudocount.use = pseudocount.use,
                mean.fxn = mean.fxn,
                fc.name = fc.name,
                base = base
        )
        de.results <- FindMarkers_UMI(
                object = data.use,
                slot = data.slot,
                counts = counts,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
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
                fc.results = fc.results,
                densify = densify,
                ...
        )
        return(de.results)
}


#' @importFrom Matrix rowMeans
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers DimReduc
#'
FindMarkers_UMI.DimReduc <- function(
        object,
        cells.1 = NULL,
        cells.2 = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = "wilcox",
        p.adjust.methods = "bonferroni",
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        mean.fxn = rowMeans,
        fc.name = NULL,
        densify = FALSE,
        ...
        
) {
        if (test.use %in% DEmethods_counts()) {
                stop("The following tests cannot be used for differential expression on a reduction as they assume a count model: ",
                     paste(DEmethods_counts(), collapse=", "))
        }
        data <- t(x = Embeddings(object = object))
        Seurat:::ValidateCellGroups(
                object = data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                min.cells.group = min.cells.group
        )
        features <- features %||% rownames(x = data)
        # reset parameters so no feature filtering is performed
        if (test.use %in% DEmethods_noprefilter()) {
                features <- rownames(x = data)
                min.diff.pct <- -Inf
                logfc.threshold <- 0
        }
        fc.results <- FoldChange(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                mean.fxn = mean.fxn,
                fc.name = fc.name
        )
        # subsample cell groups if they are too large
        if (max.cells.per.ident < Inf) {
                set.seed(seed = random.seed)
                if (length(x = cells.1) > max.cells.per.ident) {
                        cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
                }
                if (length(x = cells.2) > max.cells.per.ident) {
                        cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
                }
                if (!is.null(x = latent.vars)) {
                        latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
                }
        }
        de.results <- Seurat:::PerformDE(
                object = data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                test.use = test.use,
                p.adjust.methods = p.adjust.methods,
                verbose = verbose,
                min.cells.feature = min.cells.feature,
                latent.vars = latent.vars,
                densify = densify,
                ...
        )
        de.results <- cbind(de.results, fc.results)
        if (only.pos) {
                de.results <- de.results[de.results$avg_diff > 0, , drop = FALSE]
        }
        if (test.use %in% DEmethods_nocorrect()) {
                de.results <- de.results[order(-de.results$power, -de.results$avg_diff), ]
        } else {
                de.results <- de.results[order(de.results$p_val, -de.results$avg_diff), ]
                de.results$p_val_adj = p.adjust(
                        p = de.results$p_val,
                        method = p.adjust.methods,
                        n = nrow(x = object)
                )
        }
        return(de.results)
}

#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param mean.fxn Function to use for fold change or average difference calculation.
#' If NULL, the appropriate function will be chose according to the slot used
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame. If NULL, the fold change column will be named
#' according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data
#' slot "avg_diff".
#' @param base The base with respect to which logarithms are computed.
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers Seurat
#'
FindMarkers_UMI.Seurat <- function(
        object,
        ident.1 = NULL,
        ident.2 = NULL,
        group.by = NULL,
        subset.ident = NULL,
        assay = NULL,
        slot = 'data',
        reduction = NULL,
        features = NULL,
        logfc.threshold = 0.25,
        test.use = "wilcox",
        p.adjust.methods = "bonferroni",
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        min.cells.feature = 3,
        min.cells.group = 3,
        pseudocount.use = 1,
        mean.fxn = NULL,
        fc.name = NULL,
        base = 2,
        densify = FALSE,
        ...
) {
        if (!is.null(x = group.by)) {
                if (!is.null(x = subset.ident)) {
                        object <- subset(x = object, idents = subset.ident)
                }
                Idents(object = object) <- group.by
        }
        if (!is.null(x = assay) && !is.null(x = reduction)) {
                stop("Please only specify either assay or reduction.")
        }
        if (length(x = ident.1) == 0) {
                stop("At least 1 ident must be specified in `ident.1`")
        }
        # select which data to use
        if (is.null(x = reduction)) {
                assay <- assay %||% DefaultAssay(object = object)
                data.use <- object[[assay]]
                cellnames.use <-  colnames(x = data.use)
        } else {
                data.use <- object[[reduction]]
                cellnames.use <- rownames(data.use)
        }
        cells <- Seurat:::IdentsToCells(
                object = object,
                ident.1 = ident.1,
                ident.2 = ident.2,
                cellnames.use = cellnames.use
        )
        # fetch latent.vars
        if (!is.null(x = latent.vars)) {
                latent.vars <- FetchData(
                        object = object,
                        vars = latent.vars,
                        cells = c(cells$cells.1, cells$cells.2)
                )
        }
        # check normalization method
        norm.command <- paste0("NormalizeData.", assay)
        if (norm.command %in% Command(object = object) && is.null(x = reduction)) {
                norm.method <- Command(
                        object = object,
                        command = norm.command,
                        value = "normalization.method"
                )
                if (norm.method != "LogNormalize") {
                        mean.fxn <- function(x) {
                                return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
                        }
                }
        }
        de.results <- FindMarkers_UMI(
                object = data.use,
                slot = slot,
                cells.1 = cells$cells.1,
                cells.2 = cells$cells.2,
                features = features,
                logfc.threshold = logfc.threshold,
                p.adjust.methods = p.adjust.methods,
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
                mean.fxn = mean.fxn,
                base = base,
                fc.name = fc.name,
                densify = densify,
                ...
        )
        return(de.results)
}


#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange default
FoldChange.default <- function(
        object,
        cells.1,
        cells.2,
        mean.fxn,
        fc.name,
        features = NULL,
        ...
) {
        features <- features %||% rownames(x = object)
        # Calculate percent expressed
        thresh.min <- 0
        pct.1 <- round(
                x = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
                        length(x = cells.1),
                digits = 3
        )
        pct.2 <- round(
                x = rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
                        length(x = cells.2),
                digits = 3
        )
        # Calculate fold change
        data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
        data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
        fc <- (data.1 - data.2)
        fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
        colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
        return(fc.results)
}


#' @importFrom Matrix rowMeans
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange Assay
FoldChange.Assay <- function(
        object,
        cells.1,
        cells.2,
        features = NULL,
        slot = "data",
        pseudocount.use = 1,
        fc.name = NULL,
        mean.fxn = NULL,
        base = 2,
        ...
) {
        data <- GetAssayData(object = object, slot = slot)
        mean.fxn <- mean.fxn %||% switch(
                EXPR = slot,
                'data' = function(x) {
                        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
                },
                'scale.data' = rowMeans,
                function(x) {
                        return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
                }
        )
        # Omit the decimal value of e from the column name if base == exp(1)
        base.text <- ifelse(
                test = base == exp(1),
                yes = "",
                no = base
        )
        fc.name <- fc.name %||% ifelse(
                test = slot == "scale.data",
                yes = "avg_diff",
                no = paste0("avg_log", base.text, "FC")
        )
        FoldChange(
                object = data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = features,
                mean.fxn = mean.fxn,
                fc.name = fc.name
        )
}

#' @importFrom Matrix rowMeans
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange DimReduc
FoldChange.DimReduc <- function(
        object,
        cells.1,
        cells.2,
        features = NULL,
        slot = NULL,
        pseudocount.use = NULL,
        fc.name = NULL,
        mean.fxn = NULL,
        ...
) {
        mean.fxn <- mean.fxn %||% rowMeans
        fc.name <- fc.name %||% "avg_diff"
        data <- t(x = Embeddings(object = object))
        features <- features %||% rownames(x = data)
        # Calculate avg difference
        data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
        data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
        fc <- (data.1 - data.2)
        fc.results <- data.frame(fc)
        colnames(fc.results) <- fc.name
        return(fc.results)
}

#' @param ident.1 Identity class to calculate fold change for; pass an object of class
#' \code{phylo} or 'clustertree' to calculate fold change for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to calculate fold change for
#' @param reduction Reduction to use - will calculate average difference on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to
#' calculating fold change (see example in \code{\link{FindMarkers}})
#' @param subset.ident Subset a particular identity class prior to regrouping.
#' Only relevant if group.by is set (see example in \code{\link{FindMarkers}})
#' @param assay Assay to use in fold change calculation
#' @param slot Slot to pull data from
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param mean.fxn Function to use for fold change or average difference calculation
#' @param base The base with respect to which logarithms are computed.
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#'
#' @rdname FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange Seurat
FoldChange.Seurat <- function(
        object,
        ident.1 = NULL,
        ident.2 = NULL,
        group.by = NULL,
        subset.ident = NULL,
        assay = NULL,
        slot = 'data',
        reduction = NULL,
        features = NULL,
        pseudocount.use = 1,
        mean.fxn = NULL,
        base = 2,
        fc.name = NULL,
        ...
) {
        if (!is.null(x = group.by)) {
                if (!is.null(x = subset.ident)) {
                        object <- subset(x = object, idents = subset.ident)
                }
                Idents(object = object) <- group.by
        }
        if (!is.null(x = assay) && !is.null(x = reduction)) {
                stop("Please only specify either assay or reduction.")
        }
        # select which data to use
        if (is.null(x = reduction)) {
                assay <- assay %||% DefaultAssay(object = object)
                data.use <- object[[assay]]
                cellnames.use <-  colnames(x = data.use)
        } else {
                data.use <- object[[reduction]]
                cellnames.use <- rownames(data.use)
        }
        cells <- IdentsToCells(
                object = object,
                ident.1 = ident.1,
                ident.2 = ident.2,
                cellnames.use = cellnames.use
        )
        fc.results <- FoldChange(
                object = data.use,
                cells.1 = cells$cells.1,
                cells.2 = cells$cells.2,
                features = features,
                slot = slot,
                pseudocount.use = pseudocount.use,
                mean.fxn = mean.fxn,
                base = base,
                fc.name = fc.name
        )
        return(fc.results)
}
