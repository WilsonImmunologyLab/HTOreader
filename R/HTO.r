#' @importFrom hash hash

#'
NULL

#' A seurat style theme for ggplot2 figures (Light version)
#'
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient aes element_rect element_line element_text theme margin
#' @return A theme object
#'
#' @export

LightTheme <- function(...) {
  light.background <- element_rect(fill = 'white')
  light.background.no.border <- element_rect(fill = 'white', size = 0)
  font.margin <- 4
  black.text <- element_text(
    size = 12,
    colour = 'black',
    margin = margin(
      t = font.margin,
      r = font.margin,
      b = font.margin,
      l = font.margin
    )
  )
  black.line <- element_line(colour = 'black', size = 1)
  no.line <- element_line(size = 0)
  #   Create the light theme
  light.theme <- theme(
    #   Set background colors
    plot.background = light.background,
    panel.background = light.background,
    legend.background = light.background,
    legend.box.background = light.background.no.border,
    legend.key = light.background.no.border,
    strip.background = element_rect(fill = 'grey50', colour = NA),
    #   Set text colors
    plot.title = black.text,
    plot.subtitle = black.text,
    axis.title = black.text,
    axis.text = black.text,
    legend.title = black.text,
    legend.text = black.text,
    strip.text = black.text,
    #   Set line colors
    axis.line.x = black.line,
    axis.line.y = black.line,
    panel.grid = no.line,
    panel.grid.minor = no.line,
    #   Validate the theme
    validate = TRUE,
    #   Extra parameters
    ...
  )
  return(light.theme)
}

#' findMaximum
#'
#' find Maximum value
#'
#' @param x a data array
#'
#' @export

findMaximum <- function(x){
  temp_list <- c()
  index_list <- c()

  for (m in 1:(length(x) - 2)) {
    tmp <- which.max(x[m:(m + 2)])
    if(tmp == 2) {
      temp_list <- c(temp_list, x[m + 1])
      index_list <- c(index_list, m + 1)
    }
  }

  out <- list(index = index_list, value = temp_list)
  return(out)
}

#' findMinimum
#'
#' find Minimum value
#'
#' @param x a data array
#'
#' @export

findMinimum <- function(x){
  temp_list <- c()
  index_list <- c()

  for (m in 1:(length(x) - 2)) {
    tmp <- which.min(x[m:(m + 2)])
    if(tmp == 2) {
      temp_list <- c(temp_list, x[m + 1])
      index_list <- c(index_list, m + 1)
    }
  }

  out <- list(index = index_list, value = temp_list)
  return(out)
}


#' HTOcutoff
#'
#' function to determine cutoff for HTO signals
#' @importFrom flexmix FLXMRglm flexmix parameters
#'
#' @param object Seurat object
#' @param assay the assay name of Hashtag signals, Default is "HTO"
#' @param method normalize method, could be "CLR" or 'log'. Default is "log"
#' @param min_limit min limit of cutoff, default is 1.5 (empirical value)
#' @param re_sampling_cutoff a cutoff for re-sampling when there are too many HTOs (> 4) in a single experiment, by default we set 1.5 for CLR and 3 for log method
#'
#' @export
HTOcutoff <-  function(object = NULL, assay = 'HTO', method = 'log', min_limit = 1.5, re_sampling_cutoff = NULL){
  assay.data <- object@assays[[assay]]@counts
  features <- rownames(assay.data)
  curoff_feature <- hash()
  if(is.null(re_sampling_cutoff)) {
    if(method == 'CLR') {
      re_sampling_cutoff <- 1.5
    } else {
      re_sampling_cutoff <- 3
    }
  }

  for (feature in features) {
    # fetch data
    if(method == 'log'){
      normdata <- log1p(assay.data[feature,])
    } else {
      normdata <- object@assays[[assay]]@data[feature,]
    }

    # when HTO numbers > 4, re-sampling to make the data distribution balanced
    if (length(features) > 4) {
      re_sampled_data <- normdata[which(normdata >= 3)]
      sampling_data <- sort(normdata[which(normdata < 3)])
      index <- seq(1,length(sampling_data),by= length(features)/2)
      sampling_data <- sampling_data[index]
      re_sampled_data <- c(sampling_data, re_sampled_data)

      re_sampled_data <- as.data.frame(re_sampled_data)
    } else {
      re_sampled_data <- as.data.frame(normdata)
      colnames(re_sampled_data) <- c('re_sampled_data')
    }

    # mixture modeling
    mo1 <- FLXMRglm(family = "gaussian")
    mo2 <- FLXMRglm(family = "gaussian")
    flexfit <- flexmix(re_sampled_data ~ 1, data = re_sampled_data, k = 2, model = list(mo1, mo2))

    # get fitted parameters
    c1 <- parameters(flexfit, component=1)[[1]]
    c2 <- parameters(flexfit, component=2)[[1]]

    # determine the cutoff
    # select_min <- c1[1] + (c2[1] - c1[1]) * c1[2] / (c1[2] + c2[2])
    select_min <- c1[1] + (c2[1] - c1[1]) * sqrt(c1[2]) / (sqrt(c1[2]) + sqrt(c2[2])) # when SD1 and SD2 are very different, use the ratio of SQRT instead of that of original value can improve the results
    if(select_min < min_limit){

      density.data <- density(normdata)
      if(c1[1] > c2[1]) {
        sel_index <- intersect(which(density.data$x > c2[1]), which(density.data$x < c1[1]))
      } else {
        sel_index <- intersect(which(density.data$x < c2[1]), which(density.data$x > c1[1]))
      }
      sel_y <- density.data$y[sel_index]
      sel_x <- density.data$x[sel_index]
      select_min <- sel_x[which(sel_y == min(sel_y))]
    }
    curoff_feature[feature] <- select_min
  }

  return(curoff_feature)
}

#' PlotHTO
#'
#' function to plot histgram of each HTO signal (and their cutoff)
#' @importFrom ggplot2 geom_vline geom_point ggplot ggtitle xlab ylab
#' @importFrom hash values
#'
#' @param object Seurat object
#' @param assay the assay name of Hashtag signals, Default is "HTO"
#' @param method normalize method, could be "CLR" or 'log'. Default is "log"
#' @param cutoff users can choose to display cutoff for each HTO signal by giving a data array
#' @param xlim max limit of x axis. Default is "NULL". Can be used to trim the density plot for better visualization of positive peaks and negative peaks
#' @param plot_existing_cutoff plot existing cutoff if the method fits
#'
#' @export
PlotHTO <- function(object = NULL, assay = 'HTO', method = 'log', cutoff = NULL, xlim = NULL, plot_existing_cutoff = TRUE){
  Plots <- list()

  assay.data <- object@assays[[assay]]@counts
  features <- rownames(assay.data)

  if(!is.null(cutoff)){
    if(length(features) != length(cutoff)){
      return(message('Number of cutoff is not equal to that of features!'))
    }
  } else {
    if(plot_existing_cutoff == TRUE){
      tryCatch(
        expr = {
          if(method == object@misc[['HTOmethod']]){
            cutoff <- values(object@misc[['HTOcutoff']])
          }
        },
        error = function(e){
          # (Optional)
          # Do this if an error is caught...
        },
        warning = function(w){
          # (Optional)
          # Do this if an warning is caught...
        },
        finally = {
          # (Optional)
          # Do this at the end before quitting the tryCatch structure...
        }
      )
    }
  }

  if(!is.null(xlim)){
    if(length(features) != length(xlim)){
      warning('Number of xlim is not equal to that of features, will use the first cutoff for all HTO density plots!')
    }
  }

  i = 1
  for (feature in features) {
    # fetch data
    if(method == 'log'){
      normdata <- log1p(assay.data[feature,])
    } else {
      normdata <- object@assays[[assay]]@data[feature,]
    }
    # calculate density, x is the normalized value, y is the density
    data.dense <- data.frame(density(normdata)[c('x','y')])
    plot <- ggplot(data.dense, aes(x=x, y=y)) + geom_point() + ggtitle(feature) + xlab('Normalized Value') + ylab('Density') + LightTheme() + geom_vline(xintercept=cutoff[i], colour='red', size = 2)
    if(!is.null(xlim)){
      if(length(xlim) == length(features)){
        plot <- plot + xlim(NA, xlim[i])
      }else{
        plot <- plot + xlim(NA, xlim[1])
      }
    }

    Plots[[i]] <- plot
    i <- i + 1
  }
  return(Plots)
}

#' HTOIdAssign
#'
#' assign HTO groups to each cell according to their HTO signals
#'
#' @param data nromalized HTO data
#' @param cutoff the cutoff values for HTO signals
#'
#' @export
HTOIdAssign <- function(data, cutoff){
  features <- rownames(data)
  HTOid <- rep('Negative', dim(data)[2])
  HTOid_detail <- rep('Negative', dim(data)[2])

  # for each cell
  for (cell_index in c(1:dim(data)[2])) {
    pos_feature <- c()
    # read its expression on each feature
    for (feature in features) {
      this_value <- data[feature,cell_index]
      if(this_value >= cutoff[[feature]]){
        pos_feature <- c(pos_feature, feature)
      }
    }

    if(length(pos_feature) == 1){
      HTOid[cell_index] <- pos_feature[1]
      HTOid_detail[cell_index] <- pos_feature[1]
    } else if (length(pos_feature) > 0){
      HTOid[cell_index] <- 'Doublet'
      HTOid_detail[cell_index] <- paste(pos_feature, collapse = '_')
    }
  }

  out <- list(HTOid  = HTOid, HTOid_detail = HTOid_detail)
  return(out)
}

#' HTOClassification
#'
#' main function for HTO group Classification
#'
#' @param object Seurat object
#' @param assay the assay name of Hashtag signals, Default is "HTO"
#' @param method normalize method, could be "CLR" or 'log'. If not specified, will use both methods and return the best results. Default is NULL
#' @param specify_cutoff users can specify their own cutoffs for HTOs instead of letting program to determine them
#' @param min_limit the min limit of cutoff, 1.5 for CLR and 3 for log (empirical value)
#'
#' @export
HTOClassification <-  function(object = NULL, assay = 'HTO', method = NULL, specify_cutoff = NULL, min_limit = NULL){
  if(is.null(method)){
      normdata2 <- object@assays[[assay]]@data
      normdata1 <- log1p(object@assays[[assay]]@counts)
      features <- rownames(normdata1)
    } else {
      if(method == 'CLR'){
        normdata <- object@assays[[assay]]@data
      } else {
        normdata <- log1p(object@assays[[assay]]@counts)
      }
      features <- rownames(normdata)
    }

  # determine cutoff and cluster cells
  if(!is.null(specify_cutoff)){
    if(is.null(method)){
      return(message('Please specify the method if you would like to provide the specify_cutoff!'))
    }
    # if cutoff is given
    if(length(specify_cutoff) != length(features)){
      return(message('Provided cutoff number is not equal to HTO feature number!'))
    }
    cutoff <- hash(features, specify_cutoff)
    hto_res <- HTOIdAssign(normdata, cutoff)
  } else {
    if(is.null(method)) {
      # if method is NULL, run both methods
      cutoff1 <- HTOcutoff(object, assay = assay,  method =  'log', min_limit = 3)
      hto_res1 <- HTOIdAssign(normdata1, cutoff1)
      n1 <- length(which(hto_res1$HTOid %in% c('Doublet','Negative')))

      cutoff2 <- HTOcutoff(object, assay = assay,  method =  'CLR', min_limit = 1.5)
      hto_res2 <- HTOIdAssign(normdata2, cutoff2)
      n2 <- length(which(hto_res2$HTOid %in% c('Doublet','Negative')))

      if(n1 > n2) {
        cutoff <- cutoff2
        hto_res <- hto_res2
        method <- 'CLR'
      } else {
        cutoff <- cutoff1
        hto_res <- hto_res1
        method <- 'log'
      }

    } else {
       # if cutoff is not given
      if (is.null(min_limit)){
        if(method == 'CLR') {
          min_limit = 1.5
        } else {
          min_limit = 3
        }
      }

      cutoff <- HTOcutoff(object, assay = assay,  method =  method, min_limit = min_limit) # also output cutoff to the public environment
      hto_res <- HTOIdAssign(normdata, cutoff)
    }
  }

  # save run info
  object@misc[['HTOmethod']] <- method
  object@misc[['HTOassay']] <- assay
  object@misc[['HTOcutoff']] <- cutoff

  # info display
  Msg = paste('Normalization method is: ', method, '\n')
  cat(Msg)
  for (v in ls(cutoff)) {
    str = paste('Cutoff for ',v,'is',cutoff[[v]],'\n', sep = ' ')
    cat(str)
  }

  # annotate Seurat object
  object$HTOid <- hto_res$HTOid
  object$HTOid_detail <- hto_res$HTOid_detail

  return(object)
}


#' HybridDemultiplexing
#'
#' function for HybridDemultiplexing
#'
#' @param object Seurat object
#' @param cellhashing_label label name for cell hashing in the meta data
#' @param genotype_label label name for snp-based method in the meta data
#' @param hto_names a vector of all HTO names
#' @param c_threshold threshold of convergence score, default is 0.7.
#' @param snp_doublet_label doublet label in SNP results.
#' @param snp_unassign_label unassigned label in SNP results.
#' @param cellhashing_doublet_label doublet label in cell hashing results.
#' @param cellhashing_negative_label negative label in cell hashing results.
#'
#' @export
#'
HybridDemultiplexing <- function(object = NULL, cellhashing_label = NULL, genotype_label = NULL, hto_names = NULL, c_threshold = 0.7, snp_doublet_label = 'doublet', snp_unassign_label = 'unassigned', cellhashing_doublet_label = 'Doublet', cellhashing_negative_label = 'Negative') {
  # object <- DN
  # cellhashing_label <- 'bff_cluster'
  # genotype_label <- 'souporcell'
  # hto_names <- c("S282","S283","S284","S344","S397","S417","S421","S423")

  working_data <- object@meta.data[,c(cellhashing_label, genotype_label)]
  colnames(working_data) <- c('Cellhash','Genotype')
  working_data$Cellhash <- as.character(working_data$Cellhash)
  working_data$Genotype <- as.character(working_data$Genotype)

  waiting_list <- c()
  for (hto in hto_names) {
    if(hto %in% names(table(working_data$Cellhash))) {
      # do nothing
    } else {
      waiting_list <- c(waiting_list, hto)
    }
  }

  working_data$Cellhash_cat <- 'singlet'
  working_data$Cellhash_cat[which(working_data$Cellhash == cellhashing_doublet_label)] <- 'doublet'
  working_data$Cellhash_cat[which(working_data$Cellhash == cellhashing_negative_label)] <- 'unassigned'
  working_data$Cellhash[which(working_data$Cellhash == cellhashing_doublet_label)] <- 'doublet'
  working_data$Cellhash[which(working_data$Cellhash == cellhashing_negative_label)] <- 'unassigned'

  working_data$Genotype_cat <- 'singlet'
  working_data$Genotype_cat[which(working_data$Genotype == snp_doublet_label)] <- 'doublet'
  working_data$Genotype_cat[which(working_data$Genotype == snp_unassign_label)] <- 'unassigned'
  working_data$Genotype[which(working_data$Genotype == snp_doublet_label)] <- 'doublet'
  working_data$Genotype[which(working_data$Genotype == snp_unassign_label)] <- 'unassigned'

  # link cell hash clusters with genotype groups
  snp_singlet_index <- which(working_data$Genotype_cat == 'singlet')
  cellhash_singlet_index <- which(working_data$Cellhash_cat == 'singlet')
  singlet_set_working_data <- working_data[intersect(snp_singlet_index, cellhash_singlet_index),]

  corr_data <- table(singlet_set_working_data$Cellhash, singlet_set_working_data$Genotype)

  pair_res <- list()
  for (cellhashing_label in rownames(corr_data)) {
    cur_row_data <- corr_data[cellhashing_label,]
    cur_row_data <- sort(cur_row_data, decreasing = TRUE)
    pair_res[[cellhashing_label]] <- names(cur_row_data[1])
  }
  if(length(waiting_list) > 1) {
    stop("More than 1 hashtags missing!")
  } else if (length(waiting_list) == 1) {
    pair_res[[waiting_list[1]]] <- setdiff(colnames(corr_data), unlist(pair_res))
    msg <- paste0('The missing HTO ',waiting_list[1], ' is assigned to ',  pair_res[[waiting_list[1]]])
    message(msg)
  }

  pair_res <- unlist(pair_res)
  info_str = 'We matched: \n'
  for (index in c(1:length(pair_res))) {
    info_str = paste0(info_str, names(pair_res)[index], ' with ', pair_res[index], '\n')
  }
  cat(info_str)
  # update object
  for (index in c(1:length(pair_res))) {
    working_data[which(working_data$Genotype == pair_res[index]), 'Genotype'] <- names(pair_res)[index]
  }

  # calculate convergence score
  singlet_set_working_data <- working_data[intersect(snp_singlet_index, cellhash_singlet_index),]
  c_score <- length(which(singlet_set_working_data$Cellhash == singlet_set_working_data$Genotype)) / dim(singlet_set_working_data)[1]

  if (c_score > c_threshold) {
    working_data$hybridID <- 'unassigned'

    consisstent_index <- which(working_data$Cellhash == working_data$Genotype)
    working_data$hybridID[consisstent_index] <- working_data$Genotype[consisstent_index]
    inconsisstent_case1_index <- intersect(which(working_data$Genotype_cat == 'singlet'), which(working_data$Cellhash_cat != 'singlet'))
    working_data$hybridID[inconsisstent_case1_index] <- working_data$Genotype[inconsisstent_case1_index]

    object$hybridID <- working_data$hybridID
    a <- table(working_data$hybridID)
    singlet_rate <- 1 - (( a['doublet'] + a['unassigned']) / sum(a))
    message("Hybrid demultiplexing done! Results in hybridID")
    return(object)

  } else {
    object$hybridID <- object@meta.data[,cellhashing_label]
    warning_msg <- paste("Corrolation between your cell hashing demultiplexing and SNP-based demultiplexing is very poor. The convergance score is", c_score, ", which is much lower than normal. In this case, we can not trust your SNP-based demultiplexing results. We set your cell hashing ID to the hybrid labels. However, a low convergance score usually indicates high doublet rates or even messed genetic background among different samples. We highly recommand the users to double check the data.")
    warning(warning_msg)
    return(object)
  }
}
