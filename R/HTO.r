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
#' @param method normalize method, could be "CLR" or 'log'. Default is "CLR"
#' @param min_limit min limit of cutoff, default is 1
#'
#' @export
HTOcutoff <-  function(object = NULL, assay = 'HTO', method = 'CLR', min_limit = 1){
  assay.data <- object@assays[[assay]]@counts
  features <- rownames(assay.data)
  curoff_feature <- hash()

  for (feature in features) {
    # fetch data
    if(method == 'log'){
      normdata <- log1p(assay.data[feature,])
    } else {
      normdata <- object@assays[[assay]]@data[feature,]
    }

    normdata <- as.data.frame(normdata)
    # mixture modeling
    mo1 <- FLXMRglm(family = "gaussian")
    mo2 <- FLXMRglm(family = "gaussian")
    flexfit <- flexmix(normdata ~ 1, data = normdata, k = 2, model = list(mo1, mo2))

    # get fitted parameters
    c1 <- parameters(flexfit, component=1)[[1]]
    c2 <- parameters(flexfit, component=2)[[1]]

    # determine the cutoff
    # select_min <- c1[1] + (c2[1] - c1[1]) * c1[2] / (c1[2] + c2[2])
    select_min <- c1[1] + (c2[1] - c1[1]) * sqrt(c1[2]) / (sqrt(c1[2]) + sqrt(c2[2])) # when SD1 and SD2 are very different, use the ratio of SQRT instead of that of original value can improve the results
    if(select_min < min_limit){
      select_min <- min_limit
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
#' @param method normalize method, could be "CLR" or 'log'. Default is "CLR"
#' @param cutoff users can choose to display cutoff for each HTO signal by giving a data array
#' @param xlim max limit of x axis. Default is "NULL". Can be used to trim the density plot for better visualization of positive peaks and negative peaks
#' @param plot_existing_cutoff plot existing cutoff if the method fits
#'
#' @export
PlotHTO <- function(object = NULL, assay = 'HTO', method = 'CLR', cutoff = NULL, xlim = NULL, plot_existing_cutoff = TRUE){
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
#' @param method normalize method, could be "CLR" or 'log'. Default is "CLR"
#' @param specify_cutoff users can specify their own cutoffs for HTOs instead of letting program to determine them
#' @param min_limit the min limit of cutoff
#'
#' @export
HTOClassification <-  function(object = NULL, assay = 'HTO', method = 'CLR', specify_cutoff = NULL, min_limit = 1){
  if(method == 'CLR'){
    normdata <- object@assays[[assay]]@data
  } else {
    normdata <- log1p(object@assays[[assay]]@counts)
  }
  features <- rownames(normdata)

  # determine cutoff and cluster cells
  if(!is.null(specify_cutoff)){
    # if cutoff is given
    if(length(specify_cutoff) != length(features)){
      return(message('Provided cutoff number is not equal to HTO feature number!'))
    }
    cutoff <- hash(features, specify_cutoff)
    hto_res <- HTOIdAssign(normdata, cutoff)
  } else {
    # if cutoff is not given
    cutoff <- HTOcutoff(object, assay = assay,  method =  method, min_limit = min_limit) # also output cutoff to the public environment
    hto_res <- HTOIdAssign(normdata, cutoff)
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
