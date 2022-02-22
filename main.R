library('tidyverse')
library('RColorBrewer')
library(dplyr)

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
  sdata <- read.csv(intensity_data, header = TRUE, sep=delimiter)
  return(sdata)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  #dplyr::pivot_longer(intensity,names_to = "Samples",values_to = "Value")
  variance_value <- pca_results$sdev^2 #variance is the square of standard deviation
   V<-variance_value/ sum(variance_value) #variance by each PC
  return(V)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  new_tibble <- tibble(variance_explained = pca_ve,principal_components = as_factor(colnames(pca_results$x)),levels=colnames(pca_results$x), cumulative = cumsum(variance_explained))
  return(new_tibble)

}

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble gnerated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples
plot_pca_variance <- function(variance_tibble,fig.cap= caption_4) {
  lineplot_cols <- c('Cumulative'='red')
  barchart_cols <-c('Variance Explained'='black')
 # variance_tibble +lineplot_cols+barchart_cols 
 
   c<- ggplot(variance_tibble) + 
    geom_bar(aes(x=principal_components, y=variance_explained, fill='Variance Explained'), stat='identity', color='black') + #plot bw PC and Variances 
    geom_line(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) + #line connecting  cumulative sum of PCs
    geom_point(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) + # barchart bw  PC and cumulative sum of PCs
    scale_colour_manual(name="Cumulative",values=lineplot_cols) +
    scale_fill_manual(name="Variance Explained",values=barchart_cols) +
    labs(x='PC', y='% variance') +
    theme_bw( base_size = 11,base_family = "") +
    theme(axis.text.x = element_text(face="bold", color="#993333", size=5, angle=90))+
    ggtitle(fig.cap)
  
   return(c)
   
}

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results,fig.cap = caption_5) {
  meta <- readr::read_csv(metadata) %>% 
  dplyr::select(geo_accession, SixSubtypesClassification)
    
    labeled <- pca_results$x %>% 
      as_tibble(pca_results,rownames='geo_accession') %>% 
      left_join(meta, by='geo_accession')
    
    
    biplot <- ggplot(labeled) +
      geom_point(aes(x=PC1, y=PC2, color=SixSubtypesClassification)) +
      theme_gray()+ggtitle(fig.cap)
    return(biplot)
  }


#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_csv (str): The path to the differential expression results
#'   file we have provided
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the CSV.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
  significance_ids <- read_data(diff_exp_csv,',') %>% #reads data into dataframe
      as_tibble(rownames='probeid') %>% #rows considered are probe id 
      filter(padj < fdr_threshold) %>% #filter more than threshold value
      pull(probeid) #gives in a form of vector 
    return(significance_ids)
  
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  c <- intensity %>% 
    as_tibble(rownames='probeid') %>% 
    filter(probeid %in% sig_ids_list) %>% # probes in the list of  significant probes
    column_to_rownames(var='probeid') %>% #naming probeids as rownames
    as.matrix() #return as matrix
  return(c)
}


#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette,fig.cap= caption_8) {
  a<-heatmap(de_intensity,col=brewer.pal(num_colors,palette),main =fig.cap )
  return(a)
}
