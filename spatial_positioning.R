############### PACKAGES ###############

library(dplyr)
library(ggplot2)
library(ggthemes)
library(fpc)
library(dbscan)
library(factoextra)
library(ggpubr)
library(tidyverse)
library(spatialEco)
library(sp)
library(cowplot)
library(ptinpoly)
library(tidyr)
library(sf)
library(viridis)



############### USER INPUT ###############

# run names
run <- c("samplename")

# choose appropriate scaling factor
# scaling_factor <- 2*0.73 # 5.5 x 5.5 mm square puck
scaling_factor <- 0.73 # 3 mm circular puck

# DBSCAN parameters
eps = 50
minPts_vec = c(4,5,6,7,8,9,10,11,12,13,14,15)

# number of random cell barcode plots to generate
number_to_sample = 50

# high nUMI cut-off for SB combined filtering 
nUMI_cutoff_combined = 256 



############### LOAD DATA ###############

list_sb_matched <- list()
list_sb_whitelist <- list()
for(i in seq_along(run)){
  list_sb_whitelist[[i]] <- read.delim(paste0("/path/to/slidetags_", run[i], "/df_whitelist_", run[i], ".txt"))
  list_sb_matched[[i]] <- read.csv(paste0("/path/to/slide-tags/slidetags_", run[i], "/matching_result_", run[i], ".csv"))
}



############### SPATIAL POSITIONING ###############

for(i in seq_along(run)){
  


  ############### OUTPUT PATHS ###############
  
  ## create a directory for output if it doesn't already exist
  dir.create(paste0("/path/to/slidetags_", run[i], "/output"))
  dir.create(paste0("/path/to/slidetags_", run[i], "/output/cell_bc_plots"))
  dir.create(paste0("/path/to/slidetags_", run[i], "/output/plots"))
  
  # specify paths
  coordinates_df_output_path <- paste0("/path/to/slidetags_", run[i], "/output/coordinates_df_", run[i], ".txt")
  cell_bc_plots_path <- paste0("/path/to/slidetags_", run[i], "/output/cell_bc_plots/")
  general_plots_path <- paste0("/path/to/slidetags_", run[i], "/output/plots/")
  


  ############### DATA SETUP ###############

  # scale coordinates to microns
  list_sb_matched[[i]]$x_um <- list_sb_matched[[i]]$xcoord*scaling_factor
  list_sb_matched[[i]]$y_um <- list_sb_matched[[i]]$ycoord*scaling_factor
  
  # clean up data frame (get rid of unscaled coordinates)
  list_sb_matched[[i]] <- list_sb_matched[[i]] %>% dplyr::select(matched_beadbarcode, Illumina_barcode, x_um, y_um)
  
  # merge matching results with whitelist
  cb_sb_spatial_data <- merge(list_sb_matched[[i]], list_sb_whitelist[[i]], by.x = "Illumina_barcode", by.y = "bead_bc")
  
  # remove non-distinct rows (if there are any, there technically shouldn't be)
  cb_sb_spatial_data <- cb_sb_spatial_data %>% distinct()
  
  

  ############### PLOT SPATIAL BARCODES OF RANDOM CELLS ###############
  
  # generate list of distinct SB
  cell_bc_10x_distinct <- as.vector(cb_sb_spatial_data %>% distinct(cell_bc_10x))$cell_bc_10x
  
  random_sample_cell_bc_10x <- sample(cell_bc_10x_distinct, number_to_sample, replace = FALSE)
  
  for (j in 1:number_to_sample) {
    cell_bc <- random_sample_cell_bc_10x[j]
    cell <- cb_sb_spatial_data %>% dplyr::filter(cell_bc_10x == cell_bc)
    cell <- cell[order(cell$nUMI),]
    p <- ggplot(cell, aes(x = x_um, y = y_um, colour = nUMI)) + 
      geom_point(alpha = 0.5) + 
      theme_few() +
      ggtitle(paste(cell_bc)) + coord_fixed(ratio = 1)
    ggsave(paste(cell_bc_plots_path, cell_bc, ".jpeg", sep = ""), plot = p, width = 6, height = 4)
  }
  
  

  ############### PLOT SPATIAL BARCODE UMIS IN SPACE ###############
  
  # cell barcode - spatial barcode combinations 
  cb_sb_spatial_data <- cb_sb_spatial_data[order(cb_sb_spatial_data$nUMI), ]
  p<-ggplot(data=cb_sb_spatial_data,aes(x= x_um, y=y_um, col= nUMI))+
    geom_point(size=0.8)+
    coord_fixed()+
    theme_classic() +
    scale_color_viridis(trans = "log", option = "B") +
    labs(x = "x (um)", y = "y (um)") +
    ggtitle("Cell barcode grouped spatial barcode nUMIs")
  ggsave(paste(general_plots_path, "SB_spatial_UMI_grouped.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # collapse illumunia_barcode into matched_beadbarcode and recalculate nUMI
  cb_sb_spatial_data_collapse <- aggregate(cb_sb_spatial_data$nUMI, list(cb_sb_spatial_data$matched_beadbarcode), FUN=sum) 
  names(cb_sb_spatial_data_collapse) <- c("matched_beadbarcode", "nUMI_collapsed")
  coordinates <- cb_sb_spatial_data %>% dplyr::select(matched_beadbarcode, x_um, y_um)
  cb_sb_spatial_data_collapse <- merge(cb_sb_spatial_data_collapse, coordinates, by = "matched_beadbarcode", all = FALSE)
  cb_sb_spatial_data_collapse <- cb_sb_spatial_data_collapse[!duplicated(cb_sb_spatial_data_collapse),]
  
  # plot collapsed
  cb_sb_spatial_data_collapse <- cb_sb_spatial_data_collapse[order(cb_sb_spatial_data_collapse$nUMI), ]
  p<-ggplot(data=cb_sb_spatial_data_collapse,aes(x= x_um, y=y_um, col= nUMI_collapsed))+
    geom_point(size=0.8)+
    coord_fixed()+
    theme_classic() +
    scale_color_viridis(trans = "log", option = "B") +
    labs(x = "x (um)", y = "y (um)") +
    ggtitle("Total spatial barcode nUMIs") + 
    labs(colour = "nUMI") 
  ggsave(paste(general_plots_path, "SB_spatial_UMI_total.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  

  ############### PLOT SPATIAL BARCODES UMI DISTRIBUTION ###############
  
  # cb-sb combinations
  p <- ggplot(cb_sb_spatial_data, aes(x = " ", y = log2(nUMI))) + 
    geom_jitter(alpha = 0.2) +
    geom_violin() + 
    theme_few() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab("")
  ggsave(paste(general_plots_path, "SB_violin_grouped_log2.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # total sb nUMIs
  p <- ggplot(cb_sb_spatial_data_collapse, aes(x = " ", y = log2(nUMI_collapsed))) + 
    geom_jitter(alpha = 0.2) +
    geom_violin() + 
    theme_few() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_hline(yintercept = log2(nUMI_cutoff_combined), linetype = "dashed", colour = "red") +
    xlab("")
  ggsave(paste(general_plots_path, "SB_violin_total_log2.pdf", sep = ""), plot = p, width = 6, height = 4)
  

  
  ############### REMOVE VERY HIGH UMIS BASED ON DISTRIBUTION ABOVE ###############
  
  cb_sb_spatial_data_collapse_filtered <- cb_sb_spatial_data_collapse %>% dplyr::filter(nUMI_collapsed < nUMI_cutoff_combined)
  
  cb_sb_spatial_data_filtered_by_total <- merge(cb_sb_spatial_data_collapse_filtered, cb_sb_spatial_data, by = "matched_beadbarcode")
  cb_sb_spatial_data_filtered_by_total <- cb_sb_spatial_data_filtered_by_total %>% dplyr::select(matched_beadbarcode,
                                                                                                 nUMI_collapsed,
                                                                                                 x_um.x,
                                                                                                 y_um.x,
                                                                                                 Illumina_barcode,
                                                                                                 CB_SB,
                                                                                                 nUMI,
                                                                                                 cell_bc_10x
  )
  names(cb_sb_spatial_data_filtered_by_total) <- c("matched_beadbarcode",
                                                   "nUMI_collapsed",
                                                   "x_um",
                                                   "y_um",
                                                   "Illumina_barcode",
                                                   "CB_SB",
                                                   "nUMI",
                                                   "cell_bc_10x"
  )
  
  

  
  # aggregate proportion cells mapped across minPts parameters
  proportion_cells_mapped_vec <- NULL
  
  
  coordinates_df_rbind <- data.frame(matrix(ncol = 16, nrow = 0))
  x <- c("cell_bc", "x_um", "y_um", "number_clusters","SB_total","SB_UMI_total","SB_noise", "SB_UMI_noise","SB_top_cluster","SB_UMI_top_cluster", "proportion_SB_noise" ,"proportion_SB_UMI_noise","proportion_SB_top_cluster","proportion_SB_UMI_top_cluster","minPts","eps")
  colnames(coordinates_df_rbind) <- x
  
  

  ############### RUN DBSCAN ###############
  
  for (minPts in minPts_vec) {
    
    print(paste("minPts =", minPts))
    
    # create list of distinct cell barcodes to iterate through in for loop (new list in case filtering removes cells)
    list_of_cell_bc <- cb_sb_spatial_data_filtered_by_total %>% distinct(cell_bc_10x)
    list_of_cell_bc <- as.vector(list_of_cell_bc$cell_bc_10x)
    
    # set vectors to NULL for dbscan performance parameters (some of these are legacy and no longer needed)
    cell_bc_vec <- NULL
    centroid_x_vec <- NULL
    centroid_y_vec <- NULL
    number_clusters_vec <- NULL
    SB_total_vec <- NULL
    SB_UMI_total_vec <- NULL
    SB_noise_vec <- NULL
    SB_UMI_noise_vec <- NULL
    SB_top_cluster_vec <- NULL
    SB_UMI_top_cluster_vec <- NULL
    coordinates_df <- NULL
    coordinates_df_labeled <- NULL
    coordinates_df_labeled_filtered <- NULL
    CA1 <- NULL
    CA1_plot <- NULL
    DG1 <- NULL
    DG1_plot <- NULL
    dbscan_df <- data.frame()
    
    # run dbscan for the spatial barcodes within each cell
    for (cell_bc in list_of_cell_bc) {
      
      # output single cell
      cell <- cb_sb_spatial_data_filtered_by_total %>% filter(cell_bc_10x == cell_bc)
      
      # run dbscan
      db <- dbscan::dbscan(cell[3:4], eps = eps, minPts = minPts, weights = cell$nUMI)
      
      # collect dbscan cluster assignments
      cell$cluster <- db$cluster
      number_clusters <- max(cell$cluster)
      cell$cluster <- as.factor(cell$cluster)
      
      # collect dbscan metrics
      SB_total <- length(cell$matched_beadbarcode)
      SB_UMI_total <- sum(cell$nUMI)
      cell_noise <- cell %>% filter(cluster == 0)
      SB_noise <- length(cell_noise$matched_beadbarcode)
      SB_UMI_noise <- sum(cell_noise$nUMI)
      
      # assign cell coordinates as centroid of top cluster (later i filter out cells with > 1 cluster)
      # otherwise, assign coordinates as NA
      cell_top_cluster <- cell %>% filter(cluster == 1)
      
      if (dim(cell_top_cluster)[1] > 0) {
        
        # depends on version of sp package
        # cell_sp <- SpatialPointsDataFrame(cell_top_cluster[3:4], cell_top_cluster)
        # cell_sp_df <- as.data.frame(cell_sp)
        # centroid <- wt.centroid(st_as_sf(cell_sp), "nUMI", FALSE)
        # centroid_x <- centroid[1]
        # centroid_y <- centroid[2]
        
        cell_sp <- SpatialPointsDataFrame(cell_top_cluster[3:4], cell_top_cluster)
        cell_sp_df <- as.data.frame(cell_sp)
        centroid <- wt.centroid(st_as_sf(cell_sp), "nUMI", FALSE)
        centroid_x <- centroid$X
        centroid_y <- centroid$Y
        
        SB_top_cluster <- length(cell_top_cluster$matched_beadbarcode)
        SB_UMI_top_cluster <- sum(cell_top_cluster$nUMI)
        
        dbscan_df <- rbind(dbscan_df, cell_sp_df)
        
      } else {
        
        centroid_x <- NA
        centroid_y <- NA
        SB_top_cluster <- NA
        SB_UMI_top_cluster <- NA
        
      }
      
      # add dbscan parameters to vector
      cell_bc_vec <- append(cell_bc_vec, paste0(cell_bc))
      centroid_x_vec <- c(centroid_x_vec, paste0(centroid_x))
      centroid_y_vec <- c(centroid_y_vec, paste0(centroid_y))
      number_clusters_vec <- c(number_clusters_vec, paste0(number_clusters))
      SB_total_vec <- c(SB_total_vec, paste0(SB_total))
      SB_UMI_total_vec <- c(SB_UMI_total_vec, paste0(SB_UMI_total))
      SB_noise_vec <- c(SB_noise_vec, paste0(SB_noise))
      SB_UMI_noise_vec <- c(SB_UMI_noise_vec, paste0(SB_UMI_noise))
      SB_top_cluster_vec <- c(SB_top_cluster_vec, paste0(SB_top_cluster))
      SB_UMI_top_cluster_vec <- c(SB_UMI_top_cluster_vec, paste0(SB_UMI_top_cluster))
      
    }
    
    # aggregate dbscan information into dataframe
    coordinates_df_test <- data.frame(cell_bc_vec, centroid_x_vec, centroid_y_vec, number_clusters_vec, 
                                      SB_total_vec, SB_UMI_total_vec, SB_noise_vec, SB_UMI_noise_vec,
                                      SB_top_cluster_vec, SB_UMI_top_cluster_vec)
    
    # remove duplicates
    coordinates_df_test <- coordinates_df_test %>% unique()
    
    # rename columns 
    names(coordinates_df_test) <- c("cell_bc", "x_um", "y_um", "number_clusters", "SB_total", 
                                    "SB_UMI_total", "SB_noise", "SB_UMI_noise", "SB_top_cluster",
                                    "SB_UMI_top_cluster")
    
    # make columns numeric as appropriate
    coordinates_df_test$x_um <- as.numeric(coordinates_df_test$x_um)
    coordinates_df_test$y_um <- as.numeric(coordinates_df_test$y_um)
    coordinates_df_test$number_clusters <- as.numeric(coordinates_df_test$number_clusters)
    coordinates_df_test$SB_total <- as.numeric(coordinates_df_test$SB_total)
    coordinates_df_test$SB_UMI_total <- as.numeric(coordinates_df_test$SB_UMI_total)
    coordinates_df_test$SB_noise <- as.numeric(coordinates_df_test$SB_noise)
    coordinates_df_test$SB_UMI_noise <- as.numeric(coordinates_df_test$SB_UMI_noise)
    coordinates_df_test$SB_top_cluster <- as.numeric(coordinates_df_test$SB_top_cluster)
    coordinates_df_test$SB_UMI_top_cluster <- as.numeric(coordinates_df_test$SB_UMI_top_cluster)
    
    # replace NA with 0 in appropriate columns
    coordinates_df_test$SB_top_cluster[is.na(coordinates_df_test$SB_top_cluster)] = 0
    coordinates_df_test$SB_UMI_top_cluster[is.na(coordinates_df_test$SB_UMI_top_cluster)] = 0
    
    # create column with proportion SB noise and proportion SB UMI noise
    coordinates_df_test <- coordinates_df_test %>% 
      mutate(proportion_SB_noise = SB_noise / SB_total) %>%
      mutate(proportion_SB_UMI_noise = SB_UMI_noise / SB_UMI_total) %>%
      mutate(proportion_SB_top_cluster = SB_top_cluster / SB_total) %>%
      mutate(proportion_SB_UMI_top_cluster = SB_UMI_top_cluster / SB_UMI_total) %>%
      mutate(minPts = minPts) %>%
      mutate(eps = eps)
    
    coordinates_df_rbind <- rbind(coordinates_df_rbind, coordinates_df_test)
    
    
    
    ############### DBSCAN SUMMARY PLOT ###############
    
    p<-ggplot(coordinates_df_test, aes(x = number_clusters)) + 
      geom_histogram(colour = "black", fill = "grey", bins = 50) + 
      theme_few() +
      xlab("DBSCAN Clusters") +
      ylab("Count")
    ggsave(paste(general_plots_path, minPts, "_dbscan_clusters_histogram.png", sep = ""), plot = p, width = 6, height = 4)
    
    
    
    ############### DBSCAN SUMMARY FILE ###############
    
    # assemble summary metrics for summary output file
    eps
    minPts
    
    # number cells mapped
    spatially_mapped_top_cluster <- (coordinates_df_test %>% filter(number_clusters == 1) %>% count())$n
    spatially_mapped_any_cluster_number <- (coordinates_df_test %>% filter(number_clusters > 0) %>% count())$n
    total_cells <- dim(coordinates_df_test)[1]
    proportion_cells_mapped <- (spatially_mapped_top_cluster/total_cells)[1]
    
    # sb per cell
    median_unique_SB_per_cell <- median(coordinates_df_test$SB_total)
    mean_unique_SB_per_cell <- mean(coordinates_df_test$SB_total)
    median_SB_UMI_per_cell <- median(coordinates_df_test$SB_UMI_total)
    mean_SB_UMI_per_cell <- mean(coordinates_df_test$SB_UMI_total)
    
    # signal to noise ratios total
    median_proportion_unique_SB_top_cluster_total <- median(coordinates_df_test$proportion_SB_top_cluster)
    mean_proportion_unique_SB_top_cluster_total <- mean(coordinates_df_test$proportion_SB_top_cluster)
    median_proportion_SB_UMI_top_cluster_total <- median(coordinates_df_test$proportion_SB_UMI_top_cluster)
    mean_proportion_SB_UMI_top_cluster_total <- mean(coordinates_df_test$proportion_SB_UMI_top_cluster)
    
    # signal to noise ratios for mapped cells
    median_mapped_cells_proportion_unique_SB_top_cluster <- (coordinates_df_test %>% 
                                                               filter(number_clusters == 1) %>% 
                                                               select(proportion_SB_top_cluster) %>% 
                                                               summarise(across(everything(), median)))$proportion_SB_top_cluster
    mean_mapped_cells_proportion_unique_SB_top_cluster <- (coordinates_df_test %>% 
                                                             filter(number_clusters == 1) %>% 
                                                             select(proportion_SB_top_cluster) %>% 
                                                             summarise(across(everything(), mean)))$proportion_SB_top_cluster
    median_mapped_cells_proportion_SB_UMI_top_cluster <- (coordinates_df_test %>% 
                                                            filter(number_clusters == 1) %>% 
                                                            select(proportion_SB_UMI_top_cluster) %>% 
                                                            summarise(across(everything(), median)))$proportion_SB_UMI_top_cluster
    mean_mapped_cells_proportion_SB_UMI_top_cluster <- (coordinates_df_test %>% 
                                                          filter(number_clusters == 1) %>% 
                                                          select(proportion_SB_UMI_top_cluster) %>% 
                                                          summarise(across(everything(), mean)))$proportion_SB_UMI_top_cluster
    
    summary_metrics <- data.frame(eps, 
                                  minPts, 
                                  spatially_mapped_top_cluster,
                                  spatially_mapped_any_cluster_number,
                                  total_cells,
                                  proportion_cells_mapped,
                                  median_unique_SB_per_cell,
                                  mean_unique_SB_per_cell,
                                  median_SB_UMI_per_cell,
                                  mean_SB_UMI_per_cell,
                                  median_proportion_unique_SB_top_cluster_total,
                                  mean_proportion_unique_SB_top_cluster_total,
                                  median_proportion_SB_UMI_top_cluster_total,
                                  mean_proportion_SB_UMI_top_cluster_total,
                                  median_mapped_cells_proportion_unique_SB_top_cluster,
                                  mean_mapped_cells_proportion_unique_SB_top_cluster,
                                  median_mapped_cells_proportion_SB_UMI_top_cluster,
                                  mean_mapped_cells_proportion_SB_UMI_top_cluster
    )
    
    write.csv(summary_metrics, paste(general_plots_path, minPts, "_summary_metrics.csv", sep = ""))
    
    
    proportion_cells_mapped_vec <- c(proportion_cells_mapped_vec, paste0(proportion_cells_mapped))
    
    
    
    
  }
  


  ############### SELECT OPTIMAL MINPTS ###############
  
  proportion_cells_mapped_df <- data.frame(minPts_vec, proportion_cells_mapped_vec)
  names(proportion_cells_mapped_df) <- c("minPts", "proportion_cells_mapped")
  proportion_cells_mapped_df$minPts <- as.numeric(proportion_cells_mapped_df$minPts)
  proportion_cells_mapped_df$proportion_cells_mapped <- as.numeric(proportion_cells_mapped_df$proportion_cells_mapped)
  
  # plot proportion cells mapped
  p <- ggplot(proportion_cells_mapped_df, aes(x = minPts, y = proportion_cells_mapped)) +
    geom_point() +
    geom_line() +
    theme_few() +
    ylab("Proportion nuclei mapped") +
    xlab("minPts parameter")
  print(p)
  ggsave(paste(general_plots_path, "dbscan_sensitivity_optimization_curve.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # most sensitive dbscan parameter
  most_sensitive_minPts <- proportion_cells_mapped_df[which(proportion_cells_mapped_df$proportion_cells_mapped == max(proportion_cells_mapped_df$proportion_cells_mapped))
                                                      ,]$minPts
  
  print(most_sensitive_minPts)
  # if a tie between parameters, take the biggest one? maybe more specific? hard to tell...
  if (length(most_sensitive_minPts) > 1) {
    most_sensitive_minPts = max(most_sensitive_minPts)
  }
  
  print(most_sensitive_minPts)
  


  ############### RERUN UNDER OPTIMAL PARAMETERS ###############
    
  # parameters
  minPts = most_sensitive_minPts
  
  # create list of distinct cell barcodes 
  list_of_cell_bc <- cb_sb_spatial_data_filtered_by_total %>% distinct(cell_bc_10x)
  list_of_cell_bc <- as.vector(list_of_cell_bc$cell_bc_10x)
  
  # set vectors to NULL for dbscan performance parameters (some of these are legacy and no longer needed)
  cell_bc_vec <- NULL
  centroid_x_vec <- NULL
  centroid_y_vec <- NULL
  number_clusters_vec <- NULL
  SB_total_vec <- NULL
  SB_UMI_total_vec <- NULL
  SB_noise_vec <- NULL
  SB_UMI_noise_vec <- NULL
  SB_top_cluster_vec <- NULL
  SB_UMI_top_cluster_vec <- NULL
  coordinates_df_labeled <- NULL
  coordinates_df_labeled_filtered <- NULL
  coordinates_df <- NULL
  dbscan_df <- data.frame()
  
  for (cell_bc in list_of_cell_bc) {
    
    # output single cell
    cell <- cb_sb_spatial_data_filtered_by_total %>% filter(cell_bc_10x == cell_bc)
    
    # run dbscan
    db <- dbscan::dbscan(cell[3:4], eps = eps, minPts = minPts, weights = cell$nUMI)
    
    # collect dbscan cluster assignments
    cell$cluster <- db$cluster
    number_clusters <- max(cell$cluster)
    cell$cluster <- as.factor(cell$cluster)
    
    # collect dbscan metrics
    SB_total <- length(cell$matched_beadbarcode)
    SB_UMI_total <- sum(cell$nUMI)
    cell_noise <- cell %>% filter(cluster == 0)
    SB_noise <- length(cell_noise$matched_beadbarcode)
    SB_UMI_noise <- sum(cell_noise$nUMI)
    
    # assign cell coordinates as centroid of top cluster (later i filter out cells with > 1 cluster)
    # otherwise, assign coordinates as NA
    cell_top_cluster <- cell %>% filter(cluster == 1)
    
    if (dim(cell_top_cluster)[1] > 0) {
      
      # cell_sp <- SpatialPointsDataFrame(cell_top_cluster[3:4], cell_top_cluster)
      # cell_sp_df <- as.data.frame(cell_sp)
      # centroid <- wt.centroid(st_as_sf(cell_sp), "nUMI", FALSE)
      # centroid_x <- centroid[1]
      # centroid_y <- centroid[2]
      
      cell_sp <- SpatialPointsDataFrame(cell_top_cluster[3:4], cell_top_cluster)
      cell_sp_df <- as.data.frame(cell_sp)
      centroid <- wt.centroid(st_as_sf(cell_sp), "nUMI", FALSE)
      centroid_x <- centroid$X
      centroid_y <- centroid$Y
      
      SB_top_cluster <- length(cell_top_cluster$matched_beadbarcode)
      SB_UMI_top_cluster <- sum(cell_top_cluster$nUMI)
      
      dbscan_df <- rbind(dbscan_df, cell_sp_df)
      
    } else {
      
      centroid_x <- NA
      centroid_y <- NA
      SB_top_cluster <- NA
      SB_UMI_top_cluster <- NA
      
    }
    
    # add dbscan parameters to vector
    cell_bc_vec <- append(cell_bc_vec, paste0(cell_bc))
    centroid_x_vec <- c(centroid_x_vec, paste0(centroid_x))
    centroid_y_vec <- c(centroid_y_vec, paste0(centroid_y))
    number_clusters_vec <- c(number_clusters_vec, paste0(number_clusters))
    SB_total_vec <- c(SB_total_vec, paste0(SB_total))
    SB_UMI_total_vec <- c(SB_UMI_total_vec, paste0(SB_UMI_total))
    SB_noise_vec <- c(SB_noise_vec, paste0(SB_noise))
    SB_UMI_noise_vec <- c(SB_UMI_noise_vec, paste0(SB_UMI_noise))
    SB_top_cluster_vec <- c(SB_top_cluster_vec, paste0(SB_top_cluster))
    SB_UMI_top_cluster_vec <- c(SB_UMI_top_cluster_vec, paste0(SB_UMI_top_cluster))
    
  }
  
  # aggregate dbscan information into dataframe
  coordinates_df <- data.frame(cell_bc_vec, centroid_x_vec, centroid_y_vec, number_clusters_vec, 
                               SB_total_vec, SB_UMI_total_vec, SB_noise_vec, SB_UMI_noise_vec,
                               SB_top_cluster_vec, SB_UMI_top_cluster_vec)
  
  # remove duplicates
  coordinates_df <- coordinates_df %>% unique()
  
  # rename columns 
  names(coordinates_df) <- c("cell_bc", "x_um", "y_um", "number_clusters", "SB_total", 
                             "SB_UMI_total", "SB_noise", "SB_UMI_noise", "SB_top_cluster",
                             "SB_UMI_top_cluster")
  
  # make columns numeric as appropriate
  coordinates_df$x_um <- as.numeric(coordinates_df$x_um)
  coordinates_df$y_um <- as.numeric(coordinates_df$y_um)
  coordinates_df$number_clusters <- as.numeric(coordinates_df$number_clusters)
  coordinates_df$SB_total <- as.numeric(coordinates_df$SB_total)
  coordinates_df$SB_UMI_total <- as.numeric(coordinates_df$SB_UMI_total)
  coordinates_df$SB_noise <- as.numeric(coordinates_df$SB_noise)
  coordinates_df$SB_UMI_noise <- as.numeric(coordinates_df$SB_UMI_noise)
  coordinates_df$SB_top_cluster <- as.numeric(coordinates_df$SB_top_cluster)
  coordinates_df$SB_UMI_top_cluster <- as.numeric(coordinates_df$SB_UMI_top_cluster)
  
  # replace NA with 0 in appropriate columns
  coordinates_df$SB_top_cluster[is.na(coordinates_df$SB_top_cluster)] = 0
  coordinates_df$SB_UMI_top_cluster[is.na(coordinates_df$SB_UMI_top_cluster)] = 0
  
  # create column with proportion SB noise and proportion SB UMI noise
  coordinates_df <- coordinates_df %>% 
    mutate(proportion_SB_noise = SB_noise / SB_total) %>%
    mutate(proportion_SB_UMI_noise = SB_UMI_noise / SB_UMI_total) %>%
    mutate(proportion_SB_top_cluster = SB_top_cluster / SB_total) %>%
    mutate(proportion_SB_UMI_top_cluster = SB_UMI_top_cluster / SB_UMI_total)
  

  
  ############### PLOT DBSCAN RESULTS SUMMARY ###############
  
  # Plot cells with coordinates in space 
  p<-ggplot(coordinates_df, aes(x = x_um, y = y_um)) + 
    geom_point(alpha = 0.5) + 
    coord_fixed(ratio = 1) + 
    theme_few() +
    labs(colour = "DBSCAN Clusters") +
    ggtitle("No cluster filter")
  ggsave(paste(general_plots_path, "spatial_mapping_no_cluster_filter.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  p <- ggplot(coordinates_df %>% filter(number_clusters == 1), aes(x = x_um, y = y_um)) + 
    geom_point(alpha = 0.5) + 
    coord_fixed(ratio = 1) + 
    theme_few() +
    labs(colour = "DBSCAN Clusters") +
    ggtitle("Single DBSCAN clusters")
  ggsave(paste(general_plots_path, "spatial_mapping.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot number of clusters histogram
  p<-ggplot(coordinates_df, aes(x = number_clusters)) + 
    geom_histogram(colour = "black", fill = "grey", bins = 50) + 
    theme_few() +
    xlab("DBSCAN Clusters") +
    ylab("Count")
  ggsave(paste(general_plots_path, "dbscan_clusters_histogram.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot SB_total distribution
  p<-ggplot(coordinates_df, aes(x = SB_total)) + 
    geom_histogram(colour = "black", fill = "grey", bins = 50) + 
    theme_few() +
    xlab("Unique Spatial Barcodes per Cell") +
    ylab("Count")
  ggsave(paste(general_plots_path, "SB_unique_per_cell.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot SB_UMIs_total distribution
  p<-ggplot(coordinates_df, aes(x = SB_UMI_total)) + 
    geom_histogram(colour = "black", fill = "grey", bins = 50) + 
    theme_few() +
    xlab("Spatial Barcode UMIs per Cell") +
    ylab("Count")
  ggsave(paste(general_plots_path, "SB_UMI_per_cell.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot proportion of SB in top cluster
  p<-ggplot(coordinates_df, aes(x = "", y = proportion_SB_top_cluster)) + 
    geom_violin(colour = "black", fill = "grey") + 
    geom_jitter(alpha = 0.5) +
    theme_few() +
    xlab("") +
    ylab("Proportion Unique Spatial Barcodes in Top Cluster") + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  ggsave(paste(general_plots_path, "top_cluster_signal_to_noise_unique.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot proportion of SB UMI in top cluster
  p<-ggplot(coordinates_df, aes(x = "", y = proportion_SB_UMI_top_cluster)) + 
    geom_violin(colour = "black", fill = "grey") + 
    geom_jitter(alpha = 0.5) +
    theme_few() +
    xlab("") +
    ylab("Proportion Spatial Barcode UMIs in Top Cluster") + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  ggsave(paste(general_plots_path, "top_cluster_signal_to_noise_UMI.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot proportion of SB noise
  p<-ggplot(coordinates_df, aes(x = "", y = proportion_SB_noise)) + 
    geom_violin(colour = "black", fill = "grey") + 
    geom_jitter(alpha = 0.5) +
    theme_few() +
    xlab("") +
    ylab("Proportion Unique Spatial Barcodes Denoted Noise") + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  ggsave(paste(general_plots_path, "proportion_unique_SB_noise.pdf", sep = ""), plot = p, width = 6, height = 4)
  
  # Plot proportion of SB UMI noise
  p<-ggplot(coordinates_df, aes(x = "", y = proportion_SB_UMI_noise)) + 
    geom_violin(colour = "black", fill = "grey") + 
    geom_jitter(alpha = 0.5) +
    theme_few() +
    xlab("") +
    ylab("Proportion Spatial Barcode UMIs Denoted Noise") + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  ggsave(paste(general_plots_path, "proportion_SB_UMI_noise.pdf", sep = ""), plot = p, width = 6, height = 4)
  

  
  ############### DBSCAN SUMMARY FILE ###############
  
  # assemble summary metrics for summary output file
  eps
  minPts
  
  # number cells mapped
  spatially_mapped_top_cluster <- (coordinates_df %>% filter(number_clusters == 1) %>% count())$n
  spatially_mapped_any_cluster_number <- (coordinates_df %>% filter(number_clusters > 0) %>% count())$n
  total_cells <- dim(coordinates_df)[1]
  proportion_cells_mapped <- (spatially_mapped_top_cluster/total_cells)[1]
  
  # sb per cell
  median_unique_SB_per_cell <- median(coordinates_df$SB_total)
  mean_unique_SB_per_cell <- mean(coordinates_df$SB_total)
  median_SB_UMI_per_cell <- median(coordinates_df$SB_UMI_total)
  mean_SB_UMI_per_cell <- mean(coordinates_df$SB_UMI_total)
  
  # signal to noise ratios total
  median_proportion_unique_SB_top_cluster_total <- median(coordinates_df$proportion_SB_top_cluster)
  mean_proportion_unique_SB_top_cluster_total <- mean(coordinates_df$proportion_SB_top_cluster)
  median_proportion_SB_UMI_top_cluster_total <- median(coordinates_df$proportion_SB_UMI_top_cluster)
  mean_proportion_SB_UMI_top_cluster_total <- mean(coordinates_df$proportion_SB_UMI_top_cluster)
  
  # signal to noise ratios for mapped cells
  median_mapped_cells_proportion_unique_SB_top_cluster <- (coordinates_df %>% 
                                                             filter(number_clusters == 1) %>% 
                                                             select(proportion_SB_top_cluster) %>% 
                                                             summarise(across(everything(), median)))$proportion_SB_top_cluster
  mean_mapped_cells_proportion_unique_SB_top_cluster <- (coordinates_df %>% 
                                                           filter(number_clusters == 1) %>% 
                                                           select(proportion_SB_top_cluster) %>% 
                                                           summarise(across(everything(), mean)))$proportion_SB_top_cluster
  median_mapped_cells_proportion_SB_UMI_top_cluster <- (coordinates_df %>% 
                                                          filter(number_clusters == 1) %>% 
                                                          select(proportion_SB_UMI_top_cluster) %>% 
                                                          summarise(across(everything(), median)))$proportion_SB_UMI_top_cluster
  mean_mapped_cells_proportion_SB_UMI_top_cluster <- (coordinates_df %>% 
                                                        filter(number_clusters == 1) %>% 
                                                        select(proportion_SB_UMI_top_cluster) %>% 
                                                        summarise(across(everything(), mean)))$proportion_SB_UMI_top_cluster
  
  summary_metrics <- data.frame(eps, 
                                minPts, 
                                nUMI_cutoff_combined,
                                spatially_mapped_top_cluster,
                                spatially_mapped_any_cluster_number,
                                total_cells,
                                proportion_cells_mapped,
                                median_unique_SB_per_cell,
                                mean_unique_SB_per_cell,
                                median_SB_UMI_per_cell,
                                mean_SB_UMI_per_cell,
                                median_proportion_unique_SB_top_cluster_total,
                                mean_proportion_unique_SB_top_cluster_total,
                                median_proportion_SB_UMI_top_cluster_total,
                                mean_proportion_SB_UMI_top_cluster_total,
                                median_mapped_cells_proportion_unique_SB_top_cluster,
                                mean_mapped_cells_proportion_unique_SB_top_cluster,
                                median_mapped_cells_proportion_SB_UMI_top_cluster,
                                mean_mapped_cells_proportion_SB_UMI_top_cluster
  )
  
  write.csv(summary_metrics, paste(general_plots_path,"summary_metrics.csv"))
  
  

  ############### WRITE COORDINATES DATAFRAME ###############
  
  write.table(coordinates_df, coordinates_df_output_path)

}


