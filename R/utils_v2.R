### SpatialQC - update


## Loading packages ---------
packages <- c("Seurat", "dplyr", "ggplot2", "shadowtext", "scales", "cowplot", "data.table", "Matrix", "matrixStats", "SingleCellExperiment", "SpatialExperiment", "SpatialFeatureExperiment", "bluster", "BiocParallel", "BioQC", "coop", "fs", "stringr")

# Check and install packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install the package from CRAN
    install.packages(pkg)
    
    # Check if it's a Bioconductor package
    if (pkg %in% c("SingleCellExperiment", "SpatialExperiment", "SpatialFeatureExperiment", "bluster", "BiocParallel", "BioQC")) {
      # Install the Bioconductor package
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    }
    
    # Load the package
    library(pkg, character.only = TRUE)
  } else {
    # Load the package
    library(pkg, character.only = TRUE)
  }
}


## Define auxiliary functions ---------
Read10X_h5_gz <- function(filename, use.names = TRUE, unique.features = TRUE) {
  # Check if the file is gzipped
  if (grepl("\\.gz$", filename)) {
    # Create a temporary file for decompressed data
    temp_file <- tempfile(fileext = ".h5")
    
    # Decompress the gzipped file
    gunzip(filename, destname = temp_file, remove = FALSE)
    
    # Read the decompressed HDF5 file
    data <- Read10X_h5(temp_file, use.names = use.names, unique.features = unique.features)
    
    # Clean up: remove the temporary file
    unlink(temp_file)
    
    return(data)
  } else {
    # If not gzipped, read directly
    return(Read10X_h5(filename, use.names = use.names, unique.features = unique.features))
  }
}



# Function to calculate size in GB with logging
size_in_gb <- function(path) {
  if (length(path) == 0) {
    return(NA)  # Return NA if path is not available
  }
  
  file_size <- file.size(path)
  if (file_size < 0) {
    return(NA)
  }
  
  size_gb <- round(file_size / (1024^3), 3)  # Convert bytes to gigabytes
  return(size_gb)
}






## Define Metrics functions ---------

getGlobalFDR <- function(seu_obj = NULL, 
                         features = NULL, 
                         tx_file = 'path_to_txFile', 
                         cellSegMeta = 'path_to_cellMeta', 
                         platform = NULL) {
  # Initialize variable
  tx_df <- NULL
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      tx_df <- data.table::fread(tx_file, header = TRUE)
      data.table::setnames(tx_df, "feature_name", "target")
      
    } else if (platform == "CosMx") {
      tx_df <- data.table::fread(tx_file, header = TRUE)
      
      
    } else if (platform == "Merscope") {
      tx_df <- data.table::fread(tx_file, header = TRUE)
      tx_df$target <- tx_df$gene
      
    }
    
    # Filter and process data
    negProbes <- tx_df$target[grep('Neg*|Blank*|BLANK*', tx_df$target)]
    allGenes <- unique(tx_df[!target %in% negProbes, target])  # List of unique genes (non-control or blank probes) in panel
    
    # Create table with expression per each gene in panel
    expTableAll <- tx_df[, .(Count = .N), by = target]
    expTable <- expTableAll[expTableAll$target %in% allGenes, ]
    expNeg <- sum(expTableAll[expTableAll$target %in% unique(negProbes), ]$Count)  # Sum of all negative control or blank or unassigned barcodes (i.e. non-specific)
    
    numGenes <- length(expTable$target)
    numNeg <- length(expTableAll[expTableAll$target %in% unique(negProbes), ]$target)
    
    fdr <- (expNeg / (sum(expTable$Count) + expNeg)) * (numGenes / numNeg) * 1 / 100
    
    if (is.null(features)) {
      return(fdr)
    } else {
      fdr <- (expNeg / (sum(expTable[expTable$target %in% features, ]$Count) + expNeg)) * (numGenes / numNeg) * 1 / 100
      return(fdr)
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
    
    return()
    
  }
}  
  
  


getTxPerArea <- function(seu_obj = NULL, 
                         features = NULL,
                         platform = NULL, 
                         cellSegMeta = 'path_to_cellMeta', 
                         tx_file = NULL, 
                         exp_file = 'path_to_exp_mtx') {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      cell_meta <- data.table::fread(file.path(cellSegMeta))
      mean_tx_norm <- mean(cell_meta$total_counts / cell_meta$cell_area)
      
      # If features are specified - have to calculate differently
      if (!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # Subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[tx_df$cell_id %in% cell_meta$cell_id, ]
        tx_df <- tx_df[tx_df$feature_name %in% features, ]
        # Count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell_id')
        tx_df <- tx_df %>% group_by(cell_id, cell_area) %>% tally()
        
        mean_tx_norm <- mean(tx_df$n / tx_df$cell_area)
      }
      
    } else if (platform == "CosMx") {
      cell_meta <- data.table::fread(file.path(cellSegMeta))
      cell_meta <- data.frame(cell_meta)
      
      if ("nCount_RNA" %in% names(cell_meta) == FALSE) {
        # Creating nCount_RNA column
        nCount_RNA <- cell_meta %>%
          group_by(cell_ID) %>%
          tally() %>%
          data.frame()
        names(nCount_RNA)[2] <- "nCount_RNA"
        
        cell_meta <- merge(cell_meta, nCount_RNA, by = "cell_ID")
      }
      
      if ("Area.um2" %in% names(cell_meta) == TRUE) {
        names(cell_meta)[names(cell_meta$Area) %in% "Area.um2"] <- "Area"
      }
      
      mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area)
      
      if (!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # Subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[cell_ID != 0 & target %in% features]
        # Count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell')
        tx_df <- tx_df %>% group_by(cell, Area.um2) %>% tally()
        
        mean_tx_norm <- mean(tx_df$n / tx_df$Area.um2)
      }
      
    } else if (platform == "Merscope") {
      if (any(grepl("metadata", cellSegMeta))) {
        cell_meta <- data.table::fread(file.path(cellSegMeta)) %>% data.frame()
        #names(cell_meta)[names(cell_meta) == "V1"] <- "cell_ID"
        names(cell_meta)[1] <- "cell_ID"
        
        # Function to calculate area
        calculate_area <- function(min_x, max_x, min_y, max_y) {
          width <- max_x - min_x
          height <- max_y - min_y
          area <- width * height
          return(area)
        }
        
        # Apply function to each row (assuming all necessary values are filled)
        cell_meta$Area.um2 <- mapply(calculate_area,
                                     cell_meta$min_x,
                                     cell_meta$max_x,
                                     cell_meta$min_y,
                                     cell_meta$max_y)
        
        if ("Area.um2" %in% names(cell_meta) == TRUE) {
          names(cell_meta)[names(cell_meta) == "Area.um2"] <- "Area"
        }
        
        # Function to calculate nCount_RNA
        exp_df <- data.table::fread(exp_file)
        exp_df <- data.frame(exp_df)
        rownames(exp_df) <- exp_df[, 1]
        exp_df[, 1] <- NULL
        exp_df <- t((exp_df))
        
        nCount_RNA <- colSums(exp_df)
        nCount_RNA <- data.frame(cell_ID = colnames(exp_df), nCount_RNA = nCount_RNA)
        
        cell_meta <- merge(cell_meta, nCount_RNA, by = "cell_ID")
        
        # Calculate mean tx
        mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area)
        
        if (!is.null(features)) {
          exp_df <- exp_df[rownames(exp_df) %in% features, ]
          
          nCount_RNA <- colSums(exp_df)
          nCount_RNA <- data.frame(cell_ID = colnames(exp_df), nCount_RNA = nCount_RNA)
          cell_meta$nCount_RNA <- NULL
          cell_meta <- merge(cell_meta, nCount_RNA, by = "cell_ID")
          
          # Calculate mean tx
          mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area)
          
        }
      }
      
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(mean_tx_norm)
}
  


getTxPerNuc <- function(seu_obj = NULL, 
                        features = NULL, 
                        tx_file = NULL, 
                        platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      tx_df <- data.table::fread(tx_file)
      
      if (is.null(features)) {
        # Subset the tx file with only existing assigned cells and features
        # Remove neg control probes
        negProbes <- unique(tx_df$feature_name[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$feature_name)])
        # Number of txs in nucleus
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes])[1]
        # Number of cells
        nCells <- length(unique(tx_df$cell_id))
      } else {
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes & feature_name %in% features])[1]
      }
      
    } else if (platform == "CosMx") {
      tx_df <- data.table::fread(tx_file)
      
      if (is.null(features)) {
        # Subset the tx file with only existing assigned cells and features
        # Remove neg control probes
        negProbes <- unique(tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$target)])
        # Number of txs in nucleus
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes])
        # Number of cells
        nCells <- length(unique(tx_df$cell[tx_df$cell_ID != 0]))
      } else {
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes & target %in% features])
      }
      
    } else if (platform == "Merscope") {
      tx_df <- data.table::fread(tx_file)
      tx_df$target <- tx_df$gene
      tx_df$cell_ID <- tx_df$barcode_id
      
      if (is.null(features)) {
        # Subset the tx file with only existing assigned cells and features
        # Remove neg control probes
        negProbes <- unique(tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$target)])
        # Number of txs in nucleus
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & !target %in% negProbes])
        # Number of cells
        nCells <- length(unique(tx_df$cell_ID[tx_df$cell_ID != 0]))
      } else {
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & !target %in% negProbes & target %in% features])
      }
      
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(nTx_nuc / nCells)
}





getEntropy <- function(seu_obj = NULL, 
                       features = NULL, 
                       expMat = 'path_to_expMat', 
                       platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      exp <- exp[, -c(1:2)]
      exp <- t(exp)
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      exp <- exp[, -c(1:2)]
      exp <- t(exp)
      
    }
    
    if (is.null(features)) {
      entropy_out <- BioQC::entropy(as.matrix(exp))
    } else {
      entropy_out <- BioQC::entropy(as.matrix(exp[rownames(exp) %in% features, ]))
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(entropy_out)
}






getMeanSignalRatio <- function(seu_obj = NULL, 
                               features = NULL, 
                               platform = NULL, 
                               expMat = 'path_to_expMat') {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- t(exp)
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- t(exp)
      
    }
    
    noise <- exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ]
    exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ]
    
    if (is.null(features)) {
      MeanSignalRatio_out <- suppressWarnings(mean(log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1))))
    } else { 
      MeanSignalRatio_out <- suppressWarnings(mean(log10(rowMeans(exp[rownames(exp) %in% features, ] + .1)) - log10(rowMeans(noise + .1))))
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(MeanSignalRatio_out)
}




getNcells <- function(seu_obj = NULL, 
                      expMat = 'path_to_expMat', 
                      platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
        ncell <- ncol(exp)
        
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
        ncell <- ncol(exp)
      }
    } else if (platform == "CosMx") {
      ncell <- nrow(data.table::fread(expMat))
      
    } else if (platform == "Merscope") {
      ncell <- nrow(data.table::fread(expMat))
      
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(ncell)
}



getTxPerCell <- function(seu_obj = NULL, 
                         features = NULL, 
                         expMat = 'path_to_exprMatrix', 
                         platform) { 
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      exp <- t(exp)  ## Transposing for consistency - row = genes, column = cells
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat)) %>% data.frame()
      rownames(exp) <- exp[, 1]
      exp[, 1] <- NULL
      exp <- t(exp)  ## Transposing for consistency - row = genes, column = cells
      
    }
    
    if (is.null(features)) {
      features <- rownames(exp)
      # Remove non-specific probes 
      features <- features[!grepl('Unassigned*|NegControl*|BLANK*|Blank*|SystemControl*|NegPrb*', features)] 
    }
    
    # Calculate average number of transcripts per cell 
    mean_tx <- mean(colSums(exp[features, ]))
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(mean_tx)
}




getCellTxFraction <- function(seu_obj = NULL, 
                              features = NULL, 
                              tx_file = 'path_to_tx', 
                              platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    tx_df <- data.table::fread(tx_file, header = TRUE)
    total_tx_count <- nrow(tx_df)
    
    if (platform == 'Xenium') {
      if (is.null(features)) {
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')
      } else {
        tx_df <- tx_df[tx_df$feature_name %in% features, ]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')
      }
      
    } else if (platform == "CosMx") {
      if (is.null(features)) {
        unassigned_tx_count <- sum(tx_df$CellComp == 'None')
        total_tx_count <- nrow(tx_df)
      } else {
        tx_df <- tx_df[tx_df$target %in% features, ]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(tx_df$CellComp == 'None')
      }
      
    } else if (platform == "Merscope") {
      tx_df$target <- tx_df$gene
      if (is.null(features)) {
        unassigned_tx_count <- sum(grepl("Neg|SystemControl|Blank|UNASSIGNED", tx_df$target))
        total_tx_count <- nrow(tx_df)
      } else {
        tx_df <- tx_df[tx_df$target %in% features, ]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(grepl("Neg|SystemControl|Blank|UNASSIGNED", tx_df$target))
      }
      
    }
    
    CellTxFraction_out <- (total_tx_count - unassigned_tx_count) / total_tx_count 
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(CellTxFraction_out)
}




getMaxRatio <- function(seu_obj = NULL, 
                        features = NULL, 
                        expMat = 'path_to_expMat', 
                        platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      exp <- t(exp)  ## Transposing for consistency - row = genes, column = cells
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      exp <- t(exp)  ## Transposing for consistency - row = genes, column = cells
      
    }
    
    tx_means <- rowMeans(exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ])
    neg_probe_means <- rowMeans(exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ])
    
    if (is.null(features)) {
      MaxRatio_out <- log10(max(tx_means)) - log10(mean(neg_probe_means))
      
    } else {
      MaxRatio_out <- log10(max(tx_means[features])) - log10(mean(neg_probe_means)) 
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(MaxRatio_out)
}



getMECR <- function(seu_obj = NULL, 
                    expMat = 'path_to_expMat', 
                    platform = NULL) { 
  # This function comes from Hartman & Satija, bioRxiv, 2024
  # We are using a custom marker table. The original publication bases it on
  # scRNA-seq from matched tissue.
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## Transposing for consistency - row = genes, column = cells
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## Transposing for consistency - row = genes, column = cells
      
    }
    
    # Detect Mouse dataset
    nMouse_genes <- rownames(exp)[grepl("^[A-Z][a-z0-9]*$", rownames(exp))] %>% length() 
    nGenes_total <- rownames(exp)[!grepl("Codeword|Probe|Blank", rownames(exp))] %>% length()
    
    if (nMouse_genes / nGenes_total >= 0.3) {
      marker_df <- data.frame(
        gene = c("EPCAM", "KRT19", "KRT8", 
                 "CD3E", "CD3D", "CD8A", "NKG7",
                 "MS4A1", "CD79A",
                 "PECAM1", "CLDN5", "VWF",
                 "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
                 "PDGFRA", "DPT", "COL1A1",
                 "MYH11", "ACTG2"),
        cell_type = c("Epithelial", "Epithelial", "Epithelial",
                      "T", "T", "T", "T",
                      "B", "B",
                      "Endo", "Endo", "Endo",
                      "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                      "Fibro", "Fibro", "Fibro",
                      "Muscle", "Muscle"))
      rownames(marker_df) <- marker_df$gene
      marker_df$gene <- tools::toTitleCase((tolower(marker_df$gene)))
      rownames(marker_df) <- marker_df$gene
      
    } else {
      marker_df <- data.frame(
        gene = c("EPCAM", "KRT19", "KRT8", 
                 "CD3E", "CD3D", "CD8A", "NKG7",
                 "MS4A1", "CD79A",
                 "PECAM1", "CLDN5", "VWF",
                 "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
                 "PDGFRA", "DPT", "COL1A1",
                 "MYH11", "ACTG2"),
        cell_type = c("Epithelial", "Epithelial", "Epithelial",
                      "T", "T", "T", "T",
                      "B", "B",
                      "Endo", "Endo", "Endo",
                      "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                      "Fibro", "Fibro", "Fibro",
                      "Muscle", "Muscle"))
      rownames(marker_df) <- marker_df$gene
    }
    
    genes <- intersect(rownames(exp), rownames(marker_df))
    mtx <- as.matrix(exp[genes, ])
    
    coexp.rates <- c()
    # print(paste0("Marker count: ", length(genes)))
    if (length(genes) > 25) { 
      genes <- sample(genes, 25) 
    }
    
    for (g1 in genes) {
      for (g2 in genes) {
        if ((g1 != g2) && (g1 > g2) && (marker_df[g1, 'cell_type'] != marker_df[g2, 'cell_type'])) {
          c1 <- mtx[g1, ]
          c2 <- mtx[g2, ]
          coexp.rates <- c(
            coexp.rates,
            sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0)) # >0 too liberal of an expression threshold?
        }
      }
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == 'Xenium') {
      #print("")
      
    }
    if (platform == 'CosMx') {
      #print("")
      
    }
    if (platform == 'Merscope') {
      #print("")
      
    }
  }
  
  return(round(mean(coexp.rates), digits = 3))
}




getSparsity <- function(seu_obj = NULL, 
                        features = NULL, 
                        expMat = 'path_to_expMat', 
                        platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- t(exp)
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      exp <- t(exp)
      
    }
    
    if (is.null(features)) {
      Sparsity_out <- coop::sparsity(as.matrix(exp))
    } else {
      Sparsity_out <- coop::sparsity(as.matrix(exp[rownames(exp) %in% features, ]))
    }
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  return(Sparsity_out)
}



getComplexity <- function(seu_obj = NULL, 
                          features = NULL, 
                          expMat = 'path_to_expMat', 
                          platform = NULL) {
  
  ### Seurat Object OFF   
  if (is.null(seu_obj)) {
    
    if (platform == 'Xenium') {
      expMatDir <- dirname(expMat)
      mtx_bar_feat_path <- fs::dir_ls(expMatDir, recurse = TRUE, type = "file")
      
      if (grepl('.h5', expMat) == FALSE) {
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz$', mtx_bar_feat_path)] 
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz$', mtx_bar_feat_path)] 
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz$', mtx_bar_feat_path)]
        
        # Check if files exist
        if (length(mtx_path) == 0 || length(bar_path) == 0 || length(feat_path) == 0) {
          stop("Required files not found.")
        }
        
        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = FALSE)
        rows <- data.table::fread(file.path(feat_path), header = FALSE)
        rownames(exp) <- rows$V2  # Gene symbols
        colnames(exp) <- cols$V1   # Barcodes
      } else {
        # Handle .h5 file case
        mtx_h5_path <- expMat
        
        if (length(mtx_h5_path) == 0) {
          stop("No .h5 files found in the specified directory.")
        }
        
        exp <- Read10X_h5_gz(filename = file.path(mtx_h5_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp)
      }
    } else if (platform == "CosMx") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      # exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## Transposing for consistency - row = genes, column = cells
      
    } else if (platform == "Merscope") {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1", "sampleID", "slide", "case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## Remove first 2 columns - usually FOV and Cell_ID information
      # exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## Transposing for consistency - row = genes, column = cells
      
    }
    
    total_sum <- sum(exp)
    
    if (is.null(features)) {
      row_sums <- rowSums(exp)
      cumulative_sums <- cumsum(row_sums)
      Complexity_out <- as.integer(which(cumulative_sums >= total_sum / 2)[1])
      
    } else {
      total_sum_feat <- sum(exp[rownames(exp) %in% features, ])
      row_sums <- rowSums(as.data.frame(exp[rownames(exp) %in% features, ]))
      cumulative_sums <- cumsum(row_sums)
    }
    
    result <- which(cumulative_sums >= total_sum / 2)[1]
    #weigh <- 31
    
    ### Seurat Object ON   
  } else {  # Check if sobj is NOT NULL
    if (platform == "Xenium") {
      #print("")
      
    }
    if (platform == "CosMx") {
      #print("")
      
    }
    if (platform == "Merscope") {
      #print("")
      
    }
  }
  
  # return( as.integer(result) / (weigh * ( nrow(exp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)))) ) )  
  return(as.integer(result) / (nrow(exp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)))))
}




##### Function Skeleton 
# getTxPerNuc <- function(seu_obj = NULL, 
#                         features = NULL, 
#                         tx_file = NULL, 
#                         platform = NULL) {
#   
#   ### Seurat Object OFF   
#   if (is.null(seu_obj)) {
#     
#     if (platform == 'Xenium') {
#       
#       
#     } else if (platform == "CosMx") {
#       
#       
#     } else if (platform == "Merscope") {
#       
#       
#     }
#     
#     ### Seurat Object ON   
#   } else {  # Check if sobj is NOT NULL
#     if (platform == "Xenium") {
#       #print("")
#       
#     }
#     if (platform == "CosMx") {
#       #print("")
#       
#     }
#     if (platform == "Merscope") {
#       #print("")
#       
#     }
#   }
#   
#   return()
# }






# What is this ? ----------
# if(is.null(features)){
#   features <- rownames(seu_obj)
# } else{
#   features <- features
# }
# 
# path <- unique(seu_obj$path)
# platform <- unique(seu_obj$platform)
#----


# What is this? ----
# # Read Tx localization data
# tx_df <- readTxMeta(path, platform)
# 
# if(platform == "Xenium"){
#   tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) &
#                     overlaps_nucleus == 1 &
#                     features %in% features) %>%
#     group_by(cell_id) %>%
#     summarize(nuc_counts = n())
#   
#   
# } else if(platform == "CosMx"){
#   tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
#   tx_df <- tx_df %>%
#     filter(cell_id %in% colnames(seu_obj) &
#              CellComp == "Nuclear" &
#              target %in% features) %>%
#     group_by(cell_id) %>%
#     summarize(nuc_counts = n())
#   
# } else if(platform == "Merscope"){
#   print("Working on support")
#   
# } else{
#   print("Platform not supported")
# }
# 
# res <- data.frame(
#   sample_id = unique(seu_obj$sample_id),
#   platform = unique(seu_obj$platform),
#   value=mean(tx_df$nuc_counts)
# )
# 
# return(res)
#----




      
      
