# modify from https://github.com/ChengxiangQiu/tome_code/blob/main/help_code/help_code.R
createLineage_Knn <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5){
    print(unique(pd$stage))
    print(dim(emb))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
    pd$state = pd$Anno
    
    res = list()
    
    rep_i = 1
    
    while(rep_i < (replication_times+1)){
        
        sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))
        
        emb_sub = emb[sampling_index,]
        pd_sub = pd[sampling_index,]
        
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
        
        pre_state_min = min(table(as.vector(pd_sub1$state)))
        
        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min
            print(k_neigh)
        }
        
        if (pre_state_min < 3){
            next
        }
        
        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
        
        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
            tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
        
        tmp2 <- matrix(NA,length(state2),length(state1))
        for(i in 1:length(state2)){
            x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
            for(j in 1:length(state1)){
                tmp2[i,j] <- sum(x==state1[j])
            }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1
        
        res[[rep_i]] = tmp2
        
        rep_i = rep_i + 1
        
    }
    
    return(res)
}


#' Subway_map.step_1.UMAP_list
#' 
#' @param archrPig input ArchR object
#' @param useMatrix used matrix for UMAP
#' @param order_date order of the developmental stages
#' @param batch_group column name for batch correction
# modify from https://github.com/ChengxiangQiu/tome_code/blob/main/help_code/help_code.R
Subway_map.step_1.UMAP_list.seurat <- function(
archrPig,
useMatrix = "PeakMatrix",
order_date = c("E14","E18","E43","E55","E93","7d","14d","60d","180d"),
batch_group = "Sample",
regionsList = NULL,
subset_by = NULL,
subset_count = NULL,
useFeatureCutoff = "q50",
k_anchor = 6,
k_weight = 100,
min_dist = 0.3,
seeds = 1234,
saveSeurat = F,
plotlist = NULL,
outname = "test"
){
	library(Seurat)
	library(Signac)
	library(future.apply)
	set.seed(seeds)
	## check input
	if(length(grep("q", useFeatureCutoff)) == 0){
		useFeatureCutoff <- as.numeric(useFeatureCutoff)
	}
	#### get UMAP for adjacent time points
	order_date <- order_date %>% intersect(., unique(archrPig@cellColData$Date))
	umap_list <- list()
	system(paste("rm -rf", paste0(outname, ".used.Kweight.txt")))
	for(date_idx in 1:(length(order_date)-1)){
		print(paste("#################", order_date[date_idx],order_date[date_idx+1],sep = ":"))
		obj <- subset(archrPig@cellColData,Date %in% c(order_date[date_idx],order_date[date_idx+1])) %>% rownames %>% archrPig[., ]
		if(!is.null(subset_by) & !is.null(subset_count)){
			obj1 <- subset(obj@cellColData,Date %in% c(order_date[date_idx])) %>% rownames %>% obj[., ] %>% subset_cells_archr(., subset_by, subset_count, 1234)
			obj2 <- subset(obj@cellColData,Date %in% c(order_date[date_idx+1])) %>% rownames %>% obj[., ] %>% subset_cells_archr(., subset_by, subset_count, 1234)
			obj <- obj[c(obj1, obj2),]
		}
		#saveRDS(obj@cellColData, paste0(outname, ".metadata.rds"))
		# get peak matrix
		archrOBJ.peakMTX <- getMatrixFromProject(ArchRProj = obj, useMatrix = useMatrix, useSeqnames = NULL, binarize = FALSE)
		region_posi <- as.data.frame(archrOBJ.peakMTX@rowRanges)	#only for peak matrix
		rownames(archrOBJ.peakMTX@assays@data[[1]]) <- paste(region_posi[,1],region_posi[,2],region_posi[,3],sep="-")
		# overlap peaks
		if(!is.null(regionsList)){
			unionpeaks <- reduce(c(regionsList[[order_date[date_idx]]], regionsList[[order_date[date_idx+1]]]), min.gapwidth = 0) %>% as.data.frame
			unionpeaks <- paste(unionpeaks[,1], unionpeaks[,2], unionpeaks[,3], sep = "-")
			archrOBJ.peakMTX@assays@data[[1]] <- archrOBJ.peakMTX@assays@data[[1]][intersect(rownames(archrOBJ.peakMTX@assays@data[[1]]), unionpeaks), ]
		}
		# generate signac object
		atac.assay <- CreateChromatinAssay(counts = archrOBJ.peakMTX@assays@data[[1]], sep = c("-", "-"), min.cells = 3, min.features = 100)
		pbmc.atac <- CreateSeuratObject(counts = atac.assay, assay = "ATAC", meta.data = as.data.frame(obj@cellColData))
		rm(obj, archrOBJ.peakMTX, region_posi, atac.assay)
		gc()
		# compute LSI
		pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = useFeatureCutoff)
		pbmc.atac <- RunTFIDF(pbmc.atac)
		pbmc.atac <- RunSVD(pbmc.atac)
		# compute LSI for split data
		obj.list <- SplitObject(object = pbmc.atac, split.by = batch_group)
		for (i in 1:length(x = obj.list)) {
			obj.list[[i]] <- FindTopFeatures(obj.list[[i]], min.cutoff = useFeatureCutoff)
			obj.list[[i]] <- RunTFIDF(obj.list[[i]])
			obj.list[[i]] <- RunSVD(obj.list[[i]])
		}
		# find  integration anchors
		integration.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = rownames(obj.list[[1]]), reduction = "rlsi", dims = 2:30, k.anchor = k_anchor)
		rm(obj.list)
		gc()
		## integrate LSI embeddings 
		# to get the number of max anchors (Error in idx[i, ] <- res[[i]][[1]] :number of items to replace is not a multiple of replacement length)
		if(length(integration.anchors@object.list) == 2){
			source('~/01.sxzhang/scripts/getMaxAnchor.R')
			#lapply(list.files(path = "~/01.sxzhang/01.pig.snatac/03.archr_process/03.3rd/Seurat.Rscript-4.4.0/01.SeuratR4Maxanchors", pattern = "\\.R$", full.names = TRUE), source)
			maxanchors <- IntegrateEmbeddings.xxx(anchorset = integration.anchors, reductions = pbmc.atac[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 1:30, k.weight = 100)
			if(maxanchors < k_weight){
				k_weight = maxanchors
				system(paste("echo", paste(order_date[date_idx], order_date[date_idx+1], k_weight, sep = "\t"), ">>", paste0(outname, ".used.Kweight.txt")))
			}
		}
		integrated <- IntegrateEmbeddings(anchorset = integration.anchors, reductions = pbmc.atac[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 1:30, k.weight = k_weight)
		rm(pbmc.atac)
		gc()
		# create a new UMAP using the integrated embeddings
		integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30, n.components = 3, min.dist = min_dist, seed.use = seeds)
		if(saveSeurat){
			saveRDS(integrated, paste0(outname, ".", paste(order_date[date_idx], order_date[date_idx+1], sep = "_"), ".seurat.rds"))
		}
		# plot
		Clustering_plot(integrated,plotlist = plotlist,col_cluster = plotlist[1],ggtitle = NULL,dotsize = 0.5,p.width = 14,p.height = 8,outprefix = paste0(outname, ".", paste(order_date[date_idx],order_date[date_idx+1],sep="_")))
		# get umap
		umap_list <- c(umap_list, list(data.frame(Embeddings(object = integrated, reduction = "umap"))))
		rm(integrated)
		gc()
	}
	names(umap_list) <- order_date[1:(length(order_date)-1)]
	saveRDS(umap_list, paste0(outname, ".adjacent_UMAP_list.rds"))
}

#### running Knn to find ancestor 
Subway_map.step_2.knn_list <- function(
umap_list,
metadata,
CelltypeN = "CellType",
DateN = "Date",
order_date = c("E14","E18","E43","E55","E93","7d","14d","60d","180d"),
seeds = 1234,
outname = "test"
){
	require(FNN)
	set.seed(seeds)
	# prepare metadata
	all_metadata <- readRDS(metadata)
	all_metadata$day <- all_metadata[, DateN]
	all_metadata$Anno <- all_metadata[, CelltypeN]
	all_metadata$stage <- all_metadata$day
	order_date <- order_date %>% intersect(., unique(all_metadata$day))
	# UMAP
	umap_list <- readRDS(umap_list)
	if(length(setdiff(order_date[-length(order_date)], names(umap_list))) > 0)
		stop(paste("Date name is/are not presect in UMAP list file:", setdiff(order_date, names(umap_list))))
	# KNN
	knn_list <- list()
	for(i in 1:(length(order_date)-1)){
		umapdf <- umap_list[[order_date[i]]]
		submeta <- all_metadata[rownames(umapdf), c("day","Anno","stage")]
		submeta$day <- mapvalues(submeta$day, c(order_date[i], order_date[i+1]), c("pre", "nex"))
		res = createLineage_Knn(umapdf, submeta,  k_neigh = 5)
		knn_list <- c(knn_list, list(res))
	}
	names(knn_list) <- order_date[1:(length(order_date)-1)]
	saveRDS(knn_list, paste0(outname, "_knn_list.rds"))

	#### creating the median value of the matrix
	knn_median_list <- list()
	for(date_i in 1:(length(order_date)-1)){
		replication_times=500
		dat = knn_list[[date_i]]
		state_1 = row.names(dat[[1]]) %>% paste0(order_date[date_i+1], ":", .)
		state_2 = names(dat[[1]]) %>% paste0(order_date[date_i], ":", .)
		tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
		for(i in 1:nrow(dat[[1]])){
			for(j in 1:ncol(dat[[1]])){
				xx = NULL
				for(k in 1:replication_times){
					xx = c(xx, dat[[k]][i,j])
				}
				tmp_1[i,j] = median(xx[!is.na(xx)])
			}
		}
		tmp_1 = data.frame(tmp_1)
		row.names(tmp_1) = state_1
		names(tmp_1) = state_2
		knn_median_list <- c(knn_median_list, list(tmp_1))
	}
	names(knn_median_list) <- names(knn_list)
	saveRDS(knn_median_list, paste0(outname, "_knn_median_list.rds"))
}


#### Step3_summarize_results ####
Subway_map.step_3.knn_res <- function(
knn_median_list,
edge_cutpff = 0.2,
outname = "test"
){
	library(reshape2)
	library(scales)
	library(ggplot2)
	library(dplyr)
	
	# read in knn data
	knn_median_list <- readRDS(knn_median_list)
	order_date <- names(knn_median_list)
	
	# generate edge data
	dat = NULL
	for(i in 1:length(knn_median_list)){
		print(order_date[i])
		dat = rbind(dat, melt(as.matrix(knn_median_list[[i]])))
	}
	dat = data.frame(dat)
	names(dat) = c("nex", "pre", "prob")
	dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
	dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
	dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
	dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))
	saveRDS(dat, paste0(outname, ".edge_all.rds"))

	print(paste0("how many edges: ", nrow(dat)))
	print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
	print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.2,])))
	print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
	print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
	print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
	print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))

	#### extract edges with prob > 0.2
	x = dat[dat$prob >= edge_cutpff,]
	x = x[,c("pre","nex","prob")]
	print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))

	#### introduce 1 "dummy nodes", corresponding to Dummy at Root (as a root)
	dummy = NULL
	for(first_date in grep(order_date[1], unique(x$pre), value = T)){
		dummy = rbind(dummy, c("Root:Dummy",first_date,1))
	}
	dummy = data.frame(dummy)
	names(dummy) = c("pre","nex","prob")
	res = rbind(x, dummy)
	
	# summary
	dat_sub = res
	dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
	dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))
	sum(dat_sub$pre_cell == dat_sub$nex_cell)
	sum(dat_sub$pre_cell == dat_sub$nex_cell)/nrow(dat_sub)
	sum(dat_sub$pre_cell != dat_sub$nex_cell)
	sum(dat_sub$pre_cell != dat_sub$nex_cell)/nrow(dat_sub)

	#### summary on finally used to create the tree

	print(paste0("how many edges: ", nrow(res)))
	print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))
	
	# checking lack terms
	# the lackterms means that the terms don't have linked to any cell types of previous time point du to the filtering, but any cell type should have previous links
	lackterms <- setdiff(unique(dat_sub$pre),unique(dat_sub$nex))
	print(paste("The terms that must be added into resulting 'txt' files:", paste(lackterms,collapse = ", ")))
	
	# edge info for tree plot
	nex_list = as.vector(unique(res$nex))
	tree = NULL
	for(i in 1:length(nex_list)){
		res_sub = res[res$nex==nex_list[i],]
		if(nrow(res_sub)==1){
			tree = rbind(tree, res_sub)
		} else {
			res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
			tree = rbind(tree, res_sub[1,])
		}
	}
	tree = data.frame(tree)
	tree = tree[,c("pre","nex")]
	
	write.table(res, paste0(outname, ".edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
	write.table(tree, paste0(outname, ".edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
}


#### Step3_summarize_results ####
#' @param edge edge file generated by 'Subway_map.step_3.knn_res' function
#' @param edge_prob edge prob file generated by 'Subway_map.step_3.knn_res' function
#' @param celltype_group cell type groups (two columns, 1st are cell types and 2nd are groups in number format start from 0) (the order of this file is the same as shown in tree from top to bottom)
#' @param order_date included date (tree root has beed assigned as 'Root' in 'Subway_map.step_3.knn_res' function), default (see below)
#' @param width_para control tree width (14 for 22 cell types, 3 for 9 cell types), default '14'
#' @param height_para control tree height (23 for 22 cell types, 18 for 9 cell types), default '23'
#' @param used_cols colors for tree nodes, default (see below)
#' @param pdf_height pdf height, default '744'
#' @param pdf_width pdf width, default '1800'
#' @param outname prefix of output name, default 'test'
Subway_map.step_4.plot_tree <- function(
edge,
edge_prob,
celltype_group = NULL,
order_date = c("E14","E18","E43","E55","E93","7d","14d","60d","180d"),
dotsize = 5,
width_para = 14,
height_para = 23,
used_cols = c("#000000","#1F78B4","#E31A1C","#33A02C","#FF7F00","#A6CEE3","#FB9A99","#B2DF8A","#FDBF6F","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#00E5EE","#76EEC6","#FFD700","#FF69B4","#00008B","#008B8B","#8B0000","black","#00B2EE","#FF6347","#006400","#FFA500","#551A8B"),
pdf_height = 744,
pdf_width = 1500,
outname = "test"
){
	library(dplyr)
	library(pagedown)
	
	tree <- read.table(edge, stringsAsFactor = F, sep = "\t")
	# get cell type group
	if(is.null(celltype_group)){
		warning("file of 'celltype_group' need to be modified to show colors")
		all_celltypes <- c(tree[,1], tree[,2]) %>% as.character %>% strsplit(.,":") %>% sapply(.,function(x) x[[2]]) %>% unique %>% as.data.frame
		all_celltypes$group = 0
		write.table(all_celltypes, paste0(outname, ".celltype_group.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
	}else{
		# cell group checking
		all_celltypes <- c(tree[,1], tree[,2]) %>% as.character %>% strsplit(.,":") %>% sapply(.,function(x) x[[2]]) %>% unique %>% as.data.frame
		colnames(all_celltypes)[1] <- "V1"
		celltype_order <- unique(all_celltypes[,1])
		cellgroup_file <- read.table(celltype_group, stringsAsFactor = F, sep = "\t")
		if(length(setdiff(all_celltypes[,1], cellgroup_file[,1])) > 0){
			warning(paste("Cell types in tree are not present in input group file:", paste(setdiff(all_celltypes[,1], cellgroup_file[,1]), collapse = ", ")))
			all_celltypes$group = 0
		}else{
			all_celltypes <- merge(all_celltypes, cellgroup_file, by = "V1", all.x = T, all.y = F)
			all_celltypes[,2] <- mapvalues(all_celltypes[,2], unique(all_celltypes[,2]), 0:(length(unique(all_celltypes[,2]))-1))
			rownames(all_celltypes) <- all_celltypes[,1]
			all_celltypes <- all_celltypes[celltype_order, , drop = F]
		}
		write.table(all_celltypes, paste0(outname, ".celltype_group.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
	}
	celltype_group <- paste0(outname, ".celltype_group.txt")
	## renew py script 'create_map.py'
	py.time_point <- c(tree[,1], tree[,2]) %>% as.character %>% strsplit(.,":") %>% sapply(.,function(x) x[[1]]) %>% unique %>% intersect(c("Root", order_date),.) %>% paste(., collapse = '","') %>% paste0('time_point = ["', ., '"]')
	py.edge <- paste0('file = open("', edge, '")')
	py.edge_prob <- paste0('file = open("', edge_prob, '")')
	py.celltype_group <- paste0('file = open("', celltype_group, '")')
	py.colors <- paste0('"', used_cols, '"') %>% paste(0:(length(used_cols)-1), ., sep = ":") %>% paste(., collapse = ",") %>% paste0("color_map = {", ., "}")
	py.root_edge <- 'Root:Dummy'
	py.out_json <- paste0(outname, ".json") %>% paste0('with open("', ., '", "w") as json_file:')
	
	create_map_name <- paste0(outname, ".create_map.py")
	system(paste0("cp ~/01.sxzhang/01.pig.snatac/03.archr_process/02.2nd/11.Fig3/03.organogenesis.cellTraj/test/11.TOME_roadmap/99.test_code_seurat/TOME_script/tome_code/help_code/create_map.py ", create_map_name))
	system(paste0("sed -i '4c ", py.time_point, "' ", create_map_name))
	system(paste0("sed -i '14c ", py.edge, "' ", create_map_name))
	system(paste0("sed -i '35c ", py.edge_prob, "' ", create_map_name))
	system(paste0("sed -i '53c ", py.celltype_group, "' ", create_map_name))
	system(paste0("sed -i '86c ", py.colors, "' ", create_map_name))
	system(paste0("sed -i 's/^color_map/    color_map/' ", create_map_name))
	system(paste0("sed -i '113c ", py.out_json, "' ", create_map_name))
	system(paste0("sed -i 's/E3:Morula/", py.root_edge, "/g' ", create_map_name))
	## run 'create_map.py'
	system(paste("python", create_map_name))
	## generate 'map.html'
	map_file_name <- paste0(outname, ".map.sxz.html")
	system(paste("cp ~/01.sxzhang/01.pig.snatac/03.archr_process/02.2nd/11.Fig3/03.organogenesis.cellTraj/test/11.TOME_roadmap/99.test_code_seurat/TOME_script/map.sxz.html", map_file_name))
	system(paste0("sed -i 's/stroke-width: 5px/stroke-width: ", dotsize,"px/' ", map_file_name))
	system(paste0("sed -i '30r ", outname,".json' ", map_file_name))
	system(paste0("sed -i '30d' ", map_file_name))
	system(paste0("sed -i 's/14); });/", width_para, "); });/' ", map_file_name))  # set width, 14 for 22 cell types, 3 for 9 cell types
	system(paste0("sed -i 's/23); });/", height_para, "); });/' ", map_file_name))  # set height, 23 for 22 cell types, 18 for 9 cell types
}

#### create heatmap of each pair of adjacent stages
Subway_map.step_4.plot_heatmap <- function(
edge_all,
#celltype_group,
order_date = c("E14","E18","E43","E55","E93","7d","14d","60d","180d"),
margins = c(8, 8),
orders = F,
Key = T,
pdf_height = NULL,
pdf_width = NULL,
outname = "test"
){
	library(gplots)
	library(viridis)
	
	dat = readRDS(edge_all)
	time_point = c(dat$pre_time, dat$nex_time) %>% unique %>% intersect(order_date, .)
	
	pdf_vector <- c()
	for(kk in 1:(length(time_point)-1)){
	  print(paste0(kk,"/",(length(time_point)-1)))
	  dat_sub = dat[dat$pre_time == time_point[kk],c("pre_cell","nex_cell","prob")]
	  if(length(unique(dat_sub$pre_cell)) <= 1 | length(unique(dat_sub$nex_cell)) <= 1){
	    print(paste("Skipped: Cells of ", time_point[kk], ":", length(unique(dat_sub$pre_cell)), ", Cells of ", time_point[kk+1], ":", length(unique(dat_sub$nex_cell))))
	    next
	  }
	  df = dcast(dat_sub, nex_cell~pre_cell)
	  rownames(df) <- df[,1]; df <- df[,-1]
	  if(orders){
	    df <- order_heatmap(df)
	  }
	  # set pdf size
	  if(is.null(pdf_height) & is.null(pdf_width)){
	    pdf(paste0(outname, ".edge_heatmap.tmp", kk,".pdf"), height = 3+0.2*nrow(df), width = 3+0.2*ncol(df))
	  }else{
	    pdf(paste0(outname, ".edge_heatmap.tmp", kk,".pdf"), height = pdf_height, width = pdf_width)
	  }
	  par(mar = c(0, 0, 0, 0))
	  pdf_vector <- paste(pdf_vector, paste0(outname, ".edge_heatmap.tmp", kk,".pdf"))
	  if(Key){
	    heatmap.2(as.matrix(df), col=viridis, scale="none", Rowv = FALSE, Colv = FALSE, key = Key, keysize = 1, density.info="none", trace="none", cexRow=1, cexCol = 1, 
	      margins = margins, xlab = time_point[kk], ylab = time_point[kk+1], breaks=seq(0,1,0.01), aspect = 1)
	  }else{
	    heatmap.2(as.matrix(df), col=viridis, scale="none", Rowv = FALSE, Colv = FALSE, key = Key, keysize = 1, density.info="none", trace="none", cexRow=1, cexCol = 1, 
	      margins = margins, xlab = time_point[kk], ylab = time_point[kk+1], breaks=seq(0,1,0.01), lwid = c(0.01,1), lhei = c(0.01,1), aspect = 1)
	  }
	  dev.off()
	}
	
	system(paste("pdfunite", pdf_vector, paste0(outname, ".edge_heatmap.pdf")))
	system(paste0("rm ", pdf_vector))
}
# order heatmap to generate diagonal pattern
order_heatmap <- function(mat){
	# order/cluster rows
	dist_matrix <- dist(mat, method = "euclidean")
	hc <- hclust(dist_matrix, method = "complete")
	mat <- mat[hc$order, ]
	# order columns
	col_order <- c()
	for(i in 1:nrow(mat)){
		maxindex <- which(mat[i,] == max(mat[i,]))
		col_order <- c(col_order, maxindex)
	}
	mat <- mat[, c(unique(col_order), setdiff(1:ncol(mat), unique(col_order)))]
	# return res
	return(mat)
}
