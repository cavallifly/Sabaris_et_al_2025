##########################################################################################################
.shaman_split_and_merge_matrix <- function(mat, grid_min, grid_max, grid_res, new_res) {
    if (!(new_res / grid_res == floor(new_res / grid_res))) {
        stop(paste0("plot resolution", new_res, " is not a multiple of grid resolution ", grid_res))
    }
    expand <- seq(grid_min, grid_max, by = grid_res)
    multi <- new_res / grid_res
    m <- reshape2::melt(mat)
    m$i <- as.numeric(m[, 1])
    m$j <- as.numeric(m[, 2])
    if (identical(which(expand == 0), integer(0))) {
        # center spans the 0 need to take only a single bin
        if (multi / 2 == floor(multi / 2)) {
            stop(paste0("cannot maintain center at 0 - ", multi, " is even, center spans 0"))
        }
        center_bin <- length(expand) / 2
        new_pos_bin <- c(rep(0, (multi - 1) / 2), rep(1:(max(m$i) / 2 / multi), each = multi))[1:(center_bin - 1)]
        new_bin <- c(rev(-new_pos_bin), 0, new_pos_bin)
    } else {
        if (!(multi / 2 == floor(multi / 2))) {
            stop(paste0("cannot maintain center at 0 - ", multi, " is odd, center does not span 0"))
            return(integer(0))
        }
        center_bin_pos <- which(expand == 0)
        new_pos_bin <- c(rep(0, multi / 2), rep(1:(max(m$i) / 2 / multi), each = multi))[1:(center_bin_pos - 1)]
        new_bin <- c(rev(-new_pos_bin), new_pos_bin)
    }
    m$n_i <- new_bin[m$i]
    m$n_j <- new_bin[m$j]
    d <- plyr::ddply(m, c("n_i", "n_j"), plyr::summarise, sum = sum(value, na.rm = TRUE))
    return(reshape2::dcast(d, n_i ~ n_j)[, -1])
}
#pal = colorRampPalette(c("#00008b", "#0000ff", "white", "#ff0000", "#8b000")), 
##########################################################################################################
my_shaman_plot_feature_grid <- function(grids, range, grid_resolution, plot_resolution, type = "enrichment",
                                     fig_fn = "", fig_width = 400, fig_height = 400,
				     pal = colorRampPalette(c("blue", "#87FFFF", "white", "#FF413D", "black")),
                                     zlim = c(-1, 1),
				     text = NA, N=NA)
{
    obs <- NULL
    exp <- NULL
    for (g in grids)
    {
        if (plot_resolution != grid_resolution)
	{
            o <- .shaman_split_and_merge_matrix(g$obs, -range, range, grid_resolution, plot_resolution)
            e <- .shaman_split_and_merge_matrix(g$exp, -range, range, grid_resolution, plot_resolution)
        } else {
            o <- g$obs
            e <- g$exp
        }
	print("Grid Obs")
	print(o)
	print("Grid Exp")	
	print(e)	

        if (is.null(obs))
	{
            obs <- o
            exp <- e
        } else {
            obs <- obs + o
            exp <- exp + e
        }
    }
    if (fig_fn != "")
    {
        png(fig_fn, width = fig_width, height = fig_height)
    }
    par(mar = c(0, 0, 0, 0))
    if (type == "enrichment")
    {
	tmpObs <- obs / sum(obs, na.rm = TRUE)
	tmpExp <- exp / sum(exp, na.rm = TRUE)
	print("Obs / sum")
	print(tmpObs)
	print("Exp / sum")	
	print(tmpExp)	

        data <- as.matrix(log2((obs / sum(obs, na.rm = TRUE)) / (exp / sum(exp, na.rm = TRUE))))
	min  <- abs(min(data))
	max  <- max(data)
	max  <- max(min,max)*(1.0001)
	#	zlim <- c(-max,max)
	image(data, zlim = zlim, col = pal(1000), xaxt = "n", yaxt = "n")
    } else {
        if (type == "observed")
	{
            image(as.matrix(log10(obs / sum(obs))), zlim = zlim, col = pal(1000), xaxt = "n", yaxt = "n")
        } else {
            if (type == "expected")
	    {
                image(as.matrix(log10(exp / sum(exp))), zlim = zlim, col = pal(1000), xaxt = "n", yaxt = "n")
            } else {
                stop(paste("unknown plot type:", type))
            }
        }
    }
    abline(v = 0.5, col = "grey20")
    abline(h = 0.5, col = "grey20")
    #####
    if(!is.na(text)) {
        text(x=0.2,y=0.9,labels=sprintf("%d l %.2f:%.2f",N,-max,max), cex=2, font=2)
    }
    #####
    if (fig_fn != "") {
        dev.off()
    }
    return(list(obs = obs, exp = exp))
}

############################################
shaman_generate_feature_grid_2d <- function(iter2d, obs_track_nm, exp_track_nm=paste0(obs_track_nm, "_shuffle"), range=25000, resolution=500,min_dist=2e05, max_dist=1e06)
{
    #check tracks
    if (!gtrack.exists(obs_track_nm))
    {
  	stop(paste("Missing obs_track_nm (", obs_track_nm, ") in track db"))
    }
    if (!gtrack.exists(exp_track_nm))
    {
 	stop(paste("Missing exp_track_nm (", exp_track_nm, ") in track db"))
    }
    options(gmultitasking=FALSE)
    
    expand = seq(0-range, range, by=resolution)
    grid <- shaman_norm_contact_grid_2d(obs_track_nm, exp_track_nm, iter2d, expand, expand, regularization=0,min_dist=min_dist, max_dist=max_dist)
    return(grid)
}

############################################
shaman_norm_contact_grid_2d <- function(track, shuffled_track, iter2d, expand1, expand2, regularization=30, min_dist=2e05, max_dist=1e06)
{
    o = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
    e = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
	
    grid <- iter2d
    band = c(-max_dist, -min_dist)
	
    message(paste("found ", nrow(grid), "regions"))
    interv_set_out = gsub("+", "", paste(track, paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
    message(interv_set_out)
    if (!gintervals.exists(interv_set_out))
    {
 	giterator.intervals(track, band=band, intervals.set.out=interv_set_out)
    }
    print("Observed track(s)")
    print("A data frame of the track under analysis")
    print(interv_set_out)

    message("obs neighbors...")
    a = gintervals.neighbors(interv_set_out, grid,
  	    mindist1=min(expand1), maxdist1=max(expand1),
	    mindist2=min(expand2), maxdist2=max(expand2))
    colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))	    
    print("Paired interactions of the Observed track in the matrices of peak-pairs between minDist and maxDist")
    print(nrow(a))
    print(a)	    

    o = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE),
	    cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))
    print("Observed matrix")		
    print(o)	    

    interv_set_out = gsub("+", "", paste(shuffled_track, paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
    if (!gintervals.exists(interv_set_out))
    {
 	giterator.intervals(shuffled_track, band=band, intervals.set.out=interv_set_out)
    }
    print("Expected track(s)")
    print("A data frame of the track under analysis")    
    print(interv_set_out)

    message("exp neighbors...")
    a = gintervals.neighbors(interv_set_out, grid, mindist1=min(expand1), maxdist1=max(expand1),
   	mindist2=min(expand2), maxdist2=max(expand2))
    colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
    print("Paired interactions of the Expected (shuffled) track in the matrices of peak-pairs between minDist and maxDist")
    print(nrow(a))
    print(a)	    

    e = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE),
 	      cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))

    print("Expected matrix")
    print(e)	      
	
    o[o < regularization] = NA
    e[is.na(o) | e < regularization] = NA

    return(list(obs=o, exp=e, grid_size=nrow(grid)));
}

############################################
shaman_generate_feature_grid <- function(feature1, feature2, obs_track_nm, exp_track_nm=paste0(obs_track_nm, "_shuffle"),
  score_track_nm=paste0(obs_track_nm, "_score"), score_filter=30, range=25000, resolution=500,
  min_dist=2e05, max_dist=1e06)
{
  #check tracks
  if (!gtrack.exists(obs_track_nm)) {
        stop(paste("Missing obs_track_nm (", obs_track_nm, ") in track db"))
    }
   if (!gtrack.exists(exp_track_nm)) {
        stop(paste("Missing exp_track_nm (", exp_track_nm, ") in track db"))
   }
   options(gmultitasking=FALSE)
   expand = seq(0-range, range, by=resolution)
   filter_vtrack = NA
   if (gtrack.exists(score_track_nm)) {
        gvtrack.create("v_fg_max_score", score_track_nm, "max")
        filter_vtrack = "v_fg_max_score"
   }

   grid <- .shaman_norm_contact_grid_by_shuffled_filter(obs_track_nm, exp_track_nm, filter_vtrack, score_filter,
	unique(feature1[,c("chrom", "start", "end")]), unique(feature2[,c("chrom", "start", "end")]),
	expand, expand, min_dist=min_dist, max_dist=max_dist, regularization=0)
   if (!is.na(filter_vtrack)) {
	gvtrack.rm(filter_vtrack)
   }
   return(grid)
}

############################################
.shaman_norm_contact_grid_by_shuffled_filter <- function(track, shuffled_track, filter_vtrack, filter_value, interv1, interv2, expand1, expand2, min_dist, max_dist, regularization=30)
{
  o = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
  e = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
  message("building grid...")
  grid = plyr::ddply(interv1, c("chrom"), function(i1) {
        i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
        if (nrow(i2) ==0) {
                return(c())
        }
        g = expand.grid(1:nrow(i1), 1:nrow(i2))
        g = g[i2$start[g$Var2]-i1$start[g$Var1] > min_dist &
                i2$start[g$Var2]-i1$start[g$Var1] < max_dist,]
          grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
          colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
          grid$start1 = grid$start1 + floor((grid$end1-grid$start1)/2)
          grid$start2 = grid$start2 + floor((grid$end2-grid$start2)/2)
          grid$end1 = grid$start1+1
          grid$end2 = grid$start2+1

         return(grid)
   })[,-1]

  band = c(-max_dist, -min_dist)
  message("screening for filter")
  if (!is.na(filter_vtrack)) {
        gvtrack.iterator.2d(filter_vtrack, sshift1=min(expand1), eshift1=max(expand1), sshift2=min(expand2), eshift2=max(expand2))
        grid = gscreen(sprintf("%s > %s", filter_vtrack, filter_value), intervals=grid, iterator=grid, band=band)
  }
  message(paste("found ", nrow(grid), "regions"))
  interv_set_out = gsub("+", "", paste(track, paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
  if (!gintervals.exists(interv_set_out)) {
        giterator.intervals(track, band=band, intervals.set.out=interv_set_out)
  }
  message("obs neighbors...")
  a = gintervals.neighbors(interv_set_out, grid,
        mindist1=min(expand1), maxdist1=max(expand1),
        mindist2=min(expand2), maxdist2=max(expand2))
  colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
  o = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE),
        cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))

  interv_set_out = gsub("+", "", paste(shuffled_track, paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
  if (!gintervals.exists(interv_set_out)) {
        giterator.intervals(shuffled_track, band=band, intervals.set.out=interv_set_out)
  }
  message("exp neighbors...")
  a = gintervals.neighbors(interv_set_out, grid, mindist1=min(expand1), maxdist1=max(expand1),
        mindist2=min(expand2), maxdist2=max(expand2))
  colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
  e = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE),
        cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))

  o[o < regularization] = NA
  e[is.na(o) | e < regularization] = NA
  return(list(obs=o, exp=e, grid_size=nrow(grid)));
}