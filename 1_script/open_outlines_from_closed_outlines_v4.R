#' open_outlines_from_closed_outlines
#'
#' A function to create open outlines from closed outlines. The way it is implemented, the open outline will start at the highest, left-most coordinate, run clock-wise, and will end at the lowest, left-most coordinate.
#'
#' @param outlines_combined A Momocs Out file containing closed outlines.
#' @return Returns the combined Coo objects in a single Opn file
#'
#' @export



open_outlines_from_closed_outlines <- function(outlines_combined) {
  open_outlines_list <- list()


  dist_to_start_and_end <- function(x,y, max_y, min_y){

    DistToStart <- ((0 - x)^2 + (max_y - y)^2)^0.5
    current_outline_df[current_row, "DistToStart"] <- DistToStart

    DistToEnd <- ((0 - x)^2 + (0 - y)^2)^0.5
    current_outline_df[current_row, "DistToEnd"] <- DistToEnd

    return(current_outline_df[current_row, ])

  }

  pb <- utils::txtProgressBar(min = 0, max = length(outlines_combined), style = 3)
  for (counter in 1:length(outlines_combined)) {
    current_outline_name <- names(outlines_combined$coo[counter])

    current_outline <- outlines_combined$coo[[counter]]

    current_outline_df <- as.data.frame(current_outline)
    names(current_outline_df) <- c("X", "Y")

    current_outline_df$DistToStart <- NA
    current_outline_df$DistToEnd <- NA

    test_list <- list()
    max_y <- max(current_outline_df$Y)#+100
    # min_y <- max(current_outline_df$Y)*0.3
    for (current_row in 1:nrow(current_outline_df)) {
      so <- current_outline_df[current_row,]
      test_list[[current_row]] <- dist_to_start_and_end(x = so$X, y = so$Y, max_y = max_y)
    }
    result_df <- do.call("rbind", test_list)

    # plot(result_df[,c("X","Y")], cex = 0.5)
    # points(result_df[which(result_df$DistToStart == min(result_df$DistToStart)),c("X","Y")], col = "blue", cex = 5, add = T)
    # points(result_df[which(result_df$DistToEnd == min(result_df$DistToEnd)),c("X","Y")], col = "red", cex = 5, add = T)

    min_x_current <- result_df[which(result_df$DistToStart == min(result_df$DistToStart)),c("X","Y")]$X
    subset_min_x_y <- subset(result_df, X == min_x_current)

    # points(subset_min_x_y[which(subset_min_x_y$Y == min(subset_min_x_y$Y)),c("X","Y")], col = "green", cex = 5, add = T)


    start <- as.integer(rownames(result_df[which(result_df$DistToStart == min(result_df$DistToStart)),c("X","Y")]))
    end <- as.integer(rownames(subset_min_x_y[which(subset_min_x_y$Y == min(subset_min_x_y$Y)),c("X","Y")]))

    open_outlines_list[[current_outline_name]] <- Momocs::coo_slidegap(current_outline[c(1:end, start:nrow(current_outline)),], 
                                                                       force = T)
    
    # closed <- Momocs::Opn(current_outline[c(1:end, start:nrow(current_outline)),])
    # 
    # 
    # open_outlines_list[[current_outline_name]] <- Momocs::coo_slidegap(x, force = T)


    setTxtProgressBar(pb, counter)
  }
  close(pb)

  return(Momocs::Opn(open_outlines_list))
  
}
