###################
rm(list=ls())

library(Momocs)
library(ggplot2)
library(outlineR)


# create closed outlines from prepared images
cc <- Momocs::import_jpg(list.files(file.path("2_data", "barbs"), 
                                    full.names = T))
outline_names <- names(cc)

# import CSV
point_info <- readr::read_csv(file.path("2_data", "Tsirintoulaki_OPAR_barbed_points_edit.csv"))


### we match the general ID (i.e. "ABL_1") to the ID of the individual barbs (i.e. "ABL_1_pseudo_no_1)
general_ID_to_barb_ID <- list()
for (i in 1:length(outline_names)){
  current_split <- strsplit(outline_names[i], 
                            split = "_")[[1]]
  
  ID <- paste0(current_split[1], "_", current_split[2])  
  
  general_ID_to_barb_ID[[i]] <- data.frame(ID = ID,
                                           ID_barb = outline_names[i])
}
general_ID_to_barb_ID_df <- do.call("rbind", general_ID_to_barb_ID)

### and combine it with the initial CSV
ID_df <- dplyr::left_join(general_ID_to_barb_ID_df, point_info, 
                          by = "ID")
ID_df$ID <- as.factor(ID_df$ID)


## transform outlines into a Momocs "Out" object and appended supplementary CSV information
out <- Momocs::Out(cc,
                   fac = ID_df)
### check
Momocs::panel(out,
              fac = "ID")



# create "Opn" outlines from closed "Out"-lines

## call the custom function
source(file.path("1_script", "open_outlines_from_closed_outlines_v3.R"))

### run it
open_test <- open_outlines_from_closed_outlines(out)

### name coo-list
names(open_test$coo) <- outline_names
### append CSV again
open_test <- Opn(open_test$coo,
                 fac = ID_df)


### center and scale
out_centered <- Momocs::coo_center(open_test)
open_test <- Momocs::coo_scale(out_centered)



### check result
Momocs::pile(open_test)
Momocs::panel(open_test)



output_dir <- file.path("3_output", "barbs")
dir.create(output_dir,
           recursive = T)
# save in Momocs' Opn Coo-format
saveRDS(open_test,
        file = file.path(output_dir, "barbs_open_outlines_OpnCoo.RDS"))
# save as list of coordinates
saveRDS(open_test$coo,
        file = file.path(output_dir, "barbs_open_outlines_listOfCoords.RDS"))










