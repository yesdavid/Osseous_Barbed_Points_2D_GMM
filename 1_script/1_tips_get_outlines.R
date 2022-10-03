###################
rm(list=ls())


library(Momocs)
library(ggplot2)

# create closed outlines from prepared images
cc <- Momocs::import_jpg(list.files(file.path("2_data", "tips"),
                                    full.names = T))
outline_names <- names(cc)


# import CSV
point_info <- readr::read_csv(file.path("2_data", "Tsirintoulaki_OPAR_barbed_points_edit.csv"))
ID_df <- subset(point_info, 
                ID %in% names(cc))
ID_df$ID <- as.factor(ID_df$ID)


## transform outlines into a Momocs "Out" object and appended supplementary CSV information
out <- Momocs::Out(cc,
                   fac = ID_df)
### check
Momocs::panel(out,
              fac = "ID")


# create "Opn" outlines from closed "Out"-lines

## call the custom function
source(file.path("1_script", "open_outlines_from_closed_outlines_v4.R"))

### run it
open_test <- open_outlines_from_closed_outlines(out)

### append CSV again
open_test <- Opn(open_test$coo,
                 fac = ID_df)

### center and scale
out_centered <- Momocs::coo_center(open_test)
open_test <- Momocs::coo_scale(out_centered)

### check result
Momocs::pile(open_test)
Momocs::panel(open_test)


output_dir <- file.path("3_output", "tips")
dir.create(output_dir,
           recursive = T)
# save in Momocs' Opn Coo-format
saveRDS(open_test,
        file = file.path(output_dir, "tips_open_outlines_OpnCoo.RDS"))
# save as list of coordinates
saveRDS(open_test$coo,
        file = file.path(output_dir, "tips_open_outlines_listOfCoords.RDS"))




