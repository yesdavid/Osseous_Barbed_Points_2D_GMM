# overall barb variability/disparity

library(ggplot2)
library(readr)
library(fpc)
library(Momocs)
library(dplyr)
library(dispRity)


rm(list=ls())

point_info_raw <- readr::read_csv(file.path("2_data", "Tsirintoulaki_OPAR_barbed_points_edit.csv"),
                                  col_types = cols(Period = col_factor(levels = c("Late Pleistocene", 
                                                                                  "Early Holocene"))))

point_info_raw$Age_calBP_mean <- 
sapply(point_info_raw$Age_calBP, FUN = function(x){
  ages <- strsplit(x, split = "-")[[1]]
  mean(c(as.numeric(ages[1]), as.numeric(ages[2])))
})




point_info_raw$Latitude_zone <- factor(unlist(lapply(point_info_raw$Latitude, function(x){if(x<49){"<49°"}else if(x<53){"<53°"}else if(x<57){"<57°"}})),
                                       levels = c("<57°", "<53°", "<49°"))

# create timebins 
point_info_raw$Timebins <- 
factor(sapply(point_info_raw$Age_calBP_mean, FUN = function(x){
  if(x < 10700){
    "Boreal\n~10.7-9.3 ka calBP"
  } else if(x < 11700){
    "Preboreal\n~11.7-10.7 ka calBP"
  } else if(x < 12900){
    "Younger Dryas Complex\n~12.9-11.7 ka calBP"
  } else if(x < 14600){
    "Bølling/Allerød Complex\n~14.6-12.9 ka calBP"
  } else if(x < 18000){
    "Late Pleniglacial\n~18-14.6 ka calBP"
  } 
}
),
levels = c("Late Pleniglacial\n~18-14.6 ka calBP", 
           "Bølling/Allerød Complex\n~14.6-12.9 ka calBP", 
           "Younger Dryas Complex\n~12.9-11.7 ka calBP", 
           "Preboreal\n~11.7-10.7 ka calBP",
           "Boreal\n~10.7-9.3 ka calBP"))



C14_dates_df_list <- list()
for (i in 1:nrow(point_info_raw)){
  C14_dates_df <- point_info_raw[i,c("ID", "Age_calBP")]
  date_from_to <- strsplit(C14_dates_df$Age_calBP, "-")[[1]]
  C14_dates_df$mean_age_calBP <- round((as.integer(date_from_to[1]) + as.integer(date_from_to[2]))/2, 
                                       digits = 0)
  C14_dates_df_list[[i]] <- C14_dates_df
}
C14_dates_df <- do.call(rbind.data.frame, C14_dates_df_list)

point_info_raw <- dplyr::left_join(point_info_raw, C14_dates_df, by = "ID")

output_dir <- file.path("3_output", "barbs")


open_outlines_raw <- readRDS(file.path(output_dir, "barbs_open_outlines_OpnCoo.RDS"))
# point_info <- open_outlines$fac

barb_id_df_list <- list()
for (i in names(open_outlines_raw$coo)) {
  barb_id_df <- data.frame(ID = strsplit(i, "_pseudo_no_")[[1]][1],
                           ID_barb = i)
  barb_id_df_list[[i]] <- barb_id_df
}
barb_id_df <- do.call(rbind.data.frame, barb_id_df_list)

point_info_raw_barb_IDs <- dplyr::left_join(barb_id_df, point_info_raw, by = "ID")

point_info <- subset(point_info_raw_barb_IDs, ID_barb %in% names(open_outlines_raw$coo))

open_outlines <- Momocs::Opn(x = open_outlines_raw$coo,
                             fac = point_info)
## dfourier
open_outlines_dfourier <- Momocs::dfourier(open_outlines)           # discrete cosine transform.

## PCA
open_outlines_PCA <- Momocs::PCA(open_outlines_dfourier,
                                 fac = point_info)             # we calculate a PCA on it
# rownames(open_outlines_PCA$x) <- open_outlines_PCA$fac$ID_barb

################################################################
IDs_w_more_than_X_barbs <- Momocs::filter(open_outlines, 
                                          ID %in% names(which(table(open_outlines$ID) >= 4))  # just retain IDs with more than X samples
                                          # & !(ID %in% "BTZ_3")
                                          ) 
## dfourier
open_outlines_dfourier_IDs_w_more_than_X_barbs <- Momocs::dfourier(IDs_w_more_than_X_barbs)           # discrete cosine transform.

## PCA
open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA <- Momocs::PCA(open_outlines_dfourier_IDs_w_more_than_X_barbs,
                                 fac = point_info)             # we calculate a PCA on it
rownames(open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA$x) <- open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA$fac$ID_barb




######################################################
# disparity of each Period/TimeSlice
rownames_DATASETS <- list()
for(i in levels(open_outlines_PCA$Period)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Period == i)$ID_barb)
}
rownames_DATASETS

TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 1000)
TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

# Wilcox.test
test.dispRity(TS_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "BH")
# PERMANOVA
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "BH")


TS_names <- names(TS_disp$disparity)
disparity_df_list <- list()
for(i in TS_names){
  disparity_df_list[[i]] <- data.frame(Period = paste0(i, 
                                                       "\n(n=",nrow(TS_disp$subsets[[i]]$elements),")"),
                                       disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                       nelements = nrow(TS_disp$subsets[[i]]$elements),
                                       TS = i)
}
disparity_df_TSdiscrete_armatureOutlines_perTShard <- do.call(rbind.data.frame, disparity_df_list)

disparity_df_TSdiscrete_armatureOutlines_perTShard$Period <- factor(disparity_df_TSdiscrete_armatureOutlines_perTShard$Period)
disparity_df_TSdiscrete_armatureOutlines_perTShard$Period <- relevel(disparity_df_TSdiscrete_armatureOutlines_perTShard$Period, "Late Pleistocene\n(n=102)")

disparity_TSdiscrete_armatureOutlines_ggplot_perTShard <- 
  ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTShard,
         aes(x = Period, 
             y = disparity
         )) +
  # geom_violin() + 
  geom_boxplot(notch = F, 
               width = 0.1, 
               aes(fill = Period)) +
  theme_bw() +
  ggtitle(NULL) +
  xlab("") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        axis.text.x = element_text(#angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0.5),
        axis.title=element_text(size=14)) +
  guides(color = FALSE, 
         fill = FALSE) #+ 
  # coord_flip()

disparity_TSdiscrete_armatureOutlines_ggplot_perTShard

ggplot_perTShard_height <- 6
ggplot_perTShard_width <- 6

ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
       filename = file.path(output_dir, "Figure_6_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
       filename = file.path(output_dir, "Figure_6_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.png"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")





############################################
# disparity by latitudinal zone
rownames_DATASETS <- list()
for(i in unique(open_outlines_PCA$Latitude_zone)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Latitude_zone == i)$ID_barb)
}
rownames_DATASETS

unique_artefacts_in_lat_zones <- 
  sapply(rownames_DATASETS, 
         function(x){
           unique(unlist(lapply(x, 
                                function(x){
                                  paste0(strsplit(x, "_")[[1]][1], "_", strsplit(x, "_")[[1]][2])
                                }
           )))
         }
  )

TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 1000)
TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

# Wilcox.test
test.dispRity(TS_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "BH")
# PERMANOVA
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "BH")


TS_names <- names(TS_disp$disparity)
disparity_df_list <- list()
for(i in TS_names){
  disparity_df_list[[i]] <- data.frame(Latitude_zone = paste0(i, 
                                                              "\n(nBarbs=",nrow(TS_disp$subsets[[i]]$elements),
                                                              ",\nnArtefacts=", length(unique_artefacts_in_lat_zones[[i]]),")"),
                                       disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                       nelements = nrow(TS_disp$subsets[[i]]$elements),
                                       TS = i)
}
disparity_df_TSdiscrete_armatureOutlines_perLatitude_zone <- do.call(rbind.data.frame, disparity_df_list)


disparity_df_TSdiscrete_armatureOutlines_perLatitude_zone$Latitude_zone <- factor(disparity_df_TSdiscrete_armatureOutlines_perLatitude_zone$Latitude_zone,
                                                                       levels = (unique(disparity_df_TSdiscrete_armatureOutlines_perLatitude_zone$Latitude_zone)))


disparity_TSdiscrete_armatureOutlines_ggplot_perLatitude_zone <- 
  ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perLatitude_zone,
         aes(x = Latitude_zone, 
             y = disparity)) +
  # geom_violin(notch = F, 
  #             # width = 0.1, 
  #             aes(fill = Latitude_zone)) +
  # geom_boxplot(notch = F,
  #              width = 0.1,
  #              color = "black",
  #              fill = "white") +
  geom_boxplot(notch = F,
               aes(fill = Latitude_zone)) +
  theme_bw() +
  ggtitle(NULL) +
  xlab("Latitudinal Zone") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title=element_text(size=14)) +
  guides(color = FALSE, fill = FALSE) + 
  scale_fill_manual(values = terrain.colors(4))

disparity_TSdiscrete_armatureOutlines_ggplot_perLatitude_zone

ggplot_perTShard_height <- 6
ggplot_perTShard_width <- 6

ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perLatitude_zone,
       filename = file.path(output_dir, "Figure_8_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perLatitude_zone.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perLatitude_zone,
       filename = file.path(output_dir, "Figure_8_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perLatitude_zone.png"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")

######################################################################
######################################################################

# disparity by timebins
rownames_DATASETS <- list()
for(i in levels(open_outlines_PCA$fac$Timebins)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Timebins == i)$ID_barb)
}
rownames_DATASETS
rownames_DATASETS$'Late Pleniglacial\n~18-14.6 ka calBP' <- NULL

unique_artefacts_in_lat_zones <-
  sapply(rownames_DATASETS,
         function(x){
           unique(unlist(lapply(x,
                                function(x){
                                  paste0(strsplit(x, "_")[[1]][1], "_", strsplit(x, "_")[[1]][2])
                                }
           )))
         }
  )

TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 1000)
TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

# t.test
test.dispRity(TS_disp, 
              test = t.test, 
              comparisons = "sequential",
              correction = "BH")
# PERMANOVA
adonis_timebins <- 
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "BH")
adonis_timebins
summary(adonis_timebins)


TS_names <- names(TS_disp$disparity)
disparity_df_list <- list()
for(i in TS_names){
  disparity_df_list[[i]] <- data.frame(Timebins = paste0(i, 
                                                              "\n(nBarbs=",nrow(TS_disp$subsets[[i]]$elements),
                                                              ",\nnArtefacts=", length(unique_artefacts_in_lat_zones[[i]]),")"),
                                       disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                       nelements = nrow(TS_disp$subsets[[i]]$elements),
                                       TS = i)
}
disparity_df_TSdiscrete_armatureOutlines_perTimebins <- do.call(rbind.data.frame, disparity_df_list)


disparity_df_TSdiscrete_armatureOutlines_perTimebins$Timebins <- factor(disparity_df_TSdiscrete_armatureOutlines_perTimebins$Timebins,
                                                                                  levels = c("Late Pleniglacial\n~18-14.6 ka calBP\n(nBarbs=8,\nnArtefacts=5)",
                                                                                             "Bølling/Allerød Complex\n~14.6-12.9 ka calBP\n(nBarbs=50,\nnArtefacts=8)",
                                                                                             "Younger Dryas Complex\n~12.9-11.7 ka calBP\n(nBarbs=44,\nnArtefacts=5)",
                                                                                             "Preboreal\n~11.7-10.7 ka calBP\n(nBarbs=185,\nnArtefacts=15)",
                                                                                             "Boreal\n~10.7-9.3 ka calBP\n(nBarbs=158,\nnArtefacts=16)"))

disparity_df_TSdiscrete_armatureOutlines_perTimebins_subs <- subset(disparity_df_TSdiscrete_armatureOutlines_perTimebins, !(TS %in% "Late Pleniglacial\n~18-14.6 ka calBP"))

disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins <- 
  ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTimebins_subs,
         aes(x = (Timebins), 
             y = disparity)) +
  geom_boxplot(notch = F,
               aes(fill = Timebins)) +
  theme_bw() +
  ggtitle(NULL) +
  # xlab("Epochs") + 
  xlab("") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title=element_text(size=14)) +
  guides(color = FALSE, fill = FALSE) + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "BrBG"))

disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins

ggplot_perTimebins_height <- 6
ggplot_perTimebins_width <- 12

ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins,
       filename = file.path(output_dir, "Figure_7_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTimebins.svg"),
       width = ggplot_perTimebins_width,
       height = ggplot_perTimebins_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins,
       filename = file.path(output_dir, "Figure_7_barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTimebins.png"),
       width = ggplot_perTimebins_width,
       height = ggplot_perTimebins_height,
       units = "in")



