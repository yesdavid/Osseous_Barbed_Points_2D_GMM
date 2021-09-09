# overall variability/disparity

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

IDs_w_more_than_X_barbs <- Momocs::filter(open_outlines, 
                                          ID %in% names(which(table(open_outlines$ID) >= 3))  # just retain IDs with more than X samples
                                          # & !(ID %in% "BTZ_3")
                                          ) 



## dfourier
open_outlines_dfourier <- Momocs::dfourier(IDs_w_more_than_X_barbs)           # discrete cosine transform.

## PCA
open_outlines_PCA <- Momocs::PCA(open_outlines_dfourier,
                                 fac = point_info)             # we calculate a PCA on it
rownames(open_outlines_PCA$x) <- open_outlines_PCA$fac$ID_barb

# Momocs::plot_PCA(open_outlines_PCA,
#                  ~Site,
#                  labelgroups = F,
#                  chull = T)
# Momocs::plot_PCA(open_outlines_PCA,
#                  ~Period,
#                  labelgroups = F,
#                  chull = T)


######################################################
# disparity of each ID
rownames_DATASETS <- list()
for(i in unique(open_outlines_PCA$ID)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, ID == i)$ID_barb)
}


ID_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

ID_boot <- dispRity::boot.matrix(ID_subsets, bootstraps = 100)
ID_disp <- dispRity::dispRity(ID_boot, metric = c(sum, variances))
summary(ID_disp)

# Wilcox.test
test.dispRity(ID_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "bonferroni")
# PERMANOVA
test.dispRity(ID_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "bonferroni")


ID_names <- names(ID_disp$disparity)
disparity_df_list <- list()
for(i in ID_names){
  disparity_df_list[[i]] <- data.frame(ID_w_counts = paste0(i, 
                                                       "\n(n=",nrow(ID_disp$subsets[[i]]$elements),")"),
                                       disparity = as.vector(ID_disp$disparity[[i]][[2]]),
                                       nelements = nrow(ID_disp$subsets[[i]]$elements),
                                       ID = i)
}
disparity_df_IDdiscrete_armatureOutlines_perID <- do.call(rbind.data.frame, disparity_df_list)

disparity_IDdiscrete_armatureOutlines_ggplot_perID <- 
  ggplot(data = disparity_df_IDdiscrete_armatureOutlines_perID,
         aes(x = ID_w_counts, y = disparity,
             )) +
  # geom_violin() + 
  geom_boxplot(notch = T, width = 0.4, aes(fill = ID)) +
  theme_bw() +
  ggtitle(NULL) +
  xlab("") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title=element_text(size=14)) +
  # ggthemes::scale_fill_colorblind() +
  guides(color = FALSE, fill = FALSE) + 
  coord_flip()

disparity_IDdiscrete_armatureOutlines_ggplot_perID

ggplot_perID_height <- 20
ggplot_perID_width <- 15

ggsave(disparity_IDdiscrete_armatureOutlines_ggplot_perID,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perID.svg"),
       width = ggplot_perID_width,
       height = ggplot_perID_height,
       units = "in")
ggsave(disparity_IDdiscrete_armatureOutlines_ggplot_perID,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perID.png"),
       width = ggplot_perID_width,
       height = ggplot_perID_height,
       units = "in")

disparity_df_IDdiscrete_armatureOutlines_perID$disparity_mean_per_ID <- 0
for(i in unique(disparity_df_IDdiscrete_armatureOutlines_perID$ID)) {
  mean_disp <- mean(disparity_df_IDdiscrete_armatureOutlines_perID[which(disparity_df_IDdiscrete_armatureOutlines_perID$ID == i), "disparity"])
  disparity_df_IDdiscrete_armatureOutlines_perID[which(disparity_df_IDdiscrete_armatureOutlines_perID$ID == i), "disparity_mean_per_ID"] <- mean_disp
}

aaa <- dplyr::left_join(disparity_df_IDdiscrete_armatureOutlines_perID, 
                 distinct(point_info, ID, .keep_all = T), by = c("ID" = "ID"))


ggplot() +
  geom_polygon(data = world_clip_f,
               aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey") +
  # ggrepel::geom_label_repel(data = point_info,
  #                           aes(x = Longitude, y = Latitude,
  #                               label = ID,
  #                               fill = Period
  #                           ),
  #                           size = 3,
  #                           box.padding = 0.15,
  #                           point.padding = 0.25,
  #                           segment.color = 'grey50') +
  # geom_point(data = point_info[,c("Site", "Latitude", "Longitude")],  
  #            aes(x = Longitude, y = Latitude, alpha = 0.9), 
  #            shape = 3) +
  geom_point(data = aaa,
              aes(x = Longitude, y = Latitude,
                  color = disparity_mean_per_ID,
                  # shape = Typology
                  ),
              shape = 21,
              size = 1, #stroke = 1,
              width = 0.35, 
              height = 0.35) +
  scale_shape_identity() +
  scale_color_gradient(low = "green",
                         high = "red")
  # scale_shape_manual(values=shapes_countries) +
  # scale_shape_manual(values = shapes_types) +
  coord_quickmap() +  
  theme_classic() +  
  xlab("Longitude") +
  ylab("Latitude") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) # +

cor(x = aaa$disparity, y = aaa$Latitude)
cor.test(x = aaa$disparity, y = aaa$Latitude)

ggplot(data = aaa,
       aes(y = disparity_mean_per_ID,
           x = Latitude)) + 
  geom_point() +
  geom_smooth(method = lm)
######################################################
# disparity of each ID arranged on X-axis by 14C-date


disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP <- dplyr::left_join(disparity_df_IDdiscrete_armatureOutlines_perID, point_info, by = "ID")


disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP <- 
  ggplot(data = disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP,
         aes(x = mean_age_calBP, 
             y = disparity,
         )) +
  geom_boxplot(notch = T, 
               width = 100,
               aes(fill = Location, 
                   group = ID_w_counts)) +
  xlim(16000,9500) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_reverse(expand = c(0.01,0.01),
                  limits = c(15850,9500),
                  breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + 
  # geom_smooth(method = "gam", formula = y~s(x),
  #             span = 0.1,
  #             fullrange = F) +
  # facet_wrap(~Location) +
  ylab("Disparity (sum of variances)") +
  xlab("Years in calBP") +
  theme(legend.position = "bottom")

disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP
  
ggplot_perID_mean_age_calBP_height <- 10
ggplot_perID_mean_age_calBP_width <- 15

ggsave(disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP.svg"),
       width = ggplot_perID_mean_age_calBP_width,
       height = ggplot_perID_mean_age_calBP_height,
       units = "in")
ggsave(disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP.png"),
       width = ggplot_perID_mean_age_calBP_width,
       height = ggplot_perID_mean_age_calBP_height,
       units = "in")



NGRIP1_data <- readr::read_csv(file.path("2_data", "NGRIP1_d18O_and_dust_5cm.csv")) %>% 
  select(c(`Delta O18 (permil)`, `GICC05 age (yr b2k)`))
NGRIP2_data <- readr::read_csv(file.path("2_data", "NGRIP2_d18O_and_dust_5cm.csv")) %>% 
  select(c(`Delta O18 (permil)`, `GICC05 age (yr b2k)`))

NGRIP_data <- rbind(NGRIP1_data, NGRIP2_data)

NGRIP_data$`GICC05 age (yr b2k)` <- NGRIP_data$`GICC05 age (yr b2k)` - 50 # to fit the calBP, since they have 1950 as "0"

NGRIP_data_subset <- 
  subset(NGRIP_data, 
       `GICC05 age (yr b2k)` > min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$mean_age_calBP,
                                   na.rm = T-100) &
         `GICC05 age (yr b2k)` < max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$mean_age_calBP,
                                     na.rm = T)+100) %>% 
  select(c(`Delta O18 (permil)`, `GICC05 age (yr b2k)`))


ngrip_plot <-
  ggplot(data = NGRIP_data_subset) +
  geom_line(aes(x = `GICC05 age (yr b2k)`,
                y = `Delta O18 (permil)`,
                alpha = 0.8)) +
  geom_smooth(aes(x = `GICC05 age (yr b2k)`, 
                  y = `Delta O18 (permil)`,
                  color = "red"),
              span = 0.001) +
  xlim(16000,9500) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_reverse(expand = c(0.01,0.01),
                  limits = c(15850,9500),
                  breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + 
  ylab("Delta O18 (permil)") +
  xlab("Years BP") +
  theme(legend.position = "none")

# disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP + 
#   annotation_custom(ggplotGrob(ngrip_plot), 
#                     xmin = 900, xmax = 1900,
#                     ymin = 50, ymax = 100)


cowplot::plot_grid(
  ngrip_plot,
  disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP,
  labels = "AUTO", ncol = 1
)

# # Function factory for secondary axis transforms, from https://stackoverflow.com/a/66055331
# library(scales)


######################################################
# disparity of each Period/TimeSlice
rownames_DATASETS <- list()
for(i in unique(open_outlines_PCA$Period)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Period == i)$ID_barb)
}
rownames_DATASETS

TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 100)
TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

# Wilcox.test
test.dispRity(TS_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "bonferroni")
# PERMANOVA
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "bonferroni")


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

disparity_df_TSdiscrete_armatureOutlines_perTShard$Period <- relevel(disparity_df_TSdiscrete_armatureOutlines_perTShard$Period, "Late Pleistocene\n(n=75)")

disparity_TSdiscrete_armatureOutlines_ggplot_perTShard <- 
  ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTShard,
         aes(x = Period, 
             y = disparity
         )) +
  # geom_violin() + 
  geom_boxplot(notch = T, 
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
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.png"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")


######################################################
# disparity of each Typology
rownames_DATASETS <- list()
for(i in unique(open_outlines_PCA$Typology)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Typology == i)$ID_barb)
}
rownames_DATASETS

TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
                                       group = rownames_DATASETS)

TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 100)
TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

# Wilcox.test
test.dispRity(TS_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "bonferroni")
# PERMANOVA
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "bonferroni")


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

disparity_df_TSdiscrete_armatureOutlines_perTShard$Period <- relevel(disparity_df_TSdiscrete_armatureOutlines_perTShard$Period, "Late Pleistocene\n(n=75)")

disparity_TSdiscrete_armatureOutlines_ggplot_perTShard <- 
  ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTShard,
         aes(x = Period, 
             y = disparity
         )) +
  # geom_violin() + 
  geom_boxplot(notch = T, 
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
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.png"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
