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
rownames(open_outlines_PCA$x) <- open_outlines_PCA$fac$ID_barb

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
for(i in unique(open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA$ID)){
  rownames_DATASETS[[i]] <- as.character(subset(open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA$fac, ID == i)$ID_barb)
}


ID_subsets <- dispRity::custom.subsets(open_outlines_dfourier_IDs_w_more_than_X_barbs_PCA$x, 
                                       group = rownames_DATASETS)

ID_boot <- dispRity::boot.matrix(ID_subsets, bootstraps = 1000)
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

disparity_df_IDdiscrete_armatureOutlines_perID <- dplyr::left_join(disparity_df_IDdiscrete_armatureOutlines_perID, point_info, by = "ID")
disparity_df_IDdiscrete_armatureOutlines_perID <- disparity_df_IDdiscrete_armatureOutlines_perID[order(disparity_df_IDdiscrete_armatureOutlines_perID$Latitude_zone),]
disparity_df_IDdiscrete_armatureOutlines_perID$ID <- factor(disparity_df_IDdiscrete_armatureOutlines_perID$ID,
                                                               levels = unique(disparity_df_IDdiscrete_armatureOutlines_perID[order(disparity_df_IDdiscrete_armatureOutlines_perID$Latitude_zone),"ID"]))
disparity_df_IDdiscrete_armatureOutlines_perID$ID_w_counts <- factor(disparity_df_IDdiscrete_armatureOutlines_perID$ID_w_counts,
                                                            levels = (unique(disparity_df_IDdiscrete_armatureOutlines_perID[order(disparity_df_IDdiscrete_armatureOutlines_perID$Latitude),"ID_w_counts"])))


disparity_IDdiscrete_armatureOutlines_ggplot_perID <- 
  ggplot(data = disparity_df_IDdiscrete_armatureOutlines_perID,
         aes(x = ID_w_counts, y = disparity,
             )) +
  # geom_violin() + 
  geom_boxplot(notch = F, width = 0.4, aes(fill = Latitude_zone)) +
  theme_bw() +
  ggtitle(NULL) +
  xlab("") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title=element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  # ggthemes::scale_fill_colorblind() +
  guides(color = FALSE) + 
  coord_flip() +
  labs(fill="Latitudinal Zone") +
  scale_fill_manual(values = terrain.colors(4))

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

# disparity_df_IDdiscrete_armatureOutlines_perID$disparity_mean_per_ID <- 0
# for(i in unique(disparity_df_IDdiscrete_armatureOutlines_perID$ID)) {
#   mean_disp <- mean(disparity_df_IDdiscrete_armatureOutlines_perID[which(disparity_df_IDdiscrete_armatureOutlines_perID$ID == i), "disparity"])
#   disparity_df_IDdiscrete_armatureOutlines_perID[which(disparity_df_IDdiscrete_armatureOutlines_perID$ID == i), "disparity_mean_per_ID"] <- mean_disp
# }
# 
# aaa <- dplyr::left_join(disparity_df_IDdiscrete_armatureOutlines_perID, 
#                  distinct(point_info, ID, .keep_all = T), by = c("ID" = "ID"))
# 
# 
# ggplot() +
#   geom_polygon(data = world_clip_f,
#                aes(x = long, y = lat, group = group),
#                fill = NA, colour = "grey") +
#   # ggrepel::geom_label_repel(data = point_info,
#   #                           aes(x = Longitude, y = Latitude,
#   #                               label = ID,
#   #                               fill = Period
#   #                           ),
#   #                           size = 3,
#   #                           box.padding = 0.15,
#   #                           point.padding = 0.25,
#   #                           segment.color = 'grey50') +
#   # geom_point(data = point_info[,c("Site", "Latitude", "Longitude")],  
#   #            aes(x = Longitude, y = Latitude, alpha = 0.9), 
#   #            shape = 3) +
#   geom_point(data = aaa,
#               aes(x = Longitude, y = Latitude,
#                   color = disparity_mean_per_ID,
#                   # shape = Typology
#                   ),
#               shape = 21,
#               size = 1, #stroke = 1,
#               width = 0.35, 
#               height = 0.35) +
#   scale_shape_identity() +
#   scale_color_gradient(low = "green",
#                          high = "red") +
#   # scale_shape_manual(values=shapes_countries) +
#   # scale_shape_manual(values = shapes_types) +
#   coord_quickmap() +  
#   theme_classic() +  
#   xlab("Longitude") +
#   ylab("Latitude") + 
#   scale_y_continuous(expand = c(0,0)) + 
#   scale_x_continuous(expand = c(0,0)) +
#   theme(legend.position = "right",
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16)) # +
# 
# cor(x = aaa$disparity, y = aaa$Latitude.x)
# cor.test(x = aaa$disparity, y = aaa$Latitude.x)
# 
# ggplot(data = aaa,
#        aes(y = disparity_mean_per_ID,
#            x = Latitude.x)) + 
#   geom_point() +
#   geom_smooth(method = lm)
######################################################
# disparity of each ID arranged on X-axis by 14C-date


# disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP <- dplyr::left_join(disparity_df_IDdiscrete_armatureOutlines_perID, point_info, by = "ID")
disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP <- disparity_df_IDdiscrete_armatureOutlines_perID

disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP <- 
  ggplot(data = disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP,
         aes(x = Age_calBP_mean, 
             y = disparity)) +
  geom_boxplot(notch = F, 
               width = 100,
               aes(fill = Latitude_zone, 
                   group = ID_w_counts)) +
  scale_fill_manual(values = terrain.colors(4)) +
  labs(fill = "Latitudinal Zone") +
  xlim(max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)+250,
       min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)-250) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_reverse(expand = c(0.01,0.01),
                  limits = c(max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)+250,
                             min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)-250),
                  breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + 
  # geom_smooth(method = "gam", formula = y~s(x),
  #             # span = 0.5,
  #             fullrange = F) +
  # facet_wrap(~Location) +
  ylab("Disparity (sum of variances)") +
  xlab("14C Age in calBP") +
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


library(magrittr)
library(dplyr)
NGRIP1_data <- readr::read_csv(file.path("2_data", "NGRIP1_d18O_and_dust_5cm.csv")) %>% 
  select(c(`Delta O18 (permil)`, `GICC05 age (yr b2k)`))
NGRIP2_data <- readr::read_csv(file.path("2_data", "NGRIP2_d18O_and_dust_5cm.csv")) %>% 
  select(c(`Delta O18 (permil)`, `GICC05 age (yr b2k)`))

NGRIP_data <- rbind(NGRIP1_data, NGRIP2_data)

NGRIP_data$`GICC05 age (yr b2k)` <- NGRIP_data$`GICC05 age (yr b2k)` + 1950 # as dates are in b2k -> years before _common era_; 50 years difference becuase of BP

NGRIP_data_subset <- 
  subset(NGRIP_data, 
       `GICC05 age (yr b2k)` > min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean,
                                   na.rm = T)-250 &
         `GICC05 age (yr b2k)` < max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean,
                                     na.rm = T)+250) %>% 
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
  xlim(max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)+250,
        min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)-250) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_reverse(expand = c(0,0),
                  limits = c(max(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)+250,
                             min(disparity_df_IDdiscrete_armatureOutlines_perID_mean_age_calBP$Age_calBP_mean)-250),
                  breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() + 
  ylab("Delta O18 (permil)") +
  xlab("Years BP") +
  theme(legend.position = "none")


# disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP + 
#   annotation_custom(ggplotGrob(ngrip_plot), 
#                     xmin = 900, xmax = 1900,
#                     ymin = 50, ymax = 100)

ngrip_cowplot <- 
cowplot::plot_grid(
  ngrip_plot,
  disparity_IDdiscrete_armatureOutlines_ggplot_perID_mean_age_calBP,
  labels = "AUTO", ncol = 1,
  align = "hv"
)

ngrip_cowplot_height <- 10
ngrip_cowplot_width <- 15

ggsave(ngrip_cowplot,
       filename = file.path(output_dir, "barbs_ngrip_cowplot.svg"),
       width = ngrip_cowplot_width,
       height = ngrip_cowplot_height,
       units = "in")
ggsave(ngrip_cowplot,
       filename = file.path(output_dir, "barbs_ngrip_cowplot.png"),
       width = ngrip_cowplot_width,
       height = ngrip_cowplot_height,
       units = "in")



# # Function factory for secondary axis transforms, from https://stackoverflow.com/a/66055331
# library(scales)


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
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.png"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")


# ######################################################
# # disparity of each Typology
# rownames_DATASETS <- list()
# for(i in unique(open_outlines_PCA$Typology)){
#   rownames_DATASETS[[i]] <- as.character(subset(open_outlines_PCA$fac, Typology == i)$ID_barb)
# }
# rownames_DATASETS
# 
# TS_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
#                                        group = rownames_DATASETS)
# 
# TS_boot <- dispRity::boot.matrix(TS_subsets, bootstraps = 1000)
# TS_disp <- dispRity::dispRity(TS_boot, metric = c(sum, variances))
# summary(TS_disp)
# 
# # Wilcox.test
# test.dispRity(TS_disp, 
#               test = wilcox.test, 
#               comparisons = "pairwise",
#               correction = "bonferroni")
# # PERMANOVA
# test.dispRity(TS_disp, 
#               test = adonis.dispRity, 
#               comparisons = "pairwise",
#               correction = "bonferroni")
# 
# 
# TS_names <- names(TS_disp$disparity)
# disparity_df_list <- list()
# for(i in TS_names){
#   disparity_df_list[[i]] <- data.frame(Period = paste0(i, 
#                                                        "\n(n=",nrow(TS_disp$subsets[[i]]$elements),")"),
#                                        disparity = as.vector(TS_disp$disparity[[i]][[2]]),
#                                        nelements = nrow(TS_disp$subsets[[i]]$elements),
#                                        TS = i)
# }
# disparity_df_TSdiscrete_armatureOutlines_perTShard <- do.call(rbind.data.frame, disparity_df_list)
# 
# disparity_df_TSdiscrete_armatureOutlines_perTShard$Period <- relevel(disparity_df_TSdiscrete_armatureOutlines_perTShard$Period, "Late Pleistocene\n(n=75)")
# 
# disparity_TSdiscrete_armatureOutlines_ggplot_perTShard <- 
#   ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTShard,
#          aes(x = Period, 
#              y = disparity
#          )) +
#   # geom_violin() + 
#   geom_boxplot(notch = T, 
#                width = 0.1, 
#                aes(fill = Period)) +
#   theme_bw() +
#   ggtitle(NULL) +
#   xlab("") + 
#   ylab("Disparity (sum of variances)") +
#   theme(plot.title = element_text(hjust = 0.5, 
#                                   size = 14, 
#                                   face = "bold"),
#         axis.text=element_text(size=14), #,face="bold"
#         axis.text.x = element_text(#angle = 90, 
#           vjust = 0.5, 
#           hjust = 0.5),
#         axis.title=element_text(size=14)) +
#   guides(color = FALSE, 
#          fill = FALSE) #+ 
# # coord_flip()
# 
# disparity_TSdiscrete_armatureOutlines_ggplot_perTShard
# 
# ggplot_perTShard_height <- 6
# ggplot_perTShard_width <- 6
# 
# ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
#        filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.svg"),
#        width = ggplot_perTShard_width,
#        height = ggplot_perTShard_height,
#        units = "in")
# ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
#        filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTShard.png"),
#        width = ggplot_perTShard_width,
#        height = ggplot_perTShard_height,
#        units = "in")




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
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perLatitude_zone.svg"),
       width = ggplot_perTShard_width,
       height = ggplot_perTShard_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perLatitude_zone,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perLatitude_zone.png"),
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
  geom_violin(notch = F,
              # width = 0.1,
              aes(fill = Timebins)) +
  geom_boxplot(notch = F,
               width = 0.2,
               color = "black",
               fill = "white") +
  # geom_boxplot(notch = F,
  #              aes(fill = Timebins)) +
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
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTimebins.svg"),
       width = ggplot_perTimebins_width,
       height = ggplot_perTimebins_height,
       units = "in")
ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins,
       filename = file.path(output_dir, "barbs_disparity_IDdiscrete_armatureOutlines_ggplot_perTimebins.png"),
       width = ggplot_perTimebins_width,
       height = ggplot_perTimebins_height,
       units = "in")




ggplot(point_info_raw) +
  ggridges::geom_density_ridges(aes(x = Age_calBP_mean, y = 1)) + 
  scale_x_reverse()

ggplot(point_info_raw) +
  ggridges::geom_density_ridges(aes(x = Age_calBP_mean, y = 1), 
                                jittered_points = TRUE,
                                position = ggridges::position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', 
                                point_size = 3, 
                                point_alpha = 1, 
                                alpha = 0.7,
                                panel_scaling = T) + 
  scale_x_reverse() +
  theme_bw()

ggplot(point_info_raw) +
  # geom_histogram(aes(x = Age_calBP_mean, fill = Period), 
  #                binwidth = 100,
  #                color = "black") + 
  geom_density(aes(x = Age_calBP_mean, fill = Period), 
               alpha = 0.3) +
  geom_rug(aes(x = Age_calBP_mean, y = 0, color = Period), position = position_jitter(height = 0)) +
  scale_x_reverse() +
  theme_bw()

geom_histogram(breaks=breaks,aes(x=vector,y=..density..), position="identity") + 
  geom_density(aes(x=vector,y=..density..))

######################################################################
######################################################################

# disparity by timebin, for each lat. zone

disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat <- list()
for(current_lat_zone in unique(open_outlines_PCA$fac$Latitude_zone)){
  
  current_lat_zone_subset <- subset(open_outlines_PCA$fac, Latitude_zone == current_lat_zone)
  
  
  # disparity by timebins
  rownames_DATASETS <- list()
  for(i in unique(open_outlines_PCA$fac$Timebins)){
    if(length(as.character(subset(current_lat_zone_subset, Timebins == i)$ID_barb)) > 0){
      rownames_DATASETS[[i]] <- as.character(subset(current_lat_zone_subset, Timebins == i)$ID_barb)
      }
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
  # summary(TS_disp)
  
  # # Wilcox.test
  test.dispRity(TS_disp,
                test = wilcox.test,
                comparisons = "pairwise",
                correction = "bonferroni")
  # PERMANOVA
  adonis_timebins <-
    test.dispRity(TS_disp,
                  test = adonis.dispRity,
                  comparisons = "pairwise",
                  correction = "bonferroni")
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
                                                                          levels = rev(unique(disparity_df_TSdiscrete_armatureOutlines_perTimebins$Timebins)))
  
  
  disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat[[current_lat_zone]] <- 
    ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTimebins,
           aes(x = (Timebins), 
               y = disparity)) +
    # geom_violin(notch = F, 
    #             # width = 0.1, 
    #             aes(fill = Timebins)) +
    # geom_boxplot(notch = F,
    #              width = 0.1,
    #              color = "black",
    #              fill = "white") +
    geom_boxplot(notch = F,
                 aes(fill = Timebins)) +
    theme_bw() +
    ggtitle(NULL) +
    xlab("Epochs") + 
    ylab("Disparity (sum of variances)") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text=element_text(size=14), #,face="bold"
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
          axis.title=element_text(size=14)) +
    guides(color = FALSE, fill = FALSE) + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, name = "BrBG"))
  
  
  
  ggplot_perTimebins_height <- 6
  ggplot_perTimebins_width <- 12

  ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat,
         filename = file.path(output_dir, paste0("disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat_", current_lat_zone, ".svg")),
         width = ggplot_perTShard_width,
         height = ggplot_perTShard_height,
         units = "in")
  ggsave(disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat,
         filename = file.path(output_dir, paste0("disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat_", current_lat_zone, ".png")),
         width = ggplot_perTimebins_width,
         height = ggplot_perTimebins_height,
         units = "in")
  
  
  
}
disparity_TSdiscrete_armatureOutlines_ggplot_perTimebins_per_lat

