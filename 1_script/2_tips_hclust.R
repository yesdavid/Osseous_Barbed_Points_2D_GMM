library(Momocs)
library(ggplot2)
library(NbClust)
library(ggtree)

rm(list=ls())

point_info_raw <- readr::read_csv(file.path("2_data", "Tsirintoulaki_OPAR_barbed_points_edit.csv"),
                                  col_types = cols(Period = col_factor(levels = c("Late Pleistocene", 
                                                                                  "Early Holocene"))))

output_dir <- file.path("3_output", "tips")


open_outlines_raw <- readRDS(file.path(output_dir, "tips_open_outlines_OpnCoo.RDS"))
# point_info <- open_outlines$fac


point_info <- subset(point_info_raw, ID %in% names(open_outlines_raw$coo))

open_outlines <- Momocs::Opn(x = open_outlines_raw$coo,
                             fac = point_info)


## dfourier
open_outlines_dfourier <- Momocs::dfourier(open_outlines)           # discrete cosine transform.

## PCA
open_outlines_PCA <- Momocs::PCA(open_outlines_dfourier,
                                 fac = point_info)             # we calculate a PCA on it
# rownames(open_outlines_PCA$x) <- open_outlines_PCA$fac$ID



scree_plot <- Momocs::scree_plot(open_outlines_PCA,
                                 nax = 1:5) +  
  ggplot2::theme_bw()

scree_plot

# minimum_no_of_pcs <- ncol(open_outlines_PCA$x)
minimum_no_of_pcs <- 4

gg <- Momocs::PCcontrib(open_outlines_PCA,
                        nax = 1:5,
                        sd.r = c(-1.5,-1,0,1,1.5))
pc_contrib_plot <- gg$gg + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pc_contrib_plot





## hclust
dist <- dist(open_outlines_PCA$x[,c(1:minimum_no_of_pcs)], 
             method = "euclidean")
ward_hclust <- hclust(dist,
                      method = "ward.D2")
plot(ward_hclust)


tanged_points_PCA_NbClust_ward <- NbClust::NbClust(data = open_outlines_PCA$x[,1:minimum_no_of_pcs],
                                                   distance = "euclidean",
                                                   method = "ward.D2",
                                                   index = c("gap", "silhouette"))
TP_NbClust <- tanged_points_PCA_NbClust_ward$All.index
TP_NbClust_df <- as.data.frame(TP_NbClust)
TP_NbClust_df$NClust <- 1:nrow(TP_NbClust)+1
silhouette_plot <- ggplot(TP_NbClust_df, aes(x = NClust, y = Silhouette)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(TP_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Average silhouette value")

silhouette_plot


n_clusters <- 7 

color_palette <- cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#FF6347", "#0072B2", "#D55E00", "#CC79A7")

##### tree cut
current_treecut_df <- data.frame(ID = as.character(names(cutree(ward_hclust,
                                                             k = n_clusters))),
                              cluster = cutree(ward_hclust,
                                                      k = n_clusters))
# current_treecut$ID <- as.character(current_treecut$ID)

current_treecut <- dplyr::left_join(current_treecut_df, point_info, by = "ID")

rownames(current_treecut) <- current_treecut$ID
open_outlines_w_cluster <- Momocs::Opn(open_outlines$coo,
                                       fac = current_treecut)

mypath_svg <- file.path(output_dir, "open_outlines_panel_tips_w_ID")
svg(file=mypath_svg, width = 8, height = 8)
Momocs::panel(open_outlines_w_cluster, 
              names = T)
dev.off()

## dfourier
open_outlines_w_cluster_dfourier <- Momocs::dfourier(open_outlines_w_cluster)           # discrete cosine transform.

## PCA
open_outlines_w_cluster_PCA <- Momocs::PCA(open_outlines_w_cluster_dfourier,
                                           fac = current_treecut)             # we calculate a PCA on it
# rownames(open_outlines_w_cluster_PCA$x) <- open_outlines_w_cluster_PCA$fac$ID


open_outlines_w_cluster_PCA_df <- as.data.frame(open_outlines_w_cluster_PCA$x)
open_outlines_w_cluster_PCA_df$Cluster <- open_outlines_w_cluster_PCA$fac$cluster
open_outlines_w_cluster_PCA_df$Country <- as.factor(open_outlines_w_cluster_PCA$fac$Location)
open_outlines_w_cluster_PCA_df$ID <- as.factor(open_outlines_w_cluster_PCA$fac$ID)
open_outlines_w_cluster_PCA_df$Period <- as.factor(open_outlines_w_cluster_PCA$fac$Period)
open_outlines_w_cluster_PCA_df$Typology <- as.factor(open_outlines_w_cluster_PCA$fac$Typology)

periods <- unique(point_info_raw$Period)
shapes_periods <- seq(1,length(periods))
names(shapes_periods) <- periods

types <- unique(point_info_raw$Typology)
shapes_types <- seq(1,length(types))
names(shapes_types) <- types     

countries <- unique(point_info_raw$Location)
shapes_countries <- seq(1,length(countries))
names(shapes_countries) <- countries


#### PCA plot w cluster

Momocs::plot_PCA(open_outlines_w_cluster_PCA, 
                 ~cluster, 
                 labelpoints = T,
                 legend = T,
                 morphospace = T,
                 chull = F,
                 chullfilled = T, 
                 center_origin = F)


a <- ggplot(data = open_outlines_w_cluster_PCA_df, 
            aes(shape = Typology, 
                color = Period, 
                label = ID)) +
  geom_point(size = 5) +
  geom_text(hjust = -.3) +
  coord_fixed(ratio =1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  # scale_fill_manual(values = color_palette) +
  # scale_shape_manual(values=shapes_countries) +
  scale_shape_manual(values = shapes_types) +
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) +
  coord_cartesian(xlim = c(min(open_outlines_w_cluster_PCA_df[,c("PC1")]) + min(open_outlines_w_cluster_PCA_df[,c("PC1")])*0.1, 
                           max(open_outlines_w_cluster_PCA_df[,c("PC1")]) + max(open_outlines_w_cluster_PCA_df[,c("PC1")])*0.1),
                  ylim = c(min(open_outlines_w_cluster_PCA_df[,c("PC2")]) + min(open_outlines_w_cluster_PCA_df[,c("PC2")])*0.1, 
                           max(open_outlines_w_cluster_PCA_df[,c("PC2")]) + max(open_outlines_w_cluster_PCA_df[,c("PC2")])*0.1))


b <- a + aes(x = PC1, y = PC2) +
  xlab(paste0("PC1 (", round(open_outlines_w_cluster_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC2 (", round(open_outlines_w_cluster_PCA$eig[2]*100, digits = 0), "%)")) #+
  # guides(color = FALSE, fill = FALSE)
b

tips_PCA_pc1pc2_height <- 8
tips_PCA_pc1pc2_width <- 8
ggsave(b,
       filename = file.path(output_dir, "tips_PCA_pc1pc2.svg"),
       width = tips_PCA_pc1pc2_width,
       height = tips_PCA_pc1pc2_height,
       units = "in")
ggsave(b,
       filename = file.path(output_dir, "tips_PCA_pc1pc2.png"),
       width = tips_PCA_pc1pc2_width,
       height = tips_PCA_pc1pc2_height,
       units = "in")


# mean shapes
min_no_of_coordinates <- list()
for (cluster_index in unique(open_outlines_w_cluster$fac$cluster)){
  
  current_shapes <- Momocs::slice(open_outlines_w_cluster, cluster == cluster_index)
  
  min_no_of_coordinates[[cluster_index]] <- rep(NA, length(current_shapes))
  
  for (i in 1:length(current_shapes)){
    min_no_of_coordinates[[cluster_index]][i] <- nrow(current_shapes$coo[[i]])
  }
  
  min_no_of_coordinates[[cluster_index]] <- min(min_no_of_coordinates[[cluster_index]])
  
}

## shapes get INTERPOLATED to common number of landmarks (lowest number of landmarks per cluster)
mean_shapes_cluster <- list()
for (cluster_index in unique(open_outlines_w_cluster$fac$cluster)){
  mean_shapes_cluster[[cluster_index]] <- 
    Momocs::MSHAPES(Momocs::coo_interpolate(
      Momocs::slice(open_outlines_w_cluster, 
                    cluster == cluster_index),
      n = min_no_of_coordinates[[cluster_index]])$coo)
}

open_outlines_w_cluster_mean_shapes_cluster_out <- Momocs::Out(mean_shapes_cluster,
                                                               fac = data.frame(cluster = paste0("cluster_", c(unique(open_outlines_w_cluster$fac$cluster)))))
  
Momocs::panel(open_outlines_w_cluster_mean_shapes_cluster_out,
              # fac = "cluster",
              names = T)



## save mean shapes
mean_shapes_dir <- file.path(output_dir, "mean_shapes")
dir.create(mean_shapes_dir,
           recursive = T)

for (i in unique(open_outlines_w_cluster$fac$cluster)){
  
  mypath_png <- file.path(mean_shapes_dir, paste0("Cluster ", i, " (n=",table(current_treecut$cluster)[[i]], ")_mean_shp.png"))
  
  png(file=mypath_png,
      width = 800, height = 800, units = "px")
  
  Momocs::panel(Momocs::slice(open_outlines_w_cluster_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                main = NULL,
                col = color_palette[as.integer(i)])
  
  dev.off()
  
  mypath_svg <- file.path(mean_shapes_dir, paste0("Cluster ", i, " (n=",table(current_treecut$cluster)[[i]], ")_mean_shp.svg"))
  
  svg(file=mypath_svg, width = 8, height = 8)
  
  Momocs::panel(Momocs::slice(open_outlines_w_cluster_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                main = NULL,
                col = color_palette[as.integer(i)])
  
  dev.off()
}

# save panels for each cluster
panels_dir <- file.path(output_dir, "panels")
dir.create(panels_dir,
           recursive = T)

closed_open_outlines_w_cluster <- Momocs::Out(x = open_outlines_w_cluster$coo,
                                              fac = current_treecut)

for (i in unique(open_outlines_w_cluster$fac$cluster)){
  
  current_panel <- Momocs::slice(closed_open_outlines_w_cluster, cluster == i)
  
  mypath_png <- file.path(panels_dir, paste0("open_outlines_w_cluster_colors_cluster_", i, ".png"))
  
  png(file=mypath_png,
      width = 800, height = 800, units = "px")
  
  Momocs::panel(current_panel,
                main = NULL,
                col = rep(color_palette[as.integer(i)], length(current_panel)),
                names = T)
  
  dev.off()
  
  
  mypath_svg <- file.path(panels_dir, paste0("open_outlines_w_cluster_colors_cluster_", i, ".svg"))
  
  svg(file=mypath_svg, width = 8, height = 8)
  
  Momocs::panel(current_panel,
                main = NULL,
                col = rep(color_palette[as.integer(i)], length(current_panel)),
                names = T)
  
  dev.off()
}


### ggtree


chronozones_TP_df <- point_info[,c("ID", "Location", "Period", "Typology")]
chronozones_TP_df$Location <- as.factor(chronozones_TP_df$Location)
chronozones_TP_df$Period <- as.factor(chronozones_TP_df$Period)
chronozones_TP_df$Typology <- as.factor(chronozones_TP_df$Typology)
rownames(chronozones_TP_df) <- point_info$ID

# subset(open_outlines_w_cluster_PCA_df, Cluster == 7) # check the IDs for each cluster to find the right clades


tree_w_clusterlabels <- ggtree(ward_hclust) %<+% chronozones_TP_df + 
  # geom_tiplab(size=3,
  #             aes(color = Location)) +
  geom_tiplab(size=4, offset = 0.5,
              aes(color = Period)) +
  geom_tippoint(aes(shape=Typology, hjust = 1), size = 4) +
  scale_shape_manual(values=shapes_types) +
  geom_treescale() + 
  scale_colour_discrete(na.translate = F) + 
  #geom_text(aes(label=node)) + # can be used to determine the node numbers for the chosen clusters
  geom_cladelabel(node=37, label="Cluster 1", align=T, geom='text', offset = 2.5, color=c("#E69F00"), barsize = 2, fontsize = 6) + 
  geom_cladelabel(node=39, label="Cluster 2", align=T, geom='text', offset = 2.5, color=c("#56B4E9"), barsize = 2, fontsize = 6) + 
  geom_cladelabel(node=7, label="Cluster 3", align=T, geom='text', offset = 2.5, color=c("#009E73"), barsize = 2, fontsize = 6) +
  geom_cladelabel(node=38, label="Cluster 4", align=T, geom='text', offset = 2.5, color=c("#FF6347"), barsize = 2, fontsize = 6) +
  geom_cladelabel(node=12, label="Cluster 5", align=T, geom='text', offset = 2.5, color=c("#0072B2"), barsize = 2, fontsize = 6) + 
  geom_cladelabel(node=40, label="Cluster 6", align=T, geom='text', offset = 2.5, color=c("#D55E00"), barsize = 2, fontsize = 6) +
  geom_cladelabel(node=27, label="Cluster 7", align=T, geom='text', offset = 2.5, color=c("#CC79A7"), barsize = 2, fontsize = 6) +
  xlim_tree(46) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

tree_w_clusterlabels

tree_w_clusterlabels_height <- 12
tree_w_clusterlabels_width <- 16

ggsave(tree_w_clusterlabels,
       filename = file.path(output_dir, "tips_tree_w_clusterlabels.svg"),
       width = tree_w_clusterlabels_width,
       height = tree_w_clusterlabels_height,
       units = "in")
ggsave(tree_w_clusterlabels,
       filename = file.path(output_dir, "tips_tree_w_clusterlabels.png"),
       width = tree_w_clusterlabels_width,
       height = tree_w_clusterlabels_height,
       units = "in")

