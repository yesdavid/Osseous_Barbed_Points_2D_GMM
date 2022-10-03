# sum of variance between tips, bases, and barbs

current_elementText_size <- 16

barbs_open_outlines <- readRDS(file.path("3_output", "barbs", "barbs_open_outlines_OpnCoo.RDS"))
barbs_open_outlines$fac$type <- "Barbs"
barbs_open_outlines <- Momocs::Opn(barbs_open_outlines$coo,
                                   fac = barbs_open_outlines$fac[,c("ID_barb", "type")])
names(barbs_open_outlines$fac) <- c("ID", "type")

tips_open_outlines <- readRDS(file.path("3_output", "tips", "tips_open_outlines_OpnCoo.RDS"))
tips_open_outlines$fac$type <- "Tips"
tips_open_outlines <- Momocs::Opn(tips_open_outlines$coo,
                                   fac = tips_open_outlines$fac[,c("ID", "type")])

bases_open_outlines <- readRDS(file.path("3_output", "bases", "bases_open_outlines_OpnCoo.RDS"))
bases_open_outlines$fac$type <- "Bases"
bases_open_outlines <- Momocs::Opn(bases_open_outlines$coo,
                                   fac = bases_open_outlines$fac[,c("ID", "type")])

# combined pca?
combined_outlines <- 
Momocs::combine(barbs_open_outlines,
                tips_open_outlines,
                bases_open_outlines)

## dfourier
open_outlines_dfourier <- Momocs::dfourier(combined_outlines)           # discrete cosine transform.

## PCA
open_outlines_PCA <- Momocs::PCA(open_outlines_dfourier,
                                 fac = point_info)             # we calculate a PCA on it
rownames(open_outlines_PCA$x) <- combined_outlines$fac$ID


Momocs::plot_PCA(open_outlines_PCA,
                 ~type)

open_outlines_PCA_df <- as.data.frame(open_outlines_PCA$x)
open_outlines_PCA_df$ID <- as.factor(open_outlines_PCA$fac$ID)
open_outlines_PCA_df$type <- as.factor(open_outlines_PCA$fac$type)


a <- ggplot(data = open_outlines_PCA_df, 
            aes(fill = type#, 
                # label = ID
                )) +
  geom_point(size = 5,
             shape = 21) +
  # shadowtext::geom_shadowtext(hjust = -.3, 
  #                             size = 6,
  #                             bg.color = "white") + 
  # coord_fixed(ratio =1) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = current_elementText_size),
        axis.text = element_text(size = current_elementText_size),
        legend.text = element_text(size = current_elementText_size-2),
        legend.title = element_text(size = current_elementText_size)) +
  guides(shape = guide_legend(nrow =7,
                              title.position = c("top")),
         color = guide_legend(nrow =7,
                              title.position = c("top"))) +
  # scale_fill_manual(values = color_palette) +
  # scale_shape_manual(values=shapes_countries) +
  # scale_shape_manual(values = shapes_types) +
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) +
  coord_cartesian(xlim = c(min(open_outlines_PCA_df[,c("PC1")]) + min(open_outlines_PCA_df[,c("PC1")])*0.1, 
                           max(open_outlines_PCA_df[,c("PC1")]) + max(open_outlines_PCA_df[,c("PC1")])*0.1),
                  ylim = c(min(open_outlines_PCA_df[,c("PC2")]) + min(open_outlines_PCA_df[,c("PC2")])*0.1, 
                           max(open_outlines_PCA_df[,c("PC2")]) + max(open_outlines_PCA_df[,c("PC2")])*0.1)) + 
  labs(fill ="Type")


b <- a + aes(x = PC1, y = PC2) +
  xlab(paste0("PC1 (", round(open_outlines_PCA$eig[1]*100, digits = 1), "%)")) +
  ylab(paste0("PC2 (", round(open_outlines_PCA$eig[2]*100, digits = 1), "%)")) #+
# guides(color = FALSE, fill = FALSE)
b 


# disparity of each ID
rownames_DATASETS <- list()
for(i in unique(combined_outlines$type)){
  rownames_DATASETS[[i]] <- as.character(subset(combined_outlines$fac, type == i)$ID)
}



ID_subsets <- dispRity::custom.subsets(open_outlines_PCA$x, 
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


disparity_IDdiscrete_armatureOutlines_ggplot_perID <- 
  ggplot(data = disparity_df_IDdiscrete_armatureOutlines_perID,
         aes(x = ID_w_counts, y = disparity,
         )) +
  # geom_violin() + 
  geom_boxplot(notch = F, width = 0.4) +
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
  # coord_flip() +
  labs(fill="Latitudinal Zone") +
  scale_fill_manual(values = terrain.colors(4))

disparity_IDdiscrete_armatureOutlines_ggplot_perID
