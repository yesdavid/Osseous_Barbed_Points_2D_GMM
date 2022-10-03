library(ggplot2)
library(readr)

rm(list=ls())


# import CSV
# point_info <- readr::read_csv("Barbed_points_contextual_umlaute.csv")
point_info <- readr::read_csv(file.path("2_data", "Tsirintoulaki_OPAR_barbed_points_edit.csv"),
                              col_types = cols(Period = col_factor(levels = c("Late Pleistocene", 
                                                                              "Early Holocene"))))
# point_info$`Techno-complex` <- as.factor(point_info$`Techno-complex`)

countries <- unique(point_info$Location)
shapes_countries <- seq(1,length(countries))
names(shapes_countries) <- countries

shapes_countries_df <- as.data.frame(shapes_countries)
shapes_countries_df$Location <- row.names(shapes_countries_df)
point_info <- dplyr::left_join(point_info, shapes_countries_df, by = "Location")

techno_complex <- unique(point_info$`Techno-complex`)
shapes_typology <- seq(1,length(techno_complex))
names(shapes_typology) <- techno_complex
names(shapes_typology)[which(is.na(names(shapes_typology)))] <- "unknown"

periods <- unique(point_info$Period)
shapes_periods <- seq(1,length(periods))
names(shapes_periods) <- periods

types <- unique(point_info$Typology)
shapes_types <- seq(1,length(types))
names(shapes_types) <- types     

countries <- unique(point_info$Location)
shapes_countries <- seq(1,length(countries))
names(shapes_countries) <- countries

shapes_typology_df <- data.frame(technocomplex_shape = shapes_typology)
shapes_typology_df$`Techno-complex` <- row.names(shapes_typology_df)
point_info <- dplyr::left_join(point_info, shapes_typology_df, by = "Techno-complex")

world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
clipper_europe <- as(raster::extent(-5, 30, 41, 57.5), "SpatialPolygons")
sp::proj4string(clipper_europe) <- sp::CRS(sp::proj4string(world))
world_clip <- raster::intersect(world, clipper_europe)
world_clip_f <- ggplot2::fortify(world_clip)


### distribution map

current_distribution_map <- 
  ggplot() +
    geom_polygon(data = world_clip_f,
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "grey") +
    ggrepel::geom_label_repel(data = point_info,
                              aes(x = Longitude, y = Latitude,
                                  label = ID,
                                  fill = Period
                                  ),
                              size = 3,
                              box.padding = 0.15,
                              point.padding = 0.25,
                              segment.color = 'grey50') +
  # geom_point(data = point_info[,c("Site", "Latitude", "Longitude")],  
  #            aes(x = Longitude, y = Latitude, alpha = 0.9), 
  #            shape = 3) +
  geom_jitter(data = point_info,
              aes(x = Longitude, y = Latitude,
                  # color = Period,
                  shape = Typology),
              size = 6, #stroke = 1,
              width = 0.15, 
              height = 0.15) +
    scale_shape_identity() +
  # scale_shape_manual(values=shapes_countries) +
    scale_shape_manual(values = shapes_types) +
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
  # # theme(legend.position = "none")

plot_width_cm <- 40
plot_height_cm <- 25

ggsave(filename = file.path("3_output", "Figure_1_distribution_map.svg"),
       plot = current_distribution_map,
       device = "svg",
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm")
ggsave(filename = file.path("3_output", "Figure_1_distribution_map.png"),
       plot = current_distribution_map,
       device = "png",
       width = plot_width_cm,
       height = plot_height_cm,
       units = "cm")
