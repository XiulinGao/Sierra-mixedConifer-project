library(ggplot2)

c_data <- read.csv('data/PRISM_normals_1991-2020.csv')

p <- ggplot(c_data, aes(x = tmean_degrees_C, y = ppt_mm, label = Site)) +
  geom_point() + 
  geom_text_repel(
                  fontface = 'bold',
                  nudge_x = c(-0.5,0.5,0.5),
                  nudge_y = c(0,0,0)) +
  labs(x = 'Mean temperature (degrees C)',
       y = 'Mean precipitation (mm)') +
  theme(panel.background = element_rect(fill = NA, colour = 'NA')) +
  theme(axis.line = element_line(color='gray'),
        axis.text=element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave('figures/site_climate.png', width = 8, height = 6, units ='in')

