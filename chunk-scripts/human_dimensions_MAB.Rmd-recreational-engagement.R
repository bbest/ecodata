
com<-ecodata::engagement %>% 
  dplyr::filter(Region == "Mid-Atlantic", 
                Fishery == "Recreational")

com2<-com %>% 
  ggplot2::ggplot()+
  ggplot2::geom_point(aes(x = Eng, y = Rel, color = Rating), size = 2)+
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black")+
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(aes(x = Eng, #geom_text_repel auto-jitters text around points
                      y = Rel,
                      label = Community,
                      color = Rating), show.legend = FALSE, direction = "both", box.padding = 0.3, size = 3)+
  ggplot2::scale_color_brewer(palette = "Dark2", #Change legend labels for clarity
                     breaks = com$Rating) +
  xlim(-1,12)+
  ylim(-1,16)+
  theme(legend.position=c(0.8, 0.85), 
        legend.title = element_blank(),       
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ggplot2::xlab("Recreation Engagament Score") +
  ggplot2::ylab("Recreation Reliance Score") +
  ggplot2::ggtitle("Social Vulnerability in Top Recreational Fishing Communities")+
  #ggplot2::guides(color = FALSE) +
  ecodata::theme_ts()
  
  
  gridExtra::grid.arrange(com2, bottom = textGrob("Low <--------------------------------------------------------------------------------------------------------------------------------------> High", 
                                     x = 0.5, y = 1, gp = gpar(fontsize = 7)),
                          left = textGrob("Low <----------------------------------------------------------------------------------------> High", rot = 90,
                                   x = 1, y = 0.5, gp = gpar(fontsize = 7)))
