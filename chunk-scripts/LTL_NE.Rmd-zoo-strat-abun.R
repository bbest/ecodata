
zoo_abund <- ecodata::zoo_strat_abun %>% 
  dplyr::mutate(Value = log10(Value+1)) %>% 
  dplyr::filter(EPU %in% c("GOM", "GB"),
         !stringr::str_detect(Var, "Small|Large")) %>% 
  dplyr::group_by(Var, EPU) %>% 
  dplyr::mutate(hline = mean(Value)) 

gom_zoo<-zoo_abund %>% 
  dplyr::filter(EPU == "GOM") %>% 
  ggplot2::ggplot() +
  ggplot2::annotate("rect", fill = shade.fill, alpha = shade.alpha,
      xmin = x.shade.min , xmax = x.shade.max,
      ymin = -Inf, ymax = Inf) +
  ecodata::geom_gls(aes(x = Time, y = Value, group = Var)) +
  ggplot2::geom_line(aes(x = Time, y = Value)) +
  ggplot2::geom_point(aes(x = Time, y = Value)) +
  ggplot2::ylab(expression("Log Stratified Abundance")) +
  ggplot2::xlab(element_blank())+
  ggtitle("GOM Zooplankton abundance") +
  ggplot2::facet_wrap(Var~., ncol = 3) +
  ggplot2::scale_x_continuous(expand = c(0.01, 0.01))+
  ggplot2::geom_hline(aes(yintercept = hline),
           size = hline.size,
           alpha = hline.alpha,
           linetype = hline.lty)+
  ecodata::theme_facet() +
  ggplot2::theme(strip.text=element_text(hjust=0,
                                face = "italic"))
  
gb_zoo<-zoo_abund %>% 
  dplyr::filter(EPU == "GB") %>% 
  ggplot2::ggplot(aes(x = Time, y = Value)) +
  ggplot2::annotate("rect", fill = shade.fill, alpha = shade.alpha,
      xmin = x.shade.min , xmax = x.shade.max,
      ymin = -Inf, ymax = Inf) +
  ecodata::geom_gls() +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::ylab(expression("Log Stratified Abundance")) +
  ggplot2::xlab(element_blank())+
  ggplot2::ggtitle("GB Zooplankton abundance") +
  ggplot2::facet_wrap(Var~., ncol = 3) +
  ggplot2::scale_x_continuous(expand = c(0.01, 0.01))+
  ggplot2::geom_hline(aes(yintercept = hline),
           size = hline.size,
           alpha = hline.alpha,
           linetype = hline.lty)+
  ecodata::theme_facet() +
  ggplot2::theme(strip.text=element_text(hjust=0,
                                face = "italic"))
  
cowplot::plot_grid(gb_zoo, gom_zoo, nrow = 2)
