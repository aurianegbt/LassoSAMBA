## from  https://benjaminlouis-stat.fr/en/blog/2020-05-21-astuces-ggplot-rmarkdown/
'%+replace%' <- ggplot2::'%+replace%'

theme_riri <- function(base_size = 14){
  ggplot2::theme_bw(base_size = base_size) %+replace%
    ggplot2::theme(
      # L'ensemble de la figure
      plot.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold", margin = ggplot2::margin(0,0,5,0), hjust = 0),
      # Zone où se situe le graphique
      #panel.background = element_rect(fill='transparent'),
      plot.background = ggplot2::element_rect(fill='transparent', color=NA), 
      panel.grid.minor = ggplot2::element_line(linewidth = 0.25,
                                      linetype = 1),
      panel.border = ggplot2::element_blank(),
      # Les axes
      axis.title = ggplot2::element_text(size = ggplot2::rel(0.85), face = "bold",color="#5e5e5e"),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.70), face = "bold",color="#5e5e5e"),
      axis.line = ggplot2::element_line(color = "#5e5e5e", arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "lines"), type = "closed")),
      # La légende
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.85), face = "bold"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.70), face = "bold"),
      legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key.size = ggplot2::unit(1.5, "lines"),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      # Les étiquettes dans le cas d'un facetting
      strip.background = ggplot2::element_rect(fill = "#383838", color = "#383838"),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.85), face = "bold", color = "white", margin = ggplot2::margin(5,0,5,0))
    )
}

# The palette with grey:
cbPalette <- c("#8e6aa0","#ee6c4d","#007194","#a6c46a","#9D8F80","#FFD447")
  
  
pastelcbPalette <- c("#d0c1d7","#f8c2b4","#99e7ff","#e0e6c6","#cac2ba","#ffe591")
ggplot2::theme_set(theme_riri())
