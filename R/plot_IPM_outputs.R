#' Plots outputs variables
#' 
#' Generates a panel with plotted values of IPM outputs.
#' 
#' @param output_data Dataset containing output variables for the different phases. Must have long format and following variables: \describe{ 
#' \item{phase}{used as x value, defines which phase the output variable belongs to}
#' \item{value}{value of the output variable, used as y value in the plot}
#' \item{output}{defines which output variable 'value' refers to}} 
#' @param which_to_show Character vector containing a subset of \code{unique(output)}. Defines which output values are plotted.
#' @param y_spacing Numeric value defining what percentage of white space is placed around the limits of the y axis. Based on the difference between minimun and maximum value in a plot.
#' @param n_cols Defines how many columns constitute the panel layout.
#' @export

plot_IPM_outputs <- function(output_data, which_to_show, y_spacing = 10, n_cols = 4){
  
  plops <- list()
  
  for (i in 1:length(which_to_show)){
    min_lim <- min(output_data %>% 
                     filter(output == which_to_show[i]) %>%
                     pull(value)) - abs(diff(range(output_data %>% 
                                                     filter(output == which_to_show[i]) %>%
                                                     pull(value))))/y_spacing
    max_lim <- max(output_data %>% 
                     filter(output == which_to_show[i]) %>%
                     pull(value)) + abs(diff(range(output_data %>% 
                                                     filter(output == which_to_show[i]) %>%
                                                     pull(value))))/y_spacing
    
    plops[[i]] <- assign(paste("p", i, sep = ""), 
                         ggplot(output_data %>% filter(output == which_to_show[i]), aes(x = phase)) +
                           geom_point(aes(y = value, col = phase)) +
                           scale_color_manual(values = c(Cols)) + 
                           theme(legend.position = "none") +
                           ylim(min_lim, max_lim) + ylab("") + xlab("") +
                           ggtitle(which_to_show[i]))
  }
  gridExtra::grid.arrange(grobs = plops, ncol = n_cols)
  
}