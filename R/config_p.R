#' a theme mimic the base system
#' @param base_size, default 16
#' @param base_family, default Times
#' @export
theme_base <- function (base_size = 16, base_family = "Times") 
{
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(strip.background = element_rect(fill = "grey90", 
            colour = "grey90"), panel.border = element_rect(fill = "transparent", 
            colour = "black"), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, margin=margin(b = 10, unit = "pt")))
}