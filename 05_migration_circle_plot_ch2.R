# Install and load the circlize and scales packages if you haven't already
# install.packages("circlize")
# install.packages("scales")
library(circlize)
# citation(circlize) # This line is for getting citation info, not part of the plot script
library(scales) # Load the scales package for the alpha() function

# Population Index -> Population Label mapping
# Match order to matrix order 0, 1, 2, ... 
pop_index_label <- c("CC", "FB", "KZ", "LB", "PP")
names(pop_index_label) <- 0:4

# Define a color for each population (Ensuring FF suffix for full opacity)
# Ensure the order of colors matches the order of population labels
# pop_color_vector <- population_colors[pop_index_label] # This variable is not used

# Define a color for each population in the correct order
population_colors_ordered <- c(
  "CC" = "#414487FF",
  "FB" = "#7AD151FF",
  "KZ" = "#2A788EFF",
  "LB" = "#22A884FF",
  "PP" = "#440154FF"
)


# Create the migration rate matrix (off-diagonal elements)

migration_rates <- matrix(
  c(
    0,	0,	0.278,	0,	0, # to CC 
    0,	0,	0.270,	0,	0, # to FB 
    0.278,	0,	0,	0,	0.250, # to KZ 
    0,	0,	0.275,	0,	0, # to LB 
    0,	0,	0.278,	0,	0 #to PP 
  ),
  nrow = 5, byrow = TRUE
)
rownames(migration_rates) <- pop_index_label[as.character(0:4)]
colnames(migration_rates) <- pop_index_label[as.character(0:4)]

# New order for plotting only
pop_plot_order <- c("CC", "PP", "KZ", "LB", "FB")

# Function to create the migration circle plot
create_migration_plot <- function(migration_matrix, population_labels, population_colors, plot_order, title = "Migration between Seal Populations") {
  
  n_pops <- length(plot_order)
  circos.clear()
  circos.par(
    gap.after = c(rep(7, n_pops - 1), 5),
    cell.padding = c(0, 0, 0, 0)
  )
  
  # Initialize sectors
  circos.initialize(factors = plot_order, xlim = c(0, 1))
  
  # Draw population segments and labels
  circos.track(track.index = 1, ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                CELL_META$sector.index, facing = "clockwise", adj = c(0, 0.5))
    circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                CELL_META$xlim[2], CELL_META$ylim[2],
                col = population_colors[CELL_META$sector.index], border = "black")
  }, bg.border = "black", track.height = mm_h(5))
  
  # Link offsets
  link_offsets <- setNames(rep(0, n_pops), plot_order)
  
  # Visual parameters
  offset_increment <- 0.2
  min_arrow_length <- 0.09
  min_line_width <- 6
  max_linewidth <- 17
  max_arrow_length <- 0.5
  
  # Prepare migration data frame
  migration_df <- expand.grid(from = plot_order, to = plot_order, stringsAsFactors = FALSE)
  migration_df <- subset(migration_df, from != to)
  migration_df$rate <- mapply(function(f, t) migration_matrix[t, f], migration_df$from, migration_df$to)
  migration_df <- migration_df[order(-migration_df$rate), ]
  
  # Draw links
  for (row in seq_len(nrow(migration_df))) {
    from_pop <- migration_df$from[row]
    to_pop   <- migration_df$to[row]
    rate     <- migration_df$rate[row]
    
    if (rate > 0) {
      line_width <- pmax(min_line_width, rate / max(migration_matrix) * max_linewidth)
      base_col <- population_colors[from_pop]
      arrow_col <- if ((from_pop == "PP" && to_pop == "KZ") || (from_pop == "CC" && to_pop == "KZ")) {
        alpha(base_col, 0.5)
      } else {
        base_col
      }
      arrow_length <- pmax(min_arrow_length, rate / max(migration_matrix) * max_arrow_length)
      
      pos_from <- 0.5 + (link_offsets[[from_pop]] * offset_increment / 2)
      pos_to   <- 0.5 + (-link_offsets[[to_pop]] * offset_increment / 2)
      
      circos.link(from_pop, pos_from,
                  to_pop, pos_to,
                  lwd = line_width,
                  col = arrow_col,
                  directional = 1,
                  arr.length = arrow_length,
                  arr.width = arrow_length * 2,
                  rou = 0.80)
      
      link_offsets[[to_pop]]   <- link_offsets[[to_pop]] + 1
      link_offsets[[from_pop]] <- link_offsets[[from_pop]] + 1
    }
  }
  
  title(title, cex.main = 1.2)
  circos.clear()
}
# Create the circle plot, passing the color vector and the new plot order
create_migration_plot(migration_rates, pop_index_label, population_colors_ordered, pop_plot_order)


