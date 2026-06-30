library(tidyverse)

## Input
symbols <- read_tsv("C:/.../symbols_&_info__upd.tsv") |>
					mutate(symbols = str_trim(symbols)) |>
					arrange(class) |>
					arrange(type)

## Assess the intersection
# The idea is to create a heatmap using ggplot geom_tile(), which produces rectangle parameterised by the center of the rect and its size, https://ggplot2-book.org/individual-geoms.html
# So, it is possible to calculate rectangle for each pair of classes and color it according to the intersection (full, partial, zero).
# Additional rows and columns of rectangles will hold the labels for types and classes.
data <- symbols |> select(type, class, symbols) |>
				   cross_join(symbols |> select(type, class, symbols)) |>
				   rename(type_row = type.x, class_row = class.x, symbols_row = symbols.x, type_col = type.y, class_col = class.y, symbols_col = symbols.y) |>
				   select(type_row, type_col, class_row, class_col, symbols_row, symbols_col) |>
				   mutate(intersection = if_else(symbols_row == symbols_col, 1, 0)) |>
				   # check for partial intersection
				   rowwise() |>
				   mutate(intersection = if_else(intersection == 0 & any( (symbols_row |> str_split(pattern = ", ") |> unlist()) %in% (symbols_col |> str_split(pattern = ", ") |> unlist())), .5, intersection)) |>
				   ungroup()
# Calculate rectangles in a simple nested loop
data$x_  <- NA
data$y_ <- NA
first_col = data[1, 4] |> pull()
y_ <- 0
x_ <- 0
for (row in seq(1:nrow(data))) {
	# Check if next row is needed
	if (row > 1 & first_col == data[row, 4] |> pull()) {
		y_ <- y_ - 50
		x_ <- 0
	}
	data[row,8] <- x_
	data[row,9] <- y_ 
	# Move right
	x_ <- x_ + 50
}
x_labels <- data |> pull(class_col) |> unique()
y_labels <- data |> pull(class_row) |> unique() |> rev()
# Get labels for rows and columns
class_row <- data |> select(class_row, y_) |>
						distinct() |>
						rename(y_class = y_) |>
						mutate(x_class = -330, intersection = 0)
type_row <- data |> select(type_row, y_) |>
						distinct() |>
						rename(y_type = y_) |>
						mutate(x_type = -930, intersection = 0) |>
						group_by(type_row) |>
						mutate(h_type = abs(min(y_type) - max(y_type)) + 50) |>
						mutate(y_type = mean(y_type)) |>
						slice_head(n = 1) |>
						ungroup()
class_col <- data |> select(class_col, x_) |>
						distinct() |>
						rename(x_class = x_) |>
						mutate(y_class = -2330, intersection = 0)
type_col <- data |> select(type_col, x_) |>
						distinct() |>
						rename(x_type = x_) |>
						mutate(y_type = -2930, intersection = 0) |>
						group_by(type_col) |>
						mutate(w_type = max(x_type) - min(x_type) + 50) |>
						mutate(x_type = mean(x_type)) |>
						slice_head(n = 1) |>
						ungroup()
# Get hlines and vlines
hline_class <- class_row |> select(y_class) |>
						mutate(y2_class = y_class - 25) |>
						mutate(y_class =  y_class + 25) |>
						pivot_longer(c(y_class, y2_class), names_to = "name", values_to = "y_") |>
						select(y_) |>
						arrange()  |>
						distinct()
hline_type <- type_row |> select(y_type, h_type) |>
						mutate(y2_type = y_type - h_type/2) |>
						mutate(y_type =  y_type + h_type/2) |>
						pivot_longer(c(y_type, y2_type), names_to = "name", values_to = "y_") |>
						select(y_) |>
						arrange()  |>
						distinct()
vline_type <- type_col |> select(x_type, w_type) |>
						mutate(x2_type = x_type - w_type/2) |>
						mutate(x_type =  x_type + w_type/2) |>
						pivot_longer(c(x_type, x2_type), names_to = "name", values_to = "x_") |>
						select(x_) |>
						arrange()  |>
						distinct()
vline_class <- class_col |> select(x_class) |>
						mutate(x2_class = x_class - 25) |>
						mutate(x_class =  x_class + 25) |>
						pivot_longer(c(x_class, x2_class), names_to = "name", values_to = "x_") |>
						select(x_) |>
						arrange()  |>
						distinct()
# Plot
hmap <- ggplot() +
  			geom_tile(data = data, aes(x_, y_,  fill = intersection)) +
  			scale_fill_gradient2(low = "white", mid = "#F1828D", high = "#765D69", midpoint = .5) +
  			scale_x_continuous(breaks = seq(0, 2000, by = 50),   labels = x_labels) +
  			scale_y_continuous(breaks = seq(-2000, 0, by = 50),  labels = y_labels) +
  			# Add symbol classes and types for rows
  			geom_tile(data = class_row, aes(x_class, y_class, width = 600), color = "black", fill = "white", linewidth = .03) +
  			geom_tile(data = type_row,  aes(x_type, y_type, width = 600, height = h_type),   color = "black", fill = "white", linewidth = .25) +
  			# Add text labels for rows
  			geom_text(data = class_row, aes(label = class_row, x = x_class, y = y_class), size = 2.5) +
  			geom_text(data = type_row,  aes(label = type_row,  x = x_type,  y = y_type),  size = 2.5) +
  			# Add symbol classes and types for cols
  			geom_tile(data = class_col, aes(x_class, y_class, width = 50, height = 600), color = "black", fill = "white", linewidth = .03) +
  			geom_tile(data = type_col,  aes(x_type, y_type, height = 600, width = w_type),   color = "black", fill = "white", linewidth = .25) +
  			# Add text labels for cols
  			geom_text(data = class_col, aes(label = class_col, x = x_class, y = y_class), angle = 90, size = 2.5) +
  			geom_text(data = type_col,  aes(label = type_col, x = x_type, y = y_type), angle = 90, size = 2.5) +
  			# Add hlines for classes and types
  			geom_segment(data = hline_class, aes(x = -25, y = y_, xend = 2025, yend = y_), linewidth = .03) +
  			geom_segment(data = hline_type,  aes(x = -930, y = y_, xend = 2025, yend = y_), linewidth = .25) +
  			# Add vlines for classes and types
  			geom_segment(data = vline_type,  aes(x = x_, y = 25, xend = x_, yend = -2630), linewidth = .25) +
  			geom_segment(data = vline_class, aes(x = x_, y = 25, xend = x_, yend = -2630), linewidth = .03) +
  			theme_void() +
  			theme(legend.position="none") +
  			coord_fixed()
hmap

## Export
ggsave("C:/.../hmap.png",
			scale = 1, width = 7, height = 7, units = "in", dpi = 300, bg = "white")