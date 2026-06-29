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
# Plot
hmap <- ggplot(data, aes(x_, y_, fill = intersection)) +
  			geom_tile() +
  			scale_fill_gradient2(low = "white", mid = "#F1828D", high = "#765D69", midpoint = .5) +
  			scale_x_continuous(breaks = seq(0, 2000, by = 50),  labels = x_labels) +
  			scale_y_continuous(breaks = seq(-2000, 0, by = 50), labels = y_labels) +
  			theme_void() +
  			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  			theme(axis.text.y = element_text(size = 8)) +
  			theme(legend.position="none") +
  			coord_fixed()
hmap

## Export this draft
ggsave("C:/.../hmap.png",
			scale = 1, width = 7, height = 7, units = "in", dpi = 300, bg = "white")