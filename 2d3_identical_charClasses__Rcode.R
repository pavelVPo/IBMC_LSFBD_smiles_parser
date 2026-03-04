library(tidyverse)

## Import the data
data_raw <- read_tsv("C:/.../chars_&_symbs.tsv")

## Templates to generate dot-code to be interpreted with the cmd version of Graphviz (graphviz-14.1.2, https://graphviz.org/download/) and adjusted using Inkscape (https://inkscape.org/)
initiator_str  <- "digraph G {\nlayout = sfdp;\nsplines = ortho;\n\n"
terminator_str <- "\n\n}"

## Summarize the data
data_summ <- data_raw |> group_by(chars) |>
					 	 summarise(n_classes = n())

## Prepare the data
# data
data <- data_raw |> group_by(chars) |>
					mutate(charClass_combined = str_c(charClass, sep = ", "),
						   charClass_id_combined = str_c(charClass_id, sep = ", "),
						   charClass_description_combined = str_c(charClass_description),
						   chars_id = cur_group_id()) |>
					ungroup() |>
					distinct() |>
					select(charClass, charClass_id, chars, chars_id, symbClass_id, symbClass, symbType, charClass_description, charClass_combined, charClass_id_combined, charClass_description_combined)
# schema
data_schema <- data |> select(chars, chars_id, charClass, symbClass, symbType) |> distinct()
#|> mutate_all(~str_replace_all(., "-", "")) |> mutate_all(~str_replace_all(., "/", "")) |> mutate_all(~str_replace_all(., " ", "_"))

## Prepare the data to be visualized using D3.JS
links_all_raw <- tibble(source = rep(NA, 1000), target = rep(NA, 1000), value = rep(NA, 1000),
					start_type = rep(NA, 1000), end_type = rep(NA, 1000),
					source_color = rep(NA, 1000), target_color = rep(NA, 1000), link_color = rep(NA, 1000))
row_n <- 1
# Base nodes
base_nodes <- data_schema |> pull(symbType) |> unique()
for (i in seq(1:length(base_nodes))) {
	links <- data_schema |> filter(symbType == base_nodes[i]) |> pull(symbClass) |> unique()
	for (k in seq(1:length(links))) {
		links_all_raw[row_n, 1] <- base_nodes[i]
		links_all_raw[row_n, 2] <- links[k]
		links_all_raw[row_n, 3] <- 1
		links_all_raw[row_n, 4] <- "Symbol type"
		links_all_raw[row_n, 5] <- "Symbol class"
		links_all_raw[row_n, 6] <- "#475C7A"
		links_all_raw[row_n, 7] <- "#685D79"
		row_n <- row_n + 1
	}
}
# L0
l0_nodes <- data_schema |> pull(symbClass) |> unique()
for (i in seq(1:length(l0_nodes))) {
	links <- data_schema |> filter(symbClass == l0_nodes[i]) |> pull(charClass) |> unique()
	for (k in seq(1:length(links))) {
		links_all_raw[row_n, 1] <- l0_nodes[i]
		links_all_raw[row_n, 2] <- links[k]
		links_all_raw[row_n, 3] <- 1
		links_all_raw[row_n, 4] <- "Symbol class"
		links_all_raw[row_n, 5] <- "Character class"
		links_all_raw[row_n, 6] <- "#685D79"
		links_all_raw[row_n, 7] <- "#D8737F"
		row_n <- row_n + 1
	}
}
# L1
l1_nodes <- data_schema |> pull(charClass) |> unique()
for (i in seq(1:length(l1_nodes))) {
	links <- data_schema |> filter(charClass == l1_nodes[i]) |> pull(chars) |> unique()
	for (k in seq(1:length(links))) {
		links_all_raw[row_n, 1] <- l1_nodes[i]
		links_all_raw[row_n, 2] <- links[k]
		links_all_raw[row_n, 3] <- 1
		links_all_raw[row_n, 4] <- "Character class"
		links_all_raw[row_n, 5] <- "Chars"
		links_all_raw[row_n, 6] <- "#D8737F"
		links_all_raw[row_n, 7] <- "#FCBB6D"
		row_n <- row_n + 1
	}
}

links <- links_all_raw |> filter(!is.na(source)) |> group_by(target) |>
							mutate(value = 100*n()) |>
							ungroup() |>
							mutate(link_color = if_else(value > 100, "#8B3834", "#161212FF"))

## Export
write_csv(links, "C:/.../smiles.csv")