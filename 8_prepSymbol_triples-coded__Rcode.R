library(tidyverse)

### input
data <- read_tsv("C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/code_v2/data/chars_&_symbs.tsv") |>
					select(symbs, symbClass) |>
					# to shorten the task: change each letter to 'A|a'; each num
					mutate(symbs = str_replace_all(symbs, fixed(':[0-9],'), "?1,")) |>
					mutate(symbs = str_replace_all(symbs, fixed(':[0-9][0-9],'), "?11,")) |>
					mutate(symbs = str_replace_all(symbs, fixed(':[0-9][0-9][0-9]'), "?111")) |>
					mutate(symbs = str_replace_all(symbs, '[[:digit:]]', "1")) |>
					mutate(symbs = str_replace_all(symbs, '[[:upper:]]', "A")) |>
					mutate(symbs = str_replace_all(symbs, '[[:lower:]]', "a")) |>
					mutate(symbs = str_replace_all(symbs, '[.[^a1A,\\s]]', "?")) |>
					distinct() |>
					rowwise() |>
					mutate(symbs = symbs |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", ")) |>
					ungroup() |>
					distinct()
allowed <- read_tsv("C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/code_v2/data/theory_pairs_symbClass_rules.tsv")

### process
# get all symbols and
all_symbs <- data |> select(symbs) |>
					 	separate_longer_delim(symbs, delim = ", ") |>
						# to shorten the task: change each letter to 'A|a'; each num
						mutate(symbs = str_replace_all(symbs, fixed(':[0-9],'), "?1,")) |>
						mutate(symbs = str_replace_all(symbs, fixed(':[0-9][0-9],'), "?11,")) |>
						mutate(symbs = str_replace_all(symbs, fixed(':[0-9][0-9][0-9]'), "?111")) |>
						mutate(symbs = str_replace_all(symbs, '[[:digit:]]', "1")) |>
						mutate(symbs = str_replace_all(symbs, '[[:upper:]]', "A")) |>
						mutate(symbs = str_replace_all(symbs, '[[:lower:]]', "a")) |>
						mutate(symbs = str_replace_all(symbs, '[.[^a1A,\\s]]', "?")) |>
					 	distinct()
# unique pairs of symbols
symb_pair <- allowed |> select(left, right) |>
							distinct() |>
							rename(left_class = left, right_class = right) |>
							inner_join(data, by = c("left_class" = "symbClass")) |>
							rename(left = symbs) |>
							inner_join(data, by = c("right_class" = "symbClass")) |>
							rename(right = symbs) |>
							rowwise() |>
							mutate(classes = str_c(left_class, right_class, sep = " - ")) |>
							ungroup() |>
							select(left, right, classes) |>
							distinct() |>
							group_by(classes) |>
							separate_longer_delim(left,  delim = ", ") |>
							separate_longer_delim(right, delim = ", ") |>
							ungroup() |>
							group_by(left, right) |>
							mutate(classes = str_c(classes |> unlist() |> unique() |> sort(), collapse = ", ")) |>
							slice_head(n = 1) |>
							ungroup() |>
							distinct()
# export pairs of symbols
write_tsv(symb_pair, "C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/code_v2/data/symb_pairs_coded.tsv")

# triples of symbols
symb_triple <- symb_pair |> cross_join(all_symbs) |>
							rename(next_symb = symbs) |>
							rowwise() |>
							mutate(classes = classes |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", ")) |>
							ungroup() |>
							distinct()

# export triples of symbols
write_tsv(symb_triple, "C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/code_v2/data/symb_triples_coded.tsv")