library(tidyverse)

## Input
pairs <- read_tsv("C:/.../pairs_processed.tsv") |>
			distinct()
symbols <- read_tsv("C:/.../symbols_&_info.tsv")

## Apply rules
# atom_oar
symb__atom_oar <- symbols |> filter(class == "atom_oar") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__atom_oar <- pairs |> inner_join(symb__atom_oar, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__atom_oar, by = c("found_target" = "symbols"))) |>
					distinct()

## Check rules
# atom_oar
pct_symb_fail__atom_oar   <- 100 * nrow(data__atom_oar |> filter(result == 0)) / nrow(data__atom_oar)
class_symb_fail__atom_oar <- data__atom_oar |> filter(result == 0) |> pull(target_class) |> unique()