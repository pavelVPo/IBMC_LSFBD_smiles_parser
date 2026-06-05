library(tidyverse)

## Input
pairs <- read_tsv("C:/.../pairs_processed.tsv") |>
			distinct()
symbols <- read_tsv("C:/.../symbols_&_info.tsv")

## Apply rules
# atom_oar / atom_bar
symb__atom_obar <- symbols |> filter(class == "atom_oar" | class == "atom_bar") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__atom_obar <- pairs |> inner_join(symb__atom_obar, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__atom_obar, by = c("found_target" = "symbols"))) |>
					distinct()
# atom_oal / atom_bal
symb__atom_obal <- symbols |> filter(class == "atom_oal" | class == "atom_bal") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__atom_obal <- pairs |> inner_join(symb__atom_obal, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__atom_obal, by = c("found_target" = "symbols"))) |>
					distinct()
# atom_oal_2 / atom_bal_2
symb__atom_obal_2 <- symbols |> filter(class == "atom_oal_2" | class == "atom_bal_2") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__atom_obal_2 <- pairs |> inner_join(symb__atom_obal_2, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__atom_obal_2, by = c("found_target" = "symbols"))) |>
					distinct()
# anything
symb__any <- symbols |> filter(class == "anything" | class == "anything") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__any <- pairs |> inner_join(symb__any, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__any, by = c("found_target" = "symbols"))) |>
					distinct()
# bracket
symb__br <- symbols |> filter(class == "s_bracket" | class == "l_bracket") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__br <- pairs |> inner_join(symb__br, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__br, by = c("found_target" = "symbols"))) |>
					distinct()
# single bond
symb__sb <- symbols |> filter(class == "single_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__sb <- pairs |> inner_join(symb__sb, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__sb, by = c("found_target" = "symbols"))) |>
					distinct()
# double bond
symb__db <- symbols |> filter(class == "double_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__db <- pairs |> inner_join(symb__db, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__db, by = c("found_target" = "symbols"))) |>
					distinct()
# triple bond
symb__tb <- symbols |> filter(class == "triple_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__tb <- pairs |> inner_join(symb__tb, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__tb, by = c("found_target" = "symbols"))) |>
					distinct()
# quadruple bond
symb__qb <- symbols |> filter(class == "quadruple_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__qb <- pairs |> inner_join(symb__qb, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__qb, by = c("found_target" = "symbols"))) |>
					distinct()
# aromatic bond
symb__ab <- symbols |> filter(class == "aromatic_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__ab <- pairs |> inner_join(symb__ab, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__ab, by = c("found_target" = "symbols"))) |>
					distinct()
# aromatic bond
symb__ab <- symbols |> filter(class == "aromatic_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__ab <- pairs |> inner_join(symb__ab, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__ab, by = c("found_target" = "symbols"))) |>
					distinct()
# no bond
symb__nb <- symbols |> filter(class == "no_bond") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__nb <- pairs |> inner_join(symb__nb, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__nb, by = c("found_target" = "symbols"))) |>
					distinct()
# branch
symb__branch <- symbols |> filter(class == "bm_ibi" | class == "bm_tbi") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__branch <- pairs |> inner_join(symb__branch, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__branch, by = c("found_target" = "symbols"))) |>
					distinct()
# ring
symb__ring <- symbols |> filter(class == "bm_iri" | class == "bm_tri") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__ring <- pairs |> inner_join(symb__ring, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__ring, by = c("found_target" = "symbols"))) |>
					distinct()
# branching with the explicit bonds
symb__tibe <- symbols |> filter(class == "bm_ibe" | class == "bm_tbe_2") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__tibe <- pairs |> inner_join(symb__tibe, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__tibe, by = c("found_target" = "symbols"))) |>
					distinct()
# branching with the explicit bonds
symb__tirei <- symbols |> filter(class == "bm_ire_2" | class == "bm_ire_4" | class == "bm_iri" |
								 class == "bm_iri_3" | class == "bm_tre_2" | class == "bm_tre_4" |
								 class == "bm_tri" | class == "bm_tri_3") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__tirei <- pairs |> inner_join(symb__tirei, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__tirei, by = c("found_target" = "symbols"))) |>
					distinct()
# cis \ trans
symb__ct <- symbols |> filter(class == "l_ct" | class == "r_ct") |>
					select(symbols) |>
					separate_longer_delim(symbols, delim = ", ") |>
					distinct()
data__ct <- pairs |> inner_join(symb__ct, by = c("target" = "symbols")) |>
					bind_rows(pairs |> inner_join(symb__ct, by = c("found_target" = "symbols"))) |>
					distinct()

## Check rules
# atom_oar / atom_bar
pct_symb_fail__atom_obar   <- 100 * nrow(data__atom_obar |> filter(result == 0)) / nrow(data__atom_obar)
class_symb_fail__atom_obar <- data__atom_obar |> filter(result == 0) |> pull(target_class) |> unique()
# atom_oal / atom_bal
pct_symb_fail__atom_obal   <- 100 * nrow(data__atom_obal |> filter(result == 0)) / nrow(data__atom_obal)
class_symb_fail__atom_obal <- data__atom_obal |> filter(result == 0) |> pull(target_class) |> unique()
# atom_oal_2 / atom_bal_2
pct_symb_fail__any   <- 100 * nrow(data__any |> filter(result == 0)) / nrow(data__any)
class_symb_fail__any <- data__any |> filter(result == 0) |> pull(target_class) |> unique()
# bracket
pct_symb_fail__br   <- 100 * nrow(data__br |> filter(result == 0)) / nrow(data__br)
class_symb_fail__br <- data__br |> filter(result == 0) |> pull(target_class) |> unique()
# single bond
pct_symb_fail__sb   <- 100 * nrow(data__sb |> filter(result == 0)) / nrow(data__sb)
class_symb_fail__sb <- data__sb |> filter(result == 0) |> pull(target_class) |> unique()
# double bond
pct_symb_fail__db   <- 100 * nrow(data__db |> filter(result == 0)) / nrow(data__db)
class_symb_fail__db <- data__db |> filter(result == 0) |> pull(target_class) |> unique()
# triple bond
pct_symb_fail__tb   <- 100 * nrow(data__tb |> filter(result == 0)) / nrow(data__tb)
class_symb_fail__tb <- data__tb |> filter(result == 0) |> pull(target_class) |> unique()
# quadruple bond
pct_symb_fail__qb   <- 100 * nrow(data__qb |> filter(result == 0)) / nrow(data__qb)
class_symb_fail__qb <- data__qb |> filter(result == 0) |> pull(target_class) |> unique()
# aromatic bond
pct_symb_fail__ab   <- 100 * nrow(data__ab |> filter(result == 0)) / nrow(data__ab)
class_symb_fail__ab <- data__ab |> filter(result == 0) |> pull(target_class) |> unique()
# no bond
pct_symb_fail__nb   <- 100 * nrow(data__nb |> filter(result == 0)) / nrow(data__nb)
class_symb_fail__nb <- data__nb |> filter(result == 0) |> pull(target_class) |> unique()
# branch
pct_symb_fail__branch   <- 100 * nrow(data__branch |> filter(result == 0)) / nrow(data__branch)
class_symb_fail__branch <- data__branch |> filter(result == 0) |> pull(target_class) |> unique()
# ring
pct_symb_fail__ring   <- 100 * nrow(data__ring |> filter(result == 0)) / nrow(data__ring)
class_symb_fail__ring <- data__ring |> filter(result == 0) |> pull(target_class) |> unique()
# tibe
pct_symb_fail__tibe   <- 100 * nrow(data__tibe |> filter(result == 0)) / nrow(data__tibe)
class_symb_fail__tibe <- data__tibe |> filter(result == 0) |> pull(target_class) |> unique()
# tirei
pct_symb_fail__tirei   <- 100 * nrow(data__tirei |> filter(result == 0)) / nrow(data__tirei)
class_symb_fail__tirei <- data__tirei |> filter(result == 0) |>
							select(found_target, target_class) |>
							distinct() |>
							inner_join(symbols_ext, by = c("found_target" = "symbols"), relationship = "many-to-many") |>
							select(target_class, class) |>
							distinct() |>
							rowwise() |>
							mutate(miss = str_c(c(target_class, class), collapse = " -> ")) |>
							pull(miss) |>
							unique()
# ct
pct_symb_fail__ct   <- 100 * nrow(data__ct |> filter(result == 0)) / nrow(data__ct)
class_symb_fail__ct <- data__ct |> filter(result == 0) |>
							select(found_target, target_class) |>
							distinct() |>
							inner_join(symbols_ext, by = c("found_target" = "symbols"), relationship = "many-to-many") |>
							select(target_class, class) |>
							distinct() |>
							rowwise() |>
							mutate(miss = str_c(c(target_class, class), collapse = " -> ")) |>
							pull(miss) |>
							unique()
