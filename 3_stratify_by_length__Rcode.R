library(tidyverse)

## Input
symbols_raw <- read_tsv("C:/.../symbols_&_characters_upd.tsv")

## Process
symbols <- symbols_raw |> separate_longer_delim(symbols, delim = ", ") |>
						  mutate(symbol_id = row_number(), length = str_length(symbols))

## Stratify by the length and export symbols
for (i in seq(1:5)) {
	symbols_strata <- symbols |> filter(length == i) |> select(symbols) |> distinct()
	write_tsv(symbols_strata, str_glue("C:/.../symbols_strata_{i}.tsv"))
	symbols_strata_str <- symbols_strata |> pull(symbols) |> unique() |> str_c(collapse = "\", \"")
	symbols_strata_str <- str_c(c("[\"", symbols_strata_str, "\"]"), collapse = "")
	write_lines(symbols_strata_str, str_glue("C:/.../symbols_strata_{i}.txt"), sep = "")
}

## Re-arrange and export the whole dataset
symbols_out <- symbols |> rename(class_id = symb_class_id) |> select(type, class_id, class, symbol_id, symbols, description) |>
						  rename(symbols_id = symbol_id) |>
						  group_by(class_id) |>
						  mutate(symbols_id = str_c(symbols_id, collapse = ", "), symbols = str_c(symbols, collapse = ", ")) |>
						  ungroup()  |>
						  distinct() |>
						  write_tsv("C:/.../symbols_&_info__upd.tsv")