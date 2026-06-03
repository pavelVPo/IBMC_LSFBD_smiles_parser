library(tidyverse)

## Input
symbols <- read_tsv("C:/.../symbols_&_info.tsv")

## Prepare the data for the test
data <- symbols |> select(symbols) |>
						separate_longer_delim(symbols, delim = ", ") |>
						distinct() |>
						cross_join(symbols |> select(symbols, class) |>
											  separate_longer_delim(symbols, delim = ", ") |>
											  distinct() |>
											  group_by(symbols) |>
											  mutate(class = str_c(class |> sort(), collapse = ", "))) |>
						rename(left = symbols.x, target = symbols.y, target_class = class) |>
						rowwise() |>
						mutate(pair = str_c(left, target, collapse = "")) |>
						ungroup() |>
						select(left, target, pair, target_class)

## Export the data to test
write_tsv(data, "C:/.../pairs_test.tsv")
