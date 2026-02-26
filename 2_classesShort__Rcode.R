library(tidyverse)

### Import
# Path
path <- "C:/.../"
# Merged classes
merged_classes_raw <- read_tsv(str_glue("{path}merged_classes.tsv")) |>
						select(-description)
# Described classes
classes <- read_tsv(str_glue("{path}class_list.tsv")) |>
						select(class, description) |>
						distinct()

### Process
# Accumulate descriptions of the merged classes
merged_classes <- merged_classes_raw |> inner_join(classes, by = c("basic_class" = "class")) |>
						group_by(merged_class) |>
						mutate(description = str_c(description, collapse = " |OR| "),
						basic_class = str_c(basic_class, collapse = ", ")) |>
						slice_head(n=1) |>
						ungroup() |>
						select(merged_class, chars, basic_class, description)

# Export the updated descriptions
write_tsv(merged_classes |> mutate(), str_glue("{path}merged_classes__autoDescribed.tsv"))