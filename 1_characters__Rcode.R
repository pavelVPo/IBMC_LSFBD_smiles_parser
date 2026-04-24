library(tidyverse)

# Input
symbs_raw <- read_tsv("C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/data_v3/symbols.tsv") |>
					  rename(symb_class_id = id) |> 
					  separate_longer_delim(symbols, delim = ", ")
# Generate actual symbols from symbol folds
symb_folds <- symbs_raw |> filter(class != "bracket" & str_detect(symbols, "\\["))
for (i in seq(1:nrow(symb_folds))) {
	if (symb_folds[i, 4] |> pull() == "([-=#$:.]") {
		symbs = c("(-", "(=", "(#", "($", "(:", "(.") |> str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[-=#$:.][0-9]") {
		symbs = expand.grid(c("-", "=", "#", "$", ":", "."), c(0, 1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "%[0-9][1-9]") {
		symbs = expand.grid(c("%"), c(0,1,2,3,4,5,6,7,8,9), c(1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "%[1-9][0-9]") {
		symbs = expand.grid(c("%"), c(1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[-=#$:.]%[0-9][1-9]") {
		symbs = expand.grid(c("-", "=", "#", "$", ":", "."), c("%"), c(0,1,2,3,4,5,6,7,8,9), c(1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[-=#$:.]%[1-9][0-9]") {
		symbs = expand.grid(c("-", "=", "#", "$", ":", "."), c("%"), c(1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == ")[-=#$:.]") {
		symbs = c(")-", ")=", ")#", ")$", "):", ").") |> str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[0-9][1-9]") {
		symbs = expand.grid(c(0,1,2,3,4,5,6,7,8,9), c(1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[1-9][0-9]") {
		symbs = expand.grid(c(1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[0-9][0-9][1-9]") {
		symbs = expand.grid(c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9), c(1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[0-9][1-9][0-9]") {
		symbs = expand.grid(c(0,1,2,3,4,5,6,7,8,9), c(1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[1-9][0-9][0-9]") {
		symbs = expand.grid(c(1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@][@]") {
		symb_folds[i, 4] <- "@@"
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@]TH[1-2]") {
		symbs = expand.grid(c("@"), c("TH"), c(1,2)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@]AL[1-2]") {
		symbs = expand.grid(c("@"), c("AL"), c(1,2)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@]SP[1-3]") {
		symbs = expand.grid(c("@"), c("SP"), c(1,2,3)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@]TB[1-20]") {
		symbs = expand.grid(c("@"), c("TB"), c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[@]OH[1-30]") {
		symbs = expand.grid(c("@"), c("OH"), c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "H[2-9]") {
		symbs = expand.grid(c("H"), c(2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[+-]") {
		symb_folds[i, 4] <- "+, -"
		next
	}
	if (symb_folds[i, 4] |> pull() == "[+][+]") {
		symb_folds[i, 4] <- "++"
		next
	}
	if (symb_folds[i, 4] |> pull() == "[-][-]") {
		symb_folds[i, 4] <- "--"
		next
	}
	if (symb_folds[i, 4] |> pull() == "[+-][1-9]") {
		symbs = expand.grid(c("+", "-"), c(1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == "[+-]1[0-5]") {
		symbs = expand.grid(c("+", "-"), c("1"), c(0,1,2,3,4,5)) |>
						mutate(v = str_c(Var1, Var2, Var3, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == ":[0-9]") {
		symbs = expand.grid(c(":"), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == ":[0-9][0-9]") {
		symbs = expand.grid(c(":"), c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
	if (symb_folds[i, 4] |> pull() == ":[0-9][0-9][0-9]") {
		symbs = expand.grid(c(":"), c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,6,7,8,9)) |>
						mutate(v = str_c(Var1, Var2, Var3, Var4, collapse = ", ")) |>
						pull() |>
						unique() |>
						str_c(collapse = ", ")
		symb_folds[i, 4] <- symbs
		next
	}
}

symbs_all <- symbs_raw |> filter(class == "bracket" | !str_detect(symbols, "\\[")) |>
							bind_rows(symb_folds |> separate_longer_delim(symbols, delim = ", ")) |>
							group_by(class) |>
							mutate(symbols = symbols |> sort() |> unique() |> str_c(collapse = ", ")) |>
							slice_head(n=1) |>
							ungroup()

# Output
write_tsv(symbs_all, "C:/Users/XPS/Documents/!__tools_for_chemical_information_processing/SMILES_parser/data_v3/symbols_&_characters.tsv", quote = "all")