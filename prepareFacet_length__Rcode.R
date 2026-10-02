## Input
symbols <- read.csv2(".../SMILES_parser/data_v4/symbols.tsv", sep = "\t")[,4] |>
						paste0(collapse = ", ") |>
						strsplit(", ") |>
						unlist() |>
						unique()

## Processing
# Get symbols of length 1
symbol_one   <- symbols[ifelse(nchar(symbols) == 1, TRUE, FALSE)] |> sort()
# Get symbols of length 2
symbol_two   <- symbols[ifelse(nchar(symbols) == 2, TRUE, FALSE)] |> sort()
# Get symbols of length 3
symbol_three <- symbols[ifelse(nchar(symbols) == 3, TRUE, FALSE)] |> sort()
# Get symbols of length 4
symbol_four  <- symbols[ifelse(nchar(symbols) == 4, TRUE, FALSE)] |> sort()
# Get symbols of length 5
symbol_five  <- symbols[ifelse(nchar(symbols) == 5, TRUE, FALSE)] |> sort()

## Export
file_con         <- file(".../SMILES_parser/data_v4/facetLength.txt", "w")
symbol_one_str   <- paste0('["', paste0(symbol_one, collapse = '", "'), '"]')
symbol_two_str   <- paste0('["', paste0(symbol_two, collapse = '", "'), '"]')
symbol_three_str <- paste0('["', paste0(symbol_three, collapse = '", "'), '"]')
symbol_four_str  <- paste0('["', paste0(symbol_four, collapse = '", "'), '"]')
symbol_five_str  <- paste0('["', paste0(symbol_five, collapse = '", "'), '"]')
rslt_str <- paste0(	"1\n",
				    length(symbol_one),
				    "\n",
				    symbol_one_str,
				    "\n\n2\n",
				    length(symbol_two),
				    "\n",
				    symbol_two_str,
				    "\n\n3\n",
				    length(symbol_three),
				    "\n",
				    symbol_three_str,
				    "\n\n4\n",
				    length(symbol_four),
				    "\n",
				    symbol_four_str,
				    "\n\n5\n",
				    length(symbol_five),
				    "\n",
				    symbol_five_str)
write(rslt_str, file_con)