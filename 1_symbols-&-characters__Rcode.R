##### This code is the basis for the further analysis.
##### It will be written using base R.  

#### Functions
# Function to get the first character from word
char_1 <- function(word) {
    first_char <- strsplit(word, split = "")[[1]][1]
    return(first_char)
}
# Function to get the second character from word
char_2 <- function(word) {
    scnd_char <- strsplit(word, split = "")[[1]][2]
    return(scnd_char)
}
# Function to get the second and third characters from word
chars_2_3 <- function(word) {
    medium_char <- strsplit(word, split = "")[[1]][2:3]
    return(medium_char)
}
# Function to get the last character from word
char_last <- function(word) {
    last_char <- tail(strsplit(word, split = "")[[1]], n = 1)
    return(last_char)
}
# Function to ommit the first character from word
chars_1rest <- function(word) {
    rest_chars <- strsplit(word, split = "")[[1]][-1]
    return(rest_chars)
}
# Function to ommit the first and second character from word
chars_2rest <- function(word) {
    rest_chars <- strsplit(word, split = "")[[1]][-2:-1]
    return(rest_chars)
}
# Function to ommit the first, second and third character from word
chars_3rest <- function(word) {
    rest_chars <- strsplit(word, split = "")[[1]][-3:-1]
    return(rest_chars)
}


#### Prepare the structure to hold the list of symbols, character classes, descriptions
symb_data   <- data.frame(symbClass_id = rep(NA, 200), symbClass = rep(NA, 200), symbType = rep(NA, 200),
                    symbs = rep(NA, 200), symbClas_description = rep(NA, 200),
                    symbType_description = rep(NA, 200))
char_data   <- data.frame(charClass_id = rep(NA, 200), charClass = rep(NA, 200), chars = rep(NA, 200),
                    charClass_description = rep(NA, 200), symbClass_id = rep(NA, 200))
symbClass_cntr <- 1
charClass_cntr <- 1

#### Enumeration of the symbols and character classes allowed in SMILES
### TYPE: ATOM
symbType <- "atom__symbol"
symbType_description <- "Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string"
## 01, symbol class: Single character atom symbols of organic (from so called organic subset) aromatic atoms lacking the additional features
# Symbols
symbols <- c('b', 'c', 'n', 'o', 's', 'p') |> sort()
class <- "atom_oar"
symbClass_description <- "Single character atom symbols of organic (from so called organic subset) aromatic atoms lacking the additional features"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_atom_oar"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 02, symbol class: Single character atom symbols of organic aliphatic atoms lacking the additional features
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I') |> sort()
class <- "atom_oal"
symbClass_description <- "Single character atom symbols of organic aliphatic atoms lacking the additional features"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_atom_oal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 03, symbol class: Single character atom symbols of aromatic atoms enclosed within brackets
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('b', 'c', 'n', 'o', 's', 'p') |> sort()
class <- "atom_bar"
symbClass_description <- "Single character atom symbols of in-bracket aromatic atoms"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_atom_bar"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 04, symbol class: Single character atom symbols of bracket aliphatic atoms
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U') |> sort()
class <- "atom_bal"
symbClass_description <- "Single character atom symbols of in-bracket aliphatic atoms"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_atom_bal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 05, symbol class: Two character atom symbols of organic aliphatic atoms lacking the additional features
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('Cl', 'Br') |> sort()
class <- "atom_oal_2"
symbClass_description <- "Two character atom symbols of organic aliphatic atoms lacking the additional features"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_atom_oal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# Characters, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_atom_oal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 06, symbol class: Two character atom symbols of aromatic atoms, which should be enclosed within the square brackets
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('se', 'as', 'te') |> sort()
class <- "atom_bar_2"
symbClass_description <- "Two character atom symbols of bracket aromatic atoms"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_atom_bar"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# Characters, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_atom_bar"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 07, symbol class: Two character atom symbols of aliphatic atoms, which should be enclosed within the square brackets
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
                          'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                          'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                          'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
                          'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr') |> sort()
class <- "atom_bal_2"
symbClass_description <- "Two character atom symbols of in-bracket aliphatic atoms"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_atom_bal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# Characters, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_atom_bal"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

### 2 ANYTHING
symbType <- "anything__symbol"
symbType_description <- "Anything symbol is the way to designate any atom or basically anything"
## ANYTHING, SYMBOL
# 08, symbol class: Single character symbol of any atom or basically anything
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('*') |> sort()
class <- "anything"
symbClass_description <- "Single character symbol of any atom or basically anything"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_anything"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

## 3 BRACKETS, SYMBOLS
symbType <- "bracket__symbol"
symbType_description <- "Square bracket symbols is the SMILES way to mark the start of the atom record including its various properties and is the way to designate the end of the atom record"
# 09, symbol class: Single character square bracket symbols
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("[", "]") |> sort()
class <- "bracket"
symbClass_description <- "Single character symbol of square bracket"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
## BRACKETS, class
# Characters
# bracket, start
characters <- symbols[1]
charClass <- "s_bracket"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bracket, end
characters <- symbols[2]
charClass <- "e_bracket"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

### 4 BONDS
symbType <- "bond__symbol"
symbType_description <- "Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string"
## 4 BONDS, SYMBOLS
# 10, symbol class: Bond symbol corresponding to the single bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("-") |> sort()
class <- "single_bond"
symbClass_description <- "Single character symbol of single bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# single_bond
characters <- symbols
charClass <- "w_single_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 11, symbol class: Bond symbol corresponding to the double bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("=") |> sort()
class <- "double_bond"
symbClass_description <- "Single character symbol of double bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# double_bond
characters <- symbols
charClass <- "w_double_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 12, symbol class: Bond symbol corresponding to the triple bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("#") |> sort()
class <- "triple_bond"
symbClass_description <- "Single character symbol of triple bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# triple_bond
characters <- symbols
charClass <- "w_triple_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 13, symbol class: Bond symbol corresponding to the quadruple bond
symb__quadruple_bond <- c("$")
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("$") |> sort()
class <- "quadruple_bond"
symbClass_description <- "Single character symbol of quadruple bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# quadruple_bond
characters <- symbols
charClass <- "w_quadruple_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 14, symbol class: Bond symbol corresponding to the aromatic bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c(":") |> sort()
class <- "aromatic_bond"
symbClass_description <- "Single character symbol of aromatic bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# aromatic_bond
characters <- symbols
charClass <- "w_aromatic_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 15, symbol class: Bond symbol corresponding to the absence of the bond between the two specific atoms
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c(".") |> sort()
class <- "no_bond"
symbClass_description <- "Single character symbol of no bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# no_bond
characters <- symbols
charClass <- "w_no_bond"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

### 5 BOND MODIFIERS (MULTIPLIERS)
symbType <- "modifier__symbol"
symbType_description <- "Bond multiplying symbols is used in SMILES to extend the number of atoms, for which connections to the current atom could be written"
## 5 BOND MODIFIERS (MULTIPLIERS), SYMBOLS
# 16, symbol class: Single character bond multiplying symbols initiators of branching with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("(") |> sort()
class <- "bm_ibi"
symbClass_description <- "Single character bond multiplying symbols initiators of branching with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_bm_ibi
characters <- symbols
charClass <- "w_bm_ibi"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 17, symbol class: Single character bond multiplying symbols initiators of rings with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9') |> sort()
class <- "bm_iri"
symbClass_description <- "Single character bond multiplying symbols initiators of rings with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_bm_iri
characters <- symbols
charClass <- "w_bm_iri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 18, symbol class: Single character bond multiplying symbols terminators of branching with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c(')') |> sort()
class <- "bm_tbi"
symbClass_description <- "Single character bond multiplying symbols terminators of branching with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_bm_tbi
characters <- symbols
charClass <- "w_bm_tbi"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 19, symbol class: Single character bond multiplying symbols terminators of rings with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9') |> sort()
class <- "bm_tri"
symbClass_description <- "Single character bond multiplying symbols terminators of rings with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_bm_tri
characters <- symbols
charClass <- "w_bm_tri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 20, symbol class: Two-character bond multiplying symbols initiators of branching with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('(-', '(=', '(#', '($', '(:', '(.') |> sort()
class <- "bm_ibe_2"
symbClass_description <- "Two character bond multiplying symbols initiators of branching with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_ibe, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_ibe"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_ibe, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_bm_ibe"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 21, symbol class: Two-character bond multiplying symbols initiators of rings with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('-0', '=0', '#0', '$0', ':0', '.0',
'-1', '=1', '#1', '$1', ':1', '.1',
'-2', '=2', '#2', '$2', ':2', '.2',
'-3', '=3', '#3', '$3', ':3', '.3',
'-4', '=4', '#4', '$4', ':4', '.4',
'-5', '=5', '#5', '$5', ':5', '.5',
'-6', '=6', '#6', '$6', ':6', '.6',
'-7', '=7', '#7', '$7', ':7', '.7',
'-8', '=8', '#8', '$8', ':8', '.8',
'-9', '=9', '#9', '$9', ':9', '.9' ) |> sort()
class <- "bm_ire_2"
symbClass_description <- "Two character bond multiplying symbols initiators of rings with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_ire, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_ire"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_ire, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_bm_ire"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 22, symbol class: Three-character bond multiplying symbols initiators of rings with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' ) |> sort()
class <- "bm_iri_3"
symbClass_description <- "Three character bond multiplying symbols initiators of rings with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_iri, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_iri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_iri, rest
characters <- lapply(symbols, chars_1rest) |> unlist() |> unique() |> sort()
charClass <- "r_bm_iri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 23, symbol class: Four-character bond multiplying symbols initiators of rings with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
'-%02', '=%02', '#%02', '$%02', ':%02', '.%02',
'-%03', '=%03', '#%03', '$%03', ':%03', '.%03',
'-%04', '=%04', '#%04', '$%04', ':%04', '.%04',
'-%05', '=%05', '#%05', '$%05', ':%05', '.%05',
'-%06', '=%06', '#%06', '$%06', ':%06', '.%06',
'-%07', '=%07', '#%07', '$%07', ':%07', '.%07',
'-%08', '=%08', '#%08', '$%08', ':%08', '.%08',
'-%09', '=%09', '#%09', '$%09', ':%09', '.%09',
'-%10', '=%10', '#%10', '$%10', ':%10', '.%10',
'-%11', '=%11', '#%11', '$%11', ':%11', '.%11',
'-%12', '=%12', '#%12', '$%12', ':%12', '.%12',
'-%13', '=%13', '#%13', '$%13', ':%13', '.%13',
'-%14', '=%14', '#%14', '$%14', ':%14', '.%14',
'-%15', '=%15', '#%15', '$%15', ':%15', '.%15',
'-%16', '=%16', '#%16', '$%16', ':%16', '.%16',
'-%17', '=%17', '#%17', '$%17', ':%17', '.%17',
'-%18', '=%18', '#%18', '$%18', ':%18', '.%18',
'-%19', '=%19', '#%19', '$%19', ':%19', '.%19',
'-%20', '=%20', '#%20', '$%20', ':%20', '.%20',
'-%21', '=%21', '#%21', '$%21', ':%21', '.%21',
'-%22', '=%22', '#%22', '$%22', ':%22', '.%22',
'-%23', '=%23', '#%23', '$%23', ':%23', '.%23',
'-%24', '=%24', '#%24', '$%24', ':%24', '.%24',
'-%25', '=%25', '#%25', '$%25', ':%25', '.%25',
'-%26', '=%26', '#%26', '$%26', ':%26', '.%26',
'-%27', '=%27', '#%27', '$%27', ':%27', '.%27',
'-%28', '=%28', '#%28', '$%28', ':%28', '.%28',
'-%29', '=%29', '#%29', '$%29', ':%29', '.%29',
'-%30', '=%30', '#%30', '$%30', ':%30', '.%30',
'-%31', '=%31', '#%31', '$%31', ':%31', '.%31',
'-%32', '=%32', '#%32', '$%32', ':%32', '.%32',
'-%33', '=%33', '#%33', '$%33', ':%33', '.%33',
'-%34', '=%34', '#%34', '$%34', ':%34', '.%34',
'-%35', '=%35', '#%35', '$%35', ':%35', '.%35',
'-%36', '=%36', '#%36', '$%36', ':%36', '.%36',
'-%37', '=%37', '#%37', '$%37', ':%37', '.%37',
'-%38', '=%38', '#%38', '$%38', ':%38', '.%38',
'-%39', '=%39', '#%39', '$%39', ':%39', '.%39',
'-%40', '=%40', '#%40', '$%40', ':%40', '.%40',
'-%41', '=%41', '#%41', '$%41', ':%41', '.%41',
'-%42', '=%42', '#%42', '$%42', ':%42', '.%42',
'-%43', '=%43', '#%43', '$%43', ':%43', '.%43',
'-%44', '=%44', '#%44', '$%44', ':%44', '.%44',
'-%45', '=%45', '#%45', '$%45', ':%45', '.%45',
'-%46', '=%46', '#%46', '$%46', ':%46', '.%46',
'-%47', '=%47', '#%47', '$%47', ':%47', '.%47',
'-%48', '=%48', '#%48', '$%48', ':%48', '.%48',
'-%49', '=%49', '#%49', '$%49', ':%49', '.%49',
'-%50', '=%50', '#%50', '$%50', ':%50', '.%50',
'-%51', '=%51', '#%51', '$%51', ':%51', '.%51',
'-%52', '=%52', '#%52', '$%52', ':%52', '.%52',
'-%53', '=%53', '#%53', '$%53', ':%53', '.%53',
'-%54', '=%54', '#%54', '$%54', ':%54', '.%54',
'-%55', '=%55', '#%55', '$%55', ':%55', '.%55',
'-%56', '=%56', '#%56', '$%56', ':%56', '.%56',
'-%57', '=%57', '#%57', '$%57', ':%57', '.%57',
'-%58', '=%58', '#%58', '$%58', ':%58', '.%58',
'-%59', '=%59', '#%59', '$%59', ':%59', '.%59',
'-%60', '=%60', '#%60', '$%60', ':%60', '.%60',
'-%61', '=%61', '#%61', '$%61', ':%61', '.%61',
'-%62', '=%62', '#%62', '$%62', ':%62', '.%62',
'-%63', '=%63', '#%63', '$%63', ':%63', '.%63',
'-%64', '=%64', '#%64', '$%64', ':%64', '.%64',
'-%65', '=%65', '#%65', '$%65', ':%65', '.%65',
'-%66', '=%66', '#%66', '$%66', ':%66', '.%66',
'-%67', '=%67', '#%67', '$%67', ':%67', '.%67',
'-%68', '=%68', '#%68', '$%68', ':%68', '.%68',
'-%69', '=%69', '#%69', '$%69', ':%69', '.%69',
'-%70', '=%70', '#%70', '$%70', ':%70', '.%70',
'-%71', '=%71', '#%71', '$%71', ':%71', '.%71',
'-%72', '=%72', '#%72', '$%72', ':%72', '.%72',
'-%73', '=%73', '#%73', '$%73', ':%73', '.%73',
'-%74', '=%74', '#%74', '$%74', ':%74', '.%74',
'-%75', '=%75', '#%75', '$%75', ':%75', '.%75',
'-%76', '=%76', '#%76', '$%76', ':%76', '.%76',
'-%77', '=%77', '#%77', '$%77', ':%77', '.%77',
'-%78', '=%78', '#%78', '$%78', ':%78', '.%78',
'-%79', '=%79', '#%79', '$%79', ':%79', '.%79',
'-%80', '=%80', '#%80', '$%80', ':%80', '.%80',
'-%81', '=%81', '#%81', '$%81', ':%81', '.%81',
'-%82', '=%82', '#%82', '$%82', ':%82', '.%82',
'-%83', '=%83', '#%83', '$%83', ':%83', '.%83',
'-%84', '=%84', '#%84', '$%84', ':%84', '.%84',
'-%85', '=%85', '#%85', '$%85', ':%85', '.%85',
'-%86', '=%86', '#%86', '$%86', ':%86', '.%86',
'-%87', '=%87', '#%87', '$%87', ':%87', '.%87',
'-%88', '=%88', '#%88', '$%88', ':%88', '.%88',
'-%89', '=%89', '#%89', '$%89', ':%89', '.%89',
'-%90', '=%90', '#%90', '$%90', ':%90', '.%90',
'-%91', '=%91', '#%91', '$%91', ':%91', '.%91',
'-%92', '=%92', '#%92', '$%92', ':%92', '.%92',
'-%93', '=%93', '#%93', '$%93', ':%93', '.%93',
'-%94', '=%94', '#%94', '$%94', ':%94', '.%94',
'-%95', '=%95', '#%95', '$%95', ':%95', '.%95',
'-%96', '=%96', '#%96', '$%96', ':%96', '.%96',
'-%97', '=%97', '#%97', '$%97', ':%97', '.%97',
'-%98', '=%98', '#%98', '$%98', ':%98', '.%98',
'-%99', '=%99', '#%99', '$%99', ':%99', '.%99') |> sort()
class <- "bm_ire_4"
symbClass_description <- "Four character bond multiplying symbols initiators of rings with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_ire, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_ire_4"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_ire, next from start
characters <- lapply(symbols, char_2) |> unlist() |> unique() |> sort()
charClass <- "n_bm_ire"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_ire, rest
characters <- lapply(symbols, chars_2rest) |> unlist() |> unique() |> sort()
charClass <- "r_bm_ire"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 24, symbol class: Two-character bond multiplying symbols terminators of branching with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c(')-', ')=', ')#', ')$', '):', ').') |> sort()
class <- "bm_tbe_2"
symbClass_description <- "Two character bond multiplying symbols terminators of branching with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_tbe, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_tbe"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_tre, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_bm_tbe"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 25, symbol class: Two-character bond multiplying symbols terminators of rings with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('-0', '=0', '#0', '$0', ':0', '.0',
                '-1', '=1', '#1', '$1', ':1', '.1',
                '-2', '=2', '#2', '$2', ':2', '.2',
                '-3', '=3', '#3', '$3', ':3', '.3',
                '-4', '=4', '#4', '$4', ':4', '.4',
                '-5', '=5', '#5', '$5', ':5', '.5',
                '-6', '=6', '#6', '$6', ':6', '.6',
                '-7', '=7', '#7', '$7', ':7', '.7',
                '-8', '=8', '#8', '$8', ':8', '.8',
                '-9', '=9', '#9', '$9', ':9', '.9' ) |> sort()
class <- "bm_tre_2"
symbClass_description <- "Two character bond multiplying symbols terminators of rings with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_tbe, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_tre"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_tre, end
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_bm_tre"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 26, symbol class: Four-character bond multiplying symbols terminators of rings with explicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
'-%02', '=%02', '#%02', '$%02', ':%02', '.%02',
'-%03', '=%03', '#%03', '$%03', ':%03', '.%03',
'-%04', '=%04', '#%04', '$%04', ':%04', '.%04',
'-%05', '=%05', '#%05', '$%05', ':%05', '.%05',
'-%06', '=%06', '#%06', '$%06', ':%06', '.%06',
'-%07', '=%07', '#%07', '$%07', ':%07', '.%07',
'-%08', '=%08', '#%08', '$%08', ':%08', '.%08',
'-%09', '=%09', '#%09', '$%09', ':%09', '.%09',
'-%10', '=%10', '#%10', '$%10', ':%10', '.%10',
'-%11', '=%11', '#%11', '$%11', ':%11', '.%11',
'-%12', '=%12', '#%12', '$%12', ':%12', '.%12',
'-%13', '=%13', '#%13', '$%13', ':%13', '.%13',
'-%14', '=%14', '#%14', '$%14', ':%14', '.%14',
'-%15', '=%15', '#%15', '$%15', ':%15', '.%15',
'-%16', '=%16', '#%16', '$%16', ':%16', '.%16',
'-%17', '=%17', '#%17', '$%17', ':%17', '.%17',
'-%18', '=%18', '#%18', '$%18', ':%18', '.%18',
'-%19', '=%19', '#%19', '$%19', ':%19', '.%19',
'-%20', '=%20', '#%20', '$%20', ':%20', '.%20',
'-%21', '=%21', '#%21', '$%21', ':%21', '.%21',
'-%22', '=%22', '#%22', '$%22', ':%22', '.%22',
'-%23', '=%23', '#%23', '$%23', ':%23', '.%23',
'-%24', '=%24', '#%24', '$%24', ':%24', '.%24',
'-%25', '=%25', '#%25', '$%25', ':%25', '.%25',
'-%26', '=%26', '#%26', '$%26', ':%26', '.%26',
'-%27', '=%27', '#%27', '$%27', ':%27', '.%27',
'-%28', '=%28', '#%28', '$%28', ':%28', '.%28',
'-%29', '=%29', '#%29', '$%29', ':%29', '.%29',
'-%30', '=%30', '#%30', '$%30', ':%30', '.%30',
'-%31', '=%31', '#%31', '$%31', ':%31', '.%31',
'-%32', '=%32', '#%32', '$%32', ':%32', '.%32',
'-%33', '=%33', '#%33', '$%33', ':%33', '.%33',
'-%34', '=%34', '#%34', '$%34', ':%34', '.%34',
'-%35', '=%35', '#%35', '$%35', ':%35', '.%35',
'-%36', '=%36', '#%36', '$%36', ':%36', '.%36',
'-%37', '=%37', '#%37', '$%37', ':%37', '.%37',
'-%38', '=%38', '#%38', '$%38', ':%38', '.%38',
'-%39', '=%39', '#%39', '$%39', ':%39', '.%39',
'-%40', '=%40', '#%40', '$%40', ':%40', '.%40',
'-%41', '=%41', '#%41', '$%41', ':%41', '.%41',
'-%42', '=%42', '#%42', '$%42', ':%42', '.%42',
'-%43', '=%43', '#%43', '$%43', ':%43', '.%43',
'-%44', '=%44', '#%44', '$%44', ':%44', '.%44',
'-%45', '=%45', '#%45', '$%45', ':%45', '.%45',
'-%46', '=%46', '#%46', '$%46', ':%46', '.%46',
'-%47', '=%47', '#%47', '$%47', ':%47', '.%47',
'-%48', '=%48', '#%48', '$%48', ':%48', '.%48',
'-%49', '=%49', '#%49', '$%49', ':%49', '.%49',
'-%50', '=%50', '#%50', '$%50', ':%50', '.%50',
'-%51', '=%51', '#%51', '$%51', ':%51', '.%51',
'-%52', '=%52', '#%52', '$%52', ':%52', '.%52',
'-%53', '=%53', '#%53', '$%53', ':%53', '.%53',
'-%54', '=%54', '#%54', '$%54', ':%54', '.%54',
'-%55', '=%55', '#%55', '$%55', ':%55', '.%55',
'-%56', '=%56', '#%56', '$%56', ':%56', '.%56',
'-%57', '=%57', '#%57', '$%57', ':%57', '.%57',
'-%58', '=%58', '#%58', '$%58', ':%58', '.%58',
'-%59', '=%59', '#%59', '$%59', ':%59', '.%59',
'-%60', '=%60', '#%60', '$%60', ':%60', '.%60',
'-%61', '=%61', '#%61', '$%61', ':%61', '.%61',
'-%62', '=%62', '#%62', '$%62', ':%62', '.%62',
'-%63', '=%63', '#%63', '$%63', ':%63', '.%63',
'-%64', '=%64', '#%64', '$%64', ':%64', '.%64',
'-%65', '=%65', '#%65', '$%65', ':%65', '.%65',
'-%66', '=%66', '#%66', '$%66', ':%66', '.%66',
'-%67', '=%67', '#%67', '$%67', ':%67', '.%67',
'-%68', '=%68', '#%68', '$%68', ':%68', '.%68',
'-%69', '=%69', '#%69', '$%69', ':%69', '.%69',
'-%70', '=%70', '#%70', '$%70', ':%70', '.%70',
'-%71', '=%71', '#%71', '$%71', ':%71', '.%71',
'-%72', '=%72', '#%72', '$%72', ':%72', '.%72',
'-%73', '=%73', '#%73', '$%73', ':%73', '.%73',
'-%74', '=%74', '#%74', '$%74', ':%74', '.%74',
'-%75', '=%75', '#%75', '$%75', ':%75', '.%75',
'-%76', '=%76', '#%76', '$%76', ':%76', '.%76',
'-%77', '=%77', '#%77', '$%77', ':%77', '.%77',
'-%78', '=%78', '#%78', '$%78', ':%78', '.%78',
'-%79', '=%79', '#%79', '$%79', ':%79', '.%79',
'-%80', '=%80', '#%80', '$%80', ':%80', '.%80',
'-%81', '=%81', '#%81', '$%81', ':%81', '.%81',
'-%82', '=%82', '#%82', '$%82', ':%82', '.%82',
'-%83', '=%83', '#%83', '$%83', ':%83', '.%83',
'-%84', '=%84', '#%84', '$%84', ':%84', '.%84',
'-%85', '=%85', '#%85', '$%85', ':%85', '.%85',
'-%86', '=%86', '#%86', '$%86', ':%86', '.%86',
'-%87', '=%87', '#%87', '$%87', ':%87', '.%87',
'-%88', '=%88', '#%88', '$%88', ':%88', '.%88',
'-%89', '=%89', '#%89', '$%89', ':%89', '.%89',
'-%90', '=%90', '#%90', '$%90', ':%90', '.%90',
'-%91', '=%91', '#%91', '$%91', ':%91', '.%91',
'-%92', '=%92', '#%92', '$%92', ':%92', '.%92',
'-%93', '=%93', '#%93', '$%93', ':%93', '.%93',
'-%94', '=%94', '#%94', '$%94', ':%94', '.%94',
'-%95', '=%95', '#%95', '$%95', ':%95', '.%95',
'-%96', '=%96', '#%96', '$%96', ':%96', '.%96',
'-%97', '=%97', '#%97', '$%97', ':%97', '.%97',
'-%98', '=%98', '#%98', '$%98', ':%98', '.%98',
'-%99', '=%99', '#%99', '$%99', ':%99', '.%99') |> sort()
class <- "bm_tre_4"
symbClass_description <- "Four character bond multiplying symbols terminators of rings with explicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_tre, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_tre_4"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_tre, next from start
characters <- lapply(symbols, char_2) |> unlist() |> unique() |> sort()
charClass <- "n_bm_tre"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_tre, rest
characters <- lapply(symbols, chars_2rest) |> unlist() |> unique() |> sort()
charClass <- "r_bm_tre"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 27, symbol class: Three-character bond multiplying symbols terminators of rings with implicit bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' ) |> sort()
class <- "bm_tri_3"
symbClass_description <- "Three character bond multiplying symbols terminators of rings with implicit bond"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# bm_iri, start
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_bm_tri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# bm_iri, rest
characters <- lapply(symbols, chars_1rest) |> unlist() |> unique() |> sort()
charClass <- "r_bm_tri"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

### 6 cis_trans
symbType <- "cis_trans__symbol"
symbType_description <- "Cis / trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond"
## cis_trans, SYMBOLS
# 28, symbol class: cis_trans single character symbols on the left side of the rotary non-permissive bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("/", "\\") |> sort()
class <- "ct"
symbClass_description <- "Single character symbols describing cis / trans features"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "l_ct"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# 29, cis_trans single character symbols on the right side of the rotary non-permissive bond
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("/", "\\") |> sort()
class <- "ct"
symbClass_description <- "Single character symbols describing cis / trans features"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- characters
charClass <- "r_ct"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1


### 7 All the symb_data inside the square brackets besides the main atom symbol
symbType <- "features__symbol"
symbType_description <- "Symbols of atomic features written inside the square brackets"
## All the symb_data inside the square brackets besides the main atom symbol, SYMBOLS
# 30, symbol class: Single character isotope symbols
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9') |> sort()
class <- "isotope"
symbClass_description <- "Single character isotopic symbols"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_isotope
characters <- symbols
charClass <- "w_isotope"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 31, symbol class: Multicharacter (from 2 to 3 characters) isotope symbols
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symb__isotope_1 <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
symb__isotope_2_df <- expand.grid(symb__isotope_1, symb__isotope_1)
symb__isotope_2 <- apply(X = symb__isotope_2_df, MARGIN = 1, FUN = paste, collapse = "")
symb__isotope_3_df <- expand.grid(symb__isotope_1, symb__isotope_2)
symb__isotope_3 <- apply(X = symb__isotope_3_df, MARGIN = 1, FUN = paste, collapse = "")
symbols <- c(symb__isotope_1, symb__isotope_2, symb__isotope_3) |> unique() |> sort()
class <- "isotope_m"
symbClass_description <- "Multi character isotopic symbols"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_isotope
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_isotope"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# r_isotope
characters <- lapply(symbols, chars_1rest) |> unlist() |> unique() |> sort()
charClass <- "r_isotope"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 32, symbol class: Single character chirality symbol
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('@') |> sort()
class <- "chiral"
symbClass_description <- "Single character symbols describing chirality"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
characters <- symbols
charClass <- "w_chiral"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 33, symbol class: two-character chirality symbol
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c('@@') |> sort()
class <- "chiral_2"
symbClass_description <- "Two character symbols describing chirality"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_chiral
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_chiral"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# e_chiral
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_chiral"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 34, symbol class: Multi-character (four or five character) chirality symbols, SYMBOLS
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3", "@TB1", "@TB2", "@TB3", "@TB4",
                    "@TB5", "@TB6", "@TB7", "@TB8", "@TB9", "@TB10", "@TB11", "@TB12", "@TB13", "@TB14", "@TB15", "@TB16", "@TB17", "@TB18", "@TB19", "@TB20",
                    "@OH1", "@OH2", "@OH3", "@OH4", "@OH5", "@OH6", "@OH7", "@OH8", "@OH9", "@OH10", "@OH11", "@OH12", "@OH13", "@OH14", "@OH15", "@OH16", "@OH17", "@OH18", "@OH19", "@OH20",
                    "@OH21", "@OH22", "@OH23", "@OH24", "@OH25", "@OH26", "@OH27", "@OH28", "@OH29", "@OH30") |> sort()
class <- "chiral_m"
symbClass_description <- "Multi character symbols describing chirality"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_chiral
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_chiral_m"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# n_chiral
characters <- lapply(symbols, chars_2_3) |> unlist() |> unique() |> sort()
charClass <- "m_chiral"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# r_chiral
characters <- lapply(symbols, chars_3rest) |> unlist() |> unique() |> sort()
charClass <- "r_chiral"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 35, symbol class: Single character hydrogen symbols
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("H") |> sort()
class <- "hydro"
symbClass_description <- "Single character symbols describing hydrogens"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_hydro
characters <- symbols
charClass <- "w_hydro"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 36, symbol class: two-character hydrogen symbols
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9") |> sort()
class <- "hydro_2"
symbClass_description <- "Two character symbols describing hydrogens"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_hydro
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_hydro"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# e_hydro
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_hydro"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 37, symbol class: Single character charge symbols, SYMBOLS
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("+", "-") |> sort()
class <- "charge"
symbClass_description <- "Single character symbols describing charge"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# w_charge
characters <- symbols
charClass <- "w_charge"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 38, symbol class: Two-character charge obsolete symbols, SYMBOLS
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("++", "--") |> sort()
class <- "charge_2"
symbClass_description <- "Two character symbols describing charge"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_charge
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_charge"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# e_charge
characters <- lapply(symbols, char_last) |> unlist() |> unique() |> sort()
charClass <- "e_charge"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 39, symbol class: Multicharacter (two or three characters) charge symbols, SYMBOLS
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c("+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15",
                     "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15") |> sort()
class <- "charge_m"
symbClass_description <- "Multi character symbols describing charge"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_charge
characters <- lapply(symbols, char_1) |> unlist() |> unique() |> sort()
charClass <- "s_charge_m"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# r_charge
characters <- lapply(symbols, chars_1rest) |> unlist() |> unique() |> sort()
charClass <- "r_charge"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

# 40, symbol class: Multicharacter (from 2 to 4 characters) class symbols, SYMBOLS
# Symbols
symbClass_cntr <- symbClass_cntr + 1
symbols <- c(":[0-9]", ":[0-9][0-9]", ":[0-9][0-9][0-9]") |> sort()
class <- "class_m"
symbClass_description <- "Multi character symbols describing atom class"
symb_data[symbClass_cntr, 1] <- symbClass_cntr
symb_data[symbClass_cntr, 2] <- class
symb_data[symbClass_cntr, 3] <- symbType
symb_data[symbClass_cntr, 4] <- paste(symbols, collapse = ", ")
symb_data[symbClass_cntr, 5] <- symbClass_description
symb_data[symbClass_cntr, 6] <- symbType_description
# Characters
# s_class
characters <- c(":") |> sort()
charClass <- "s_class"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1
# r_class
characters <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9') |> sort()
charClass <- "r_class"
char_data[charClass_cntr, 1] <- charClass_cntr
char_data[charClass_cntr, 2] <- charClass
char_data[charClass_cntr, 3] <- paste(characters, collapse = ", ")
char_data[charClass_cntr, 5] <- symbClass_cntr
charClass_cntr <- charClass_cntr + 1

#### Process the data
# Process
char_data <- char_data[!is.na(char_data$charClass_id), ]
symb_data <- symb_data[!is.na(symb_data$symbClass_id), ] 
chars__symbols <- char_data |> merge(symb_data, by = "symbClass_id")
chars__symbols$charClass_description <- paste0(rep("Characters describing ", length(chars__symbols[, 9])), chars__symbols[, 9])
# Update char_data using information from chars__symbols
char_data <- chars__symbols[, 1:5]

#### Export the results
write.table(char_data, "C:/.../character_classes.tsv", sep = "\t", row.names = FALSE)
write.table(symb_data, "C:/.../symbol_classes.tsv", sep = "\t", row.names = FALSE)
write.table(chars__symbols, "C:/.../chars_&_symbs.tsv", sep = "\t", row.names = FALSE)