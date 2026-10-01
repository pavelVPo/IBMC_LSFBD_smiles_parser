##### This code is the basis for the further analysis and construction.
#### Functions
# Function to get the first character from word
char_frst <- function(word) {
    first_char <- strsplit(word, split = "")[[1]][1]
    return(first_char)
}
# Function to get the second and third characters from word
chars_x_y <- function(word, x, y) {
    chars <- strsplit(word, split = "")[[1]][x:y]
    return(chars)
}

#### Prepare the initial data according to the README.md
path <- ".../SMILES_parser/data_v4/"
data   <- data.frame(record_id = seq(1:42),
                     type = c(rep("atom", 9), rep("square", 2), rep("bond", 6), rep("modifier", 12), rep("cis_trans", 2), rep("property", 11)),
                     class = c("atom_oar", "atom_oal", "atom_bar", "atom_bal", "atom_oal_2", "atom_bar_2", "atom_bal_2", "anything", "ambient", "s_square", "e_square", "single_bond", "double_bond", "triple_bond", "quadruple_bond", "no_bond", "bm_ibi", "bm_iri", "bm_tbi", "bm_tri", "bm_ibe", "bm_ibe_2", "bm_ire_2", "bm_iri_3", "bm_ire_4", "bm_tbe_2", "bm_tre_2", "bm_tre_4", "bm_tri_3", "l_ct", "r_ct", "isotope", "isotope_m", "chiral", "chiral_2", "chiral_m", "hydro", "hydro_2", "charge", "charge_2", "charge_m", "class"),
                     symbols_draft = c("b, c, n, o, s, p", "B, C, N, O, S, P, F, I", "b, c, n, o, s, p", "H, B, C, N, O, F, P, S, K, V, Y, I, W, U", "Cl, Br", "se, as, te", "He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr", "*", "", "[", "]", "-", "=", "#", "$", ":", ".", "(", "0, 1, 2, 3, 4, 5, 6, 7, 8, 9", ")", "0, 1, 2, 3, 4, 5, 6, 7, 8, 9", "([-=#$:.]", "[-=#$:.][0:9]", "%[0:9][1:9], %[1:9][0:9]", "[-=#$:.]%[0:9][1:9], [-=#$:.]%[1:9][0:9]", ")[-=#$:.]", "[-=#$:.][0:9]", "[-=#$:.]%[0:9][1:9], [-=#$:.]%[1:9][0:9]", "%[0:9][1:9], %[1:9][0:9]", "[/\\]", "[/\\]", "1, 2, 3, 4, 5, 6, 7, 8, 9", "[0:9][1:9], [1:9][0:9], [0:9][0:9][1:9], [0:9][1:9][0:9], [1:9][0:9][0:9]", "@", "@@", "@TH[1:2], @AL[1:2], @SP[1:3], @TB[1:20], @OH[1:30]", "H", "H[2:9]", "+, -", "++, --", "[+-][1:9], [+-]1[0:5]", ":[0:9], :[0:9][0:9], :[0:9][0:9][0:9]"),
                     type_description = c(rep("node of the molecule's graph, i.e. an atom OR smth having similar behavior according to the SMILES rules", 9),
                                          rep("mark of the start/end of the record for an atom having some property (-ies)", 2),
                                          rep("edge of the molecule's graph, i.e. chemical bond OR smth having similar behavior according to the SMILES rules", 6),
                                          rep("bond extension for the linear notation", 12),
                                          rep("Cis/trans", 2),
                                          rep("properties of an atom", 11)),
                     class_description = c("Single character atom symbols of organic (from so called organic subset) aromatic atoms lacking the additional grammatical requirements and features",
                                            "Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features",
                                            "Single character atom symbols of aromatic atoms enclosed within brackets",
                                            "Single character atom symbols of bracket aliphatic atoms",
                                            "Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features",
                                            "Two character atom symbols of in-bracket aromatic atoms",
                                            "Two character atom symbols of in-bracket aliphatic atoms",
                                            "Single character symbol of any atom or basically anything",
                                            "Zero character (pseudo)symbol trailing start and end of the SMILES string",
                                            "Single character symbol to start the record for atom having property (-ies)",
                                            "Single character symbol to end the record for atom having property (-ies)",
                                            "Single character bond symbol corresponding to the single bond",
                                            "Single character bond symbol corresponding to the double bond",
                                            "Single character bond symbol corresponding to the triple bond",
                                            "Single character bond symbol corresponding to the quadruple bond",
                                            "Single character bond symbol corresponding to the aromatic bond",
                                            "Single character bond symbol corresponding to the absence of the bond",
                                            "Single character bond multiplying symbols initiators of branching with implicit bond",
                                            "Single character bond multiplying symbols initiators of rings with implicit bond",
                                            "Single character bond multiplying symbols terminators of branching with implicit bond",
                                            "Single character bond multiplying symbols terminators of rings with implicit bond",
                                            "Two-character bond multiplying symbols initiators of branching with explicit bond",
                                            "Two-character bond multiplying symbols initiators of rings with explicit bond",
                                            "Three-character bond multiplying symbols initiators of rings with implicit bond",
                                            "Four-character bond multiplying symbols initiators of rings with explicit bond",
                                            "Two-character bond multiplying symbols terminators of branching with explicit bond",
                                            "Two-character bond multiplying symbols terminators of rings with explicit bond",
                                            "Four-character bond multiplying symbols terminators of rings with explicit bond",
                                            "Three-character bond multiplying symbols terminators of rings with implicit bond",
                                            "Cis/trans single character symbols on the left side of the rotary non-permissive bond",
                                            "Cis/trans symbols on the right side of the rotary non-permissive bond",
                                            "Single character isotope symbols",
                                            "Multicharacter (from 2 to 3 characters) isotope symbols",
                                            "Single character chirality symbol",
                                            "Two-character chirality symbol",
                                            "Multicharacter (four or five character) chirality symbols",
                                            "Single character hydrogen symbol",
                                            "Two-character hydrogen symbols",
                                            "Single character charge symbols",
                                            "Two-character charge obsolete symbols",
                                            "Multicharacter (two or three characters) charge symbols",
                                            "Multicharacter (from 2 to 4 characters) atom class symbols")
                     )

#### Extend the list of symbols
for (i in seq(1:nrow(data))) {
    if (data[i,4] == "([-=#$:.]" | data[i,4] == ")[-=#$:.]") {
        data[i,4] <- paste0(char_frst(data[i,4]), chars_x_y(data[i,4], 3, 8)) |> paste0(collapse = ", ")
    } else if (data[i,4] == "[-=#$:.][0:9]") {
        data[i,4] <- expand.grid(chars_x_y(data[i,4], 2, 7), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "") |> paste0(collapse = ", ")
    } else if (data[i,4] == "%[0:9][1:9], %[1:9][0:9]") {
        draft_1   <- expand.grid("%", (seq(0:9)-1), seq(1:9) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_2   <- expand.grid("%", seq(1:9), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft     <- c(draft_1, draft_2) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    } else if (data[i,4] == "[-=#$:.]%[0:9][1:9], [-=#$:.]%[1:9][0:9]") {
        draft_1   <- expand.grid(chars_x_y(data[i,4], 2, 7), "%", (seq(0:9)-1), seq(1:9) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_2   <- expand.grid(chars_x_y(data[i,4], 2, 7), "%", seq(1:9), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft     <- c(draft_1, draft_2) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    } else if (data[i,4] == "[/\\]") {
        data[i,4] <- "/, \\"
    } else if (data[i,4] == "[0:9][1:9], [1:9][0:9], [0:9][0:9][1:9], [0:9][1:9][0:9], [1:9][0:9][0:9]") {
        draft_1   <- expand.grid((seq(0:9)-1), seq(1:9) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_2   <- expand.grid(seq(1:9), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_3   <- expand.grid((seq(0:9)-1), (seq(0:9)-1), seq(1:9) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_4   <- expand.grid((seq(0:9)-1), seq(1:9), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_5   <- expand.grid(seq(1:9), (seq(0:9)-1), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft     <- c(draft_1, draft_2, draft_3, draft_4, draft_5) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    } else if(data[i,4] == "@TH[1:2], @AL[1:2], @SP[1:3], @TB[1:20], @OH[1:30]") {
        draft_1   <- paste0("@TH", seq(1:2)) |> paste0(collapse = ", ")
        draft_2   <- paste0("@AL", seq(1:2)) |> paste0(collapse = ", ")
        draft_3   <- paste0("@SP", seq(1:3)) |> paste0(collapse = ", ")
        draft_4   <- paste0("@TB", seq(1:20)) |> paste0(collapse = ", ")
        draft_5   <- paste0("@OH", seq(1:2)) |> paste0(collapse = ", ")
        draft     <- c(draft_1, draft_2, draft_3, draft_4, draft_5) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    } else if(data[i,4] == "H[2:9]") {
        data[i,4] <- paste0("H", (seq(1:8)+1)) |> paste0(collapse = ", ")
    } else if(data[i,4] == "[+-][1:9], [+-]1[0:5]") {
        draft_1   <- expand.grid(c("+", "-"), seq(1:9) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_2   <- expand.grid(c("+", "-"), "1", (seq(1:6)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft     <- c(draft_1, draft_2) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    } else if(data[i,4] == ":[0:9], :[0:9][0:9], :[0:9][0:9][0:9]") {
        draft_1   <- paste0(":", (seq(0:9)-1)) |> paste0(collapse = ", ")
        draft_2   <- expand.grid(":", (seq(0:9)-1), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft_3   <- expand.grid(":", (seq(0:9)-1), (seq(0:9)-1), (seq(0:9)-1) |> as.character()) |> apply(1, paste0, collapse = "")
        draft     <- c(draft_1, draft_2, draft_3) |> unique() |> sort() |> paste0(collapse = ", ")
        data[i,4] <- draft
    }
}

#### Export the results
write.table(data, paste0(path, "symbols.tsv"), sep = "\t", row.names = FALSE)