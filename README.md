# SMILES Parser (work in progress)

For the SMILES (Simplified Molecular Input Line Entry System) reference, please, SEE:

-   <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>

-   <http://opensmiles.org/opensmiles.html>

-   <https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System>

In short:

-   SMILES is a language, which is used to describe a molecule in the terms of and according to the rules of chemical valence and mathemathical graph theories.

-   SMILES string is the linear notation of a spanning tree of a graph representing molecule.

-   SMILES string includes atom symbols, bonds between them and some characteristics of these entities.

-   Also, according to the Wikipedia, "From the view point of a formal language theory, SMILES is a word. A SMILES is parsable with a context-free parser".

-   To build such a parser it will be useful to understand, what are the symbols and corresponding characters allowed in SMILES and what are their allowed combinations.

## Analysis of characters allowed in SMILES

### Atom symbols and corresponding characters

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two categories by their length:

-   Symbols consisting of the single character

-   Symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two categories by their grammatic requirements:

-   Symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   Symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on their nature and surroundings:

-   Symbols of the aromatic atoms

-   Symbols of the aliphatic atoms

Thus, the following type of atom symbols probably could be found in SMILES:

-   Single character atom symbols of organic aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could designated as **watom_oar**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

-   Single character atom symbols of organic aliphatic atoms:

| B, C, N, O, S, P, F, I

Corresponding characters could designated as **watom_oal**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

-   Single character atom symbols of bracket aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could designated as **watom_bar**,where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic.

-   Single character atom symbols of bracket aliphatic atoms:

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could designated as **watom_bal**,where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic.

-   Two character atom symbols of organic aromatic atoms:

| Such symbols do not exist

-   Two character atom symbols of organic aliphatic atoms:

| Cl, Br

Corresponding characters could designated as **satom_oal** & **eatom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic.

-   Two character atom symbols of bracket aromatic atoms:

| se, as, te

Corresponding characters could designated as **satom_bar** & **eatom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

-   Two character atom symbols of bracket aliphatic atoms:

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could designated as **satom_bal** & **eatom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

#### Question

So, "SMILES is parsable with a context-free parser", which in practice probably means that given the set of SMILES strings, each of them could be correctly parsed without considering any other string.

Maybe it is also possible to guess the atom and its properties expressed as the type of atom symbol in the single SMILES string by a single character, without considering any other character in this SMILES string?

To be honest, the answer is definately not, but the details matter.

Trying to answer this question by checking the intersection between the sets of characters describing atom symbols of different types, it is possible to get some insights into what also should be considered, besides the one character at time, while parsing SMILES string.

This check up is quite easy to do using R (<https://cran.rstudio.org/bin/windows/>), Tidyverse (<https://tidyverse.org/>) and ggupset (<https://github.com/const-ae/ggupset>):

``` {#char_intersect .R .illustration}
library(tidyverse)
library(ggupset)
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
# each charcter corresponds to the whole symbol
watom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
# each charcter corresponds to the whole symbol
watom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
# each charcter corresponds to the whole symbol
watom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
# each charcter corresponds to the whole symbol
watom_bal <- symb__atom_bal
#05, two character atom symbol of organic aromatic atoms
# symb_atom_oar IS EMPTY
# 06, two character atom symbol of organic aliphatic atoms
symb__atom_oal <- c('Cl', 'Br')
# Get the first character for each atom symbol
satom_oal <- map(symb__atom_oal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_oal <- map(symb__atom_oal, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 07, two character atom symbol of bracket aromatic atom
symb__atom_bar <- c('se', 'as', 'te')
# Get the first character for each atom symbol
satom_bar <- map(symb__atom_bar, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bar <- map(symb__atom_bar, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 08, two character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar',
                    'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb',
                    'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
                    'Te', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt',
                    'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
                    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl',
                    'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np',
                    'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
satom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()

## Analyze characters constituting atom symbols from different classes
atom_symbol_characters <- tibble(
    character = c(watom_oar,
                  watom_oal,
                  watom_bar,
                  watom_bal,
                  satom_oal,
                  eatom_oal,
                  satom_bar,
                  eatom_bar,
                  satom_bal,
                  eatom_bal),
    class = c(rep("watom_oar", length(watom_oar)),
              rep("watom_oal", length(watom_oal)),
              rep("watom_bar", length(watom_bar)),
              rep("watom_bal", length(watom_bal)),
              rep("satom_oal", length(satom_oal)),
              rep("eatom_oal", length(eatom_oal)),
              rep("satom_bar", length(satom_bar)),
              rep("eatom_bar", length(eatom_bar)),
              rep("satom_bal", length(satom_bal)),
              rep("eatom_bal", length(eatom_bal)))
    ) |>
    group_by(character) |>
    summarise(class = list(class))
                    
# UpSet plot for classes
plot_class <- ggplot(atom_symbol_characters, aes(x=class)) +
                geom_bar() +
                scale_x_upset() +
                theme_minimal()
plot_class
```

#### Result

![**Figure 1.** Intersections of the character sets describing different atom symbol types.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_atomchars.png)

As it can be clearly seen from the **Figure 1:** yes, it will not be possible to guess the atom and its type operating single atom symbol character at time while parsing SMILES string, since there are intersections between the different sets of characters. However, on the positive side: intersections are not uniform:

| e.g., satom_bal, watom_bal & watom_oal sets have the **largest** intersection while intersection between eatom_bal, watom_bar, watom_oarm, satom_bar andeatom_bar sets is the **smallest**

Probably, this observations could be used latter to provide the future SMILES parser with some efficiency.

### Bond symbols and corresponding characters

#### What are they?

....
