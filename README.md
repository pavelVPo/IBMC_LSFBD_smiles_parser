# SMILES Parser (work in progress)

For the SMILES (Simplified Molecular Input Line Entry System) reference, please, SEE:

-   <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>

-   <http://opensmiles.org/opensmiles.html>

-   <https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System>

In short:

-   SMILES is a language, which is used to describe a molecule in the terms of and according to the rules of chemical valence and mathematical graph theories.

-   SMILES string is the linear notation of a spanning tree of a graph representing molecule.

-   SMILES string includes atom symbols, bonds between them and some characteristics of these entities.

-   Also, according to the Wikipedia, "From the view point of a formal language theory, SMILES is a word. A SMILES is parsable with a context-free parser".

-   To build such a parser it will be useful to understand, what are the symbols and corresponding characters allowed in SMILES and what are their allowed combinations.

## Enumeration of symbols and characters allowed in SMILES

| *This part may require some further adjustments and corrections, but should be OK in general.*

### Atom symbols and corresponding characters

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two *categories* by their length:

-   Symbols consisting of the single character

-   Symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two categories by their grammatical requirements:

-   Symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   Symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on their nature and surroundings:

-   Symbols of the aromatic atoms

-   Symbols of the aliphatic atoms

##### Thus, the following ***classes*** of atom symbols probably could be found in SMILES:

-   Single character atom symbols of organic (from so called *organic* subset) aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could be designated as **watom_oar**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

-   Single character atom symbols of organic aliphatic atoms:

| B, C, N, O, S, P, F, I

Corresponding characters could be designated as **watom_oal**,where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

-   Single character atom symbols of bracket aromatic atoms:

| b, c, n, o, s, p

Corresponding characters could be designated as **watom_bar**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic. As it can be seen, this class contains the same symbols as watom_oar, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the watom_bar.

-   Single character atom symbols of bracket aliphatic atoms:

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **watom_bal**,where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic. As it can be seen, this class contains the same symbols as watom_oal, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the watom_bar.

-   Two character atom symbols of organic aliphatic atoms:

| Cl, Br

Corresponding characters could be designated as **satom_oal** & **eatom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic.

-   Two character atom symbols of bracket aromatic atoms:

| se, as, te

Corresponding characters could be designated as **satom_bar** & **eatom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

-   Two character atom symbols of bracket aliphatic atoms:

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **satom_bal** & **eatom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

Also, **\*** is an allowed symbol in SMILES, which corresponds to the any atom symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms:

-   Single character atom symbol of any atom:

| \*

Corresponding character could be designated as **anything**, since it could mean basically anything.

#### Question

So, "SMILES is parsable with a context-free parser", which in practice means that given the set of SMILES strings, each of them could be correctly parsed without considering any other string.

Maybe it is also possible to guess the atom and its properties expressed as the class of an atom symbol in the single SMILES string by a single character, without considering any other character in this SMILES string?

To be honest, the answer is definitely not, but the details matter.

Trying to answer this question by checking the intersection between the sets of characters describing atom symbols of different classes, it is possible to get some insights into what also should be considered, besides the one character at time, while parsing SMILES string.

This check up is quite easy to do using R (<https://cran.rstudio.org/bin/windows/>), Tidyverse (<https://tidyverse.org/>) and ggupset (<https://github.com/const-ae/ggupset>):

``` {#char_intersect__atomSymbols .R .illustration}
library(tidyverse)
library(ggupset)

### Enumeration of the character classes allowed in SMILES
## ATOMS
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
watom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
watom_bal <- symb__atom_bal
#06, two character atom symbol of organic aliphatic atoms
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
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
satom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()
## ANYTHING
# 09, single character symbol of anything
symb__anything <- c('*')
anything <- symb__anything
## SQUARE BRACKETS
# 10, single character marking the start of the atom record
bracket_start <- c('[')
# 11, single character marking the end of the atom record
bracket_end <- c(']')

## Analyze characters constituting atom symbols from different classes
# Prepare the data
atom_symbol_characters <- tibble(character = c( watom_oar,
watom_oal,
watom_bar,
watom_bal,
satom_oal,
eatom_oal,
satom_bar,
eatom_bar,
satom_bal,
eatom_bal,
anything,
bracket_start,
bracket_end),
class = c( rep("watom_oar", length(watom_oar)),
rep("watom_oal", length(watom_oal)),
rep("watom_bar", length(watom_bar)),
rep("watom_bal", length(watom_bal)),
rep("satom_oal", length(satom_oal)),
rep("eatom_oal", length(eatom_oal)),
rep("satom_bar", length(satom_bar)),
rep("eatom_bar", length(eatom_bar)),
rep("satom_bal", length(satom_bal)),
rep("eatom_bal", length(eatom_bal)),
rep("anything", length(anything)),
rep("bracket_start", length(bracket_start)),
rep("bracket_end", length(bracket_end)) )) |>
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

![**Figure 1.** Intersections of the character sets describing different atom symbol classes.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_atomClasses.png)

**Figure 1.** Intersections of the character sets describing different atom symbol classes.

As it can be clearly seen from the **Figure 1:** yes, it will not be possible to guess the atom and its class operating single atom symbol character at time while parsing SMILES string, since there are intersections between the different sets of characters. However, on the positive side: intersections are not uniform:

| e.g.,
| satom_bal, watom_bal & watom_oal sets have the **largest** intersection, while intersection between eatom_bal, watom_bar, watom_oarm, satom_bar and eatom_bar sets is the **smallest**

Probably, this observations could be used latter to provide the future SMILES parser with some efficiency.

### Square bracket symbols and corresponding characters

#### What are they?

Square bracket symbol **[\*\* is the SMILES way to mark the start of the atom record including its various properties and \*\*]** is the way to designate the end of the atom record.

| [, ]

Corresponding characters could be designated as **bracket_start & bracket_end**.

### Bond symbols and corresponding characters

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, five of which corresponds to the certain type of chemical bond:**

-   Single character bond symbol corresponding to the single bond:

| -

This single character symbol could be and typically is omitted, since by default all the atoms, which symbols are written side by side in SMILES string, are presumed to be connected by this type of bond. Corresponding character will be designated as **single_bond**.

-   Single character bond symbol corresponding to the double bond:

| =

Corresponding character will be designated as **double_bond**.

-   Single character bond symbol corresponding to the triple bond:

| \#

Corresponding character will be designated as **triple_bond**.

-   Single character symbol corresponding to the quadruple bond:

| \$

Corresponding character will be designated as **quadruple_bond**.

-   Single character symbol corresponding to the aromatic bond:

| :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications. Corresponding character will be designated as **aromatic_bond_obsolete**.

-   Single character symbol corresponding to the absence of the bond between the two specific atoms:

| .

As it was said earlier, by default each atom in SMILES string is considered to be connected with its immediate neighbors via the single bond (**-**). Thus, symbol corresponding to the negation of the bond is needed sometimes, and here it is: **.** Corresponding character will be designated as **no_bond**.

### Bond modifying (multiplying) symbols and corresponding characters

#### What are they?

Bond multiplying symbols is used in SMILES to extend the number of atoms, for which connections to the current atom could be written using linear notation (SMILES).

Bond multiplying symbols allowed in SMILES could be divided into two categories by their length:

-   Symbols consisting of the single character

-   Symbols consisting of two characters

-   Symbols consisting of three characters

-   Symbols consisting of four characters

Bond multiplying symbols allowed in SMILES could be divided into two categories according to their role in completing the task:

-   Symbols initiators

-   Symbols terminators

Bond multiplying symbols allowed in SMILES could be divided into two categories by their task:

-   Symbols used to indicate simple additional bond for the current atom (branch)

-   Symbols used to indicate additional bond, which allows for the cycle (ring) to be formed

Bond multiplying symbols allowed in SMILES could be divided into two categories according to the state and quality of the additional bond:

-   Symbols explicitly including additional bond (any bond could be added)

-   Symbols implicitly including additional bond (only single bond could be added)

##### Thus, the following classes of bond multiplying symbols probably could be found in SMILES:

-   Single character bond multiplying symbols initiators of branching with implicit bond:

| (

Corresponding characters could be designated as **wbm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols initiators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wbm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols terminators of branching with implicit bond:

| )

Corresponding characters could be designated as **wbm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Single character bond multiplying symbols terminators of rings with implicit bond:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wbm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

-   Two-character bond multiplying symbols initiators of branching with explicit bond:

| (-, (=, (#, (\$, (:, (.

Corresponding characters could be designated as **sbm_ibe** & **ebm_ibe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

-   Two-character bond multiplying symbols initiators of rings with explicit bond:

| -0, =0, #0, \$0, :0, .0,
| -1, =1, #1, \$1, :1, .1,
| -2, =2, #2, \$2, :2, .2,
| -3, =3, #3, \$3, :3, .3,
| -4, =4, #4, \$4, :4, .4,
| -5, =5, #5, \$5, :5, .5,
| -6, =6, #6, \$6, :6, .6,
| -7, =7, #7, \$7, :7, .7,
| -8, =8, #8, \$8, :8, .8,
| -9, =9, #9, \$9, :9, .9

Corresponding characters could be designated as **sbm_ire_2** & **ebm_ire_2**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **2** - for two-character.

-   Three-character bond multiplying symbols initiators of rings with implicit bond:

| %01, %02, %03, %04, %05, %06, %07, %08, %09, %10,
| %11, %12, %13, %14, %15, %16, %17, %18, %19, %20,
| %21, %22, %23, %24, %25, %26, %27, %28, %29, %30,
| %31, %32, %33, %34, %35, %36, %37, %38, %39, %40,
| %41, %42, %43, %44, %45, %46, %47, %48, %49, %50,
| %51, %52, %53, %54, %55, %56, %57, %58, %59, %60,
| %61, %62, %63, %64, %65, %66, %67, %68, %69, %70,
| %71, %72, %73, %74, %75, %76, %77, %78, %79, %80,
| %81, %82, %83, %84, %85, %86, %87, %88, %89, %90,
| %91, %92, %93, %94, %95, %96, %97, %98, %99

Corresponding characters could be designated as **sbm_iri** & **rbm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and last suffix **i** - for implicit.

-   Four-character bond multiplying symbols initiators of rings with explicit bond:

| -%0[1-9], =%0[1-9], #%0[1-9], \$%0[1-9], :%0[1-9], .%0[1-9],
| -%[1-9][0-9], =%[1-9][0-9], #%[1-9][0-9], \$%[1-9][0-9], :%[1-9][0-9], .%[1-9][0-9]

Corresponding characters could be designated as **sbm_ire_4, nbm_ire_4** & **rbm_ire_4**,where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **4** - for four-character.

-   Two-character bond multiplying symbols terminators of branching with explicit bond:

| )-, )=, )#, )\$, ):, ).

Corresponding characters could be designated as **sbm_tbe** & **ebm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

-   Two-character bond multiplying symbols terminators of rings with explicit bond:

| -0, =0, #0, \$0, :0, .0,
| -1, =1, #1, \$1, :1, .1,
| -2, =2, #2, \$2, :2, .2,
| -3, =3, #3, \$3, :3, .3,
| -4, =4, #4, \$4, :4, .4,
| -5, =5, #5, \$5, :5, .5,
| -6, =6, #6, \$6, :6, .6,
| -7, =7, #7, \$7, :7, .7,
| -8, =8, #8, \$8, :8, .8,
| -9, =9, #9, \$9, :9, .9

Corresponding characters could be designated as **sbm_tre_2** & **ebm_tre_2**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **2** - for two-character.

-   Four-character bond multiplying symbols terminators of rings with explicit bond:

| -%0[1-9], =%0[1-9], #%0[1-9], \$%0[1-9], :%0[1-9], .%0[1-9],
| -%[1-9][0-9], =%[1-9][0-9], #%[1-9][0-9], \$%[1-9][0-9], :%[1-9][0-9], .%[1-9][0-9]

Corresponding characters could be designated as **sbm_tre_4, nbm_tre_4** & **rbm_tre_4**,where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit and last suffix **4** - for four-character.

-   Three-character bond multiplying symbols terminators of rings with implicit bond:

| %01, %02, %03, %04, %05, %06, %07, %08, %09, %10,
| %11, %12, %13, %14, %15, %16, %17, %18, %19, %20,
| %21, %22, %23, %24, %25, %26, %27, %28, %29, %30,
| %31, %32, %33, %34, %35, %36, %37, %38, %39, %40,
| %41, %42, %43, %44, %45, %46, %47, %48, %49, %50,
| %51, %52, %53, %54, %55, %56, %57, %58, %59, %60,
| %61, %62, %63, %64, %65, %66, %67, %68, %69, %70,
| %71, %72, %73, %74, %75, %76, %77, %78, %79, %80,
| %81, %82, %83, %84, %85, %86, %87, %88, %89, %90,
| %91, %92, %93, %94, %95, %96, %97, %98, %99

Corresponding characters could be designated as **sbm_tri** & **rbm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, and last suffix i - for implicit.

#### Sort of the same question

From simple enumeration its quite clear that it will not possible to discriminate between the different classes of bond modifying (multiplying) symbols operating one character at time, also it will be not possible to discriminate between the some classes of bond multipliers and bonds while checking the first character of the symbol.

However, what are the intersections between the character sets describing different classes of bond multipliers?

Here is the code:

``` {#char_intersect__bmSymbols .R .illustration}
library(tidyverse)
library(ggupset)

### Enumeration of the character classes allowed in SMILES
##  Bond modifiers (multipliers)
# 18, single character bond multiplying symbols initiators of branch with implicit bond
symb__wbm_ibi <- c('(')
wbm_ibi <- symb__wbm_ibi
# 19, single character bond multiplying symbols terminators of branch with implicit bond
symb__wbm_tbi <- c(')')
wbm_tbi <- symb__wbm_tbi
# 20, single character bond multiplying symbols initiators of rings with implicit bond
symb__wbm_iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_iri <- symb__wbm_iri
# 21, single character bond multiplying symbols terminators of rings with implicit bond
symb__wbm_tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_tri <- symb__wbm_tri
# 22, two-character bond multiplying symbols initiators of branching with explicit bond
symb__bm_ibe <- c('(-', '(=', '(#', '($', '(:', '(.')
sbm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 23, Two-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_2  <- c('-0', '=0', '#0', '$0', ':0', '.0',
'-1', '=1', '#1', '$1', ':1', '.1',
'-2', '=2', '#2', '$2', ':2', '.2',
'-3', '=3', '#3', '$3', ':3', '.3',
'-4', '=4', '#4', '$4', ':4', '.4',
'-5', '=5', '#5', '$5', ':5', '.5',
'-6', '=6', '#6', '$6', ':6', '.6',
'-7', '=7', '#7', '$7', ':7', '.7',
'-8', '=8', '#8', '$8', ':8', '.8',
'-9', '=9', '#9', '$9', ':9', '.9' )
sbm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 24, four-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_4  <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
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
'-%99', '=%99', '#%99', '$%99', ':%99', '.%99')
sbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
nbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 25, Three-character bond multiplying symbols initiators of rings with implicit bond
symb__bm_iri <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' )
sbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 26, two-character bond multiplying symbols terminators of branch with explicit bond
symb__tbe <- c(')-', ')=', ')#', ')$', '):', ').')
sbm_tbe   <- map(symb__tbe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tbe   <- map(symb__tbe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 27, Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_2  <- symb__bm_ire_2
sbm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 28, Four-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_4  <- symb__bm_ire_4
sbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
nbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 29, three-character bond multiplying symbols terminators of rings with implicit bond
symb__bm_tri <- symb__bm_iri
sbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()

## Analyze characters constituting atom symbols from different classes
# Prepare the data
bondMod_symbol_characters <- tibble(character = c( wbm_ibi,
wbm_tbi,
wbm_iri,
sbm_ibe,
ebm_ibe,
sbm_ire_2,
ebm_ire_2,
sbm_ire_4,
nbm_ire_4,
rbm_ire_4,
sbm_iri,
rbm_iri,
sbm_tbe,
ebm_tbe,
sbm_tre_2,
ebm_tre_2,
sbm_tre_4,
nbm_tre_4,
rbm_tre_4,
sbm_tri,
rbm_tri),
class = c( rep("wbm_ibi", length(wbm_ibi)),
rep("wbm_tbi", length(wbm_tbi)),
rep("wbm_iri", length(wbm_iri)),
rep("sbm_ibe", length(sbm_ibe)),
rep("ebm_ibe", length(ebm_ibe)),
rep("sbm_ire_2", length(sbm_ire_2)),
rep("ebm_ire_2", length(ebm_ire_2)),
rep("sbm_ire_4", length(sbm_ire_4)),
rep("nbm_ire_4", length(nbm_ire_4)),
rep("rbm_ire_4", length(rbm_ire_4)),
rep("sbm_iri", length(sbm_iri)),
rep("rbm_iri", length(rbm_iri)),
rep("sbm_tbe", length(sbm_tbe)),
rep("ebm_tbe", length(ebm_tbe)),
rep("sbm_tre_2", length(sbm_tre_2)),
rep("ebm_tre_2", length(ebm_tre_2)),
rep("sbm_tre_4", length(sbm_tre_4)),
rep("nbm_tre_4", length(nbm_tre_4)),
rep("rbm_tre_4", length(rbm_tre_4)),
rep("sbm_tri", length(sbm_tri)),
rep("rbm_tri", length(rbm_tri)) )) |>
group_by(character) |>
summarise(class = list(class))                    
# UpSet plot for classes
plot_class <- ggplot(bondMod_symbol_characters, aes(x=class)) +
                geom_bar() +
                scale_x_upset() +
                theme_minimal()
plot_class
```

![**Figure 2.** Intersections of the character sets describing different bond multiplying symbols classes.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_bmodifiersClasses.png)

**Figure 2.** Intersections of the character sets describing different bond multiplying symbols classes.

As it can be clearly seen from the **Figure 2:** the sets of characters have large intersections in this case, it will be useful later while creating rules for parsing SMILES. At the moment these intersection suggest that the characters describing symbols from different classes were enumerated well: the largest intersection is between the set of characters describing end of the ring symbols or whole rings. As it should be ([0-9]).

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two categories of cis/trans symbols are allowed in SMILES:

-   Cis/trans single character symbols on the left side of the rotary non-permissive bond:

| /, \\

Corresponding single character characters could be designated as **lct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

-   Cis/trans symbols on the right side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **rct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\F)=C/F are equivalent, since

> The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.
>
> *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****
>
> **\*** <http://opensmiles.org/opensmiles.html>

#### Question

Just of curiosity, how many SMILES strings do contain cis/trans symbols, which could be misinterpreted by the *several SMILES interpreters*?

It is quite easy to approximate the answer to this question using

-   ChEMBL data (Zdrazil, Barbara. "Fifteen years of ChEMBL and its role in cheminformatics and drug discovery." *Journal of Cheminformatics* 17.1 (2025): 1-9.)

-   R ([https://cran.rstudio.org/bin/windows/](#0){style="font-size: 11pt;"})

-   Tidyverse ([https://tidyverse.org/](#0){style="font-size: 11pt;"})

-   DBI (<https://cran.r-project.org/web/packages/DBI>)

-   RMariaDB (<https://cran.r-project.org/web/packages/RMariaDB/index.html>)

Here is the code:

``` r
library(RMariaDB)
library(DBI)
library(tidyverse)

# 30, Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__lct <- c("/", "\\\\")
lct <- symb__lct
# 31, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__rct <- c("/", "\\\\")
rct <- symb__rct

## Patterns to search for
# Parsing SMILES is a task, which is harder than it may appears on the first glance
# Thus, at this stage the regexps will be used, which allows to extact substring containing only the first pair of cis/trans symbols
# ^ - stands for the start of the string
# [^\\\\/]* - means thath from the start of the string and up to the next meaningful part of regexp there should not be matches with the characters of cis/trans symbols
# The last meaningful (and variable) part of the regexps stands for the one of the variants of cis/trans symbols usage from http://opensmiles.org/opensmiles.html
# for example, in pattern_baseOne, [^\\(]/[^\\\\/]*=[^\\\\/]*[aA-zZ]/. matches F/C=C/F and does not match C(/F)=C/F
# for example, in pattern_hardOne, *C\\(/.\\)=C/. matches C(/F)=C/F and does not match F/C=C/F
pattern_baseOne <- "^[^\\\\/]*[^\\(]/[Cc]=[Cc]/."
pattern_baseTwo <- "^[^\\\\/]*[^\\(]\\\\.[Cc]=[Cc]/."
pattern_hardOne <- "^[^\\\\/]*[Cc]\\(/.\\)=[Cc]/."
pattern_hardTwo <- "^[^\\\\/]*[Cc]\\(\\\\.\\)=[Cc]/."
str_extract('C(/F)=C/F', pattern_baseOne)
str_extract('F/C=C/F', pattern_baseOne)
str_extract('C(/F)=C/F', pattern_hardOne)
str_extract('F/C=C/F', pattern_hardOne)
## Connect to DB
mysql_password = '***'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_36',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)
## Extract SMILES
cs__query <- dbSendQuery(con, 'SELECT canonical_smiles FROM compound_structures')
cs_smiles <- dbFetch(cs__query) |> distinct()
dbClearResult(cs__query)
# Close the connection
dbDisconnect(con)
## Check patterns against SMILES
cs_smiles_checked <- cs_smiles |> rowwise() |>
                    mutate(
                        pattern_baseOne = str_extract(canonical_smiles, pattern_baseOne),
                        pattern_baseTwo = str_extract(canonical_smiles, pattern_baseTwo),
                        pattern_hardOne = str_extract(canonical_smiles, pattern_hardOne),
                        pattern_hardTwo = str_extract(canonical_smiles, pattern_hardTwo)
                    ) |>
                    ungroup()
## Filter and Count
cs_smiles_matched <- cs_smiles_checked |> filter(if_any(starts_with("pattern"), ~ !is.na(.)))                       # 49 673  records matched one of the patterns
cs_smiles_matched_hard <- cs_smiles_matched |> filter(if_any(starts_with("pattern_hard"), ~ !is.na(.)))             # 0       records matched one of the hard patterns
cs_smiles_matched_base <- cs_smiles_matched |> filter(if_any(starts_with("pattern_base"), ~ !is.na(.)))             # 49 673  records matched one of the base patterns
```

The first pairs of cis/trans symbols in ChEMBL data (canonical smiles) do not include the hard cases similar to the ones described in <http://opensmiles.org/opensmiles.html>

Probably, it is safe to say that the percentage of such cases should be quite low at least in the curated databases.

From this, it is possible to assume that reliability of the SMILES as a form of representation of chemical structures comes not only from the basic rules of this language, but also from the standards of its usage adopted in the community.

Still, ability to parse SMILES using basic rules are essential to maintain this status.

### All the symbols and corresponding characters inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating speciffic isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

-   Single character isotope symbols:

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **wisotope**, where prefix **w** stands for the whole symbol.

-   Multicharacter (from 2 to 3 characters) isotope symbols:

| [0-9][0-9], [0-9][0-9][0-9]

Corresponding characters could be designated as **sisotope_m & risotope_m**, where prefix **s** stands for the start and prefix r stands for the rest of the symbol and suffix **m** stands for the multi.

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

-   Single character chirality symbols:

| \@

Corresponding characters could be designated as **wchiral**, where prefix **w** stands for the whole symbol.

-   Two-character chirality symbols:

| \@\@

Corresponding characters could be designated as **schiral & echiral**, where prefix **e** stands for the start and prefix **r** stands for the end of the symbol.

-   Multicharacter (from three to four characters) chirality symbols:

| \@, T, H, A, L, S, P, B, O, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **schiral_m, nchiral_m & rchiral_m**, where prefix **s** stands for the start, **n** stands for the next from the start, prefix **r** stands for the rest of the symbol and suffix **m** stands for the multi.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 categories by their length:

-   Single character hydrogen symbols:

| H

Corresponding characters could be designated as **whydrogen**, where prefix **w** stands for the whole symbol.

-   Two-character hydrogen symbols:

| H0, H1, H2, H3, H4, H5, H6, H7, H8, H9

Corresponding characters could be designated as **shydrogen & ehydrogen**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 categories by their length:

-   Single character charge symbols:

| +, -

Corresponding characters could be designated as **wcharge**, where prefix **w** stands for the whole symbol.

-   Two-character charge obsolete symbols:

| ++, - -

Corresponding characters could be designated as **scharge_obsolete & echarge_obsolete**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

-   Multicharacter (two or three characters) charge symbols:

| +1, +2, +3, +4, +5, +6, +7, +8, +9, +10, +11, +12, +13, +14, +15,
|  -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15

Corresponding characters could be designated as **scharge_m & rcharge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi.

##### Class symbols

Class symbols designate the class of the atom, which is thing defined by the author of SMILES string.

Class symbols allowed in SMILES could be divided into 3 categories by their length:

-   Multicharacter (from 2 to 4 characters) class symbols:

| :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]

Corresponding characters could be designated as **sclass & rclass**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

As it can be seen without the further analysis, the aforementioned in-brackets categories of symbols are highly interconnected, but still it is quite easy to discriminate between the different symbols inside the square brackets, since they appear in fixed order.

### Outline

At this point it is possible to summarize the current results in the code, which could be used to represent all the mentioned here and allowed in SMILES symbol classes as sets of characters and assess the intersection between them.

``` r
library(tidyverse)
library(ggupset)

### Enumeration of the character classes allowed in SMILES
## ATOMS
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
watom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
watom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
watom_bal <- symb__atom_bal
# 05, two character atom symbol of organic aliphatic atoms
symb__atom_oal <- c('Cl', 'Br')
# Get the first character for each atom symbol
satom_oal <- map(symb__atom_oal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_oal <- map(symb__atom_oal, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 06, two character atom symbol of bracket aromatic atom
symb__atom_bar <- c('se', 'as', 'te')
# Get the first character for each atom symbol
satom_bar <- map(symb__atom_bar, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bar <- map(symb__atom_bar, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 07, two character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
                          'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                          'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                          'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
                          'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
satom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
eatom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()
## ANYTHING
# 08, single character symbol of anything
symb__anything <- c('*')
anything <- symb__anything
## SQUARE BRACKETS
# 09, single character marking the start of the atom record
bracket_start <- c('[')
# 10, single character marking the end of the atom record
bracket_end <- c(']')
## BONDS
# 11, single character symbol corresponding to the single bond
symb__single_bond <- c("-")
single_bond <- symb__single_bond
# 12, single character symbol corresponding to the double bond
symb__double_bond <- c("=")
double_bond <- symb__double_bond
# 13, single character symbol corresponding to the triple bond
symb__triple_bond <- c("#")
triple_bond <- symb__triple_bond
# 14, single character symbol corresponding to the quadruple bond
symb__quadruple_bond <- c("$")
quadruple_bond <- symb__quadruple_bond
# 15, obsolete single character symbol corresponding to the aromatic bond
symb__aromatic_bond_obsolete <- c(":")
aromatic_bond_obsolete <- symb__aromatic_bond_obsolete
# 16, single character symbol corresponding to the abscence of bond
symb__no_bond <- c(".")
no_bond <- symb__no_bond
## BOND modifiers (multipliers)
# 17, single character bond multiplying symbols initiators of branch with implicit bond
symb__wbm_ibi <- c('(')
wbm_ibi <- symb__wbm_ibi
# 18, single character bond multiplying symbols terminators of branch with implicit bond
symb__wbm_tbi <- c(')')
wbm_tbi <- symb__wbm_tbi
# 19, single character bond multiplying symbols initiators of rings with implicit bond
symb__wbm_iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_iri <- symb__wbm_iri
# 20, single character bond multiplying symbols terminators of rings with implicit bond
symb__wbm_tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
wbm_tri <- symb__wbm_tri
# 21, two-character bond multiplying symbols initiators of branching with explicit bond
symb__bm_ibe <- c('(-', '(=', '(#', '($', '(:', '(.')
sbm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 22, Two-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_2  <- c('-0', '=0', '#0', '$0', ':0', '.0',
'-1', '=1', '#1', '$1', ':1', '.1',
'-2', '=2', '#2', '$2', ':2', '.2',
'-3', '=3', '#3', '$3', ':3', '.3',
'-4', '=4', '#4', '$4', ':4', '.4',
'-5', '=5', '#5', '$5', ':5', '.5',
'-6', '=6', '#6', '$6', ':6', '.6',
'-7', '=7', '#7', '$7', ':7', '.7',
'-8', '=8', '#8', '$8', ':8', '.8',
'-9', '=9', '#9', '$9', ':9', '.9' )
sbm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 23, four-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_4  <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
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
'-%99', '=%99', '#%99', '$%99', ':%99', '.%99')
sbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
nbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 24, Three-character bond multiplying symbols initiators of rings with implicit bond
symb__bm_iri <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' )
sbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_iri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 25, two-character bond multiplying symbols terminators of branch with explicit bond
symb__tbe <- c(')-', ')=', ')#', ')$', '):', ').')
sbm_tbe   <- map(symb__tbe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tbe   <- map(symb__tbe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 26, Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_2  <- symb__bm_ire_2
sbm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ebm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 27, Four-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_4  <- symb__bm_ire_4
sbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
nbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rbm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 28, three-character bond multiplying symbols terminators of rings with implicit bond
symb__bm_tri <- symb__bm_iri
sbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rbm_tri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 29, Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__lct <- c("/", "\\\\")
lct <- symb__lct
# 30, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__rct <- c("/", "\\\\")
rct <- symb__rct
## All the things inside the square brackets besides the main atom symbol
# 31, single character isotope symbols:
symb__isotope_1 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
wisotope <- symb__isotope_1
# 32, multicharacter isotope symbols:
symb__isotope_2 <- map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__isotope_3 <- map(symb__isotope_2, \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
sisotope_m <- map(symb__isotope_3, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
risotope_m <- map(symb__isotope_3, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 33, single character chirality symbol:
symb__chiral_1 <- c('@')
wchiral <- symb__chiral_1
# 34, two-character chirality symbol:
symb__chiral_2 <- c('@@')
schiral_2 <- map(symb__chiral_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
echiral_2 <- map(symb__chiral_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 35, Multicharacter (from three to four characters) chirality symbols
symb__chir_m <- c("@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3", "@TB1", "@TB2", "@TB3", "@TB4",
                    "@TB5", "@TB6", "@TB7", "@TB8", "@TB9", "@TB10", "@TB11", "@TB12", "@TB13", "@TB14", "@TB15", "@TB16", "@TB17", "@TB18", "@TB19", "@TB20",
                    "@OH1", "@OH2", "@OH3", "@OH4", "@OH5", "@OH6", "@OH7", "@OH8", "@OH9", "@OH10", "@OH11", "@OH12", "@OH13", "@OH14", "@OH15", "@OH16", "@OH17", "@OH18", "@OH19", "@OH20",
                    "@OH21", "@OH22", "@OH23", "@OH24", "@OH25", "@OH26", "@OH27", "@OH28", "@OH29", "@OH30")
schiral_m <- map(symb__chir_m, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
nchiral_m <- map(symb__chir_m, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
rchiral_m <- map(symb__chir_m, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 36, single character hydrogen symbol:
symb__hydrogen_1 <- c("H")
whydrogen <- symb__hydrogen_1
# 37, two-character hydrogen symbol:
symb__hydrogen_2 <- c("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9")
shydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
ehydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 38, single character charge symbol:
symb__charge_1 <- c("+", "-")
wcharge <- symb__charge_1
# 39, two-character charge symbol, obsolete:
symb__charge_2 <- c("++", "--")
scharge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
echarge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 40, multicharacter charge symbol:
symb__charge_m <- c("+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15",
                     "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15")
scharge_m <- map(symb__charge_m, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rcharge_m <- map(symb__charge_m, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 41, multicharacter class symbols:
symb__class_2 <- map(c(":"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_3 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_4 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(symb__isotope_2, \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class <- c(symb__class_2, symb__class_3, symb__class_4)
sclass <- map(symb__class, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
rclass <- map(symb__class, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()

### Analyze characters from different symbol classes
character_classes <- tibble(character = c(watom_oar,
watom_oal,
watom_bar,
watom_bal,
satom_oal,
eatom_oal,
satom_bar,
eatom_bar,
satom_bal,
eatom_bal,
anything,
bracket_start,
bracket_end,
single_bond,
double_bond,
triple_bond,
quadruple_bond,
aromatic_bond_obsolete,
no_bond,
wbm_ibi,
wbm_tbi,
wbm_iri,
wbm_tri,
sbm_ibe,
ebm_ibe,
sbm_ire_2,
ebm_ire_2,
sbm_ire_4,
nbm_ire_4,
rbm_ire_4,
sbm_iri,
rbm_iri,
sbm_tbe,
ebm_tbe,
sbm_tre_2,
ebm_tre_2,
sbm_tre_4,
nbm_tre_4,
rbm_tre_4,
sbm_tri,
rbm_tri,
lct,
rct,
wisotope,
sisotope_m,
risotope_m,
wchiral,
schiral_2,
echiral_2,
schiral_m,
nchiral_m,
rchiral_m,
whydrogen,
shydrogen,
ehydrogen,
wcharge,
scharge_obsolete,
echarge_obsolete,
scharge_m,
rcharge_m,
sclass,
rclass),
class = c(rep("atom", length(watom_oar)),
rep("atom", length(watom_oal)),
rep("atom", length(watom_bar)),
rep("atom", length(watom_bal)),
rep("atom", length(satom_oal)),
rep("atom", length(eatom_oal)),
rep("atom", length(satom_bar)),
rep("atom", length(eatom_bar)),
rep("atom", length(satom_bal)),
rep("atom", length(eatom_bal)),
rep("anything", length(anything)),
rep("bracket", length(bracket_start)),
rep("bracket", length(bracket_end)),
rep("bond", length(single_bond)),
rep("bond", length(double_bond)),
rep("bond", length(triple_bond)),
rep("bond", length(quadruple_bond)),
rep("bond", length(aromatic_bond_obsolete)),
rep("bond", length(no_bond)),
rep("bond modifier", length(wbm_ibi)),
rep("bond modifier", length(wbm_tbi)),
rep("bond modifier", length(wbm_iri)),
rep("bond modifier", length(wbm_tri)),
rep("bond modifier", length(sbm_ibe)),
rep("bond modifier", length(ebm_ibe)),
rep("bond modifier", length(sbm_ire_2)),
rep("bond modifier", length(ebm_ire_2)),
rep("bond modifier", length(sbm_ire_4)),
rep("bond modifier", length(nbm_ire_4)),
rep("bond modifier", length(rbm_ire_4)),
rep("bond modifier", length(sbm_iri)),
rep("bond modifier", length(rbm_iri)),
rep("bond modifier", length(sbm_tbe)),
rep("bond modifier", length(ebm_tbe)),
rep("bond modifier", length(sbm_tre_2)),
rep("bond modifier", length(ebm_tre_2)),
rep("bond modifier", length(sbm_tre_4)),
rep("bond modifier", length(nbm_tre_4)),
rep("bond modifier", length(rbm_tre_4)),
rep("bond modifier", length(sbm_tri)),
rep("bond modifier", length(rbm_tri)),
rep("cis/trans", length(lct)),
rep("cis/trans", length(rct)),
rep("atom_property", length(wisotope)),
rep("atom_property", length(sisotope_m)),
rep("atom_property", length(risotope_m)),
rep("atom_property", length(wchiral)),
rep("atom_property", length(schiral_2)),
rep("atom_property", length(echiral_2)),
rep("atom_property", length(schiral_m)),
rep("atom_property", length(nchiral_m)),
rep("atom_property", length(rchiral_m)),
rep("atom_property", length(whydrogen)),
rep("atom_property", length(shydrogen)),
rep("atom_property", length(ehydrogen)),
rep("atom_property", length(wcharge)),
rep("atom_property", length(scharge_obsolete)),
rep("atom_property", length(echarge_obsolete)),
rep("atom_property", length(scharge_m)),
rep("atom_property", length(rcharge_m)),
rep("atom_property", length(sclass)),
rep("atom_property", length(rclass)) )
) |>
distinct() |>
group_by(character) |>
summarise(class = list(class))
# UpSet plot for classes
plot_class <- ggplot(character_classes, aes(x=class)) +
                geom_bar() +
                scale_x_upset() +
                theme_minimal()
plot_class
```

![**Figure 3.** Intersections between the character sets describing different groups of symbol classes.](https://github.com/pavelVPo/IBMC_LSFBD_smiles_parser/blob/main/ggupset_charClasses.png)

**Figure 3.** Intersections between the character sets describing different groups of symbol classes.

The results from Figure 3 will be considered in the further development.

## Pairs of characters in SMILES and in SMILES from ChEMBL (for example)

In short:

It is possible to read SMILES string character by character. In this case, each new character is the point where some questions should be answered and some decisions should be taken:

-   Is this character allowed in SMILES?

-   Is this character allowed after the previous character?

-   What does this character could mean in the current context (with the given previous characters)?

-   What are the questions, which should be addressed by the next character given this one and previous ones (in the updated context)?

Thus, it will be useful to prepare the list of all possible pairs of character classes and characters and check viability of their pairs in theory and in the real-world data to direct the further process.

### Possible pairs of characters in SMILES

To check the viability of the pairs of character classes, i.e. to answer the question:

| is it allowed in SMILES for the particular class of characters from those described earlier to be followed by the other?

The simple set of rules could be deduced from the existing SMILES description (<http://opensmiles.org/opensmiles.html>) and common sense applied.

The code to do that is provided below, it should be noted that this is a draft, which is intended to be improved after the empirical study. Also, the results depend on the order of application of rules, since the *case_when* from dplyr is used, which works this way (<https://dplyr.tidyverse.org/reference/case_when.html>). The final result will take into account the probability for each of the rules to be useful given the distribution of the character classes and characters (latter) in the available SMILES string, so, probably, the whole thing will be re-worked.

``` r
library(tidyverse)

### Enumeration of the character classes allowed in SMILES
## ATOMS
# 01, single character atom symbol of organic aromatic atom
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
w_atom_oar <- symb__atom_oar
# 02, single character atom symbol of organic aliphatic atom
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
w_atom_oal <- symb__atom_oal
# 03, single character atom symbol of bracket aromatic atom
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
w_atom_bar <- symb__atom_bar
# 04, single character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
w_atom_bal <- symb__atom_bal
# 05, two character atom symbol of organic aliphatic atoms
symb__atom_oal <- c('Cl', 'Br')
# Get the first character for each atom symbol
s_atom_oal <- map(symb__atom_oal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
e_atom_oal <- map(symb__atom_oal, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 06, two character atom symbol of bracket aromatic atom
symb__atom_bar <- c('se', 'as', 'te')
# Get the first character for each atom symbol
s_atom_bar <- map(symb__atom_bar, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
e_atom_bar <- map(symb__atom_bar, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 07, two character atom symbol of bracket aliphatic atom
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
# Get the first character for each atom symbol
s_atom_bal <- map(symb__atom_bal, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
# Get the last character for each atom symbol
e_atom_bal <- map(symb__atom_bal, \(x) str_sub(x, -1)) |> unlist() |> unique()
## ANYTHING
# 08, single character symbol of anything
symb__anything <- c('*')
anything <- symb__anything
## SQUARE BRACKETS
# 09, single character marking the start of the atom record
bracket_start <- c('[')
# 10, single character marking the end of the atom record
bracket_end <- c(']')
## BONDS
# 11, single character symbol corresponding to the single bond
symb__single_bond <- c("-")
single_bond <- symb__single_bond
# 12, single character symbol corresponding to the double bond
symb__double_bond <- c("=")
double_bond <- symb__double_bond
# 13, single character symbol corresponding to the triple bond
symb__triple_bond <- c("#")
triple_bond <- symb__triple_bond
# 14, single character symbol corresponding to the quadruple bond
symb__quadruple_bond <- c("$")
quadruple_bond <- symb__quadruple_bond
# 15, obsolete single character symbol corresponding to the aromatic bond
symb__aromatic_bond_obsolete <- c(":")
aromatic_bond_obsolete <- symb__aromatic_bond_obsolete
# 16, single character symbol corresponding to the abscence of bond
symb__no_bond <- c(".")
no_bond <- symb__no_bond
## BOND modifiers (multipliers)
# 17, single character bond multiplying symbols initiators of branch with implicit bond
symb__w_bm_ibi <- c('(')
w_bm_ibi <- symb__w_bm_ibi
# 18, single character bond multiplying symbols terminators of branch with implicit bond
symb__w_bm_tbi <- c(')')
w_bm_tbi <- symb__w_bm_tbi
# 19, single character bond multiplying symbols initiators of rings with implicit bond
symb__w_bm_iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
w_bm_iri <- symb__w_bm_iri
# 20, single character bond multiplying symbols terminators of rings with implicit bond
symb__w_bm_tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
w_bm_tri <- symb__w_bm_tri
# 21, two-character bond multiplying symbols initiators of branching with explicit bond
symb__bm_ibe <- c('(-', '(=', '(#', '($', '(:', '(.')
s_bm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_bm_ibe <- map(symb__bm_ibe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 22, Two-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_2  <- c('-0', '=0', '#0', '$0', ':0', '.0',
'-1', '=1', '#1', '$1', ':1', '.1',
'-2', '=2', '#2', '$2', ':2', '.2',
'-3', '=3', '#3', '$3', ':3', '.3',
'-4', '=4', '#4', '$4', ':4', '.4',
'-5', '=5', '#5', '$5', ':5', '.5',
'-6', '=6', '#6', '$6', ':6', '.6',
'-7', '=7', '#7', '$7', ':7', '.7',
'-8', '=8', '#8', '$8', ':8', '.8',
'-9', '=9', '#9', '$9', ':9', '.9' )
s_bm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_bm_ire_2  <- map(symb__bm_ire_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 23, four-character bond multiplying symbols initiators of rings with explicit bond
symb__bm_ire_4  <- c('-%01', '=%01', '#%01', '$%01', ':%01', '.%01',
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
'-%99', '=%99', '#%99', '$%99', ':%99', '.%99')
s_bm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
n_bm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
r_bm_ire_4  <- map(symb__bm_ire_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 24, Three-character bond multiplying symbols initiators of rings with implicit bond
symb__bm_iri <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' )
s_bm_iri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
r_bm_iri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 25, two-character bond multiplying symbols terminators of branch with explicit bond
symb__tbe <- c(')-', ')=', ')#', ')$', '):', ').')
s_bm_tbe   <- map(symb__tbe, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_bm_tbe   <- map(symb__tbe, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 26, Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_2  <- symb__bm_ire_2
s_bm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_bm_tre_2 <- map(symb__bm_tre_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 27, Four-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_4  <- symb__bm_ire_4
s_bm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
n_bm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
r_bm_tre_4 <- map(symb__bm_tre_4, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 28, three-character bond multiplying symbols terminators of rings with implicit bond
symb__bm_tri <- symb__bm_iri
s_bm_tri <- map(symb__bm_iri, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
r_bm_tri <- map(symb__bm_iri, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 29, Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__lct <- c("/", "\\\\")
lct <- symb__lct
# 30, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__rct <- c("/", "\\\\")
rct <- symb__rct
## All the things inside the square brackets besides the main atom symbol
# 31, single character isotope symbols:
symb__isotope_1 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
w_isotope <- symb__isotope_1
# 32, multicharacter isotope symbols:
symb__isotope_2 <- map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__isotope_3 <- map(symb__isotope_2, \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
s_isotope_m <- map(symb__isotope_3, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
r_isotope_m <- map(symb__isotope_3, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 33, single character chirality symbol:
symb__chiral_1 <- c('@')
w_chiral <- symb__chiral_1
# 34, two-character chirality symbol:
symb__chiral_2 <- c('@@')
s_chiral_2 <- map(symb__chiral_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_chiral_2 <- map(symb__chiral_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 35, Multicharacter (from three to four characters) chirality symbols
symb__chir_m <- c("@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3", "@TB1", "@TB2", "@TB3", "@TB4",
"@TB5", "@TB6", "@TB7", "@TB8", "@TB9", "@TB10", "@TB11", "@TB12", "@TB13", "@TB14", "@TB15", "@TB16", "@TB17", "@TB18", "@TB19", "@TB20",
"@OH1", "@OH2", "@OH3", "@OH4", "@OH5", "@OH6", "@OH7", "@OH8", "@OH9", "@OH10", "@OH11", "@OH12", "@OH13", "@OH14", "@OH15", "@OH16", "@OH17", "@OH18", "@OH19", "@OH20",
"@OH21", "@OH22", "@OH23", "@OH24", "@OH25", "@OH26", "@OH27", "@OH28", "@OH29", "@OH30")
s_chiral_m <- map(symb__chir_m, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
n_chiral_m <- map(symb__chir_m, \(x) str_sub(x, 2, 2)) |> unlist() |> unique()
r_chiral_m <- map(symb__chir_m, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 36, single character hydrogen symbol:
symb__hydrogen_1 <- c("H")
w_hydrogen <- symb__hydrogen_1
# 37, two-character hydrogen symbol:
symb__hydrogen_2 <- c("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9")
s_hydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_hydrogen <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 38, single character charge symbol:
symb__charge_1 <- c("+", "-")
w_charge <- symb__charge_1
# 39, two-character charge symbol, obsolete:
symb__charge_2 <- c("++", "--")
s_charge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
e_charge_obsolete <- map(symb__hydrogen_2, \(x) str_sub(x, -1)) |> unlist() |> unique()
# 40, three-character charge symbol:
symb__charge_m <- c("+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15",
"-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15")
s_charge_m <- map(symb__charge_m, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
r_charge_m <- map(symb__charge_m, \(x) str_sub(x, 3)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()
# 41, multicharacter class symbols:
symb__class_2 <- map(c(":"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_3 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class_4 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(symb__isotope_2, \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__class <- c(symb__class_2, symb__class_3, symb__class_4)
s_class <- map(symb__class, \(x) str_sub(x, 1, 1)) |> unlist() |> unique()
r_class <- map(symb__class, \(x) str_sub(x, 2)) |> unlist() |> unique() |> map(\(x) str_split(x, "")) |> unlist() |> unique()

### Assess all the theoretically possible pairs of symbol classes, which could be found in SMILES strings
available_classes <- c("w_atom_oar", "w_atom_oal", "w_atom_bar", "w_atom_bal", "s_atom_oal", "e_atom_oal", "s_atom_bar", "e_atom_bar", "s_atom_bal", "e_atom_bal",
"anything", "bracket_start", "bracket_end", "single_bond", "double_bond", "triple_bond", "quadruple_bond", "aromatic_bond_obsolete", "no_bond", "w_bm_ibi",
"w_bm_tbi", "w_bm_iri", "w_bm_tri", "s_bm_ibe", "e_bm_ibe", "s_bm_ire_2", "e_bm_ire_2", "s_bm_ire_4", "n_bm_ire_4", "r_bm_ire_4", "s_bm_iri", "r_bm_iri", "s_bm_tbe", "e_bm_tbe",
"s_bm_tre_2", "e_bm_tre_2", "s_bm_tre_4", "n_bm_tre_4", "r_bm_tre_4", "s_bm_tri", "r_bm_tri", "lct", "rct", "w_isotope", "s_isotope_m", "r_isotope_m", "w_chiral", "s_chiral_2",
"e_chiral_2", "s_chiral_m", "n_chiral_m", "r_chiral_m", "w_hydrogen", "s_hydrogen", "e_hydrogen", "w_charge", "s_charge_obsolete", "e_charge_obsolete", "s_charge_m",
"r_charge_m", "s_class", "r_class")
available_classes_description <- c("Single character atom symbols of organic aromatic atoms",
"Single character atom symbols of organic aliphatic atoms",
"Single character atom symbols of bracket aromatic atoms",
"Single character atom symbols of bracket aliphatic atoms",
"Two character atom symbols of organic aliphatic atoms, start",
"Two character atom symbols of organic aliphatic atoms, end",
"Two character atom symbols of bracket aromatic atoms, start",
"Two character atom symbols of bracket aromatic atoms, end",
"Two character atom symbols of bracket aliphatic atoms, start",
"Two character atom symbols of bracket aliphatic atoms, end",
"Anything",
"Square bracket, start",
"Square bracket, end",
"Single bond",
"Double bond",
"Triple bond",
"Quadruple bond",
"Aromatic bond, obsolete",
"No bond",
"Single character bond multiplying symbols initiators of branching with implicit bond",
"Single character bond multiplying symbols terminators of branching with implicit bond",
"Single character bond multiplying symbols initiators of rings with implicit bond",
"Single character bond multiplying symbols terminators of rings with implicit bond",
"Two-character bond multiplying symbols initiators of branching with explicit bond, start",
"Two-character bond multiplying symbols initiators of branching with explicit bond, end",
"Two-character bond multiplying symbols initiators of rings with explicit bond, start",
"Two-character bond multiplying symbols initiators of rings with explicit bond, end",
"Four-character bond multiplying symbols initiators of rings with explicit bond, start",
"Four-character bond multiplying symbols initiators of rings with explicit bond, next from start",
"Four-character bond multiplying symbols initiators of rings with explicit bond, rest",
"Three-character bond multiplying symbols initiators of rings with implicit bond, start",
"Three-character bond multiplying symbols initiators of rings with implicit bond, rest",
"Two-character bond multiplying symbols terminators of branching with explicit bond, start",
"Two-character bond multiplying symbols terminators of branching with explicit bond, end",
"Two-character bond multiplying symbols terminators of rings with explicit bond, start",
"Two-character bond multiplying symbols terminators of rings with explicit bond, end",
"Four-character bond multiplying symbols terminators of rings with explicit bond, start",
"Four-character bond multiplying symbols terminators of rings with explicit bond, next from start",
"Four-character bond multiplying symbols terminators of rings with explicit bond, rest",
"Three-character bond multiplying symbols terminators of rings with implicit bond, start",
"Three-character bond multiplying symbols terminators of rings with implicit bond, rest",
"cis/trans on the left side",
"cis/trans on the right side",
"Single character isotope symbols",
"Multicharacter isotope symbols, start",
"Multicharacter isotope symbols, end",
"Single character chirality symbols",
"Two-character chirality symbols, start",
"Two-character chirality symbols, end",
"Multicharacter chirality symbols, start",
"Multicharacter chirality symbols, next from start",
"Multicharacter chirality symbols, rest",
"Single character hydrogen symbols",
"Two-character hydrogen symbols, start",
"Two-character hydrogen symbols, end",
"Single character charge symbols",
"Two-character charge obsolete symbols, start",
"Two-character charge obsolete symbols, end",
"Multicharacter charge symbols, start",
"Multicharacter charge symbols, rest",
"Multicharacter class symbols, start",
"Multicharacter class symbols, rest")

class_description <- tibble(class = available_classes, description = available_classes_description)
# Get the classes' combinations
class_combs <- tibble(class_pair = map(available_classes, \(x) map(available_classes, \(y) paste(x,y, collapse = "", sep = "-"))) |>
unlist()) |>
distinct() |>
separate_wider_delim(class_pair, names =c("left_class", "right_class"), delim = "-")
# Add descriptions to the combinations
combs_described_raw <- class_combs |> inner_join(class_description, by = c("left_class" = "class")) |>
rename(left_description = description) |>
inner_join(class_description, by = c("right_class" = "class")) |>
rename(right_description = description) |>
mutate(allowed = NA, why = NA) |>
select(allowed, left_class, right_class, left_description, right_description)

### There are 3844 combinations
# Categorize combinations of character classes describing two- and multi-character symbols
# Here come the classes for two-character symbols
s_vec <- available_classes[startsWith(available_classes, "s_")]  |> unique()
e_vec <- available_classes[startsWith(available_classes, "e_")]  |> unique()
n_vec <- available_classes[startsWith(available_classes, "n_")] |> unique()
r_vec <- available_classes[startsWith(available_classes, "r_")] |> unique()
w_vec <- available_classes[startsWith(available_classes, "w_")] |> unique()
# Here come the bodies for the classes for two-character symbols
e_vec_body <- available_classes[startsWith(available_classes, "e_")] |> str_replace("e_", "") |> unique()
n_vec_body <- available_classes[startsWith(available_classes, "n_")] |> str_replace("n_", "") |> unique()
r_vec_body <- available_classes[startsWith(available_classes, "r_")] |> str_replace("r_", "") |> unique()
# Here come the classes for atoms
orgsubset_start_vec <- c(available_classes[startsWith(available_classes, "w_atom_o")], available_classes[startsWith(available_classes, "s_atom_o")]) |> unique()
singleorg_vec <- c(available_classes[startsWith(available_classes, "w_atom_o")])
singlebracket_vec <- c(available_classes[startsWith(available_classes, "w_atom_b")])
orgend_vec <- c(available_classes[startsWith(available_classes, "e_atom_o")])
bracket_vec <- c(available_classes[str_detect(available_classes, "atom_b")])
bond_vec <- c(available_classes[str_detect(available_classes, "_bond")])
bondmod_start_vec <- c(available_classes[startsWith(available_classes, "w_bm")], available_classes[startsWith(available_classes, "s_bm")])
bondmod_notstart_vec <- c(available_classes[startsWith(available_classes, "n_bm")], available_classes[startsWith(available_classes, "r_bm")], available_classes[startsWith(available_classes, "e_bm")])
internals_vec <- c("w_isotope", "s_isotope_m", "r_isotope_m", "w_chiral", "s_chiral_2", "e_chiral_2", "s_chiral_m", "n_chiral_m", "r_chiral_m", "w_hydrogen", "s_hydrogen", "e_hydrogen", "w_charge", "s_charge_obsolete", "e_charge_obsolete", "s_charge_m", "r_charge_m", "s_class", "r_class")
chiralitystart_vec <- available_classes[str_detect(available_classes, "[ws]_chiral")]
hstart_vec <- available_classes[str_detect(available_classes, "[ws]_hydrogen")]
chargestart_vec <- available_classes[str_detect(available_classes, "[ws]_charge")]
classstart_vec <- available_classes[startsWith(available_classes, "s_class")]
bracketatom_end_vec <- c("e_atom_bal", "e_atom_bar")
bracketatom_start_vec <- c("s_atom_bal", "s_atom_bar")


# Formulate and apply the rules of combinations' validity
combs_described_m <- combs_described_raw |>
mutate(left_noprefix = str_replace(left_class, "^._", "")) |>
mutate(right_noprefix = str_replace(right_class, "^._", "")) |>
rowwise() |>
mutate( allowed = case_when(
## multicharacter, allowed
# 1: s->e, same body
left_noprefix == right_noprefix & left_class %in% s_vec & right_class %in% e_vec ~ "1%correct multichar",
# 2: s->r, same body
left_noprefix == right_noprefix & ! left_noprefix %in% n_vec_body & left_class %in% s_vec & right_class %in% r_vec ~ "1%correct multichar",
# 3: s->n, same body
left_noprefix == right_noprefix & left_noprefix %in% n_vec_body & left_class %in% s_vec & right_class %in% n_vec ~ "1%correct multichar",
# 4: n->r, same body
left_noprefix == right_noprefix & left_noprefix %in% n_vec_body & left_class %in% n_vec & right_class %in% r_vec ~ "1%correct multichar",
# 4: r->r, same body
left_noprefix == right_noprefix & left_class %in% r_vec & right_class %in% r_vec ~ "1%correct multichar",
## multicharacter, not allowed
# 1: s->e, s->r, s->n, n->r not the same body
left_noprefix != right_noprefix & left_class %in% s_vec & right_class %in% e_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & ! left_noprefix %in% n_vec_body & left_class %in% s_vec & right_class %in% r_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & ! left_noprefix %in% n_vec_body & left_class %in% s_vec & right_class %in% r_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & left_noprefix %in% n_vec_body & left_class %in% n_vec & right_class %in% r_vec ~ "0%wrong multichar",
# 2: same body, wrong sequence
# - s ->(n)-> r
left_noprefix == right_noprefix & left_noprefix %in% n_vec_body & left_class %in% s_vec & right_class %in% r_vec ~ "0%wrong multichar",
# - e->s
left_noprefix == right_noprefix & left_class %in% e_vec & right_class %in% s_vec ~ "0%wrong multichar",
# - r->s
left_noprefix == right_noprefix & left_class %in% r_vec & right_class %in% s_vec ~ "0%wrong multichar",
# - n->s
left_noprefix == right_noprefix & left_class %in% n_vec & right_class %in% s_vec ~ "0%wrong multichar",
# 3: different bodies, start/intermediate to new body
left_noprefix != right_noprefix & (left_class %in% s_vec | left_class %in% n_vec) ~ "0%wrong multichar",
# 4: * bodies, one single character position to another
left_class %in% s_vec & right_class %in% s_vec ~ "0%wrong multichar",
left_class %in% n_vec & right_class %in% n_vec ~ "0%wrong multichar",
left_class %in% e_vec & right_class %in% e_vec ~ "0%wrong multichar",
# 5: * bodies, end not to start
left_class %in% e_vec & ! right_class %in% s_vec ~ "0%wrong multichar",
# 6: different bodies, one multichar to another
left_noprefix != right_noprefix & left_class %in% s_vec & right_class %in% s_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & left_class %in% e_vec & right_class %in% e_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & left_class %in% n_vec & right_class %in% n_vec ~ "0%wrong multichar",
left_noprefix != right_noprefix & left_class %in% r_vec & right_class %in% r_vec ~ "0%wrong multichar",
# 7: same bodies, one s/n to w
left_noprefix == right_noprefix & (left_class %in% s_vec | left_class %in% n_vec) & right_class %in% w_vec ~ "0%wrong multichar",

## atoms, watom_oar and watom_oal, allowed
# 1: to the start of organic subset
left_class %in% singleorg_vec & right_class %in% orgsubset_start_vec ~ "1%correct atom-atom",
# 2: to the `anything`
left_class %in% singleorg_vec & right_class == "anything" ~ "1%correct atom-anything",
# 3: to the opening bracket
left_class %in% singleorg_vec & right_class == "bracket_start" ~ "1%correct atom-atom",
# 4: to the * bond
left_class %in% singleorg_vec & right_class %in% bond_vec ~ "1%correct atom-bond",
# 5: to the start of bond modifiers
left_class %in% singleorg_vec & right_class %in% bondmod_start_vec ~ "1%correct atom-bondMod",
# 6: to cis/trans
left_class %in% singleorg_vec & right_class == "lct" ~ "1%correct atom-cis/trans",
left_class %in% singleorg_vec & right_class == "rct" ~ "1%correct atom-cis/trans",
## atoms, w_atom_oar, w_atom_oal, s_atom_oal, not allowed
# 1: to the end of organic atoms
left_class %in% orgsubset_start_vec & right_class %in% orgend_vec ~ "0%wrong atom-atom",
# 2: to any point of bracket atoms
left_class %in% orgsubset_start_vec & right_class %in% bracket_vec ~ "0%wrong atom-atom",
# 3: to the closing bracket
left_class %in% singleorg_vec & right_class == "bracket_end" ~ "0%wrong atom-atom",
# 4: to the not start of bond modifiers
left_class %in% singleorg_vec & right_class %in% bondmod_notstart_vec ~ "0%wrong atom-bondMod",
# 5: to the brackets' internals
left_class %in% singleorg_vec & right_class %in% internals_vec ~ "0%wrong atom-property",

## atoms, s_atom_oal & e_atom_oal, allowed
# 1: end to the start of organic subset
left_class == "e_atom_oal" & right_class %in% orgsubset_start_vec ~ "1%correct atom-atom",
# 2: to the `anything`
left_class == "e_atom_oal" & right_class == "anything" ~ "1%correct atom-anything",
# 3: to the opening bracket
left_class == "e_atom_oal" & right_class == "bracket_start" ~ "1%correct atom-atom",
# 4: to the * bond
left_class == "e_atom_oal" & right_class %in% bond_vec ~ "1%correct atom-bond",
# 5: to the start of bond modifiers
left_class == "e_atom_oal" & right_class %in% bondmod_start_vec ~ "1%correct atom-bondMod",
# 6: to cis/trans
left_class == "e_atom_oal" & right_class == "lct" ~ "1%correct atom-cis/trans",
left_class == "e_atom_oal" & right_class == "rct" ~ "1%correct atom-cis/trans",
# 7: fron start to finish
left_class == "s_atom_oal" & right_class == "e_atom_oal" ~ "1%correct in-atom",
## atoms, s_atom_oal & e_atom_oal, not allowed
# 1: to the end of organic atoms
left_class == "e_atom_oal" & right_class %in% orgend_vec ~ "0%wrong atom-atom",
# 2: to any point of bracket atoms
left_class == "e_atom_oal" & right_class %in% bracket_vec ~ "0%wrong atom-atom",
# 3: to the closing bracket
left_class == "e_atom_oal" & right_class == "bracket_end" ~ "0%wrong atom-atom",
# 4: to the not start of bond modifiers
left_class == "e_atom_oal" & right_class %in% bondmod_notstart_vec ~ "0%wrong atom-bondMod",
# 5: to the brackets' internals
left_class == "e_atom_oal" & right_class %in% internals_vec ~ "0%wrong atom-property",
# 6: from start not to the end
left_class == "s_atom_oal" & right_class != "e_atom_oal" ~ "1%wrong in-atom",

## atoms, watom_bar and watom_bal, allowed
left_class %in% singlebracket_vec & right_class == "bracket_end" ~ "1%correct bracket",
left_class %in% singlebracket_vec & right_class %in% chiralitystart_vec ~ "1%correct chirality",
left_class %in% singlebracket_vec & right_class %in% hstart_vec ~ "1%correct Hs",
left_class %in% singlebracket_vec & right_class %in% chargestart_vec ~ "1%correct charge",
left_class %in% singlebracket_vec & right_class %in% classstart_vec ~ "1%correct class",
## atoms, watom_bar and watom_bal, not allowed
# * else:
left_class %in% singlebracket_vec ~ "0%wrong in-bracket placement",
## atoms, satom_bar & eatom_bar, satom_bal & eatom_bal, allowed
left_class %in% bracketatom_end_vec  & right_class == "bracket_end" ~ "1%correct bracket",
left_class %in% bracketatom_end_vec & right_class %in% chiralitystart_vec ~ "1%correct chirality",
left_class %in% bracketatom_end_vec & right_class %in% hstart_vec ~ "1%correct Hs",
left_class %in% bracketatom_end_vec & right_class %in% chargestart_vec ~ "1%correct charge",
left_class %in% bracketatom_end_vec & right_class %in% classstart_vec ~ "1%correct class",
left_class %in% bracketatom_start_vec & right_class %in% bracketatom_end_vec ~ "1%correct in-atom",
## atoms, satom_bar & eatom_bar, satom_bal & eatom_bal, not allowed
left_class %in% bracketatom_start_vec ~ "0%wrong in-bracket placement",
left_class %in% bracketatom_end_vec ~ "0%wrong in-bracket placement",

## Anything, allowed
# 1: to the start of organic subset
left_class == "anything" & right_class %in% orgsubset_start_vec ~ "1%correct atom-atom",
# 2: to the `anything`
left_class == "anything" & right_class == "anything" ~ "1%correct atom-anything",
# 3: to the opening bracket
left_class == "anything" & right_class == "bracket_start" ~ "1%correct atom-atom",
# 4: to the * bond
left_class == "anything" & right_class %in% bond_vec ~ "1%correct atom-bond",
# 5: to the start of bond modifiers
left_class == "anything" & right_class %in% bondmod_start_vec ~ "1%correct atom-bondMod",
# 6: to cis/trans
left_class == "anything" & right_class == "lct" ~ "1%correct atom-cis/trans",
left_class == "anything" & right_class == "rct" ~ "1%correct atom-cis/trans",
# 7: as if in brackets
left_class == "anything" & right_class == "bracket_end" ~ "1%correct bracket",
left_class == "anything" & right_class %in% chiralitystart_vec ~ "1%correct chirality",
left_class == "anything" & right_class %in% hstart_vec ~ "1%correct Hs",
left_class == "anything" & right_class %in% chargestart_vec ~ "1%correct charge",
left_class == "anything" & right_class %in% classstart_vec ~ "1%correct class",
## Anything, not allowed
left_class == "anything" ~ "0%wrong placement",

## Square brackets, allowed
left_class == "bracket_start" & right_class %in% c("w_isotope", "s_isotope_m") ~ "1%correct bracket",
left_class == "bracket_start" & right_class %in% c("s_atom_bal", "s_atom_bar", "w_atom_bal", "w_atom_bar") ~ "1%correct bracket",
left_class == "bracket_end" & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "e_atom_oal", "anything", "bracket_start",
"single_bond", "double_bond", "triple_bond", "quadruple_bond", "aromatic_bond_obsolete",
"no_bond", "w_bm_ibi", "w_bm_tbi", "w_bm_iri", "w_bm_tri", "s_bm_ibe", "s_bm_ire_2", "s_bm_ire_4",
"s_bm_iri", "s_bm_tbe", "s_bm_tre_2", "s_bm_tre_4", "s_bm_tri", "lct", "rct") ~ "1%correct bracket",
## Square brackets, not allowed
left_class == "bracket_start"  ~ "1%wrong bracket",
left_class == "bracket_end"  ~ "1%wrong bracket",

## Bonds, allowed
left_class %in% bond_vec & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "anything", "bracket_start",
"w_bm_ibi", "w_bm_tbi", "w_bm_iri", "w_bm_tri", "s_bm_ibe", "s_bm_ire_2", "s_bm_ire_4",
"s_bm_iri", "s_bm_tbe", "s_bm_tre_2", "s_bm_tre_4", "s_bm_tri", "lct", "rct")  ~ "1%correct bracket",
## Bonds, not allowed
left_class %in% bond_vec ~ "0%wrong bond",

## End of Bond modifying (multiplying), branch iniators, allowed
left_class %in% c("w_bm_ibi", "e_bm_ibe") & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "anything", "bracket_start", "w_bm_ibi", "s_bm_ibe","lct", "rct")  ~ "1%correct branch initiator",
## Bond modifying (multiplying), branch iniator, not allowed
left_class %in% c("w_bm_ibi", "e_bm_ibe") ~ "0%wrong branch initiator",
## End of bond multiplying symbols initiators of rings, allowed
left_class %in% c("w_bm_iri", "e_bm_ire_2", "r_bm_iri", "r_bm_ire_4") & right_class %in% c("w_atom_oar", "w_atom_oal", "w_atom_bar", "w_atom_bal", "s_atom_oal", "anything", "bracket_start",
"single_bond", "double_bond", "triple_bond", "quadruple_bond", "aromatic_bond_obsolete", "no_bond",
"w_bm_ibi", "w_bm_tbi", "w_bm_iri", "w_bm_tri", "s_bm_ibe", "s_bm_ire_2",
"s_bm_ire_4", "s_bm_iri", "s_bm_tbe", "s_bm_tre_2", "s_bm_tre_4", "s_bm_tri", "lct", "rct")  ~ "1%correct ring initiator",
## End of bond modifying (multiplying), iniators of rings, not allowed
left_class %in% c("w_bm_iri", "e_bm_ire_2", "r_bm_iri", "r_bm_ire_4") ~ "0%wrong ring initiator",
## End of Bond modifying (multiplying), branch terminators, allowed
left_class %in% c("w_bm_tbi", "e_bm_tbe") & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "s_atom_bar", "s_atom_bal", "anything", "bracket_start",
"single_bond", "double_bond", "triple_bond", "quadruple_bond", "aromatic_bond_obsolete", "no_bond",
"w_bm_ibi", "w_bm_tbi", "w_bm_iri", "w_bm_tri", "s_bm_ibe", 
"s_bm_ire_2", "s_bm_ire_4", "s_bm_iri", "s_bm_tbe", "s_bm_tre_2", "s_bm_tre_4", "s_bm_tri", "lct", "rct") ~ "1% correct branch terminator",
## End of Bond modifying (multiplying), branch terminators, not allowed
left_class %in% c("w_bm_tbi", "e_bm_tbe") ~ "0%wrong branch terminator",
## End of bond multiplying symbols terminators of rings, allowed
left_class %in% c("w_bm_tri", "e_bm_tre_2", "r_bm_tri", "r_bm_tre_4") & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "anything", "bracket_start",
"single_bond", "double_bond", "triple_bond", "quadruple_bond",
"aromatic_bond_obsolete", "no_bond", "w_bm_ibi", "w_bm_tbi", "w_bm_iri", "w_bm_tri",
"s_bm_ibe", "s_bm_ire_2", "s_bm_ire_4", "s_bm_iri","s_bm_tbe", "s_bm_tre_2",
"s_bm_tre_4", "s_bm_tri", "lct", "rct") ~ "1%correct ring terminator",
## End of bond multiplying symbols terminators of rings, not allowed
left_class %in% c("w_bm_tri", "e_bm_tre_2", "r_bm_tri", "r_bm_tre_4") ~ "0%wrong branch terminator",

## Cis/Trans
left_class %in% c("lct", "rct") & right_class %in% c("w_atom_oar", "w_atom_oal", "s_atom_oal", "anything", "bracket_start") ~ "1%correct cis/trans marker",
left_class %in% c("lct", "rct") ~ "0%wrong cis/trans marker",

## Inner bracket things
# Isotope
left_class %in% c("w_isotope", "r_isotope_m") & right_class %in% c("w_atom_bar", "w_atom_bal", "s_atom_bar", "s_atom_bal") ~ "1%correct isotope",
left_class %in% c("w_isotope", "r_isotope_m") ~ "0%wrong isotope",
# Chirality
left_class %in% c("w_chiral", "e_chiral_2", "r_chiral_m") & right_class %in% c("w_hydrogen", "s_hydrogen", "w_charge", "s_charge_obsolete", "s_charge_m", "s_class") ~ "1%correct chirality",
left_class %in% c("w_chiral", "e_chiral_2", "r_chiral_m") ~ "0%wrong chirality",
# Hydrogen
left_class %in% c("w_hydrogen", "e_hydrogen") & right_class %in% c("w_charge", "s_charge_obsolete", "s_charge_m", "s_class") ~ "1%correct hydrogen",
left_class %in% c("w_hydrogen", "e_hydrogen") ~ "0%wrong hydrogen",
# Charge
left_class %in% c("w_charge", "e_charge_obsolete", "r_charge_m") & right_class %in% c("s_class") ~ "1%correct charge",
left_class %in% c("w_charge", "e_charge_obsolete", "r_charge_m") ~ "0%wrong charge",
# Class
left_class %in% c("r_class") & right_class == "bracket_end" ~ "1%correct class",
left_class %in% c("r_class") ~ "0%wrong class",


.default = as.character(allowed)
) ) |>
ungroup() |>
separate_wider_delim(allowed, names = c('allowed', 'why'), delim = "%") |>
arrange(desc(allowed))     
```

As the result, the set of rules was formulated to categorize the pairs of character classes as allowed and not allowed in SMILES, according to this rules:

| Only 652 of 3844 possible combinations of the described character classes are in line with the existing SMILES general specifications

Here is the table containing pairs of character classes allowed in SMILES according to the formulated rules:

| left_class | right_class | left_description | right_description | why |
|---------------|---------------|---------------|---------------|---------------|
| w_atom_oar | w_atom_oar | Single character atom symbols of organic aromatic atoms | Single character atom symbols of organic aromatic atoms | correct atom-atom |
| w_atom_oar | w_atom_oal | Single character atom symbols of organic aromatic atoms | Single character atom symbols of organic aliphatic atoms | correct atom-atom |
| w_atom_oar | s_atom_oal | Single character atom symbols of organic aromatic atoms | Two character atom symbols of organic aliphatic atoms, start | correct atom-atom |
| w_atom_oar | anything | Single character atom symbols of organic aromatic atoms | Anything | correct atom-anything |
| w_atom_oar | bracket_start | Single character atom symbols of organic aromatic atoms | Square bracket, start | correct atom-atom |
| w_atom_oar | single_bond | Single character atom symbols of organic aromatic atoms | Single bond | correct atom-bond |
| w_atom_oar | double_bond | Single character atom symbols of organic aromatic atoms | Double bond | correct atom-bond |
| w_atom_oar | triple_bond | Single character atom symbols of organic aromatic atoms | Triple bond | correct atom-bond |
| w_atom_oar | quadruple_bond | Single character atom symbols of organic aromatic atoms | Quadruple bond | correct atom-bond |
| w_atom_oar | aromatic_bond_obsolete | Single character atom symbols of organic aromatic atoms | Aromatic bond, obsolete | correct atom-bond |
| w_atom_oar | no_bond | Single character atom symbols of organic aromatic atoms | No bond | correct atom-bond |
| w_atom_oar | w_bm_ibi | Single character atom symbols of organic aromatic atoms | Single character bond multiplying symbols initiators of branching with implicit bond | correct atom-bondMod |
| w_atom_oar | w_bm_tbi | Single character atom symbols of organic aromatic atoms | Single character bond multiplying symbols terminators of branching with implicit bond | correct atom-bondMod |
| w_atom_oar | w_bm_iri | Single character atom symbols of organic aromatic atoms | Single character bond multiplying symbols initiators of rings with implicit bond | correct atom-bondMod |
| w_atom_oar | w_bm_tri | Single character atom symbols of organic aromatic atoms | Single character bond multiplying symbols terminators of rings with implicit bond | correct atom-bondMod |
| w_atom_oar | s_bm_ibe | Single character atom symbols of organic aromatic atoms | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_ire_2 | Single character atom symbols of organic aromatic atoms | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_ire_4 | Single character atom symbols of organic aromatic atoms | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_iri | Single character atom symbols of organic aromatic atoms | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_tbe | Single character atom symbols of organic aromatic atoms | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_tre_2 | Single character atom symbols of organic aromatic atoms | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_tre_4 | Single character atom symbols of organic aromatic atoms | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oar | s_bm_tri | Single character atom symbols of organic aromatic atoms | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct atom-bondMod |
| w_atom_oar | lct | Single character atom symbols of organic aromatic atoms | cis/trans on the left side | correct atom-cis/trans |
| w_atom_oar | rct | Single character atom symbols of organic aromatic atoms | cis/trans on the right side | correct atom-cis/trans |
| w_atom_oal | w_atom_oar | Single character atom symbols of organic aliphatic atoms | Single character atom symbols of organic aromatic atoms | correct atom-atom |
| w_atom_oal | w_atom_oal | Single character atom symbols of organic aliphatic atoms | Single character atom symbols of organic aliphatic atoms | correct atom-atom |
| w_atom_oal | s_atom_oal | Single character atom symbols of organic aliphatic atoms | Two character atom symbols of organic aliphatic atoms, start | correct atom-atom |
| w_atom_oal | anything | Single character atom symbols of organic aliphatic atoms | Anything | correct atom-anything |
| w_atom_oal | bracket_start | Single character atom symbols of organic aliphatic atoms | Square bracket, start | correct atom-atom |
| w_atom_oal | single_bond | Single character atom symbols of organic aliphatic atoms | Single bond | correct atom-bond |
| w_atom_oal | double_bond | Single character atom symbols of organic aliphatic atoms | Double bond | correct atom-bond |
| w_atom_oal | triple_bond | Single character atom symbols of organic aliphatic atoms | Triple bond | correct atom-bond |
| w_atom_oal | quadruple_bond | Single character atom symbols of organic aliphatic atoms | Quadruple bond | correct atom-bond |
| w_atom_oal | aromatic_bond_obsolete | Single character atom symbols of organic aliphatic atoms | Aromatic bond, obsolete | correct atom-bond |
| w_atom_oal | no_bond | Single character atom symbols of organic aliphatic atoms | No bond | correct atom-bond |
| w_atom_oal | w_bm_ibi | Single character atom symbols of organic aliphatic atoms | Single character bond multiplying symbols initiators of branching with implicit bond | correct atom-bondMod |
| w_atom_oal | w_bm_tbi | Single character atom symbols of organic aliphatic atoms | Single character bond multiplying symbols terminators of branching with implicit bond | correct atom-bondMod |
| w_atom_oal | w_bm_iri | Single character atom symbols of organic aliphatic atoms | Single character bond multiplying symbols initiators of rings with implicit bond | correct atom-bondMod |
| w_atom_oal | w_bm_tri | Single character atom symbols of organic aliphatic atoms | Single character bond multiplying symbols terminators of rings with implicit bond | correct atom-bondMod |
| w_atom_oal | s_bm_ibe | Single character atom symbols of organic aliphatic atoms | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_ire_2 | Single character atom symbols of organic aliphatic atoms | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_ire_4 | Single character atom symbols of organic aliphatic atoms | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_iri | Single character atom symbols of organic aliphatic atoms | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_tbe | Single character atom symbols of organic aliphatic atoms | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_tre_2 | Single character atom symbols of organic aliphatic atoms | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_tre_4 | Single character atom symbols of organic aliphatic atoms | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| w_atom_oal | s_bm_tri | Single character atom symbols of organic aliphatic atoms | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct atom-bondMod |
| w_atom_oal | lct | Single character atom symbols of organic aliphatic atoms | cis/trans on the left side | correct atom-cis/trans |
| w_atom_oal | rct | Single character atom symbols of organic aliphatic atoms | cis/trans on the right side | correct atom-cis/trans |
| w_atom_bar | bracket_end | Single character atom symbols of bracket aromatic atoms | Square bracket, end | correct bracket |
| w_atom_bar | w_chiral | Single character atom symbols of bracket aromatic atoms | Single character chirality symbols | correct chirality |
| w_atom_bar | s_chiral_2 | Single character atom symbols of bracket aromatic atoms | Two-character chirality symbols, start | correct chirality |
| w_atom_bar | s_chiral_m | Single character atom symbols of bracket aromatic atoms | Multicharacter chirality symbols, start | correct chirality |
| w_atom_bar | w_hydrogen | Single character atom symbols of bracket aromatic atoms | Single character hydrogen symbols | correct Hs |
| w_atom_bar | s_hydrogen | Single character atom symbols of bracket aromatic atoms | Two-character hydrogen symbols, start | correct Hs |
| w_atom_bar | w_charge | Single character atom symbols of bracket aromatic atoms | Single character charge symbols | correct charge |
| w_atom_bar | s_charge_obsolete | Single character atom symbols of bracket aromatic atoms | Two-character charge obsolete symbols, start | correct charge |
| w_atom_bar | s_charge_m | Single character atom symbols of bracket aromatic atoms | Multicharacter charge symbols, start | correct charge |
| w_atom_bar | s_class | Single character atom symbols of bracket aromatic atoms | Multicharacter class symbols, start | correct class |
| w_atom_bal | bracket_end | Single character atom symbols of bracket aliphatic atoms | Square bracket, end | correct bracket |
| w_atom_bal | w_chiral | Single character atom symbols of bracket aliphatic atoms | Single character chirality symbols | correct chirality |
| w_atom_bal | s_chiral_2 | Single character atom symbols of bracket aliphatic atoms | Two-character chirality symbols, start | correct chirality |
| w_atom_bal | s_chiral_m | Single character atom symbols of bracket aliphatic atoms | Multicharacter chirality symbols, start | correct chirality |
| w_atom_bal | w_hydrogen | Single character atom symbols of bracket aliphatic atoms | Single character hydrogen symbols | correct Hs |
| w_atom_bal | s_hydrogen | Single character atom symbols of bracket aliphatic atoms | Two-character hydrogen symbols, start | correct Hs |
| w_atom_bal | w_charge | Single character atom symbols of bracket aliphatic atoms | Single character charge symbols | correct charge |
| w_atom_bal | s_charge_obsolete | Single character atom symbols of bracket aliphatic atoms | Two-character charge obsolete symbols, start | correct charge |
| w_atom_bal | s_charge_m | Single character atom symbols of bracket aliphatic atoms | Multicharacter charge symbols, start | correct charge |
| w_atom_bal | s_class | Single character atom symbols of bracket aliphatic atoms | Multicharacter class symbols, start | correct class |
| s_atom_oal | e_atom_oal | Two character atom symbols of organic aliphatic atoms, start | Two character atom symbols of organic aliphatic atoms, end | correct multichar |
| e_atom_oal | s_bm_ibe | Two character atom symbols of organic aliphatic atoms, end | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_ire_2 | Two character atom symbols of organic aliphatic atoms, end | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_ire_4 | Two character atom symbols of organic aliphatic atoms, end | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_iri | Two character atom symbols of organic aliphatic atoms, end | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_tbe | Two character atom symbols of organic aliphatic atoms, end | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_tre_2 | Two character atom symbols of organic aliphatic atoms, end | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_tre_4 | Two character atom symbols of organic aliphatic atoms, end | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| e_atom_oal | s_bm_tri | Two character atom symbols of organic aliphatic atoms, end | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct atom-bondMod |
| s_atom_bar | e_atom_bar | Two character atom symbols of bracket aromatic atoms, start | Two character atom symbols of bracket aromatic atoms, end | correct multichar |
| e_atom_bar | s_chiral_2 | Two character atom symbols of bracket aromatic atoms, end | Two-character chirality symbols, start | correct chirality |
| e_atom_bar | s_chiral_m | Two character atom symbols of bracket aromatic atoms, end | Multicharacter chirality symbols, start | correct chirality |
| e_atom_bar | s_hydrogen | Two character atom symbols of bracket aromatic atoms, end | Two-character hydrogen symbols, start | correct Hs |
| e_atom_bar | s_charge_obsolete | Two character atom symbols of bracket aromatic atoms, end | Two-character charge obsolete symbols, start | correct charge |
| e_atom_bar | s_charge_m | Two character atom symbols of bracket aromatic atoms, end | Multicharacter charge symbols, start | correct charge |
| e_atom_bar | s_class | Two character atom symbols of bracket aromatic atoms, end | Multicharacter class symbols, start | correct class |
| s_atom_bal | e_atom_bal | Two character atom symbols of bracket aliphatic atoms, start | Two character atom symbols of bracket aliphatic atoms, end | correct multichar |
| e_atom_bal | s_chiral_2 | Two character atom symbols of bracket aliphatic atoms, end | Two-character chirality symbols, start | correct chirality |
| e_atom_bal | s_chiral_m | Two character atom symbols of bracket aliphatic atoms, end | Multicharacter chirality symbols, start | correct chirality |
| e_atom_bal | s_hydrogen | Two character atom symbols of bracket aliphatic atoms, end | Two-character hydrogen symbols, start | correct Hs |
| e_atom_bal | s_charge_obsolete | Two character atom symbols of bracket aliphatic atoms, end | Two-character charge obsolete symbols, start | correct charge |
| e_atom_bal | s_charge_m | Two character atom symbols of bracket aliphatic atoms, end | Multicharacter charge symbols, start | correct charge |
| e_atom_bal | s_class | Two character atom symbols of bracket aliphatic atoms, end | Multicharacter class symbols, start | correct class |
| anything | w_atom_oar | Anything | Single character atom symbols of organic aromatic atoms | correct atom-atom |
| anything | w_atom_oal | Anything | Single character atom symbols of organic aliphatic atoms | correct atom-atom |
| anything | s_atom_oal | Anything | Two character atom symbols of organic aliphatic atoms, start | correct atom-atom |
| anything | anything | Anything | Anything | correct atom-anything |
| anything | bracket_start | Anything | Square bracket, start | correct atom-atom |
| anything | bracket_end | Anything | Square bracket, end | correct bracket |
| anything | single_bond | Anything | Single bond | correct atom-bond |
| anything | double_bond | Anything | Double bond | correct atom-bond |
| anything | triple_bond | Anything | Triple bond | correct atom-bond |
| anything | quadruple_bond | Anything | Quadruple bond | correct atom-bond |
| anything | aromatic_bond_obsolete | Anything | Aromatic bond, obsolete | correct atom-bond |
| anything | no_bond | Anything | No bond | correct atom-bond |
| anything | w_bm_ibi | Anything | Single character bond multiplying symbols initiators of branching with implicit bond | correct atom-bondMod |
| anything | w_bm_tbi | Anything | Single character bond multiplying symbols terminators of branching with implicit bond | correct atom-bondMod |
| anything | w_bm_iri | Anything | Single character bond multiplying symbols initiators of rings with implicit bond | correct atom-bondMod |
| anything | w_bm_tri | Anything | Single character bond multiplying symbols terminators of rings with implicit bond | correct atom-bondMod |
| anything | s_bm_ibe | Anything | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct atom-bondMod |
| anything | s_bm_ire_2 | Anything | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| anything | s_bm_ire_4 | Anything | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct atom-bondMod |
| anything | s_bm_iri | Anything | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct atom-bondMod |
| anything | s_bm_tbe | Anything | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct atom-bondMod |
| anything | s_bm_tre_2 | Anything | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| anything | s_bm_tre_4 | Anything | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct atom-bondMod |
| anything | s_bm_tri | Anything | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct atom-bondMod |
| anything | lct | Anything | cis/trans on the left side | correct atom-cis/trans |
| anything | rct | Anything | cis/trans on the right side | correct atom-cis/trans |
| anything | w_chiral | Anything | Single character chirality symbols | correct chirality |
| anything | s_chiral_2 | Anything | Two-character chirality symbols, start | correct chirality |
| anything | s_chiral_m | Anything | Multicharacter chirality symbols, start | correct chirality |
| anything | w_hydrogen | Anything | Single character hydrogen symbols | correct Hs |
| anything | s_hydrogen | Anything | Two-character hydrogen symbols, start | correct Hs |
| anything | w_charge | Anything | Single character charge symbols | correct charge |
| anything | s_charge_obsolete | Anything | Two-character charge obsolete symbols, start | correct charge |
| anything | s_charge_m | Anything | Multicharacter charge symbols, start | correct charge |
| anything | s_class | Anything | Multicharacter class symbols, start | correct class |
| bracket_start | w_atom_oar | Square bracket, start | Single character atom symbols of organic aromatic atoms | wrong bracket |
| bracket_start | w_atom_oal | Square bracket, start | Single character atom symbols of organic aliphatic atoms | wrong bracket |
| bracket_start | w_atom_bar | Square bracket, start | Single character atom symbols of bracket aromatic atoms | correct bracket |
| bracket_start | w_atom_bal | Square bracket, start | Single character atom symbols of bracket aliphatic atoms | correct bracket |
| bracket_start | s_atom_oal | Square bracket, start | Two character atom symbols of organic aliphatic atoms, start | wrong bracket |
| bracket_start | e_atom_oal | Square bracket, start | Two character atom symbols of organic aliphatic atoms, end | wrong bracket |
| bracket_start | s_atom_bar | Square bracket, start | Two character atom symbols of bracket aromatic atoms, start | correct bracket |
| bracket_start | e_atom_bar | Square bracket, start | Two character atom symbols of bracket aromatic atoms, end | wrong bracket |
| bracket_start | s_atom_bal | Square bracket, start | Two character atom symbols of bracket aliphatic atoms, start | correct bracket |
| bracket_start | e_atom_bal | Square bracket, start | Two character atom symbols of bracket aliphatic atoms, end | wrong bracket |
| bracket_start | anything | Square bracket, start | Anything | wrong bracket |
| bracket_start | bracket_start | Square bracket, start | Square bracket, start | wrong bracket |
| bracket_start | bracket_end | Square bracket, start | Square bracket, end | wrong bracket |
| bracket_start | single_bond | Square bracket, start | Single bond | wrong bracket |
| bracket_start | double_bond | Square bracket, start | Double bond | wrong bracket |
| bracket_start | triple_bond | Square bracket, start | Triple bond | wrong bracket |
| bracket_start | quadruple_bond | Square bracket, start | Quadruple bond | wrong bracket |
| bracket_start | aromatic_bond_obsolete | Square bracket, start | Aromatic bond, obsolete | wrong bracket |
| bracket_start | no_bond | Square bracket, start | No bond | wrong bracket |
| bracket_start | w_bm_ibi | Square bracket, start | Single character bond multiplying symbols initiators of branching with implicit bond | wrong bracket |
| bracket_start | w_bm_tbi | Square bracket, start | Single character bond multiplying symbols terminators of branching with implicit bond | wrong bracket |
| bracket_start | w_bm_iri | Square bracket, start | Single character bond multiplying symbols initiators of rings with implicit bond | wrong bracket |
| bracket_start | w_bm_tri | Square bracket, start | Single character bond multiplying symbols terminators of rings with implicit bond | wrong bracket |
| bracket_start | s_bm_ibe | Square bracket, start | Two-character bond multiplying symbols initiators of branching with explicit bond, start | wrong bracket |
| bracket_start | e_bm_ibe | Square bracket, start | Two-character bond multiplying symbols initiators of branching with explicit bond, end | wrong bracket |
| bracket_start | s_bm_ire_2 | Square bracket, start | Two-character bond multiplying symbols initiators of rings with explicit bond, start | wrong bracket |
| bracket_start | e_bm_ire_2 | Square bracket, start | Two-character bond multiplying symbols initiators of rings with explicit bond, end | wrong bracket |
| bracket_start | s_bm_ire_4 | Square bracket, start | Four-character bond multiplying symbols initiators of rings with explicit bond, start | wrong bracket |
| bracket_start | n_bm_ire_4 | Square bracket, start | Four-character bond multiplying symbols initiators of rings with explicit bond, next from start | wrong bracket |
| bracket_start | r_bm_ire_4 | Square bracket, start | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | wrong bracket |
| bracket_start | s_bm_iri | Square bracket, start | Three-character bond multiplying symbols initiators of rings with implicit bond, start | wrong bracket |
| bracket_start | r_bm_iri | Square bracket, start | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | wrong bracket |
| bracket_start | s_bm_tbe | Square bracket, start | Two-character bond multiplying symbols terminators of branching with explicit bond, start | wrong bracket |
| bracket_start | e_bm_tbe | Square bracket, start | Two-character bond multiplying symbols terminators of branching with explicit bond, end | wrong bracket |
| bracket_start | s_bm_tre_2 | Square bracket, start | Two-character bond multiplying symbols terminators of rings with explicit bond, start | wrong bracket |
| bracket_start | e_bm_tre_2 | Square bracket, start | Two-character bond multiplying symbols terminators of rings with explicit bond, end | wrong bracket |
| bracket_start | s_bm_tre_4 | Square bracket, start | Four-character bond multiplying symbols terminators of rings with explicit bond, start | wrong bracket |
| bracket_start | n_bm_tre_4 | Square bracket, start | Four-character bond multiplying symbols terminators of rings with explicit bond, next from start | wrong bracket |
| bracket_start | r_bm_tre_4 | Square bracket, start | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | wrong bracket |
| bracket_start | s_bm_tri | Square bracket, start | Three-character bond multiplying symbols terminators of rings with implicit bond, start | wrong bracket |
| bracket_start | r_bm_tri | Square bracket, start | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | wrong bracket |
| bracket_start | lct | Square bracket, start | cis/trans on the left side | wrong bracket |
| bracket_start | rct | Square bracket, start | cis/trans on the right side | wrong bracket |
| bracket_start | w_isotope | Square bracket, start | Single character isotope symbols | correct bracket |
| bracket_start | s_isotope_m | Square bracket, start | Multicharacter isotope symbols, start | correct bracket |
| bracket_start | r_isotope_m | Square bracket, start | Multicharacter isotope symbols, end | wrong bracket |
| bracket_start | w_chiral | Square bracket, start | Single character chirality symbols | wrong bracket |
| bracket_start | s_chiral_2 | Square bracket, start | Two-character chirality symbols, start | wrong bracket |
| bracket_start | e_chiral_2 | Square bracket, start | Two-character chirality symbols, end | wrong bracket |
| bracket_start | s_chiral_m | Square bracket, start | Multicharacter chirality symbols, start | wrong bracket |
| bracket_start | n_chiral_m | Square bracket, start | Multicharacter chirality symbols, next from start | wrong bracket |
| bracket_start | r_chiral_m | Square bracket, start | Multicharacter chirality symbols, rest | wrong bracket |
| bracket_start | w_hydrogen | Square bracket, start | Single character hydrogen symbols | wrong bracket |
| bracket_start | s_hydrogen | Square bracket, start | Two-character hydrogen symbols, start | wrong bracket |
| bracket_start | e_hydrogen | Square bracket, start | Two-character hydrogen symbols, end | wrong bracket |
| bracket_start | w_charge | Square bracket, start | Single character charge symbols | wrong bracket |
| bracket_start | s_charge_obsolete | Square bracket, start | Two-character charge obsolete symbols, start | wrong bracket |
| bracket_start | e_charge_obsolete | Square bracket, start | Two-character charge obsolete symbols, end | wrong bracket |
| bracket_start | s_charge_m | Square bracket, start | Multicharacter charge symbols, start | wrong bracket |
| bracket_start | r_charge_m | Square bracket, start | Multicharacter charge symbols, rest | wrong bracket |
| bracket_start | s_class | Square bracket, start | Multicharacter class symbols, start | wrong bracket |
| bracket_start | r_class | Square bracket, start | Multicharacter class symbols, rest | wrong bracket |
| bracket_end | w_atom_oar | Square bracket, end | Single character atom symbols of organic aromatic atoms | correct bracket |
| bracket_end | w_atom_oal | Square bracket, end | Single character atom symbols of organic aliphatic atoms | correct bracket |
| bracket_end | w_atom_bar | Square bracket, end | Single character atom symbols of bracket aromatic atoms | wrong bracket |
| bracket_end | w_atom_bal | Square bracket, end | Single character atom symbols of bracket aliphatic atoms | wrong bracket |
| bracket_end | s_atom_oal | Square bracket, end | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| bracket_end | e_atom_oal | Square bracket, end | Two character atom symbols of organic aliphatic atoms, end | correct bracket |
| bracket_end | s_atom_bar | Square bracket, end | Two character atom symbols of bracket aromatic atoms, start | wrong bracket |
| bracket_end | e_atom_bar | Square bracket, end | Two character atom symbols of bracket aromatic atoms, end | wrong bracket |
| bracket_end | s_atom_bal | Square bracket, end | Two character atom symbols of bracket aliphatic atoms, start | wrong bracket |
| bracket_end | e_atom_bal | Square bracket, end | Two character atom symbols of bracket aliphatic atoms, end | wrong bracket |
| bracket_end | anything | Square bracket, end | Anything | correct bracket |
| bracket_end | bracket_start | Square bracket, end | Square bracket, start | correct bracket |
| bracket_end | bracket_end | Square bracket, end | Square bracket, end | wrong bracket |
| bracket_end | single_bond | Square bracket, end | Single bond | correct bracket |
| bracket_end | double_bond | Square bracket, end | Double bond | correct bracket |
| bracket_end | triple_bond | Square bracket, end | Triple bond | correct bracket |
| bracket_end | quadruple_bond | Square bracket, end | Quadruple bond | correct bracket |
| bracket_end | aromatic_bond_obsolete | Square bracket, end | Aromatic bond, obsolete | correct bracket |
| bracket_end | no_bond | Square bracket, end | No bond | correct bracket |
| bracket_end | w_bm_ibi | Square bracket, end | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| bracket_end | w_bm_tbi | Square bracket, end | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| bracket_end | w_bm_iri | Square bracket, end | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| bracket_end | w_bm_tri | Square bracket, end | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| bracket_end | s_bm_ibe | Square bracket, end | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| bracket_end | e_bm_ibe | Square bracket, end | Two-character bond multiplying symbols initiators of branching with explicit bond, end | wrong bracket |
| bracket_end | s_bm_ire_2 | Square bracket, end | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| bracket_end | e_bm_ire_2 | Square bracket, end | Two-character bond multiplying symbols initiators of rings with explicit bond, end | wrong bracket |
| bracket_end | s_bm_ire_4 | Square bracket, end | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| bracket_end | n_bm_ire_4 | Square bracket, end | Four-character bond multiplying symbols initiators of rings with explicit bond, next from start | wrong bracket |
| bracket_end | r_bm_ire_4 | Square bracket, end | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | wrong bracket |
| bracket_end | s_bm_iri | Square bracket, end | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| bracket_end | r_bm_iri | Square bracket, end | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | wrong bracket |
| bracket_end | s_bm_tbe | Square bracket, end | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| bracket_end | e_bm_tbe | Square bracket, end | Two-character bond multiplying symbols terminators of branching with explicit bond, end | wrong bracket |
| bracket_end | s_bm_tre_2 | Square bracket, end | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| bracket_end | e_bm_tre_2 | Square bracket, end | Two-character bond multiplying symbols terminators of rings with explicit bond, end | wrong bracket |
| bracket_end | s_bm_tre_4 | Square bracket, end | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| bracket_end | n_bm_tre_4 | Square bracket, end | Four-character bond multiplying symbols terminators of rings with explicit bond, next from start | wrong bracket |
| bracket_end | r_bm_tre_4 | Square bracket, end | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | wrong bracket |
| bracket_end | s_bm_tri | Square bracket, end | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| bracket_end | r_bm_tri | Square bracket, end | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | wrong bracket |
| bracket_end | lct | Square bracket, end | cis/trans on the left side | correct bracket |
| bracket_end | rct | Square bracket, end | cis/trans on the right side | correct bracket |
| bracket_end | w_isotope | Square bracket, end | Single character isotope symbols | wrong bracket |
| bracket_end | s_isotope_m | Square bracket, end | Multicharacter isotope symbols, start | wrong bracket |
| bracket_end | r_isotope_m | Square bracket, end | Multicharacter isotope symbols, end | wrong bracket |
| bracket_end | w_chiral | Square bracket, end | Single character chirality symbols | wrong bracket |
| bracket_end | s_chiral_2 | Square bracket, end | Two-character chirality symbols, start | wrong bracket |
| bracket_end | e_chiral_2 | Square bracket, end | Two-character chirality symbols, end | wrong bracket |
| bracket_end | s_chiral_m | Square bracket, end | Multicharacter chirality symbols, start | wrong bracket |
| bracket_end | n_chiral_m | Square bracket, end | Multicharacter chirality symbols, next from start | wrong bracket |
| bracket_end | r_chiral_m | Square bracket, end | Multicharacter chirality symbols, rest | wrong bracket |
| bracket_end | w_hydrogen | Square bracket, end | Single character hydrogen symbols | wrong bracket |
| bracket_end | s_hydrogen | Square bracket, end | Two-character hydrogen symbols, start | wrong bracket |
| bracket_end | e_hydrogen | Square bracket, end | Two-character hydrogen symbols, end | wrong bracket |
| bracket_end | w_charge | Square bracket, end | Single character charge symbols | wrong bracket |
| bracket_end | s_charge_obsolete | Square bracket, end | Two-character charge obsolete symbols, start | wrong bracket |
| bracket_end | e_charge_obsolete | Square bracket, end | Two-character charge obsolete symbols, end | wrong bracket |
| bracket_end | s_charge_m | Square bracket, end | Multicharacter charge symbols, start | wrong bracket |
| bracket_end | r_charge_m | Square bracket, end | Multicharacter charge symbols, rest | wrong bracket |
| bracket_end | s_class | Square bracket, end | Multicharacter class symbols, start | wrong bracket |
| bracket_end | r_class | Square bracket, end | Multicharacter class symbols, rest | wrong bracket |
| single_bond | w_atom_oar | Single bond | Single character atom symbols of organic aromatic atoms | correct bracket |
| single_bond | w_atom_oal | Single bond | Single character atom symbols of organic aliphatic atoms | correct bracket |
| single_bond | s_atom_oal | Single bond | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| single_bond | anything | Single bond | Anything | correct bracket |
| single_bond | bracket_start | Single bond | Square bracket, start | correct bracket |
| single_bond | w_bm_ibi | Single bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| single_bond | w_bm_tbi | Single bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| single_bond | w_bm_iri | Single bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| single_bond | w_bm_tri | Single bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| single_bond | s_bm_ibe | Single bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| single_bond | s_bm_ire_2 | Single bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| single_bond | s_bm_ire_4 | Single bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| single_bond | s_bm_iri | Single bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| single_bond | s_bm_tbe | Single bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| single_bond | s_bm_tre_2 | Single bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| single_bond | s_bm_tre_4 | Single bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| single_bond | s_bm_tri | Single bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| single_bond | lct | Single bond | cis/trans on the left side | correct bracket |
| single_bond | rct | Single bond | cis/trans on the right side | correct bracket |
| double_bond | w_atom_oar | Double bond | Single character atom symbols of organic aromatic atoms | correct bracket |
| double_bond | w_atom_oal | Double bond | Single character atom symbols of organic aliphatic atoms | correct bracket |
| double_bond | s_atom_oal | Double bond | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| double_bond | anything | Double bond | Anything | correct bracket |
| double_bond | bracket_start | Double bond | Square bracket, start | correct bracket |
| double_bond | w_bm_ibi | Double bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| double_bond | w_bm_tbi | Double bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| double_bond | w_bm_iri | Double bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| double_bond | w_bm_tri | Double bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| double_bond | s_bm_ibe | Double bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| double_bond | s_bm_ire_2 | Double bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| double_bond | s_bm_ire_4 | Double bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| double_bond | s_bm_iri | Double bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| double_bond | s_bm_tbe | Double bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| double_bond | s_bm_tre_2 | Double bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| double_bond | s_bm_tre_4 | Double bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| double_bond | s_bm_tri | Double bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| double_bond | lct | Double bond | cis/trans on the left side | correct bracket |
| double_bond | rct | Double bond | cis/trans on the right side | correct bracket |
| triple_bond | w_atom_oar | Triple bond | Single character atom symbols of organic aromatic atoms | correct bracket |
| triple_bond | w_atom_oal | Triple bond | Single character atom symbols of organic aliphatic atoms | correct bracket |
| triple_bond | s_atom_oal | Triple bond | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| triple_bond | anything | Triple bond | Anything | correct bracket |
| triple_bond | bracket_start | Triple bond | Square bracket, start | correct bracket |
| triple_bond | w_bm_ibi | Triple bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| triple_bond | w_bm_tbi | Triple bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| triple_bond | w_bm_iri | Triple bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| triple_bond | w_bm_tri | Triple bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| triple_bond | s_bm_ibe | Triple bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| triple_bond | s_bm_ire_2 | Triple bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| triple_bond | s_bm_ire_4 | Triple bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| triple_bond | s_bm_iri | Triple bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| triple_bond | s_bm_tbe | Triple bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| triple_bond | s_bm_tre_2 | Triple bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| triple_bond | s_bm_tre_4 | Triple bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| triple_bond | s_bm_tri | Triple bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| triple_bond | lct | Triple bond | cis/trans on the left side | correct bracket |
| triple_bond | rct | Triple bond | cis/trans on the right side | correct bracket |
| quadruple_bond | w_atom_oar | Quadruple bond | Single character atom symbols of organic aromatic atoms | correct bracket |
| quadruple_bond | w_atom_oal | Quadruple bond | Single character atom symbols of organic aliphatic atoms | correct bracket |
| quadruple_bond | s_atom_oal | Quadruple bond | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| quadruple_bond | anything | Quadruple bond | Anything | correct bracket |
| quadruple_bond | bracket_start | Quadruple bond | Square bracket, start | correct bracket |
| quadruple_bond | w_bm_ibi | Quadruple bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| quadruple_bond | w_bm_tbi | Quadruple bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| quadruple_bond | w_bm_iri | Quadruple bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| quadruple_bond | w_bm_tri | Quadruple bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| quadruple_bond | s_bm_ibe | Quadruple bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_ire_2 | Quadruple bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_ire_4 | Quadruple bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_iri | Quadruple bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| quadruple_bond | s_bm_tbe | Quadruple bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_tre_2 | Quadruple bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_tre_4 | Quadruple bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| quadruple_bond | s_bm_tri | Quadruple bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| quadruple_bond | lct | Quadruple bond | cis/trans on the left side | correct bracket |
| quadruple_bond | rct | Quadruple bond | cis/trans on the right side | correct bracket |
| aromatic_bond_obsolete | w_atom_oar | Aromatic bond, obsolete | Single character atom symbols of organic aromatic atoms | correct bracket |
| aromatic_bond_obsolete | w_atom_oal | Aromatic bond, obsolete | Single character atom symbols of organic aliphatic atoms | correct bracket |
| aromatic_bond_obsolete | s_atom_oal | Aromatic bond, obsolete | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| aromatic_bond_obsolete | anything | Aromatic bond, obsolete | Anything | correct bracket |
| aromatic_bond_obsolete | bracket_start | Aromatic bond, obsolete | Square bracket, start | correct bracket |
| aromatic_bond_obsolete | w_bm_ibi | Aromatic bond, obsolete | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| aromatic_bond_obsolete | w_bm_tbi | Aromatic bond, obsolete | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| aromatic_bond_obsolete | w_bm_iri | Aromatic bond, obsolete | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| aromatic_bond_obsolete | w_bm_tri | Aromatic bond, obsolete | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| aromatic_bond_obsolete | s_bm_ibe | Aromatic bond, obsolete | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_ire_2 | Aromatic bond, obsolete | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_ire_4 | Aromatic bond, obsolete | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_iri | Aromatic bond, obsolete | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_tbe | Aromatic bond, obsolete | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_tre_2 | Aromatic bond, obsolete | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_tre_4 | Aromatic bond, obsolete | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| aromatic_bond_obsolete | s_bm_tri | Aromatic bond, obsolete | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| aromatic_bond_obsolete | lct | Aromatic bond, obsolete | cis/trans on the left side | correct bracket |
| aromatic_bond_obsolete | rct | Aromatic bond, obsolete | cis/trans on the right side | correct bracket |
| no_bond | w_atom_oar | No bond | Single character atom symbols of organic aromatic atoms | correct bracket |
| no_bond | w_atom_oal | No bond | Single character atom symbols of organic aliphatic atoms | correct bracket |
| no_bond | s_atom_oal | No bond | Two character atom symbols of organic aliphatic atoms, start | correct bracket |
| no_bond | anything | No bond | Anything | correct bracket |
| no_bond | bracket_start | No bond | Square bracket, start | correct bracket |
| no_bond | w_bm_ibi | No bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct bracket |
| no_bond | w_bm_tbi | No bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct bracket |
| no_bond | w_bm_iri | No bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct bracket |
| no_bond | w_bm_tri | No bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct bracket |
| no_bond | s_bm_ibe | No bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct bracket |
| no_bond | s_bm_ire_2 | No bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| no_bond | s_bm_ire_4 | No bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct bracket |
| no_bond | s_bm_iri | No bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct bracket |
| no_bond | s_bm_tbe | No bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct bracket |
| no_bond | s_bm_tre_2 | No bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| no_bond | s_bm_tre_4 | No bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct bracket |
| no_bond | s_bm_tri | No bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct bracket |
| no_bond | lct | No bond | cis/trans on the left side | correct bracket |
| no_bond | rct | No bond | cis/trans on the right side | correct bracket |
| w_bm_ibi | w_atom_oar | Single character bond multiplying symbols initiators of branching with implicit bond | Single character atom symbols of organic aromatic atoms | correct branch initiator |
| w_bm_ibi | w_atom_oal | Single character bond multiplying symbols initiators of branching with implicit bond | Single character atom symbols of organic aliphatic atoms | correct branch initiator |
| w_bm_ibi | s_atom_oal | Single character bond multiplying symbols initiators of branching with implicit bond | Two character atom symbols of organic aliphatic atoms, start | correct branch initiator |
| w_bm_ibi | anything | Single character bond multiplying symbols initiators of branching with implicit bond | Anything | correct branch initiator |
| w_bm_ibi | bracket_start | Single character bond multiplying symbols initiators of branching with implicit bond | Square bracket, start | correct branch initiator |
| w_bm_ibi | w_bm_ibi | Single character bond multiplying symbols initiators of branching with implicit bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct branch initiator |
| w_bm_ibi | s_bm_ibe | Single character bond multiplying symbols initiators of branching with implicit bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct branch initiator |
| w_bm_ibi | lct | Single character bond multiplying symbols initiators of branching with implicit bond | cis/trans on the left side | correct branch initiator |
| w_bm_ibi | rct | Single character bond multiplying symbols initiators of branching with implicit bond | cis/trans on the right side | correct branch initiator |
| w_bm_tbi | w_atom_oar | Single character bond multiplying symbols terminators of branching with implicit bond | Single character atom symbols of organic aromatic atoms |  correct branch terminator |
| w_bm_tbi | w_atom_oal | Single character bond multiplying symbols terminators of branching with implicit bond | Single character atom symbols of organic aliphatic atoms |  correct branch terminator |
| w_bm_tbi | s_atom_oal | Single character bond multiplying symbols terminators of branching with implicit bond | Two character atom symbols of organic aliphatic atoms, start |  correct branch terminator |
| w_bm_tbi | s_atom_bar | Single character bond multiplying symbols terminators of branching with implicit bond | Two character atom symbols of bracket aromatic atoms, start |  correct branch terminator |
| w_bm_tbi | s_atom_bal | Single character bond multiplying symbols terminators of branching with implicit bond | Two character atom symbols of bracket aliphatic atoms, start |  correct branch terminator |
| w_bm_tbi | anything | Single character bond multiplying symbols terminators of branching with implicit bond | Anything |  correct branch terminator |
| w_bm_tbi | bracket_start | Single character bond multiplying symbols terminators of branching with implicit bond | Square bracket, start |  correct branch terminator |
| w_bm_tbi | single_bond | Single character bond multiplying symbols terminators of branching with implicit bond | Single bond |  correct branch terminator |
| w_bm_tbi | double_bond | Single character bond multiplying symbols terminators of branching with implicit bond | Double bond |  correct branch terminator |
| w_bm_tbi | triple_bond | Single character bond multiplying symbols terminators of branching with implicit bond | Triple bond |  correct branch terminator |
| w_bm_tbi | quadruple_bond | Single character bond multiplying symbols terminators of branching with implicit bond | Quadruple bond |  correct branch terminator |
| w_bm_tbi | aromatic_bond_obsolete | Single character bond multiplying symbols terminators of branching with implicit bond | Aromatic bond, obsolete |  correct branch terminator |
| w_bm_tbi | no_bond | Single character bond multiplying symbols terminators of branching with implicit bond | No bond |  correct branch terminator |
| w_bm_tbi | w_bm_ibi | Single character bond multiplying symbols terminators of branching with implicit bond | Single character bond multiplying symbols initiators of branching with implicit bond |  correct branch terminator |
| w_bm_tbi | w_bm_tbi | Single character bond multiplying symbols terminators of branching with implicit bond | Single character bond multiplying symbols terminators of branching with implicit bond |  correct branch terminator |
| w_bm_tbi | w_bm_iri | Single character bond multiplying symbols terminators of branching with implicit bond | Single character bond multiplying symbols initiators of rings with implicit bond |  correct branch terminator |
| w_bm_tbi | w_bm_tri | Single character bond multiplying symbols terminators of branching with implicit bond | Single character bond multiplying symbols terminators of rings with implicit bond |  correct branch terminator |
| w_bm_tbi | s_bm_ibe | Single character bond multiplying symbols terminators of branching with implicit bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_ire_2 | Single character bond multiplying symbols terminators of branching with implicit bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_ire_4 | Single character bond multiplying symbols terminators of branching with implicit bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_iri | Single character bond multiplying symbols terminators of branching with implicit bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_tbe | Single character bond multiplying symbols terminators of branching with implicit bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_tre_2 | Single character bond multiplying symbols terminators of branching with implicit bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_tre_4 | Single character bond multiplying symbols terminators of branching with implicit bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start |  correct branch terminator |
| w_bm_tbi | s_bm_tri | Single character bond multiplying symbols terminators of branching with implicit bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start |  correct branch terminator |
| w_bm_tbi | lct | Single character bond multiplying symbols terminators of branching with implicit bond | cis/trans on the left side |  correct branch terminator |
| w_bm_tbi | rct | Single character bond multiplying symbols terminators of branching with implicit bond | cis/trans on the right side |  correct branch terminator |
| w_bm_iri | w_atom_oar | Single character bond multiplying symbols initiators of rings with implicit bond | Single character atom symbols of organic aromatic atoms | correct ring initiator |
| w_bm_iri | w_atom_oal | Single character bond multiplying symbols initiators of rings with implicit bond | Single character atom symbols of organic aliphatic atoms | correct ring initiator |
| w_bm_iri | w_atom_bar | Single character bond multiplying symbols initiators of rings with implicit bond | Single character atom symbols of bracket aromatic atoms | correct ring initiator |
| w_bm_iri | w_atom_bal | Single character bond multiplying symbols initiators of rings with implicit bond | Single character atom symbols of bracket aliphatic atoms | correct ring initiator |
| w_bm_iri | s_atom_oal | Single character bond multiplying symbols initiators of rings with implicit bond | Two character atom symbols of organic aliphatic atoms, start | correct ring initiator |
| w_bm_iri | anything | Single character bond multiplying symbols initiators of rings with implicit bond | Anything | correct ring initiator |
| w_bm_iri | bracket_start | Single character bond multiplying symbols initiators of rings with implicit bond | Square bracket, start | correct ring initiator |
| w_bm_iri | single_bond | Single character bond multiplying symbols initiators of rings with implicit bond | Single bond | correct ring initiator |
| w_bm_iri | double_bond | Single character bond multiplying symbols initiators of rings with implicit bond | Double bond | correct ring initiator |
| w_bm_iri | triple_bond | Single character bond multiplying symbols initiators of rings with implicit bond | Triple bond | correct ring initiator |
| w_bm_iri | quadruple_bond | Single character bond multiplying symbols initiators of rings with implicit bond | Quadruple bond | correct ring initiator |
| w_bm_iri | aromatic_bond_obsolete | Single character bond multiplying symbols initiators of rings with implicit bond | Aromatic bond, obsolete | correct ring initiator |
| w_bm_iri | no_bond | Single character bond multiplying symbols initiators of rings with implicit bond | No bond | correct ring initiator |
| w_bm_iri | w_bm_ibi | Single character bond multiplying symbols initiators of rings with implicit bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring initiator |
| w_bm_iri | w_bm_tbi | Single character bond multiplying symbols initiators of rings with implicit bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring initiator |
| w_bm_iri | w_bm_iri | Single character bond multiplying symbols initiators of rings with implicit bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring initiator |
| w_bm_iri | w_bm_tri | Single character bond multiplying symbols initiators of rings with implicit bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring initiator |
| w_bm_iri | s_bm_ibe | Single character bond multiplying symbols initiators of rings with implicit bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_ire_2 | Single character bond multiplying symbols initiators of rings with implicit bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_ire_4 | Single character bond multiplying symbols initiators of rings with implicit bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_iri | Single character bond multiplying symbols initiators of rings with implicit bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_tbe | Single character bond multiplying symbols initiators of rings with implicit bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_tre_2 | Single character bond multiplying symbols initiators of rings with implicit bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_tre_4 | Single character bond multiplying symbols initiators of rings with implicit bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| w_bm_iri | s_bm_tri | Single character bond multiplying symbols initiators of rings with implicit bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring initiator |
| w_bm_iri | lct | Single character bond multiplying symbols initiators of rings with implicit bond | cis/trans on the left side | correct ring initiator |
| w_bm_iri | rct | Single character bond multiplying symbols initiators of rings with implicit bond | cis/trans on the right side | correct ring initiator |
| w_bm_tri | w_atom_oar | Single character bond multiplying symbols terminators of rings with implicit bond | Single character atom symbols of organic aromatic atoms | correct ring terminator |
| w_bm_tri | w_atom_oal | Single character bond multiplying symbols terminators of rings with implicit bond | Single character atom symbols of organic aliphatic atoms | correct ring terminator |
| w_bm_tri | s_atom_oal | Single character bond multiplying symbols terminators of rings with implicit bond | Two character atom symbols of organic aliphatic atoms, start | correct ring terminator |
| w_bm_tri | anything | Single character bond multiplying symbols terminators of rings with implicit bond | Anything | correct ring terminator |
| w_bm_tri | bracket_start | Single character bond multiplying symbols terminators of rings with implicit bond | Square bracket, start | correct ring terminator |
| w_bm_tri | single_bond | Single character bond multiplying symbols terminators of rings with implicit bond | Single bond | correct ring terminator |
| w_bm_tri | double_bond | Single character bond multiplying symbols terminators of rings with implicit bond | Double bond | correct ring terminator |
| w_bm_tri | triple_bond | Single character bond multiplying symbols terminators of rings with implicit bond | Triple bond | correct ring terminator |
| w_bm_tri | quadruple_bond | Single character bond multiplying symbols terminators of rings with implicit bond | Quadruple bond | correct ring terminator |
| w_bm_tri | aromatic_bond_obsolete | Single character bond multiplying symbols terminators of rings with implicit bond | Aromatic bond, obsolete | correct ring terminator |
| w_bm_tri | no_bond | Single character bond multiplying symbols terminators of rings with implicit bond | No bond | correct ring terminator |
| w_bm_tri | w_bm_ibi | Single character bond multiplying symbols terminators of rings with implicit bond | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring terminator |
| w_bm_tri | w_bm_tbi | Single character bond multiplying symbols terminators of rings with implicit bond | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring terminator |
| w_bm_tri | w_bm_iri | Single character bond multiplying symbols terminators of rings with implicit bond | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring terminator |
| w_bm_tri | w_bm_tri | Single character bond multiplying symbols terminators of rings with implicit bond | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring terminator |
| w_bm_tri | s_bm_ibe | Single character bond multiplying symbols terminators of rings with implicit bond | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_ire_2 | Single character bond multiplying symbols terminators of rings with implicit bond | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_ire_4 | Single character bond multiplying symbols terminators of rings with implicit bond | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_iri | Single character bond multiplying symbols terminators of rings with implicit bond | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_tbe | Single character bond multiplying symbols terminators of rings with implicit bond | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_tre_2 | Single character bond multiplying symbols terminators of rings with implicit bond | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_tre_4 | Single character bond multiplying symbols terminators of rings with implicit bond | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| w_bm_tri | s_bm_tri | Single character bond multiplying symbols terminators of rings with implicit bond | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring terminator |
| w_bm_tri | lct | Single character bond multiplying symbols terminators of rings with implicit bond | cis/trans on the left side | correct ring terminator |
| w_bm_tri | rct | Single character bond multiplying symbols terminators of rings with implicit bond | cis/trans on the right side | correct ring terminator |
| s_bm_ibe | e_bm_ibe | Two-character bond multiplying symbols initiators of branching with explicit bond, start | Two-character bond multiplying symbols initiators of branching with explicit bond, end | correct multichar |
| e_bm_ibe | s_atom_oal | Two-character bond multiplying symbols initiators of branching with explicit bond, end | Two character atom symbols of organic aliphatic atoms, start | correct branch initiator |
| s_bm_ire_2 | e_bm_ire_2 | Two-character bond multiplying symbols initiators of rings with explicit bond, start | Two-character bond multiplying symbols initiators of rings with explicit bond, end | correct multichar |
| e_bm_ire_2 | s_atom_oal | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Two character atom symbols of organic aliphatic atoms, start | correct ring initiator |
| e_bm_ire_2 | s_bm_ibe | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_ire_4 | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_iri | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_tbe | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_tre_2 | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_tre_4 | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| e_bm_ire_2 | s_bm_tri | Two-character bond multiplying symbols initiators of rings with explicit bond, end | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring initiator |
| s_bm_ire_4 | n_bm_ire_4 | Four-character bond multiplying symbols initiators of rings with explicit bond, start | Four-character bond multiplying symbols initiators of rings with explicit bond, next from start | correct multichar |
| n_bm_ire_4 | r_bm_ire_4 | Four-character bond multiplying symbols initiators of rings with explicit bond, next from start | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | correct multichar |
| r_bm_ire_4 | w_atom_oar | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character atom symbols of organic aromatic atoms | correct ring initiator |
| r_bm_ire_4 | w_atom_oal | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character atom symbols of organic aliphatic atoms | correct ring initiator |
| r_bm_ire_4 | w_atom_bar | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character atom symbols of bracket aromatic atoms | correct ring initiator |
| r_bm_ire_4 | w_atom_bal | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character atom symbols of bracket aliphatic atoms | correct ring initiator |
| r_bm_ire_4 | s_atom_oal | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Two character atom symbols of organic aliphatic atoms, start | correct ring initiator |
| r_bm_ire_4 | anything | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Anything | correct ring initiator |
| r_bm_ire_4 | bracket_start | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Square bracket, start | correct ring initiator |
| r_bm_ire_4 | single_bond | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single bond | correct ring initiator |
| r_bm_ire_4 | double_bond | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Double bond | correct ring initiator |
| r_bm_ire_4 | triple_bond | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Triple bond | correct ring initiator |
| r_bm_ire_4 | quadruple_bond | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Quadruple bond | correct ring initiator |
| r_bm_ire_4 | aromatic_bond_obsolete | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Aromatic bond, obsolete | correct ring initiator |
| r_bm_ire_4 | no_bond | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | No bond | correct ring initiator |
| r_bm_ire_4 | w_bm_ibi | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring initiator |
| r_bm_ire_4 | w_bm_tbi | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring initiator |
| r_bm_ire_4 | w_bm_iri | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring initiator |
| r_bm_ire_4 | w_bm_tri | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring initiator |
| r_bm_ire_4 | s_bm_ibe | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring initiator |
| r_bm_ire_4 | s_bm_ire_2 | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| r_bm_ire_4 | r_bm_ire_4 | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | correct multichar |
| r_bm_ire_4 | s_bm_iri | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring initiator |
| r_bm_ire_4 | s_bm_tbe | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring initiator |
| r_bm_ire_4 | s_bm_tre_2 | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| r_bm_ire_4 | s_bm_tre_4 | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| r_bm_ire_4 | s_bm_tri | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring initiator |
| r_bm_ire_4 | lct | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | cis/trans on the left side | correct ring initiator |
| r_bm_ire_4 | rct | Four-character bond multiplying symbols initiators of rings with explicit bond, rest | cis/trans on the right side | correct ring initiator |
| s_bm_iri | r_bm_iri | Three-character bond multiplying symbols initiators of rings with implicit bond, start | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | correct multichar |
| r_bm_iri | w_atom_oar | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character atom symbols of organic aromatic atoms | correct ring initiator |
| r_bm_iri | w_atom_oal | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character atom symbols of organic aliphatic atoms | correct ring initiator |
| r_bm_iri | w_atom_bar | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character atom symbols of bracket aromatic atoms | correct ring initiator |
| r_bm_iri | w_atom_bal | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character atom symbols of bracket aliphatic atoms | correct ring initiator |
| r_bm_iri | s_atom_oal | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Two character atom symbols of organic aliphatic atoms, start | correct ring initiator |
| r_bm_iri | anything | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Anything | correct ring initiator |
| r_bm_iri | bracket_start | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Square bracket, start | correct ring initiator |
| r_bm_iri | single_bond | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single bond | correct ring initiator |
| r_bm_iri | double_bond | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Double bond | correct ring initiator |
| r_bm_iri | triple_bond | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Triple bond | correct ring initiator |
| r_bm_iri | quadruple_bond | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Quadruple bond | correct ring initiator |
| r_bm_iri | aromatic_bond_obsolete | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Aromatic bond, obsolete | correct ring initiator |
| r_bm_iri | no_bond | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | No bond | correct ring initiator |
| r_bm_iri | w_bm_ibi | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring initiator |
| r_bm_iri | w_bm_tbi | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring initiator |
| r_bm_iri | w_bm_iri | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring initiator |
| r_bm_iri | w_bm_tri | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring initiator |
| r_bm_iri | s_bm_ibe | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring initiator |
| r_bm_iri | s_bm_ire_2 | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| r_bm_iri | s_bm_ire_4 | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring initiator |
| r_bm_iri | r_bm_iri | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | correct multichar |
| r_bm_iri | s_bm_tbe | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring initiator |
| r_bm_iri | s_bm_tre_2 | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| r_bm_iri | s_bm_tre_4 | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring initiator |
| r_bm_iri | s_bm_tri | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring initiator |
| r_bm_iri | lct | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | cis/trans on the left side | correct ring initiator |
| r_bm_iri | rct | Three-character bond multiplying symbols initiators of rings with implicit bond, rest | cis/trans on the right side | correct ring initiator |
| s_bm_tbe | e_bm_tbe | Two-character bond multiplying symbols terminators of branching with explicit bond, start | Two-character bond multiplying symbols terminators of branching with explicit bond, end | correct multichar |
| e_bm_tbe | s_atom_oal | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two character atom symbols of organic aliphatic atoms, start |  correct branch terminator |
| e_bm_tbe | s_atom_bar | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two character atom symbols of bracket aromatic atoms, start |  correct branch terminator |
| e_bm_tbe | s_atom_bal | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two character atom symbols of bracket aliphatic atoms, start |  correct branch terminator |
| e_bm_tbe | s_bm_ibe | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two-character bond multiplying symbols initiators of branching with explicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_ire_2 | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two-character bond multiplying symbols initiators of rings with explicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_ire_4 | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Four-character bond multiplying symbols initiators of rings with explicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_iri | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Three-character bond multiplying symbols initiators of rings with implicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_tre_2 | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Two-character bond multiplying symbols terminators of rings with explicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_tre_4 | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Four-character bond multiplying symbols terminators of rings with explicit bond, start |  correct branch terminator |
| e_bm_tbe | s_bm_tri | Two-character bond multiplying symbols terminators of branching with explicit bond, end | Three-character bond multiplying symbols terminators of rings with implicit bond, start |  correct branch terminator |
| s_bm_tre_2 | e_bm_tre_2 | Two-character bond multiplying symbols terminators of rings with explicit bond, start | Two-character bond multiplying symbols terminators of rings with explicit bond, end | correct multichar |
| e_bm_tre_2 | s_atom_oal | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Two character atom symbols of organic aliphatic atoms, start | correct ring terminator |
| e_bm_tre_2 | s_bm_ibe | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_ire_2 | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_ire_4 | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_iri | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_tbe | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_tre_4 | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| e_bm_tre_2 | s_bm_tri | Two-character bond multiplying symbols terminators of rings with explicit bond, end | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring terminator |
| s_bm_tre_4 | n_bm_tre_4 | Four-character bond multiplying symbols terminators of rings with explicit bond, start | Four-character bond multiplying symbols terminators of rings with explicit bond, next from start | correct multichar |
| n_bm_tre_4 | r_bm_tre_4 | Four-character bond multiplying symbols terminators of rings with explicit bond, next from start | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | correct multichar |
| r_bm_tre_4 | w_atom_oar | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character atom symbols of organic aromatic atoms | correct ring terminator |
| r_bm_tre_4 | w_atom_oal | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character atom symbols of organic aliphatic atoms | correct ring terminator |
| r_bm_tre_4 | s_atom_oal | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Two character atom symbols of organic aliphatic atoms, start | correct ring terminator |
| r_bm_tre_4 | anything | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Anything | correct ring terminator |
| r_bm_tre_4 | bracket_start | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Square bracket, start | correct ring terminator |
| r_bm_tre_4 | single_bond | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single bond | correct ring terminator |
| r_bm_tre_4 | double_bond | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Double bond | correct ring terminator |
| r_bm_tre_4 | triple_bond | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Triple bond | correct ring terminator |
| r_bm_tre_4 | quadruple_bond | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Quadruple bond | correct ring terminator |
| r_bm_tre_4 | aromatic_bond_obsolete | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Aromatic bond, obsolete | correct ring terminator |
| r_bm_tre_4 | no_bond | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | No bond | correct ring terminator |
| r_bm_tre_4 | w_bm_ibi | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring terminator |
| r_bm_tre_4 | w_bm_tbi | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring terminator |
| r_bm_tre_4 | w_bm_iri | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring terminator |
| r_bm_tre_4 | w_bm_tri | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring terminator |
| r_bm_tre_4 | s_bm_ibe | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring terminator |
| r_bm_tre_4 | s_bm_ire_2 | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| r_bm_tre_4 | s_bm_ire_4 | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| r_bm_tre_4 | s_bm_iri | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring terminator |
| r_bm_tre_4 | s_bm_tbe | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring terminator |
| r_bm_tre_4 | s_bm_tre_2 | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| r_bm_tre_4 | r_bm_tre_4 | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | correct multichar |
| r_bm_tre_4 | s_bm_tri | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | Three-character bond multiplying symbols terminators of rings with implicit bond, start | correct ring terminator |
| r_bm_tre_4 | lct | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | cis/trans on the left side | correct ring terminator |
| r_bm_tre_4 | rct | Four-character bond multiplying symbols terminators of rings with explicit bond, rest | cis/trans on the right side | correct ring terminator |
| s_bm_tri | r_bm_tri | Three-character bond multiplying symbols terminators of rings with implicit bond, start | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | correct multichar |
| r_bm_tri | w_atom_oar | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character atom symbols of organic aromatic atoms | correct ring terminator |
| r_bm_tri | w_atom_oal | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character atom symbols of organic aliphatic atoms | correct ring terminator |
| r_bm_tri | s_atom_oal | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Two character atom symbols of organic aliphatic atoms, start | correct ring terminator |
| r_bm_tri | anything | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Anything | correct ring terminator |
| r_bm_tri | bracket_start | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Square bracket, start | correct ring terminator |
| r_bm_tri | single_bond | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single bond | correct ring terminator |
| r_bm_tri | double_bond | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Double bond | correct ring terminator |
| r_bm_tri | triple_bond | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Triple bond | correct ring terminator |
| r_bm_tri | quadruple_bond | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Quadruple bond | correct ring terminator |
| r_bm_tri | aromatic_bond_obsolete | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Aromatic bond, obsolete | correct ring terminator |
| r_bm_tri | no_bond | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | No bond | correct ring terminator |
| r_bm_tri | w_bm_ibi | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character bond multiplying symbols initiators of branching with implicit bond | correct ring terminator |
| r_bm_tri | w_bm_tbi | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character bond multiplying symbols terminators of branching with implicit bond | correct ring terminator |
| r_bm_tri | w_bm_iri | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character bond multiplying symbols initiators of rings with implicit bond | correct ring terminator |
| r_bm_tri | w_bm_tri | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Single character bond multiplying symbols terminators of rings with implicit bond | correct ring terminator |
| r_bm_tri | s_bm_ibe | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Two-character bond multiplying symbols initiators of branching with explicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_ire_2 | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Two-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_ire_4 | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Four-character bond multiplying symbols initiators of rings with explicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_iri | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Three-character bond multiplying symbols initiators of rings with implicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_tbe | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Two-character bond multiplying symbols terminators of branching with explicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_tre_2 | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Two-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| r_bm_tri | s_bm_tre_4 | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Four-character bond multiplying symbols terminators of rings with explicit bond, start | correct ring terminator |
| r_bm_tri | r_bm_tri | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | correct multichar |
| r_bm_tri | lct | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | cis/trans on the left side | correct ring terminator |
| r_bm_tri | rct | Three-character bond multiplying symbols terminators of rings with implicit bond, rest | cis/trans on the right side | correct ring terminator |
| lct | w_atom_oar | cis/trans on the left side | Single character atom symbols of organic aromatic atoms | correct cis/trans marker |
| lct | w_atom_oal | cis/trans on the left side | Single character atom symbols of organic aliphatic atoms | correct cis/trans marker |
| lct | s_atom_oal | cis/trans on the left side | Two character atom symbols of organic aliphatic atoms, start | correct cis/trans marker |
| lct | anything | cis/trans on the left side | Anything | correct cis/trans marker |
| lct | bracket_start | cis/trans on the left side | Square bracket, start | correct cis/trans marker |
| rct | w_atom_oar | cis/trans on the right side | Single character atom symbols of organic aromatic atoms | correct cis/trans marker |
| rct | w_atom_oal | cis/trans on the right side | Single character atom symbols of organic aliphatic atoms | correct cis/trans marker |
| rct | s_atom_oal | cis/trans on the right side | Two character atom symbols of organic aliphatic atoms, start | correct cis/trans marker |
| rct | anything | cis/trans on the right side | Anything | correct cis/trans marker |
| rct | bracket_start | cis/trans on the right side | Square bracket, start | correct cis/trans marker |
| w_isotope | w_atom_bar | Single character isotope symbols | Single character atom symbols of bracket aromatic atoms | correct isotope |
| w_isotope | w_atom_bal | Single character isotope symbols | Single character atom symbols of bracket aliphatic atoms | correct isotope |
| w_isotope | s_atom_bar | Single character isotope symbols | Two character atom symbols of bracket aromatic atoms, start | correct isotope |
| w_isotope | s_atom_bal | Single character isotope symbols | Two character atom symbols of bracket aliphatic atoms, start | correct isotope |
| s_isotope_m | r_isotope_m | Multicharacter isotope symbols, start | Multicharacter isotope symbols, end | correct multichar |
| r_isotope_m | w_atom_bar | Multicharacter isotope symbols, end | Single character atom symbols of bracket aromatic atoms | correct isotope |
| r_isotope_m | w_atom_bal | Multicharacter isotope symbols, end | Single character atom symbols of bracket aliphatic atoms | correct isotope |
| r_isotope_m | s_atom_bar | Multicharacter isotope symbols, end | Two character atom symbols of bracket aromatic atoms, start | correct isotope |
| r_isotope_m | s_atom_bal | Multicharacter isotope symbols, end | Two character atom symbols of bracket aliphatic atoms, start | correct isotope |
| r_isotope_m | r_isotope_m | Multicharacter isotope symbols, end | Multicharacter isotope symbols, end | correct multichar |
| w_chiral | w_hydrogen | Single character chirality symbols | Single character hydrogen symbols | correct chirality |
| w_chiral | s_hydrogen | Single character chirality symbols | Two-character hydrogen symbols, start | correct chirality |
| w_chiral | w_charge | Single character chirality symbols | Single character charge symbols | correct chirality |
| w_chiral | s_charge_obsolete | Single character chirality symbols | Two-character charge obsolete symbols, start | correct chirality |
| w_chiral | s_charge_m | Single character chirality symbols | Multicharacter charge symbols, start | correct chirality |
| w_chiral | s_class | Single character chirality symbols | Multicharacter class symbols, start | correct chirality |
| s_chiral_2 | e_chiral_2 | Two-character chirality symbols, start | Two-character chirality symbols, end | correct multichar |
| e_chiral_2 | s_hydrogen | Two-character chirality symbols, end | Two-character hydrogen symbols, start | correct chirality |
| e_chiral_2 | s_charge_obsolete | Two-character chirality symbols, end | Two-character charge obsolete symbols, start | correct chirality |
| e_chiral_2 | s_charge_m | Two-character chirality symbols, end | Multicharacter charge symbols, start | correct chirality |
| e_chiral_2 | s_class | Two-character chirality symbols, end | Multicharacter class symbols, start | correct chirality |
| s_chiral_m | n_chiral_m | Multicharacter chirality symbols, start | Multicharacter chirality symbols, next from start | correct multichar |
| n_chiral_m | r_chiral_m | Multicharacter chirality symbols, next from start | Multicharacter chirality symbols, rest | correct multichar |
| r_chiral_m | r_chiral_m | Multicharacter chirality symbols, rest | Multicharacter chirality symbols, rest | correct multichar |
| r_chiral_m | w_hydrogen | Multicharacter chirality symbols, rest | Single character hydrogen symbols | correct chirality |
| r_chiral_m | s_hydrogen | Multicharacter chirality symbols, rest | Two-character hydrogen symbols, start | correct chirality |
| r_chiral_m | w_charge | Multicharacter chirality symbols, rest | Single character charge symbols | correct chirality |
| r_chiral_m | s_charge_obsolete | Multicharacter chirality symbols, rest | Two-character charge obsolete symbols, start | correct chirality |
| r_chiral_m | s_charge_m | Multicharacter chirality symbols, rest | Multicharacter charge symbols, start | correct chirality |
| r_chiral_m | s_class | Multicharacter chirality symbols, rest | Multicharacter class symbols, start | correct chirality |
| w_hydrogen | w_charge | Single character hydrogen symbols | Single character charge symbols | correct hydrogen |
| w_hydrogen | s_charge_obsolete | Single character hydrogen symbols | Two-character charge obsolete symbols, start | correct hydrogen |
| w_hydrogen | s_charge_m | Single character hydrogen symbols | Multicharacter charge symbols, start | correct hydrogen |
| w_hydrogen | s_class | Single character hydrogen symbols | Multicharacter class symbols, start | correct hydrogen |
| s_hydrogen | e_hydrogen | Two-character hydrogen symbols, start | Two-character hydrogen symbols, end | correct multichar |
| e_hydrogen | s_charge_obsolete | Two-character hydrogen symbols, end | Two-character charge obsolete symbols, start | correct hydrogen |
| e_hydrogen | s_charge_m | Two-character hydrogen symbols, end | Multicharacter charge symbols, start | correct hydrogen |
| e_hydrogen | s_class | Two-character hydrogen symbols, end | Multicharacter class symbols, start | correct hydrogen |
| w_charge | s_class | Single character charge symbols | Multicharacter class symbols, start | correct charge |
| s_charge_obsolete | e_charge_obsolete | Two-character charge obsolete symbols, start | Two-character charge obsolete symbols, end | correct multichar |
| e_charge_obsolete | s_class | Two-character charge obsolete symbols, end | Multicharacter class symbols, start | correct charge |
| s_charge_m | r_charge_m | Multicharacter charge symbols, start | Multicharacter charge symbols, rest | correct multichar |
| r_charge_m | r_charge_m | Multicharacter charge symbols, rest | Multicharacter charge symbols, rest | correct multichar |
| r_charge_m | s_class | Multicharacter charge symbols, rest | Multicharacter class symbols, start | correct charge |
| s_class | r_class | Multicharacter class symbols, start | Multicharacter class symbols, rest | correct multichar |
| r_class | bracket_end | Multicharacter class symbols, rest | Square bracket, end | correct class |
| r_class | r_class | Multicharacter class symbols, rest | Multicharacter class symbols, rest | correct multichar |

: **Table 1.** Pairs of character classes allowed in SMILES.

On the next step pairs of character classes from **Table 1** will be identified in SMILES strings from ChEMBL, using the results the set of rules will probably be updated.

### Possible pairs of characters in SMILES from ChEMBL

In short, to check the possible pairs of characters occurring in SMILES from ChEMBL at this stage it is possible to

1.  get all the pairs of characters from ChEMBL SMILES

2.  determine the belonging to the particular class for each character in each pair

3.  make the table containing pairs of character classes **occurring** in ChEMBL SMILES

4.  join this table with the pairs of character classes **allowed** in SMILES (Table 1)

5.  make corrections in the Table 1 if needed

6.  think how to discriminate between the different classes containing the same characters

The code for the steps 1-3 is provided bellow.

``` {#characterPairs_in_CheMBL_SMILES .R}
### Extract SMILES
cs__query <- dbSendQuery(con, 'SELECT canonical_smiles FROM compound_structures')
cs_smiles <- dbFetch(cs__query) |> distinct()
dbClearResult(cs__query)
# Close the connection
dbDisconnect(con)

### Pairs of characters from ChEMBL SMILES, for example SEE: https://stackoverflow.com/questions/71147234/extract-all-two-character-combinations-from-a-string
## Get the data
cs_smiles_chars <- cs_smiles |> rowwise() |>
								mutate(pair_list = substring(canonical_smiles, 1:(nchar(canonical_smiles) - 1), 2:nchar(canonical_smiles)) |> list()) |>
								ungroup()
# Create vector of unique pairs
occuring_pairs <- cs_smiles_chars |> pull(pair_list) |> unlist() |> unique() # only 871 really occuring pairs
## Assign the distinct characters to classes
# Tibble to store the results, its size is sufficient to store all the possible outcomes for each pair
occuring_pairs_classes <- tibble( left_char = rep(NA, length(occuring_pairs)*length(available_classes)^2), right_char = rep(NA, length(occuring_pairs)*length(available_classes)^2),
									left_class = rep(NA, length(occuring_pairs)*length(available_classes)^2), right_class = rep(NA, length(occuring_pairs)*length(available_classes)^2),
									left_description = rep(NA, length(occuring_pairs)*length(available_classes)^2), right_description = rep(NA, length(occuring_pairs)*length(available_classes)^2) )
row <- 0
for (main_i in seq(1:length(occuring_pairs))) {
	for (i in seq(1:length(available_classes))) {
		for (k in seq(1:length(available_classes))) {
			row <- row + 1
			# Check if characters of the current pair are present in the current classes
			if (occuring_pairs[main_i] |> substring(1,1) %in% (available_classes[i] |> as.name() |> eval()) & occuring_pairs[main_i] |> substring(2,2) %in% (available_classes[k] |> as.name() |> eval()) ) {
				occuring_pairs_classes[row, 1] <- occuring_pairs[main_i] |> substring(1,1)
				occuring_pairs_classes[row, 2] <- occuring_pairs[main_i] |> substring(2,2)
				occuring_pairs_classes[row, 3] <- available_classes[i]
				occuring_pairs_classes[row, 4] <- available_classes[k]
				occuring_pairs_classes[row, 5] <- available_classes_description[i]
				occuring_pairs_classes[row, 6] <- available_classes_description[k]
			}
		}
	}
}
# Leave only the cases
occuring_characters <- occuring_pairs_classes |> filter(!is.na(left_char))
occuring_classes <- occuring_pairs_classes |> filter(!is.na(left_char)) |> select(-left_char, -right_char) |> distinct()
occuring_pairs <- occuring_pairs_classes |> filter(!is.na(left_char)) |> select(left_char, right_char) |> distinct()
```

At this stage, **828** distinct pairs of characters were identified in ChEMBL SMILES, these pairs potentially could belong to the one or more of **2182** classes, total amount of records (left_character, right_character, left_class, right_class) is **46707** without applying any rules, except the initial one: **if character occurs in the class, it represents this class**.

From this, it is possible to conclude that further adjustment and latter prioritization of parsing procedure is needed.

