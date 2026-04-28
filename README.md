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

-   Also, according to the Wikipedia, "In terms of a graph-based computational procedure, SMILES is a string obtained by printing the symbol nodes encountered in a depth-first tree traversal of a chemical graph. The chemical graph is first trimmed to remove hydrogen atoms and cycles are broken to turn it into a spanning tree."

-   Thus, via direct left to right parsing of the SMILES, it is possible to obtain readily useful representation of the chemical structure.

## General SMILES parsing strategy

The main idea is as follows:

1.  Computer program accepts the SMILES string, i.e., sequence of characters constituting symbols having chemical meaning.

2.  Computer program process this string from left to right one character at time.

    **What is meant by "computer program process"?**

    | - Computer program has default state.
    | - Every time computer program encounters new (next) character, state of the computer program changes accordingly (taking into account program's current state and what character it encounters).
    | - At each step computer program takes some action to produce an output.

3.  Computer program produces an output.

    **What is "output"?**

    | Data structure appropriate for the further computer processing and filled with the chemical data encoded by the input SMILES string.

    To get an insight into what kinds of state switching will be needed and possible for the program, it will be useful to check, which pairs of characters are possible in SMILES.

    | Topic of parsers is well developed, please, see the [Grune, D., & Jacobs, C. J. (2008). Introduction to parsing. In *Parsing techniques: A practical guide* (pp. 61-102). New York, NY: Springer New York.] and documentation and theory associated with the widely accepted parser generator software for the reference.

    To construct such a parser understanding of SMILES is needed, this understanding is about knowing symbols and characters, which could be present in the SMILES string; and rules of their arrangement. Thus, at the first stage all symbols and characters allowed in SMILES will be enumrated and classified.

## Enumeration and classification of the symbols and characters allowed in SMILES

| *This part may require some further adjustments and corrections, but should be OK in general.*

The following types (upper level of classification) of SMILES symbols (sequence of characters having specific meaning) could be enumerated:

1.  Atom symbols
2.  Symbol of any atom or basically anything
3.  Square bracket symbols
4.  Bond symbols
5.  Bond modifying (multiplying) symbols
6.  Cis/Trans symbols
7.  All the symbols inside the square brackets besides the main atom symbol (isotope symbols, chirality symbols, hydrogen symbols, charge symbols, atom class symbols)

Using faceted classification scheme (<https://en.wikipedia.org/wiki/Faceted_classification>) the following aspects meaningful for the parsing task could be used to describe symbols in SMILES further:

-   number of characters constituting the symbol

-   specific grammatical requirements (symbols of some atoms could only be valid when they are enclosed within the square brackets, symbols of other atoms do not require such enclosing)

-   aromaticity (for atoms only)

-   whether symbol marks the start or the end of something in case of symbols, which go in pairs (initiator (left), terminator (right))

-   whether symbol means the single additional bond (branch) or initiation / termination of the cycle (ring) - only for the bond multiplying symbols

-   whether symbol includes bonds explicitly or implicitly - only for the bond multiplying symbols

Using the information above, it is possible to

-   construct symbol classes using meaningful combinations of facets mentioned above

-   construct corresponding classes of characters providing some convenience for parsing

-   assess their intersections and frequency in the available data, which will be useful while selecting particular parsing approach

And then, select particular parsing approach and set of rules within it and set of technologies for implementation to hopefully finally come up with the pretty normal SMILES parser.

### Atom symbols type and corresponding character classes

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. atom, in the SMILES string.

Atom symbols allowed in SMILES could be divided into two facets by their length:

-   symbols consisting of the single character

-   symbols consisting of two characters

Atom symbols allowed in SMILES could be divided into two facets by their grammatical requirements:

-   symbols, which could be written as is, corresponding atoms belong to the so called organic subset

-   symbols, which could be written only in the square brackets, so called bracket atoms and atoms from organic subset on condition that they have additional properties (charge, etc.)

Atom symbols allowed in SMILES could be divided into two categories depending on the nature of their bonding:

-   symbols of the aromatic atoms

-   symbols of the aliphatic atoms

Thus, the following classes of atom symbols allowed in SMILES could be enumerated and labeled:

1.  Single character atom symbols of organic (from so called *organic* subset) aromatic atoms lacking the additional grammatical requirements and features (**atom_oar**):

| b, c, n, o, s, p

Corresponding characters could be designated as distinct character class, **w_atom_oar**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **ar** - for aromatic.

2.  Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal**):

| B, C, N, O, S, P, F, I

Corresponding characters could be designated as distinct character class, **w_atom_oal**, where prefix **w** stands for the whole symbol, suffix **o** - for organic and suffix **al** - for aliphatic.

3.  Single character atom symbols of aromatic atoms enclosed within brackets (**atom_bar**):

| b, c, n, o, s, p

Corresponding characters could be designated as **w_atom_bar**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **ar** - for aromatic. As it can be seen, this class contains the same symbols as w_atom_oar, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

4.  Single character atom symbols of bracket aliphatic atoms (**atom_bal**):

| H, B, C, N, O, F, P, S, K, V, Y, I, W, U

Corresponding characters could be designated as **w_atom_bal**, where prefix **w** stands for the whole symbol, suffix **b** - for bracket and suffix **al** - for aliphatic. As it can be seen, this class contains the same symbols as w_atom_oal, they could be distinguished only using surrounding symbols: if atom has additional properties, its symbol should be put into the square brackets and, thus, belongs to the w_atom_bar.

5.  Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal_2**):

| Cl, Br

Corresponding character classes could be designated as **s_atom_oal** & **e_atom_oal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **o** - for organic and suffix **al** - for aliphatic. Further division of the characters describing symbol into two classes could be useful if the resulting parser will operate one character at time.

6.  Two character atom symbols of in-bracket aromatic atoms (**atom_bar_2**):

| se, as, te

Corresponding characters could be designated as **s_atom_bar** & **e_atom_bar**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **ar** - for aromatic.

7.  Two character atom symbols of in-bracket aliphatic atoms (**atom_bal_2**):

| He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

Corresponding characters could be designated as **s_atom_bal** & **e_atom_bal**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol; suffix **b** - for bracket and suffix **al** - for aliphatic.

### Symbol of anything and corresponding character class

**\*** is an allowed symbol in SMILES, it corresponds to any (atom) symbol and behaves similar to the single character atom symbols of organic aliphatic and aromatic atoms :

8.  Single character symbol of any atom or basically **anything**:

| \*

Corresponding character could be designated as **w_any**, since it could mean basically any character or symbol.

### Square bracket symbols and corresponding character classes

#### What are they?

9.  Square bracket symbols **[, ]** is the SMILES way to mark the start of the atom record including its various properties and is the way to designate the end of the atom record.

| [, ]

Corresponding character classes could be designated as **s_bracket** & **e_bracket**.

### Bond symbols and corresponding character classes

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, all of them are single character and do not have other peculiar aspects, five of them correspond to the conventional type of chemical bond:**

10. Single character bond symbol corresponding to the single bond (**single_bond**):

| -

This single character symbol could be and typically is omitted, since by default all the atoms, which symbols are written side by side in SMILES string, are presumed to be connected by this type of bond. Corresponding character class will be designated as **w_single_bond**.

11. Single character bond symbol corresponding to the double bond (**double_bond**):

| =

Corresponding character class will be designated as **w_double_bond**.

12. Single character bond symbol corresponding to the triple bond (**triple_bond**):

| \#

Corresponding character class will be designated as **w_triple_bond**.

13. Single character symbol corresponding to the quadruple bond (**quadruple_bond**):

| \$

Corresponding character class will be designated as **w_quadruple_bond**.

14. Single character symbol corresponding to the aromatic bond (**aromatic_bond**):

| :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications. Corresponding character class will be designated as **w_aromatic_bond**.

15. Single character symbol corresponding to the absence of the bond between the two specific atoms (**no_bond**):

| .

### Bond modifying (multiplying) symbols and corresponding characters

#### What are they?

Bond modifying (multiplying) symbols, i.e. **modifiers**, are used in SMILES to extend the number of atoms, for which connections to the current atom could be written using linear notation (SMILES).

Bond multiplying symbols allowed in SMILES could be divided into four facets by their length:

-   Single character symbols

-   Two-character symbols

-   Three-character symbols

-   Four-character symbols

Bond multiplying symbols allowed in SMILES go in pairs and thus could be divided into two facets according to their role in completing the task of the symbols pair:

-   Symbols initiators

-   Symbols terminators

Bond multiplying symbols allowed in SMILES could be divided into two facets by their task:

-   Symbols used to indicate simple additional bond for the current atom (branch)

-   Symbols used to indicate additional bond, which allows for the cycle (ring) to be formed

Bond multiplying symbols allowed in SMILES could be divided into two categories according to whether additional bond is explicitly written or assumed:

-   Symbols explicitly including additional bond (any bond could be added)

-   Symbols implicitly including additional bond (only single bond could be added)

##### Thus, the following classes of bond multiplying symbols could be found in SMILES:

16. Single character bond multiplying symbols initiators of branching with implicit bond (**bm_ibi**):

| (

Corresponding character class could be designated as **w_bm_ibi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, and second suffix **i** - for implicit.

17. Single character bond multiplying symbols initiators of rings with implicit bond (**bm_iri**):

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_bm_iri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, and second suffix **i** - for implicit.

18. Single character bond multiplying symbols terminators of branching with implicit bond (**bm_tbi**):

| )

Corresponding characters could be designated as **w_bm_tbi**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and second suffix **i** - for implicit.

19. Single character bond multiplying symbols terminators of rings with implicit bond (**bm_tri**):

| 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding characters could be designated as **w_bm_tri**, where prefix **w** stands for the whole symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, suffix **b** - for branching, and suffix **i** - for implicit.

20. Two-character bond multiplying symbols initiators of branching with explicit bond (**bm_ibe**):

| ([-=#\$:.]

Corresponding character classes could be designated as **s_bm_ibe** & **e_bm_ibe** where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **b** - for branching, second suffix **e** - for explicit.

21. Two-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_2**):

| [-=#\$:.][0-9]

Corresponding characters could be designated as **s_bm_ire** & **e_bm_ire**,where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit.

22. Three-character bond multiplying symbols initiators of rings with implicit bond (**bm_iri_3**):

| %[0-9][1-9], %[1-9][0-9]

Corresponding character classes could be designated as **s_bm_iri** & **r_bm_iri**, where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, suffix **r** - for ring, next suffix **i** - for implicit.

23. Four-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_4**):

| [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding character classes could be designated as **s_bm_ire_4, n_bm_ire** & **r_bm_ire**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **i** - for initiator, second suffix **r** - for ring, suffix **e** - for explicit; **4** in s_bm_ire_4 is used to differentiate this class from the s_bm_ire corresponding to the bm_ire_2, n_bm_ire and r_bm_ire are already unique.

24. Two-character bond multiplying symbols terminators of branching with explicit bond (**bm_tbe_2**):

| )[-=#\$:.]

Corresponding character classes could be designated as **s_bm_tbe** & **e_bm_tbe**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **b** - for branch, and last suffix **e** - for explicit.

25. Two-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_2**):

|  [-=#\$:.][0-9]

Corresponding character classes could be designated as **s_bm_tre** & **e_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **e** stands for the end of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit.

26. Four-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_4**):

| [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

Corresponding characters could be designated as **s_bm_tre, n_bm_tre** & **r_bm_tre**, where prefix **s** stands for the start of the symbol, prefix **n** stands for the next from start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, suffix **e** - for explicit. **4** in s_bm_tre_4 is used to differentiate this class from the s_bm_tre corresponding to the bm_tre_2, n_bm_tre and r_bm_tre are already unique.

27. Three-character bond multiplying symbols terminators of rings with implicit bond (**bm_tri_3**):

| %[0-9][1-9], %[1-9][0-9]

Corresponding characters could be designated as **s_bm_tri** & **r_bm_tri**,where prefix **s** stands for the start of the symbol, prefix **r** stands for the rest of the symbol, **bm** - bond modifying (multiplying), suffix **t** - for terminator, second suffix **r** - for ring, next suffix **i** - for implicit.

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two facets of cis/trans symbols are allowed in SMILES:

28. Cis/trans single character symbols on the left side of the rotary non-permissive bond (**l_ct**):

| /, \\

Corresponding single character characters could be designated as **l_ct**, where prefix **l** stands for the left side; **ct** - for cis/trans.

29. Cis/trans symbols on the right side of the rotary non-permissive bond:

| /, \\

Corresponding characters could be designated as **r_ct**, where prefix **r** stands for the right side; **ct** - for cis/trans.

The logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\\F)=C/F are equivalent, since

| The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.
| *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****
| **\*** <http://opensmiles.org/opensmiles.html>

### All the symbols and corresponding character classes inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating specific isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

30. Single character isotope symbols (**isotope**):

| 1, 2, 3, 4, 5, 6, 7, 8, 9

Corresponding character class could be designated as **w_isotope**, where prefix **w** stands for the whole symbol.

31. Multicharacter (from 2 to 3 characters) isotope symbols (**isotope_m**):

| [0-9][1-9], [1-9][0-9], [0-9][0-9][1-9], [0-9][1-9][0-9], [1-9][0-9][0-9]

Corresponding characters could be designated as **s_isotope & r_isotope**, where prefix **s** stands for the start and prefix **r** stands for the rest of the symbol.

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

32. Single character chirality symbol (**chiral**):

| \@

Corresponding character class could be designated as **w_chiral**, where prefix **w** stands for the whole symbol.

33. Two-character chirality symbols (**chiral_2**):

| [\@][\@]

Corresponding character classes could be designated as **s_chiral & e_chiral**, where prefix **s** stands for the start and prefix **e** stands for the end of the symbol.

34. Multicharacter (four or five character) chirality symbols (**chiral_m**):

| [\@]TH[1-2], [\@]AL[1-2], [\@]SP[1-3], [\@]TB[1-20], [\@]OH[1-30]

Corresponding character classes could be designated as **s_chiral_m, m_chiral & r_chiral**, where prefix **s** stands for the start, *prefix* **m** stands for the medium (two characters), prefix **r** stands for the rest of the symbol and *suffix* **m** stands for the multi, where it is needed.

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 facets by their length:

35. Single character hydrogen symbol (**hydro**):

| H

Corresponding character classes could be designated as **w_hydro**, where prefix **w** stands for the whole symbol.

36. Two-character hydrogen symbols (**hydro_2**):

| H[2-9]

Corresponding character classes could be designated as **s_hydro & ehydro**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

##### Charge symbols

Charge symbols are used to describe the charge of this atom (**charge**).

Charge symbols allowed in SMILES could be divided into 2 facets by their length:

37. Single character charge symbols (**charge**):

| [+-]

Corresponding character classes could be designated as **w_charge**, where prefix **w** stands for the whole symbol.

38. Two-character charge obsolete symbols (**charge_2**):

| [+][+], [-][-]

Corresponding characters could be designated as **s_charge & e_charge**, where prefix **s** stands for the start and **e** stands for the end of the symbol.

39. Multicharacter (two or three characters) charge symbols (**charge_m**):

| [+-][1-9], [+-]1[0-5]

Corresponding characters could be designated as **s_charge_m & r_charge_m**, where prefix **s** stands for the start and **r** stands for the rest of the symbol and suffix **m** stands for the multi where it is needed.

##### Class symbols

Class symbols designate user-defined class of the atom.

Class symbols allowed in SMILES may have variable length, but there is no point to divide them into facets based on this aspect, so there is only:

40. Multicharacter (from 2 to 4 characters) class symbols (**class**):

| :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]

Corresponding character classes could be designated as **s_class & r_class**, where prefix **s** stands for the start and **r** stands for the rest of the symbol.

Information on symbols is summarized in **symbols.tsv**, all symbols are provided in **symbols_all.tsv**.

## Q1: are described symbols unique, i.e. is it possible to identify each SMILES symbols based only on characters constituting it?

No, as it can be seen from **symbols.tsv**. For example,

-   terminators and initiators of branching could be and often are the same by design

-   symbols of in-bracket atoms and bracket-free atoms could be the same (organic subset)

Also, several classes of symbols could be described or are described partially by the patterns [0-9] and [1-9], which makes parsing without consideration of the environment questionable.

In the previous version the attempt was taken to divide the whole symbols into the smaller subsets of characters, the result is as follows: this strategy does not pay off, character classes probably could be useful to construct the SMILES strings, but they provide no clear benefits for parsing:

> [!NOTE]
>
> It seems to be easier (and fast enough) to read the characters one by one until the longest possible sequence describing symbol is gathered (5 characters, I guess) and decide on the actual symbol afterwards, considering matches in **current sub-string** and **state** deduced from the previous symbols and length of the remaining SMILES string.
> 

## General SMILES parsing strategy revised

1.  Computer program initializes with the SMILES string having default **state** and empty **accumulator** of characters and empty **result**.

2.  Every time computer program encounters new (next) character it accumulates this character.

3.  Every time accumulator reaches its limits (longest symbol in SMILES or end of the string) its content is being evaluated, state changes accordingly, accumulator's content is being trimmed from left to right to delete all the characters evaluated as the whole symbol at this step.

4.  Every time state changes computer program takes some action to build up an output.

5.  When the end of the string is reached, computer program outputs the result.

> [!NOTE]
>
> The following text will be rewritten accordingly.
>
