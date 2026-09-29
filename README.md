# SMILES Parser

**The whole thing besides the symbols' classification (maybe) is in progress and will be modified substantially.** **update from 29.09.2026: symbols' classification is slightly reconsidered according to the current practical attempt to develop parser, and some other minor tweaks.**

For the SMILES (Simplified Molecular Input Line Entry System) reference, please, SEE:

-   <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>

-   <http://opensmiles.org/opensmiles.html>

-   <https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System>

In short:

-   SMILES is a language, which is used to describe a molecule in terms of and according to the rules of chemical valence and mathematical graph theory.

-   SMILES string is the linear notation of spanning tree of the graph representing molecule.

-   SMILES string includes atom symbols, bonds between them and some characteristics of these entities.

-   Also, according to the Wikipedia, "From the view point of a formal language theory, SMILES is a word. A SMILES is parsable with a context-free parser".

-   To build such a parser it will be useful to understand, what are the symbols and corresponding characters allowed in SMILES and what are their allowed combinations.

-   Also, according to the Wikipedia, "In terms of a graph-based computational procedure, SMILES is a string obtained by printing the symbol nodes encountered in a depth-first tree traversal of a chemical graph. The chemical graph is first trimmed to remove hydrogen atoms and cycles are broken to turn it into a spanning tree."

-   Thus, via direct left to right parsing of the SMILES, it is possible to obtain readily useful representation of the chemical structure.

## General SMILES parsing strategy

1.  Computer program accepts the SMILES string, i.e., sequence of characters constituting symbols having chemical meaning.

2.  Computer program process this string from left to right in chunks equal in length to the longest symbol possible in SMILES.

    **What is meant by "computer program process"?**

    | - Computer program has default state.
    | - Every time computer program encounters new (next) chunk and symbol, state of the computer program changes accordingly (taking into account program's current state and what symbol it encounters).
    | - At each step computer program takes some actions to produce an output.

3.  Computer program produces an output, i.e. graph representing chemical structure in such a way that this structure could be processed in computer furhter.

    **- - -**

    To get an insight into what kinds of state switching will be needed and possible for the program, it will be useful to check, which pairs of characters are possible in SMILES.

    | Topic of parsers is well developed, please, see the [Grune, D., & Jacobs, C. J. (2008). Introduction to parsing. In *Parsing techniques: A practical guide* (pp. 61-102). New York, NY: Springer New York.] and documentation and theory associated with the widely accepted parser generator software for the reference.

    To construct such a parser understanding of SMILES is needed, this understanding is about knowing symbols and characters, which could be present in the SMILES string; and rules of their arrangement. Thus, at the first stage all symbols and characters allowed in SMILES will be enumerated and classified.

    ## Enumeration and classification of the symbols and characters allowed in SMILES

| *This part may require some further adjustments and corrections, but should be OK in general.*

The following types (upper level of classification) of SMILES symbols (sequence of characters having specific meaning) could be enumerated:

1.  Atom symbols including special ones
2.  Square bracket symbols
3.  Bond symbols including special one
4.  Bond modifying (multiplying) symbols
5.  Cis/Trans symbols
6.  All the symbols inside the square brackets besides the main atom symbol, i.e. properties (isotope symbols, chirality symbols, hydrogen symbols, charge symbols, atom class symbols)

Using faceted classification scheme (<https://en.wikipedia.org/wiki/Faceted_classification>) the following aspects meaningful for the parsing task could be used to describe symbols in SMILES further:

-   number of characters constituting the symbol

-   specific grammatical requirements (symbols of some atoms could only be valid when they are enclosed within the square brackets, symbols of other atoms do not require such enclosing)

-   aromaticity (for atoms only now and also for bonds in the earlier versions of SMILES, compatibility)

-   whether symbol marks the start or the end of something in case of symbols, which go in pairs (initiator (left), terminator (right))

-   whether symbol means the single additional bond (branch) or initiation / termination of the cycle (ring) - only for the bond multiplying symbols

-   whether symbol includes bonds explicitly or implicitly - only for the bond multiplying symbols

Using the information above, it is possible to

-   construct symbol classes using meaningful combinations of facets mentioned above

-   construct corresponding classes of characters providing some convenience for parsing -\> deprecated, now it seems that level of distinct characters is too much for the task of SMILES parsing

-   assess intersections and frequency of distinct symbols, types and classes in the available data, which will be useful while selecting particular parsing tactics

And then, select particular parsing approach and set of rules within it and set of technologies for implementation to hopefully finally come up with the pretty normal SMILES parser.

### Atom symbol type and corresponding symbol classes

#### What are they?

Atom symbol is the way to designate the node of the molecular graph, i.e. an atom, in the SMILES string.

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

> b, c, n, o, s, p

2.  Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal**):

> B, C, N, O, S, P, F, I

3.  Single character atom symbols of aromatic atoms enclosed within brackets (**atom_bar**):

> b, c, n, o, s, p

4.  Single character atom symbols of bracket aliphatic atoms (**atom_bal**):

> H, B, C, N, O, F, P, S, K, V, Y, I, W, U

5.  Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features (**atom_oal_2**):

> Cl, Br

6.  Two character atom symbols of in-bracket aromatic atoms (**atom_bar_2**):

> se, as, te

7.  Two character atom symbols of in-bracket aliphatic atoms (**atom_bal_2**):

> He, Li, Be, Ne, Na, Mg, Al, Si, Cl, Ar, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,Xe, Cs, Ba, Hf, Ta, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr

8.  Single character symbol of any atom or basically **anything**:

> \*

9.  Zero character (pseudo)symbol trailing start and end of the SMILES string (**ambient**)

> ""

### Square bracket symbol type and corresponding symbol classes

#### What are they?

Square bracket symbols **[, ]** is the SMILES way to mark the start of the record for atom having some property (-ies) and is the way to designate the end for such a record.

10. Single character symbol to start the record for atom having property (-ies) (**s_bracket**)

> [

11. Single character symbol to end the record for atom having property (-ies) (**e_bracket**)

> ]

### Bond symbols and corresponding character classes

#### What are they?

Bond symbol is the way to designate the edge of the molecular graph, i.e. chemical bond, in the SMILES string.

**There are six bond symbols allowed in SMILES, all of them are single character symbols and do not have other peculiar aspects, five of them correspond to the conventional type of chemical bond:**

12. Single character bond symbol corresponding to the single bond (**single_bond**):

> \-

13. Single character bond symbol corresponding to the double bond (**double_bond**):

> =

14. Single character bond symbol corresponding to the triple bond (**triple_bond**):

> \#

15. Single character symbol corresponding to the quadruple bond (**quadruple_bond**):

> \$

16. Single character symbol corresponding to the aromatic bond (**aromatic_bond**):

> :

It should be noted that this symbol (**:**) is deprecated and typically omitted. Aromaticity is rather described using atom symbols: **C** - aliphatic carbon, **c** - aromatic carbon; thus, bond between the **c** and **c** is considered aromatic without additional indications.

17. Single character symbol corresponding to the absence of the bond between the two specific atoms (**no_bond**):

> .

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

18. Single character bond multiplying symbols initiators of branching with implicit bond (**bm_ibi**):

> (

19. Single character bond multiplying symbols initiators of rings with implicit bond (**bm_iri**):

> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

20. Single character bond multiplying symbols terminators of branching with implicit bond (**bm_tbi**):

> )

21. Single character bond multiplying symbols terminators of rings with implicit bond (**bm_tri**):

> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

22. Two-character bond multiplying symbols initiators of branching with explicit bond (**bm_ibe**):

> ([-=#\$:.]

23. Two-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_2**):

> [-=#\$:.][0-9]

24. Three-character bond multiplying symbols initiators of rings with implicit bond (**bm_iri_3**):

> \%[0-9][1-9], %[1-9][0-9]

25. Four-character bond multiplying symbols initiators of rings with explicit bond (**bm_ire_4**):

> [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

26. Two-character bond multiplying symbols terminators of branching with explicit bond (**bm_tbe_2**):

> )[-=#\$:.]

27. Two-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_2**):

> [-=#\$:.][0-9]

28. Four-character bond multiplying symbols terminators of rings with explicit bond (**bm_tre_4**):

> [-=#\$:.]%[0-9][1-9], [-=#\$:.]%[1-9][0-9]

29. Three-character bond multiplying symbols terminators of rings with implicit bond (**bm_tri_3**):

> \%[0-9][1-9], %[1-9][0-9]

### Cis/Trans symbols and corresponding characters

#### What are they?

Cis/trans symbols is the way to designate the position of the nodes of the molecular graph, i.e. atoms, relative to the rotary non-permissive bond (=, #, \$).

Cis/trans symbols should always be paired, i.e. atoms on each side of the bond should have their own cis/trans symbol or such symbols should be omitted on each side of the bond. Thus, two facets of cis/trans symbols are allowed in SMILES:

30. Cis/trans single character symbols on the left side of the rotary non-permissive bond (**l_ct**):

> /, \\

31. Cis/trans symbols on the right side of the rotary non-permissive bond:

> /, \\

Cis/trans symbols could be seen as bond modifiers, however, there is a rationale for them to constitute the distinct type: other mean to provide the spatial information is via atoms' properties, not bond  modifying symbols.
Thus, distinct type, to avoid confusion between the properties and modifiers. 

Chemical logic behind these symbols is outstandingly well described in <http://opensmiles.org/opensmiles.html> including the fact that such combinations of these symbols as in F/C=C/F and C(\\F)=C/F are equivalent, since

> The "visual interpretation" of the "up-ness" or "down-ness" of each single bond is **relative to the carbon atom**, not the double bond, so the sense of the symbol changes when the fluorine atom moved from the left to the right side of the alkene carbon atom.

> *Note: This point was not well documented in earlier SMILES specifications, and several SMILES interpreters are known to interpret the `'/'` and `'\'` symbols incorrectly.**\****

> **\*** <http://opensmiles.org/opensmiles.html>

### All the symbols and corresponding character classes inside the square brackets besides the main atom symbol

#### What are they?

Symbols and corresponding characters inside the square brackets besides the main atom symbol describe the main bracket atom in terms of its mass number indicating specific isotope, chiral status, number of explicit hydrogens, charge and class assigned by the author of the particular SMILES string. It should be noted that any atom symbol could be found in the square brackets and any atom symbol should be put in the square brackets if corresponding atom has aforementioned properties.

These symbols will be categorized only by the length, this is sufficient for the purpose, since these symbols have the strict order of placement inside the brackets.

##### Isotope symbols

Isotope symbols are the symbols describing mass number of the specific atom.

Isotope symbols allowed in SMILES could be divided into 3 categories by their length:

32. Single character isotope symbols (**isotope**):

> 1, 2, 3, 4, 5, 6, 7, 8, 9

33. Multicharacter (from 2 to 3 characters) isotope symbols (**isotope_m**):

> [0-9][1-9], [1-9][0-9], [0-9][0-9][1-9], [0-9][1-9][0-9], [1-9][0-9][0-9]

##### Chirality symbols

Chirality symbols are used to show that an atom is a stereocenter.

Chirality symbols allowed in SMILES could be divided into 5 categories by their length:

34. Single character chirality symbol (**chiral**):

> \@

35. Two-character chirality symbols (**chiral_2**):

> [\@][\@]

36. Multicharacter (four or five character) chirality symbols (**chiral_m**):

> [\@]TH[1-2], [\@]AL[1-2], [\@]SP[1-3], [\@]TB[1-20], [\@]OH[1-30]

##### Hydrogen symbols

Hydrogen symbols are used to designate the number of explicit hydrogens of this atom.

Hydrogen symbols allowed in SMILES could be divided into 2 facets by their length:

37. Single character hydrogen symbol (**hydro**):

> H

38. Two-character hydrogen symbols (**hydro_2**):

> H[2-9]

##### Charge symbols

Charge symbols are used to describe the charge of this atom (**charge**).

Charge symbols allowed in SMILES could be divided into 2 facets by their length:

39. Single character charge symbols (**charge**):

> [+-]

40. Two-character charge obsolete symbols (**charge_2**):

> [+][+], [-][-]

41. Multicharacter (two or three characters) charge symbols (**charge_m**):

> [+-][1-9], [+-]1[0-5]

##### Class symbols

Class symbols designate user-defined class of the atom.

Class symbols allowed in SMILES may have variable length, but there is no point to divide them into facets based on this aspect, so there is only:

42. Multicharacter (from 2 to 4 characters) atom class symbol type (**class**):

> :[0-9], :[0-9][0-9], :[0-9][0-9][0-9]
