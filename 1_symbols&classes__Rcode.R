library(tidyverse)

##### Enumeration of the symbols and character classes allowed in SMILES

#### Symbols

### 1 ATOMS, SYMBOLS
## 01, superclass: Single character atom symbols of organic (from so called organic subset) aromatic atoms lacking the additional grammatical requirements and features
symb__atom_oar <- c('b', 'c', 'n', 'o', 's', 'p')
## 02, superclass: Single character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features
symb__atom_oal <- c('B', 'C', 'N', 'O', 'S', 'P', 'F', 'I')
## 03, superclass: Single character atom symbols of aromatic atoms enclosed within brackets
symb__atom_bar <- c('b', 'c', 'n', 'o', 's', 'p')
## 04, superclass: Single character atom symbols of bracket aliphatic atoms
symb__atom_bal <- c('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U')
## 05, superclass: Two character atom symbols of organic aliphatic atoms lacking the additional grammatical requirements and features
symb__atom_oal <- c('Cl', 'Br')
## 06, superclass: Two character atom symbols of aromatic atoms, which should be enclosed within the square brackets
symb__atom_bar <- c('se', 'as', 'te')
## 07, superclass: Two character atom symbols of aliphatic atoms, which should be enclosed within the square brackets
symb__atom_bal <- c('He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn',
                          'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                          'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
                          'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am',
                          'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')
### ATOMS, CLASSES
e_atom_bal <- c("e", "i", "a", "g", "l", "r", "c", "n", "o", "u", "s", "b", "h", "d", "f", "t", "v", "m", "y", "p", "k")
e_atom_bar <- c("e", "s")
e_atom_oal <- c("l", "r")
s_atom_bal <- c("H", "L", "B", "N", "M", "A", "S", "C", "T", "F", "Z", "G", "K", "R", "P", "I", "X", "O", "D", "E", "Y")
s_atom_bar <- c("s", "a", "t")
s_atom_oal <- c("C", "B")
w_atom_bal <- c("H", "B", "C", "N", "O", "F", "P", "S", "K", "V", "Y", "I", "W", "U")
w_atom_bar <- c("b", "c", "n", "o", "s", "p")
w_atom_oal <- c("B", "C", "N", "O", "S", "P", "F", "I")
w_atom_oar <- c("b", "c", "n", "o", "s", "p")

### 2 ANYTHING, SYMBOL
## 08, superclass: Single character symbol of any atom or basically anything
symb__anything <- c("*")
### ANYTHING, CLASS
anything <- c("*")

### 3 BRACKETS, SYMBOLS
## 09, superclass: Single character square bracket symbols
symb__bracket <- c("[", "]")
### BRACKETS, CLASSES
bracket_start <- c("[")
bracket_end <- c("]")

### 4 BONDS, SYMBOLS
## 10, superclass: Single character bond symbol corresponding to the single bond
symb__single_bond <- c("-")
## 11, superclass: Single character bond symbol corresponding to the double bond
symb__double_bond <- c("=")
## 12, superclass: Single character bond symbol corresponding to the triple bond
symb__triple_bond <- c("#")
## 13, superclass: Single character symbol corresponding to the quadruple bond
symb__quadruple_bond <- c("$")
## 14, superclass: Single character symbol corresponding to the aromatic bond
symb__aromatic_bond_obsolete <- c(":")
## 15, superclass: Single character symbol corresponding to the absence of the bond between the two specific atoms
symb__no_bond <- c(".")
### BONDS, CLASSES
single_bond <- c("-")
double_bond <- c("=")
triple_bond <- c("#")
quadruple_bond <- c("$")
aromatic_bond_obsolete <- c(":")
no_bond <- c(".")


### 5 BOND MODIFIERS (MULTIPLIERS), SYMBOLS
## 16, superclass: Single character bond multiplying symbols initiators of branching with implicit bond
symb__wbm_ibi <- c('(')
## 17, superclass: Single character bond multiplying symbols initiators of rings with implicit bond
symb__wbm_iri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
## 18, superclass: Single character bond multiplying symbols terminators of branching with implicit bond
symb__wbm_tbi <- c(')')
## 19, superclass: Single character bond multiplying symbols terminators of rings with implicit bond
symb__wbm_tri <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
## 20, superclass: Two-character bond multiplying symbols initiators of branching with explicit bond
symb__bm_ibe <- c('(-', '(=', '(#', '($', '(:', '(.')
## 21, superclass: Two-character bond multiplying symbols initiators of rings with explicit bond
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
## 22, superclass: Three-character bond multiplying symbols initiators of rings with implicit bond
symb__bm_iri <- c( '%01', '%02', '%03', '%04', '%05', '%06', '%07', '%08', '%09', '%10', '%11', '%12', '%13', '%14', '%15', '%16', '%17', '%18', '%19', '%20', 
'%21', '%22', '%23', '%24', '%25', '%26', '%27', '%28', '%29', '%30',
'%31', '%32', '%33', '%34', '%35', '%36', '%37', '%38', '%39', '%40', 
'%41', '%42', '%43', '%44', '%45', '%46', '%47', '%48', '%49', '%50',
'%51', '%52', '%53', '%54', '%55', '%56', '%57', '%58', '%59', '%60',
'%61', '%62', '%63', '%64', '%65', '%66', '%67', '%68', '%69', '%70',
'%71', '%72', '%73', '%74', '%75', '%76', '%77', '%78', '%79', '%80',
'%81', '%82', '%83', '%84', '%85', '%86', '%87', '%88', '%89', '%90',
'%91', '%92', '%93', '%94', '%95', '%96', '%97', '%98', '%99' )
## 23, superclass: Four-character bond multiplying symbols initiators of rings with explicit bond
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
## 24, superclass: Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tbe <- c(')-', ')=', ')#', ')$', '):', ').')
## 25, superclass: Two-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_2  <- symb__bm_ire_2
## 26, superclass: Four-character bond multiplying symbols terminators of rings with explicit bond
symb__bm_tre_4  <- symb__bm_ire_4
## 27, superclass: Three-character bond multiplying symbols terminators of rings with implicit bond
symb__bm_tri <- symb__bm_iri
### BOND MODIFIERS (MULTIPLIERS), CLASSES
e_bm_ibe <- c("-", "=", "#", "$", ":", ".")
e_bm_ire <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
e_bm_tbe <- c("-", "=", "#", "$", ":", ".")
e_bm_tre <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
n_bm_ire <- c("%")
n_bm_tre <- c("%")
r_bm_ire <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
r_bm_iri <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
r_bm_tre <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
r_bm_tri <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
s_bm_ibe <- c("(")
s_bm_ire <- c("-", "=", "#", "$", ":", ".")
s_bm_iri <- c("%")
s_bm_tbe <- c(")")
s_bm_tre <- c("-", "=", "#", "$", ":", ".")
s_bm_tri <- c("%")
w_bm_ibi <- c("(")
w_bm_iri <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
w_bm_tbi <- c(")")
w_bm_tri <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

### 6 CIS/TRANS, SYMBOLS
## 28, superclass: Cis/trans single character symbols on the left side of the rotary non-permissive bond
symb__l_ct <- c("/", "\\")
## 29, Cis/trans single character symbols on the right side of the rotary non-permissive bond
symb__r_ct <- c("/", "\\")
### CIS/TRANS, CLASSES
l_ct <- c("\\", "/")
r_ct <- c("\\", "/")

### 7 All the things inside the square brackets besides the main atom symbol, SYMBOLS
## 30, superclass: Single character isotope symbols, SYMBOLS
symb__isotope_1 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
## 31, superclass: Multicharacter (from 2 to 3 characters) isotope symbols, SYMBOLS
symb__isotope_2 <- map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__isotope_3 <- map(symb__isotope_2, \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
## 32, superclass: Single character chirality symbol, SYMBOLS
symb__chiral_1 <- c('@')
## 33, superclass: two-character chirality symbol, SYMBOLS
symb__chiral_2 <- c('@@')
## 34, superclass: Multi-character (four or five character) chirality symbols, SYMBOLS
symb__chir_m <- c("@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3", "@TB1", "@TB2", "@TB3", "@TB4",
                    "@TB5", "@TB6", "@TB7", "@TB8", "@TB9", "@TB10", "@TB11", "@TB12", "@TB13", "@TB14", "@TB15", "@TB16", "@TB17", "@TB18", "@TB19", "@TB20",
                    "@OH1", "@OH2", "@OH3", "@OH4", "@OH5", "@OH6", "@OH7", "@OH8", "@OH9", "@OH10", "@OH11", "@OH12", "@OH13", "@OH14", "@OH15", "@OH16", "@OH17", "@OH18", "@OH19", "@OH20",
                    "@OH21", "@OH22", "@OH23", "@OH24", "@OH25", "@OH26", "@OH27", "@OH28", "@OH29", "@OH30")
## 35, superclass: Single character hydrogen symbols, SYMBOLS
symb__hydrogen_1 <- c("H")
## 36, superclass: hydrogen symbols, SYMBOLS
symb__hydrogen_2 <- c("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9")
## 37, superclass: Single character charge symbols, SYMBOLS
symb__charge_1 <- c("+", "-")
## 38, superclass: Two-character charge obsolete symbols, SYMBOLS
symb__charge_2 <- c("++", "--")
# 39, superclass: Multicharacter (two or three characters) charge symbols, SYMBOLS
symb__charge_m <- c("+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15",
                     "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15")
# 40, superclass: Multicharacter (from 2 to 4 characters) class symbols, SYMBOLS
symb__a_class_2 <- map(c(":"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__a_class_3 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__a_class_4 <- map(c(":0", ":1", ":2", ":3", ":4", ":5", ":6", ":7", ":8", ":9"), \(x) map(symb__isotope_2, \(y) paste(x,y, collapse = "", sep = ""))) |> unlist()
symb__a_class <- c(symb__a_class_2, symb__a_class_3, symb__a_class_4) |> unique()
### All the things inside the square brackets besides the main atom symbol, CLASSES
w_isotope <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
s_isotope <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
r_isotope <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
w_chiral <- c("@")
s_chiral <- c("@")
e_chiral <- c("@")
n_chiral <- c("T", "A", "S", "O")
r_chiral <- c("H", "1", "2", "L", "P", "3", "B", "4", "5", "6", "7", "8", "9", "0")
w_hydro <- c("H")
s_hydro <- c("H")
e_hydro <- c("H")
w_charge <- c("+", "-")
s_charge <- c("+", "-")
e_charge <- c("+", "-")
r_charge <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
s_aclass <- c(":")
r_aclass <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")