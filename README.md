A useful tools to search a list of peptides in a fasta file (a list of proteins, or a proteome)

# Description
    
This tool consists one module:

- `PeptideSearch`: this tool searches for peptides in a fasta file including proteins, it can find Exact Match or Marches with up to One Mismatch 

# Installation
 
## Normal installation

```bash
pip install PeptideSearch
```

## Development installation

```bash
git clone https://github.com/khaledianehdieh/PeptideSearch.git
```

## Usage

```
python3
>>> from PeptideSearch import PeptideSearch as PS 
>>> P=PS.PeptideSearch(Peptides, Fasta_File) #Peptides is a .txt file that has one peptide per line, and Fasta_File is a fasta file containing the proteins

#you might use the function as you need in three different ways:
#1- Find Exact Matches
>>> df_EM, pList= P.ExactMatch()  # returns Exact Matches (df_EM) and a list of peptides (pList) that didn't find any matche

#2- Find One MisMatches
>>> peptides= P.read_Peptides()
df_OM, NotFoundList= P.OneMismatch(peptides) #returns One  Mis Matches (df_OM) for and list of peptides (NotFoundList) that didn't find any matche

#3- Combine one and two, first find a list of exacxt matches and then look for one mismatch for the peptides that we diddnt find any match
#Also saves the result in a CSV file
>>> df_all, NotFounPeptides = P.MatchFinder()   #returns All Matches (df_all) for and a list of peptides (NotFounPeptides) that couldn't find any match

#4- Find a peptide in one sequence
result,index=PS.PeptideSearch.SequenceSearch(peptide, sequence) # result=0 means exact match, result= 1 means one mismatch,  result= -1 means no match, and index is the start location of the peptide in the sequence.

```
