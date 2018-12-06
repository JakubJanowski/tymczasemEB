# tymczasemEB

A project realised for Elements of Bioinformatics subject. 
The program computes an optimal matching of 2 DNA sequences based either on edit distance or similiarity.
The sequences and edit distance or similiarity criteria are given in .txt files in /data directory.

The edit distance and similiarity criteria are 5x5 matrices with row and column values corresponding to sequence letter change or indel ('_') as shown:

   A  C  G  T  _
A  .  .  .  .  .
C  .  .  .  .  .
G  .  .  .  .  .
T  .  .  .  .  .
_  .  .  .  .  .

It is also possible to translate given DNA sequence to RNA sequence or aminoacid sequence.
