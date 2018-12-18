using Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Utilities;

namespace Matcher {
    class RNAMatcher: IMatcher {
        private readonly byte[] u;  // sequence1 of aminoacids
        private readonly byte[] w;  // sequence2 of aminoacids
        private readonly int[,] d;  // distanceMatrix
        private readonly int[,] s;  // similiarityMatrix
        private readonly string RNASequence1;
        private readonly string RNASequence2;
        private readonly Sequence[] uSequences;
        private readonly Sequence[] wSequences;

        private int[,] D;           // edit distances matrix
        private int[,] S;           // similiarities matrix

        private static readonly Dictionary<string, AminoAcids> aminoAcidCodons = new Dictionary<string, AminoAcids>() {
                {"UUU", AminoAcids.Phe},
                {"UUC", AminoAcids.Phe},
                {"UUA", AminoAcids.Leu},
                {"UUG", AminoAcids.Leu},
                {"UCU", AminoAcids.Ser},
                {"UCC", AminoAcids.Ser},
                {"UCA", AminoAcids.Ser},
                {"UCG", AminoAcids.Ser},
                {"UAU", AminoAcids.Tyr},
                {"UAC", AminoAcids.Tyr},
                {"UAA", AminoAcids.STOP},
                {"UAG", AminoAcids.STOP},
                {"UGU", AminoAcids.Cys},
                {"UGC", AminoAcids.Cys},
                {"UGA", AminoAcids.STOP},
                {"UGG", AminoAcids.Trp},
                {"CUU", AminoAcids.Leu},
                {"CUC", AminoAcids.Leu},
                {"CUA", AminoAcids.Leu},
                {"CUG", AminoAcids.Leu},
                {"CCU", AminoAcids.Pro},
                {"CCC", AminoAcids.Pro},
                {"CCA", AminoAcids.Pro},
                {"CCG", AminoAcids.Pro},
                {"CAU", AminoAcids.His},
                {"CAC", AminoAcids.His},
                {"CAA", AminoAcids.Gln},
                {"CAG", AminoAcids.Gln},
                {"CGU", AminoAcids.Arg},
                {"CGC", AminoAcids.Arg},
                {"CGA", AminoAcids.Arg},
                {"CGG", AminoAcids.Arg},
                {"AUU", AminoAcids.Ile},
                {"AUC", AminoAcids.Ile},
                {"AUA", AminoAcids.Ile},
                {"AUG", AminoAcids.Met},
                {"ACU", AminoAcids.Thr},
                {"ACC", AminoAcids.Thr},
                {"ACA", AminoAcids.Thr},
                {"ACG", AminoAcids.Thr},
                {"AAU", AminoAcids.Asn},
                {"AAC", AminoAcids.Asn},
                {"AAA", AminoAcids.Lys},
                {"AAG", AminoAcids.Lys},
                {"AGU", AminoAcids.Ser},
                {"AGC", AminoAcids.Ser},
                {"AGA", AminoAcids.Arg},
                {"AGG", AminoAcids.Arg},
                {"GUU", AminoAcids.Val},
                {"GUC", AminoAcids.Val},
                {"GUA", AminoAcids.Val},
                {"GUG", AminoAcids.Val},
                {"GCU", AminoAcids.Ala},
                {"GCC", AminoAcids.Ala},
                {"GCA", AminoAcids.Ala},
                {"GCG", AminoAcids.Ala},
                {"GAU", AminoAcids.Asp},
                {"GAC", AminoAcids.Asp},
                {"GAA", AminoAcids.Glu},
                {"GAG", AminoAcids.Glu},
                {"GGU", AminoAcids.Gly},
                {"GGC", AminoAcids.Gly},
                {"GGA", AminoAcids.Gly},
                {"GGG", AminoAcids.Gly}
        };

        private enum MatchingType { Minimal, Maximal }

        private enum AminoAcids {
            START,
            Met = START,
            Phe,
            Leu,
            Ser,
            Tyr,
            Cys,
            Trp,
            Pro,
            His,
            Gln,
            Arg,
            Ile,
            Thr,
            Asn,
            Lys,
            Val,
            Ala,
            Asp,
            Glu,
            Gly,
            STOP,
            None = STOP
        }

        private delegate bool Compare(int a, int b);

        public bool Verbose { get; set; } = false;

        public RNAMatcher(string sequence1, string sequence2, int[,] distanceMatrix, int[,] similiarityMatrix) {
            u = sequence1.SplitInParts(3).Select(s => CodonToByte(s)).ToArray();
            w = sequence2.SplitInParts(3).Select(s => CodonToByte(s)).ToArray();
            d = distanceMatrix;
            s = similiarityMatrix;
            RNASequence1 = sequence1;
            RNASequence2 = sequence2;
            uSequences = u.FindAminoAcidSequences((byte)AminoAcids.START, (byte)AminoAcids.STOP).ToArray();
            wSequences = w.FindAminoAcidSequences((byte)AminoAcids.START, (byte)AminoAcids.STOP).ToArray();
        }

        public int ComputeEditDistance(out string[] matching) {
            bool found = false;
            int bestResult = int.MaxValue;
            Sequence uBestSubsequence = null;
            Sequence wBestSubsequence = null;
            int[,] bestD = null;
            matching = null;
            foreach (Sequence uSequence in uSequences) {
                foreach (Sequence wSequence in wSequences) {
                    int result = ComputeEditDistance(uSequence.data, wSequence.data);
                    if (result < bestResult) {
                        found = true;
                        bestResult = result;
                        uBestSubsequence = uSequence;
                        wBestSubsequence = wSequence;
                        bestD = D.Clone() as int[,];
                    }
                }
            }

            if (found)
                matching = GetMatching(bestD, d, MatchingType.Minimal, uBestSubsequence, wBestSubsequence);

            return bestResult;
        }

        public int ComputeSimiliarity(out string[] matching) {
            bool found = false;
            int bestResult = int.MinValue;
            Sequence uBestSubsequence = null;
            Sequence wBestSubsequence = null;
            int[,] bestS = null;
            matching = null;

            foreach (Sequence uSequence in uSequences) {
                foreach (Sequence wSequence in wSequences) {
                    int result = ComputeSimiliarity(uSequence.data, wSequence.data);
                    if (result > bestResult) {
                        found = true;
                        bestResult = result;
                        uBestSubsequence = uSequence;
                        wBestSubsequence = wSequence;
                        bestS = S.Clone() as int[,];
                    }
                }
            }

            if (found)
                matching = GetMatching(bestS, s, MatchingType.Maximal, uBestSubsequence, wBestSubsequence);

            return bestResult;
        }

        public int ComputeLocalSimiliarity(out string[] matching) {
            bool found = false;
            int bestResult = int.MinValue;
            Field bestField = null;
            Sequence uBestSubsequence = null;
            Sequence wBestSubsequence = null;
            int[,] bestS = null;
            matching = null;

            foreach (Sequence uSequence in uSequences) {
                foreach (Sequence wSequence in wSequences) {
                    int result = ComputeLocalSimiliarity(uSequence.data, wSequence.data, out Field field);
                    if (result > bestResult) {
                        found = true;
                        bestResult = result;
                        uBestSubsequence = uSequence;
                        wBestSubsequence = wSequence;
                        bestS = S.Clone() as int[,];
                        bestField = field.Clone();
                    }
                }
            }

            if (found)
                matching = GetMatching(bestS, s, MatchingType.Maximal, uBestSubsequence, wBestSubsequence, true, bestField);

            return bestResult;
        }

        private int ComputeEditDistance(byte[] u, byte[] w) {
            int n = u.Length;
            int m = w.Length;
            D = new int[n + 1, m + 1];
            D[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    D[0, j] += d[CodonToByte("___"), w[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    D[i, 0] += d[u[k], CodonToByte("___")];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    D[i, j] = MathUtils.Min(
                        D[i - 1, j - 1] + d[u[i - 1], w[j - 1]],
                        D[i, j - 1] + d[CodonToByte("___"), w[j - 1]],
                        D[i - 1, j] + d[u[i - 1], CodonToByte("___")]);

            if (Verbose)
                PrintMatrix(D, u, w);

            return D[n, m];
        }

        private int ComputeSimiliarity(byte[] u, byte[] w) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    S[0, j] += s[CodonToByte("___"), w[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    S[i, 0] += s[u[k], CodonToByte("___")];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]],
                            S[i, j - 1] + s[CodonToByte("___"), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], CodonToByte("___")]);

            if (Verbose)
                PrintMatrix(S, u, w);

            return S[n, m];
        }

        private int ComputeLocalSimiliarity(byte[] u, byte[] w, out Field field) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(0,
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]],
                            S[i, j - 1] + s[CodonToByte("___"), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], CodonToByte("___")]);

            if (Verbose)
                PrintMatrix(S, u, w);

            return MathUtils.Max(S, out field);
        }

        private string[] GetMatching(int[,] matrix, int[,] costMatrix, MatchingType matchingType, Sequence uSequence, Sequence wSequence, bool stopAtNonPositive = false, Field startField = null) {
            Field field = startField ?? new Field(uSequence.data.Length, wSequence.data.Length);
            Field lastField = field.Clone();

            int stopSeq1 = (field.x + uSequence.index) * 3; // 3 chars for amino
            int stopSeq2 = (field.y + wSequence.index) * 3;


            StringBuilder outSequence1 = new StringBuilder();
            StringBuilder outSequence2 = new StringBuilder();

            while ((field.x != 0 || field.y != 0) && (!stopAtNonPositive || matrix[field.x, field.y] > 0)) {
                field = GetNextField(matrix, costMatrix, matchingType, field, uSequence.data, wSequence.data);

                if (field.x == lastField.x) {
                    stopSeq2 -= 3;
                    outSequence1.Insert(0, "___");
                    outSequence2.Insert(0, RNASequence2.Substring(stopSeq2, 3));
                }
                else if (field.y == lastField.y) {
                    stopSeq1 -= 3;
                    outSequence1.Insert(0, RNASequence1.Substring(stopSeq1, 3));
                    outSequence2.Insert(0, "___");
                }
                else {
                    stopSeq1 -= 3;
                    stopSeq2 -= 3;
                    outSequence1.Insert(0, RNASequence1.Substring(stopSeq1, 3));
                    outSequence2.Insert(0, RNASequence2.Substring(stopSeq2, 3));
                }

                lastField.Copy(field);
            }

            return new string[] { outSequence1.ToString(), outSequence2.ToString() };
        }

        private static Field GetNextField(int[,] matrix, int[,] costMatrix, MatchingType matchingType, Field startField, byte[] u, byte[] w) {
            Field nextField = null;
            if (startField.x == 0) {
                nextField = startField;
                nextField.y--;
            }
            else if (startField.y == 0) {
                nextField = startField;
                nextField.x--;
            }
            else {
                int value = matrix[startField.x, startField.y];
                int nextValue;
                int bestNextValue;
                Compare compare;

                switch (matchingType) {
                    case MatchingType.Maximal:
                        bestNextValue = int.MinValue;
                        compare = delegate (int a, int b) { return a > b; };
                        break;
                    case MatchingType.Minimal:
                    default:
                        bestNextValue = int.MaxValue;
                        compare = delegate (int a, int b) { return a < b; };
                        break;
                }

                Field field = new Field(startField.x - 1, startField.y - 1);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[u[field.x], w[field.y]] + nextValue == value && compare(nextValue, bestNextValue)) {   // check if the step is legal
                    bestNextValue = nextValue;
                    nextField = field.Clone();
                }
                field = new Field(startField.x, startField.y - 1);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[CodonToByte("___"), w[field.y]] + nextValue == value && compare(nextValue, bestNextValue)) {
                    bestNextValue = nextValue;
                    nextField = field.Clone();
                }
                field = new Field(startField.x - 1, startField.y);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[u[field.x], CodonToByte("___")] + nextValue == value && compare(nextValue, bestNextValue)) {
                    bestNextValue = nextValue;
                    nextField = field;
                }
            }

            return nextField;
        }

        private static byte CodonToByte(string codon) {
            return aminoAcidCodons.TryGetValue(codon, out AminoAcids aminoAcid) ? (byte)aminoAcid : (byte)AminoAcids.None;
        }

        private static char ByteToAmino(byte b) {
            return (char)(b + 'A');
        }

        private void PrintMatrix(int[,] M, byte[] u, byte[] w) {
            int n = u.Length;
            int m = w.Length;
            int padLenght = M.Cast<int>().Select(v => v.ToString().Length).Max() + 1;
            
            Console.Write("  " + "".PadLeft(padLenght));
            for (int i = 0; i < n; i++)
                Console.Write(ByteToAmino(u[i]).ToString().PadLeft(padLenght));
            Console.WriteLine();
            Console.Write("  ");

            for (int j = 0; j <= m; j++) {
                if (j > 0)
                    Console.Write(ByteToAmino(w[j - 1]) + " ");
                for (int i = 0; i <= n; i++)
                    Console.Write(M[i, j].ToString().PadLeft(padLenght));
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
