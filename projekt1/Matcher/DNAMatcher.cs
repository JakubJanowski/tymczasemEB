using Models;
using System;
using System.Linq;
using System.Text;
using Utilities;

namespace Matcher {
    public class DNAMatcher: IMatcher {
        private readonly int n;     // sequence1 length
        private readonly int m;     // sequence2 length
        private readonly byte[] u;  // sequence1
        private readonly byte[] w;  // sequence2
        private readonly int[,] d;  // distanceMatrix
        private readonly int[,] s;  // similiarityMatrix

        private int[,] D;           // edit distances matrix
        private int[,] S;           // similiarities matrix

        private enum MatchingType { Minimal, Maximal }

        private delegate bool Compare(int a, int b);

        public bool Verbose { get; set; } = false;

        public DNAMatcher(string sequence1, string sequence2, int[,] distanceMatrix, int[,] similiarityMatrix) {
            u = sequence1.Select(c => DNAToByte(c)).ToArray();
            w = sequence2.Select(c => DNAToByte(c)).ToArray();
            n = u.Length;
            m = w.Length;
            d = distanceMatrix;
            s = similiarityMatrix;
        }

        public int ComputeEditDistance(out string[] matching) {
            D = new int[n + 1, m + 1];
            D[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    D[0, j] += d[DNAToByte('_'), w[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    D[i, 0] += d[u[k], DNAToByte('_')];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    D[i, j] = MathUtils.Min(
                        D[i - 1, j - 1] + d[u[i - 1], w[j - 1]],
                        D[i, j - 1] + d[DNAToByte('_'), w[j - 1]],
                        D[i - 1, j] + d[u[i - 1], DNAToByte('_')]);

            if (Verbose)
                PrintMatrix(D);

            matching = GetMatching(D, d, MatchingType.Minimal);

            return D[n, m];
        }

        public int ComputeSimiliarity(out string[] matching) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++)
                for (int k = 0; k < j; k++)
                    S[0, j] += s[DNAToByte('_'), w[k]];

            for (int i = 1; i <= n; i++)
                for (int k = 0; k < i; k++)
                    S[i, 0] += s[u[k], DNAToByte('_')];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]],
                            S[i, j - 1] + s[DNAToByte('_'), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], DNAToByte('_')]);

            if (Verbose)
                PrintMatrix(S);

            matching = GetMatching(S, s, MatchingType.Maximal);

            return S[n, m];
        }



        public int ComputeLocalSimiliarity(out string[] matching) {
            int n = u.Length;
            int m = w.Length;
            S = new int[n + 1, m + 1];

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= m; j++)
                    S[i, j] = MathUtils.Max(0,
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]],
                            S[i, j - 1] + s[DNAToByte('_'), w[j - 1]],
                            S[i - 1, j] + s[u[i - 1], DNAToByte('_')]);

            if (Verbose)
                PrintMatrix(S);
            
            int max = MathUtils.Max(S, out Field field);

            matching = GetMatching(S, s, MatchingType.Maximal, true, field);

            return max;
        }

        private string[] GetMatching(int[,] matrix, int[,] costMatrix, MatchingType matchingType, bool stopAtNonPositive = false, Field startField = null) {
            Field field = startField ?? new Field(n, m);
            Field lastField = field.Clone();

            StringBuilder outSequence1 = new StringBuilder();
            StringBuilder outSequence2 = new StringBuilder();

            while ((field.x != 0 || field.y != 0) && (!stopAtNonPositive || matrix[field.x, field.y] > 0)) {
                field = GetNextField(matrix, costMatrix, matchingType, field);

                if (field.x == lastField.x) {
                    outSequence1.Insert(0, '_');
                    outSequence2.Insert(0, ByteToDNA(w[field.y]));    // z nierówności trójkąta będzie działać
                }
                else if (field.y == lastField.y) {
                    outSequence1.Insert(0, ByteToDNA(u[field.x]));
                    outSequence2.Insert(0, '_');
                }
                else {
                    outSequence1.Insert(0, ByteToDNA(u[field.x]));
                    outSequence2.Insert(0, ByteToDNA(w[field.y]));
                }

                lastField.Copy(field);
            }

            return new string[] { outSequence1.ToString(), outSequence2.ToString() };
        }

        private Field GetNextField(int[,] matrix, int[,] costMatrix, MatchingType matchingType, Field startField) {
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
                if (costMatrix[DNAToByte('_'), w[field.y]] + nextValue == value && compare(nextValue, bestNextValue)) {
                    bestNextValue = nextValue;
                    nextField = field.Clone();
                }
                field = new Field(startField.x - 1, startField.y);
                nextValue = matrix[field.x, field.y];
                if (costMatrix[u[field.x], DNAToByte('_')] + nextValue == value && compare(nextValue, bestNextValue)) {
                    bestNextValue = nextValue;
                    nextField = field;
                }
            }

            return nextField;
        }

        private static byte DNAToByte(char c) {
            switch (c) {
                case 'A':
                    return 0;
                case 'C':
                    return 1;
                case 'G':
                    return 2;
                case 'T':
                    return 3;
            }

            return 4;
        }

        private static char ByteToDNA(byte b) {
            switch (b) {
                case 0:
                    return 'A';
                case 1:
                    return 'C';
                case 2:
                    return 'G';
                case 3:
                    return 'T';
            }
            return '_';
        }

        private void PrintMatrix(int[,] M) {
            int padLenght = M.Cast<int>().Select(v => v.ToString().Length).Max() + 1;

            Console.Write("  " + "".PadLeft(padLenght));
            for (int i = 0; i < n; i++)
                Console.Write(ByteToDNA(u[i]).ToString().PadLeft(padLenght));
            Console.WriteLine();
            Console.Write("  ");

            for (int j = 0; j <= m; j++) {
                if (j > 0)
                    Console.Write(ByteToDNA(w[j - 1]) + " ");
                for (int i = 0; i <= n; i++)
                    Console.Write(M[i, j].ToString().PadLeft(padLenght));
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
