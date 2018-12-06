using System;
using System.Linq;
using System.Text;

namespace DNA {
    public class DNAMatcher {
        private readonly int n;     // sequence1 length
        private readonly int m;     // sequence2 length
        private readonly byte[] u;  // sequence1
        private readonly byte[] w;  // sequence2
        private readonly int[,] d;  // distanceMatrix
        private readonly int[,] s;  // similiarityMatrix

        private int[,] D;           // edit distances matrix
        private int[,] S;           // similiarities matrix

        private enum MatchingType { Minimal, Maximal }

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

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    D[0, j] += d[DNAToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    D[i, 0] += d[u[k], DNAToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    D[i, j] = Math.Min(
                        D[i - 1, j - 1] + d[u[i - 1], w[j - 1]], Math.Min(
                        D[i, j - 1] + d[DNAToByte('_'), w[j - 1]],
                        D[i - 1, j] + d[u[i - 1], DNAToByte('_')]));
                }
            }

            PrintMatrix(D);

            matching = GetMatching(D, MatchingType.Minimal);

            return D[n, m];
        }

        public int ComputeSimiliarity(out string[] matching) {
            int n = u.Length;
            int m = w.Length;
            int[,] S = new int[n + 1, m + 1];
            S[0, 0] = 0;

            for (int j = 1; j <= m; j++) {
                for (int k = 0; k < j; k++) {
                    S[0, j] += s[DNAToByte('_'), w[k]];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int k = 0; k < i; k++) {
                    S[i, 0] += s[u[k], DNAToByte('_')];
                }
            }

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    S[i, j] = Math.Max(
                        S[i - 1, j - 1] + s[u[i - 1], w[j - 1]], Math.Max(
                        S[i, j - 1] + s[DNAToByte('_'), w[j - 1]],
                        S[i - 1, j] + s[u[i - 1], DNAToByte('_')]));
                }
            }

            PrintMatrix(S);

            matching = GetMatching(S, MatchingType.Maximal);

            return S[n, m];
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

            Console.Write("".PadLeft(2 * padLenght - 1));
            for (int i = 0; i < n; i++)
                Console.Write(ByteToDNA(u[i]).ToString().PadLeft(padLenght));
            Console.WriteLine();
            Console.Write("".PadLeft(padLenght - 1));

            for (int j = 0; j <= m; j++) {
                if (j > 0)
                    Console.Write(ByteToDNA(w[j - 1]) + " ");
                for (int i = 0; i <= n; i++) {
                    Console.Write(M[i, j].ToString().PadLeft(padLenght));
                }
                Console.WriteLine();
            }
        }


        private string[] GetMatching(int[,] matrix, MatchingType matchingType) {
            int x = n;
            int y = m;
            int lastX = x;
            int lastY = y;

            StringBuilder sequence1 = new StringBuilder();
            StringBuilder sequence2 = new StringBuilder();

            while (x != 0 && y != 0) {
                GetNextField(matrix, matchingType, ref x, ref y);

                if (x == lastX) {
                    sequence1.Insert(0, '_');
                    sequence2.Insert(0, ByteToDNA(w[y]));    // z nierówności trójkąta będzie działać
                }
                else if (y == lastY) {
                    sequence1.Insert(0, ByteToDNA(u[x]));
                    sequence2.Insert(0, '_');
                }
                else {
                    sequence1.Insert(0, ByteToDNA(u[x]));
                    sequence2.Insert(0, ByteToDNA(w[y]));
                }

                lastX = x;
                lastY = y;
            }

            return new string[] { sequence1.ToString(), sequence2.ToString() };
        }

        private static void GetNextField(int[,] matrix, MatchingType matchingType, ref int x, ref int y) {
            if (x == 0) {
                y--;
            }
            else if (y == 0) {
                x--;
            }
            else {
                int val = matrix[x, y];
                int next;

                switch (matchingType) {
                    case MatchingType.Maximal:
                        next = Math.Max(matrix[x - 1, y - 1], Math.Max(matrix[x, y - 1], matrix[x - 1, y]));
                        break;
                    case MatchingType.Minimal:
                    default:
                        next = Math.Min(matrix[x - 1, y - 1], Math.Min(matrix[x, y - 1], matrix[x - 1, y]));
                        break;
                }

                if (matrix[x - 1, y - 1] == next) {
                    x--;
                    y--;
                }
                else if (matrix[x, y - 1] == next) {
                    y--;
                }
                else if (matrix[x - 1, y] == next) {
                    x--;
                }
            }
        }
    }
}
