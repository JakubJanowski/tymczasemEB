using DNA;
using System;
using System.IO;
using System.Linq;
using System.Text;

namespace projekt1 {
    class Program {
        private const string dataDirectory = @"..\..\..\data\set1\";

        static void Main(string[] args) {
            string sequence1 = File.ReadAllText(dataDirectory + "sequence1.txt");
            string sequence2 = File.ReadAllText(dataDirectory + "sequence2.txt");
            int[,] distanceMatrix = GetMatrix(dataDirectory + "distanceMatrix.txt");
            int[,] similiarityMatrix = GetMatrix(dataDirectory + "similiarityMatrix.txt");

            DNAMatcher matcher = new DNAMatcher(sequence1, sequence2, distanceMatrix, similiarityMatrix);
            
            int editDistance = matcher.ComputeEditDistance(out string[] matching);
            PrintMatching(matching);
            Console.WriteLine($"Optimal edit distance: {editDistance}\n");

            int similiarity = matcher.ComputeSimiliarity(out matching);
            PrintMatching(matching);
            Console.WriteLine($"Optimal global similiarity: {similiarity}\n");

            Console.ReadKey();
        }

        private static int[,] GetMatrix(string filepath) {
            string[] matrixLines = File.ReadAllLines(filepath);
            int[,] matrix = new int[5, 5];

            for (int y = 0; y < 5; y++) {
                int[] row = matrixLines[y].Split(new char[] { ' ' }).Select(s => int.Parse(s)).ToArray();
                for (int x = 0; x < 5; x++)
                    matrix[y, x] = row[x];
            }
            return matrix;
        }

        private static void PrintMatching(string[] matching) {
            Console.WriteLine();
            Console.WriteLine(matching[0]);
            Console.WriteLine(matching[1]);
            Console.WriteLine();
        }
    }
}
