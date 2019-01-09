using Models;
using System.Linq;

namespace Utilities {
    static class MathUtils {
        public static int Min(params int[] values) {
            return Enumerable.Min(values);
        }

        public static int Max(params int[] values) {
            return Enumerable.Max(values);
        }

        public static int Max(int[,] array, out Field field) {
            int n = array.GetLength(0);
            int m = array.GetLength(1);
            int max = 0;
            field = new Field(0, 0);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    if (max < array[i, j]) {
                        max = array[i, j];
                        field.x = i;
                        field.y = j;
                    }
                }
            }

            return max;
        }
    }
}
