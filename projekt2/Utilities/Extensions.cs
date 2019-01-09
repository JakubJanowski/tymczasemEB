using Models;
using System;
using System.Collections.Generic;

namespace Utilities {
    static class Extensions {
        public static IEnumerable<string> SplitInParts(this string text, int partLength) {
            if (text == null)
                throw new ArgumentNullException("text");
            if (partLength <= 0)
                throw new ArgumentException("Part length has to be positive.", "partLength");

            for (var i = 0; i < text.Length; i += partLength)
                yield return text.Substring(i, Math.Min(partLength, text.Length - i));
        }

        public static IEnumerable<Sequence> FindAminoAcidSequences(this byte[] rna, byte start, byte stop) {
            int startIdx = 0;
            int stopIdx = 0;
            while (true) {
                startIdx = Array.IndexOf(rna, start, stopIdx);
                if (startIdx == -1)
                    break;
                startIdx++;
                stopIdx = Array.IndexOf(rna, stop, startIdx);
                if (stopIdx == -1)
                    break;

                yield return new Sequence(rna.SubArray(startIdx, stopIdx - startIdx), startIdx);
            }
        }

        public static T[] SubArray<T>(this T[] data, int index, int length)
        {
            T[] result = new T[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }
    }
}
