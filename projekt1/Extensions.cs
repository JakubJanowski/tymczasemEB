using System;
using System.Collections.Generic;

namespace DNA {
    static class Extensions {
        public static IEnumerable<string> SplitInParts(this string text, int partLength) {
            if (text is null)
                throw new ArgumentNullException("text");
            if (partLength <= 0)
                throw new ArgumentException("Part length has to be positive.", "partLength");

            for (var i = 0; i < text.Length; i += partLength)
                yield return text.Substring(i, Math.Min(partLength, text.Length - i));
        }
    }
}
