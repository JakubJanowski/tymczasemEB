namespace Models {
    internal class Sequence {
        public byte[] data;
        public int index;

        public Sequence(byte[] data, int index) {
            this.data = data;
            this.index = index;
        }
    }
}