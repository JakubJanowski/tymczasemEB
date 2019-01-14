namespace Models {
    public class Field {
        public int x;
        public int y;

        public Field(int x, int y) {
            this.x = x;
            this.y = y;
        }

        public Field(Field field) {
            Copy(field);
        }

        public Field Clone() {
            return new Field(this);
        }

        public void Copy(Field field) {
            x = field.x;
            y = field.y;
        }
    }
}
