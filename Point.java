public class Point {
        float x, y;

        public Point() {
                x = 0;
                y = 0;
        }

        public Point (float x, float y) {
                this.x  = x;
                this.y = y;
        }

        public Point(Point a) {
                x = a.x;
                y = a.y;
        }

        public void made_zero() {
                x = 0;
                y = 0;
        }

}
