public class Match {
        Point p, q;
        int ai, bi;
        float distance;

    public Match() {
        p = new Point();
        q = new Point();
        ai = 0;
        bi = 0;
        distance = 0;
    }

    public void made_zero() {
                p = new Point();
                q = new Point();
                p.made_zero();
                q.made_zero();
                ai = 0;
                bi = 0;
                distance = 0;
        }

}
