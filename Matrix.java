import java.util.Random;

import static java.lang.Math.abs;
import static java.lang.Math.random;

public class Matrix {
        int rows, cols;
        double[][] data;

    public Matrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.data =  new double[rows][cols];
        for(int i = 0; i < rows; ++i) {
          data[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
          data[i][j] = 0;
    }
   }
    }

    public Matrix() {
        this(1,1);
    }

    Matrix make_identity_homography()
    {
        Matrix H = make_matrix(3,3);
        H.data[0][0] = 1;
        H.data[1][1] = 1;
        H.data[2][2] = 1;
        return H;
    }


  public double[][] solve(double[][] mat, double[][] constants)
        {

            int n  = mat.length;
            //inverse of matrix mat[][]
            double inverted_mat[][] = invert(mat);

            //Multiplication of mat inverse and constants
            double result[][] = new double[n][1];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        result[i][j] = result[i][j] + inverted_mat[i][k] * constants[k][j];
                    }
                }
            }
    return result;
        }



// Method to carry out the partial-pivoting Gaussian
// elimination.  Here index[] stores pivoting order.

    public static void gaussian(double[][] a, int[] index)
    {
        int n = index.length;
        double c[] = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            double c1 = 0;
            for (int j=0; j<n; ++j)
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i)
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1)
                {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = (float) pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }

    public static double[][] invert(double[][] a)
    {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i)
            b[i][i] = 1;

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                            -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }

    Matrix make_translation_homography(float dx, float dy)
    {
        Matrix H = make_identity_homography();
        H.data[0][2] = dx;
        H.data[1][2] = dy;
        return H;
    }

    Matrix make_matrix(int rows, int cols)
    {
        Matrix m = new Matrix(rows, cols);
        return m;
    }

    Matrix copy_matrix()
    {
        int i,j;
        Matrix c = make_matrix(rows, cols);
        for(i = 0; i < rows; ++i){
            for(j = 0; j < cols; ++j){
                c.data[i][j] = data[i][j];
            }
        }
        return c;
    }

    Matrix augment_matrix()
    {
        Matrix c = make_matrix(rows, cols*2);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                c.data[i][j] = data[i][j];
            }
        }
        for(int j = 0; j < rows; ++j){
            c.data[j][j+cols] = 1;
        }
        return c;
    }

    Matrix make_identity(int rows, int cols)
    {
        Matrix m = make_matrix(rows, cols);
        for(int i = 0; i < rows && i < cols; ++i){
            m.data[i][i] = 1;
        }
        return m;
    }

    Matrix matrix_mult_matrix(Matrix b)
    {
        assert(cols == b.rows);
        Matrix p = make_matrix(rows, b.cols);
        for(int i = 0; i < p.rows; ++i){
            for(int j = 0; j < p.cols; ++j){
                p.data[i][j] = 0;
                for(int k = 0; k < cols; ++k){
                    p.data[i][j] += data[i][k]*b.data[k][j];
                }
            }
        }
        return p;
    }

    Matrix matrix_sub_matrix(Matrix b)
    {
        assert(cols == b.cols);
        assert(rows == b.rows);
        int i, j;
        Matrix p = make_matrix(rows, cols);
        for(i = 0; i < p.rows; ++i){
            for(j = 0; j < p.cols; ++j){
                p.data[i][j] = data[i][j] - b.data[i][j];
            }
        }
        return p;
    }

    Matrix transpose_matrix()
    {
        Matrix t = new Matrix();
        t.rows = cols;
        t.cols = rows;
        t.data =  new double[t.rows][t.cols];
        int i, j;
        for(i = 0; i < t.rows; ++i){
            t.data[i] = new double[t.cols];
            for(j = 0; j < t.cols; ++j){
                t.data[i][j] = data[j][i];
            }
        }
        return t;
    }

    void scale_matrix(double s)
    {
        int i, j;
        for(i = 0; i < rows; ++i){
            for(j =0 ; j < cols; ++j){
                data[i][j] *= s;
            }
        }
    }

    double[] matrix_mult_vector(double[] v)
    {
        double[] p =  new double[rows];
        for(int i = 0; i < rows; i++){
            p[i] = 0;
        }
        int i, j;
        for(i = 0; i < rows; ++i){
            for(j = 0; j < cols; ++j){
                p[i] += data[i][j]*v[j];
            }
        }
        return p;
    }

    void print_matrix()
    {
        int i, j;
        System.out.print(" __");
        for(j = 0; j < 16*cols-1; ++j) {
            System.out.print(" ");
        }
        System.out.print("__ \n");

        System.out.print("|  ");
        for(j = 0; j < 16*cols-1; ++j)
            System.out.print(" ");
        System.out.print("  |\n");

        for(i = 0; i < rows; ++i){
            System.out.print("|  ");
            for(j = 0; j < cols; ++j){
                System.out.print(data[i][j]);
                System.out.print("        ");
            }
            System.out.print(" |\n");
        }
        System.out.print("|__");
        for(j = 0; j < 16*cols-1; ++j) System.out.print(" ");
        System.out.print("__|\n");
    }

    double[] LUP_solve(Matrix U, int[] p, double[] b)
    {
        int i, j;
        double[] c = new double[rows];
        for(i = 0; i < rows; ++i){
            int pi = p[i];
            c[i] = b[pi];
            for(j = 0; j < i; ++ j){
                c[i] -= data[i][j]*c[j];
            }
        }
        for(i = U.rows-1; i >= 0; --i){
            for(j = i+1; j < U.cols; ++j){
                c[i] -= U.data[i][j]*c[j];
            }
            c[i] /= U.data[i][i];
        }
        return c;
    }

    Matrix matrix_invert()
    {
        //print_matrix(m);
        Matrix none = new Matrix();
        if(rows != cols){
            System.out.print("Matrix not square\n");
            return none;
        }
        Matrix c = augment_matrix();
        //print_matrix(c);

        int i, j, k;
        for(k = 0; k < c.rows; ++k){
            double p = 0.;
            int index = -1;
            for(i = k; i < c.rows; ++i){
                double val = abs(c.data[i][k]);
                if(val > p){
                    p = val;
                    index = i;
                }
            }
            if(index == -1){
                System.out.print("Can't do it, sorry!\n");
                return none;
            }

            double[] swap = c.data[index];
            c.data[index] = c.data[k];
            c.data[k] = swap;

            double val = c.data[k][k];
            c.data[k][k] = 1;
            for(j = k+1; j < c.cols; ++j){
                c.data[k][j] /= val;
            }
            for(i = k+1; i < c.rows; ++i){
                float s = (float) -c.data[i][k];
                c.data[i][k] = 0;
                for(j = k+1; j < c.cols; ++j){
                    c.data[i][j] +=  s*c.data[k][j];
                }
            }
        }
        for(k = c.rows-1; k > 0; --k){
            for(i = 0; i < k; ++i){
                double s = -c.data[i][k];
                c.data[i][k] = 0;
                for(j = k+1; j < c.cols; ++j){
                    c.data[i][j] += s*c.data[k][j];
                }
            }
        }
        //print_matrix(c);
        Matrix inv = make_matrix(rows, cols);
        for(i = 0; i < rows; ++i){
            for(j = 0; j < cols; ++j){
                inv.data[i][j] = c.data[i][j+cols];
            }
        }
        //print_matrix(inv);
        return inv;
    }

    int[] in_place_LUP()
    {
        int[] pivot = new int[rows];
        if(rows != cols){
            System.out.print("Matrix not square\n");
            return pivot;
        }

        for(int k = 0; k < rows; ++k) pivot[k] = k;
        for(int k = 0; k < rows; ++k){
            double p = 0.;
            int index = -1;
            for(int i = k; i < rows; ++i){
                double val = abs(data[i][k]);
                if(val > p){
                    p = val;
                    index = i;
                }
            }
            if(index == -1){
                System.out.print("Matrix is singular\n");
                return pivot;
            }

            int swapi = pivot[k];
            pivot[k] = pivot[index];
            pivot[index] = swapi;

            double[] swap = data[index];
            data[index] = data[k];
            data[k] = swap;

            for(int i = k+1; i < rows; ++i){
                data[i][k] = data[i][k]/data[k][k];
                for(int j = k+1; j < cols; ++j){
                    data[i][j] -= data[i][k] * data[k][j];
                }
            }
        }
        return pivot;
    }

    Matrix random_matrix(int rows, int cols)
    {
        Matrix m = make_matrix(rows, cols);
        int i, j;
        for(i = 0; i < rows; ++i){
            for(j = 0; j < cols; ++j){
                m.data[i][j] = random() %100 - 50;
            }
        }
        return m;
    }

    double[] sle_solve(double[] b)
    {
        int[] p = this.in_place_LUP();
        return this.LUP_solve(this, p, b);
    }

    Matrix solve_system(Matrix b)
    {
        Matrix none = new Matrix();
        Matrix Mt = this.transpose_matrix();
        Matrix MtM = Mt.matrix_mult_matrix(this);
        Matrix MtMinv = MtM.matrix_invert();
        if(MtMinv.data == null) return none;
        Matrix Mdag = MtMinv.matrix_mult_matrix(Mt);
        Matrix a = Mdag.matrix_mult_matrix(b);

        return a;
    }

    void test_matrix()
    {
        int i;
        for(i = 0; i < 100; ++i){
            int s = (int) (random()%4 + 3);
            Matrix m = random_matrix(s, s);
            Matrix inv = m.matrix_invert();
            Matrix res = inv.matrix_mult_matrix(m);
            res.print_matrix();
        }
    }
}
