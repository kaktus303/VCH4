import java.rmi.MarshalException;
import java.text.DecimalFormat;
import java.util.Arrays;

public class Matrix {
    private double[][] matrix;
    private int size;


    public Matrix(int userSize)
    {
     matrix = new double[userSize][userSize];
     size = userSize;
    }
    public Matrix(int userSize, int flag)
    {
        matrix = new double[userSize][userSize];
        size = userSize;
        if(flag == 2)
        {
            for(int i = 0;i<userSize;++i)
            {
                for(int j = 0;j<userSize;++j)
                {
                    matrix[i][j] = 1.0 / (i+j+1);
                }
            }
        }
        if(flag == 1)
        {
            for(int i = 0;i<userSize;++i)
            {
                for(int j = 0;j<userSize;++j)
                {
                    if(i == j) matrix[i][j] = 1;
                    else matrix[i][j] = 0;
                }
            }
        }
        if(flag == 0)
        {
            for(int i = 0;i<userSize;++i)
            {
                for(int j = 0;j<userSize;++j)
                {
                    matrix[i][j] = 0;
                }
            }
        }
        if(flag == 3)
        {
            for(int i = 0;i<userSize;++i)
            {
                for(int j = 0;j<userSize;++j)
                {
                    if(i>=j)
                    {
                        matrix[i][j] =  userSize - i;

                    }
                    else {
                        matrix[i][j] = userSize - j;
                    }
                }
            }
        }
    }
    public int getSize() {
        return size;
    }

    public void setElement(int row, int col, double newElement)
    {
        matrix[row][col] = newElement;
    }
    public double getElemnt(int row, int col)
    {
        return matrix[row][col];
    }
    public Matrix sumMatrix(Matrix A, Matrix B)
    {
        Matrix sumMatrix = new Matrix(A.getSize());
        for(int i = 0;i<A.getSize();++i)
        {
            for(int j = 0;j<A.getSize();++j)
            {
               sumMatrix.setElement(i,j,A.getElemnt(i,j) + B.getElemnt(i,j));
            }
        }
        return sumMatrix;
    }
    public Matrix multMatrix(Matrix B)
    {
        double newElemnet = 0;
        Matrix multMatrix = new Matrix(this.getSize());
        for(int i = 0;i<this.getSize();++i)
        {
            for(int j = 0;j<this.getSize();++j)
            {
                for(int k = 0;k<this.getSize();++k) {
                    newElemnet += this.getElemnt(i,k)*B.getElemnt(k,j);
                }
                multMatrix.setElement(i, j, newElemnet);
                newElemnet = 0;
            }
        }
        return multMatrix;

    }
    public void collect(double epsilon)
    {
        for(int i = 0;i <size;++i)
        {
            for(int j = 0;j<size;++j)
            {
                if(matrix[i][j]<=epsilon)
                {
                    matrix[i][j] = 0;
                }
            }
        }
    }
    public Vector trace()
    {
        Vector trace = new Vector(size);
        for(int i = 0;i<size;++i)
        {
            trace.setElement(i,this.getElemnt(i,i));
        }
        return trace;
    }
    public Matrix transposition()
    {
        Matrix A = new Matrix(size,0);
        for(int i = 0;i<this.getSize();++i) {
            for (int j = 0; j < this.getSize(); ++j)
            {
                A.setElement(i,j,this.getElemnt(j,i));
            }
        }
        return A;
    }
    public void equal(Matrix A)
    {
        for(int i = 0;i<size;++i)
        {
            for(int j= 0;j<size;++j)
            {
                this.setElement(i,j,A.getElemnt(i,j));
            }
        }
    }

    public Vector multOnVector(Vector b)
    {
        double sum =0;
     Vector answer = new Vector(b.getSize());
        for(int i = 0;i<this.getSize();++i)
        {
            for(int j = 0;j<this.getSize();++j)
            {
               sum += this.getElemnt(i,j)*b.getElement(j);
            }
            answer.setElement(i,sum);
            sum = 0;
        }
        return answer;
    }
    static void getCofactor(double mat[][], double temp[][],
                            int p, int q, int n)
    {
        int i = 0, j = 0;

        // Looping for each element
        // of the matrix
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                // Copying into temporary matrix
                // only those element which are
                // not in given row and column
                if (row != p && col != q) {
                    temp[i][j++] = mat[row][col];
                    // Row is filled, so increase
                    // row index and reset col index
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    /* Recursive function for finding determinant
    of matrix. n is current dimension of mat[][]. */
    static double determinantOfMatrix(double mat[][], int n,int N)
    {
        int D = 0; // Initialize result

        // Base case : if matrix
        // contains single element
        if (n == 1)
            return mat[0][0];

        // To store cofactors
        double temp[][] = new double[N][N];

        // To store sign multiplier
        int sign = 1;

        // Iterate for each element of first row
        for (int f = 0; f < n; f++) {
            // Getting Cofactor of mat[0][f]
            getCofactor(mat, temp, 0, f, n);
            D += sign * mat[0][f]
                    * determinantOfMatrix(temp, n - 1,N);

            // terms are to be added
            // with alternate sign
            sign = -sign;
        }

        return D;
    }
    public double det()
    {
        return determinantOfMatrix(this.matrix,size,size);
    }
    public double determinant(Matrix matrix1, int size) {
        double det = 0;
        if (size == 2) {

            det = (matrix1.getElemnt(0, 0) * matrix1.getElemnt(1, 1)) - (matrix1.getElemnt(0, 1) * matrix1.getElemnt(1, 0));
            return det;
        } else if (size == 1) {

            return matrix1.getElemnt(0, 0);
        } else {
            Matrix f = new Matrix(size - 1, 1);
            int a, b;
            for (int i = 0; i < size; ++i) {
                a = 0;
                for (int j = 1; j < size; ++j) {
                    b = 0;
                    for (int k = 0; k < size; ++k) {
                        if (k != i) {
                            f.setElement(a, b, matrix1.getElemnt(j, k));

                            b++;
                        }
                    }
                    a++;
                }
                det += Math.pow(-1, i) * matrix1.getElemnt(0, i) * determinant(f, size - 1);
            }
            System.out.println(4);
            return det;
        }

    }
    public double delta()
    {
        double max=0, sum = 0;
        for(int i = 0;i<this.getSize();++i)
        {
            sum =0;
            for(int j = 0;j<this.getSize();++j)
            {
                if(i!=j)
                {
                    sum+=Math.abs(this.getElemnt(i,j));
                }
            }
            if(max< sum) max = sum;
        }
        return max;
    }
    @Override
    public String toString() {
        DecimalFormat df = new DecimalFormat("0.000000000000000E0");
        String str = new String();
        str = "{\n";
        for(int i = 0;i<size;++i)
        {
            for(int j = 0;j<size;++j)
            {
                str = str + " " +df.format(matrix[i][j]) + " ";
            }
            str += "\n";
        }
        str += "}\n";
        return str;
    }
}
