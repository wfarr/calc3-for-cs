import java.util.Random;

public class Matrix
{
    //instance data for matrices
    double[][] entry;
    int length;
    
    //constructor for matrix object
    public Matrix(int m, int n) {
        entry = new double[m][n];
        length = entry.length;
    }
    
    //a second matrix constructor
    public Matrix (double[] array) {
        entry = new double[array.length][1];
        for (int i = 0; i < array.length; i++)
            entry[i][0] = array[i];
        length = entry.length;
    }
    
    //fills a matrix to correctly encode the y0 stream    
    public Matrix fillMatrixZero() {
        int length = this.length;
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++)
                if (i >= j) {
                    if (i == j) {
                        this.entry[i][j] = 1;
                    }
                    if (i == j + 2) {
                        this.entry[i][j] = 1; }
                    if (i == j + 3) {
                        this.entry[i][j] = 1; }
                } else {
                    this.entry[i][j] = 0;
                }
        return this;
    }
    
    //fills a matrix to encode the y1 stream
    public Matrix fillMatrixOne() {
        int length = this.length;
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i >= j) {
                    if (i == j) {
                        this.entry[i][j] = 1;
                    }
                    if (i == j + 1) {
                        this.entry[i][j] = 1;
                    }
                    if(i == j + 3) {
                        this.entry[i][j] = 1;
                    }
                } else {
                    this.entry[i][j] = 0;
                }
            }
        }
        return this;
    }

    //randomly generates an input stream to encode
    public static Matrix generateInput() {
        int i = 0;
        Random generator = new Random();
        Matrix x = new Matrix(Math.abs(generator.nextInt(10) + 4), 1);
        for (i = 0; i < x.length; i++) {
            x.entry[i][0] = generator.nextInt(2);
        }
        x.entry[i - 1][0] = 0;
        x.entry[i - 2][0] = 0;
        x.entry[i - 3][0] = 0;
        return x;
    }
    
    //randomly generates an output stream to decode    
    public static Matrix generateOutput() {
        int i = 0;
        Random generator = new Random();
        Matrix x = new Matrix(generator.nextInt(10) + 4, 1);
        for (i = 0; i < x.length; i++) {
            x.entry[i][0] = generator.nextInt(2);
        }
        return x;
    }
    
    //multiples two matrices together
    public Matrix matrixMultiply(Matrix B) {
        int k,i,j;
        Matrix C = new Matrix(this.length, B.entry[0].length);
        for (k = 0; k < B.entry[0].length; k++) {
            for (i = 0; i < this.length; i++) {
                for (j = 0; j < B.length; j++) {
                    C.entry[i][k] = ((this.entry[i][j] * B.entry[j][k]) + C.entry[i][k]) % 2;
                }
            }
        }
        return C;
    }

    //Jacobi iteration for binary matrices
    public Matrix jacobiIteration(Matrix A, int maxit) {
        int k, i, j;
        
        Matrix result = new Matrix (A.length, maxit);
        for (k = 0; k < maxit; k++) {
            for (i = 0; i < this.length; i++) {
                double sum = 0;
                for (j = 0; j < A.length; j++) {
                    if (i != j) {
                        if (k == 0) {
                            sum = 0;
                        } else {
                            sum = ((A.entry[i][j] * result.entry[j][k - 1]) + sum) % 2;
                        }
                    }
                }
                result.entry[i][k] = (1 / A.entry[i][i]) * ((Math.abs(this.entry[i][0] - sum))) % 2;
            }
        } 
        return result;
    }
    
    //Gauss-Seidel Iteration for binary matrices
    public Matrix gaussseidelIteration(Matrix A, int maxit) {
        int k, i, j;
        
        Matrix result = new Matrix (A.length, maxit);
        for (k = 0; k < maxit; k++) {
            for (i = 0; i < this.length; i++) {
                double sum = 0;
                for (j = 0; j <= i-1; j++) {
                    sum = ((A.entry[i][j] * result.entry[j][k]) + sum) % 2;    
                }
                for (j = i + 1; j < A.length; j++) {    
                    if (k == 0) {
                        sum = (A.entry[i][j] * 0 + sum) % 2;
                    } else {
                        sum = ((A.entry[i][j] * result.entry[j][k - 1]) + sum) % 2;
                    }
                    result.entry[i][k] = (1 / A.entry[i][i]) * ((Math.abs(this.entry[i][0] - sum))) % 2;
                }
            }
        }
        return result;
    }
            
    //Method to print out matrices to the console
    public String toString() {
        for (int i = 0; i < this.length; i++) {
            for(int j = 0; j < entry[i].length; j++) {
                System.out.print("\t" + entry[i][j]);
            }
            System.out.print("\n");
        }
        return (this.length + " by " + this.entry[0].length + "\n\n");
    }
}