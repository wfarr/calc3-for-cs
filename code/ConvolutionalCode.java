import java.util.Random;

public class ConvolutionalCode
{
    public static void main (String[] args)
    {

        Matrix input, output, Y0, Y1, A, B;
    
        //Encode specific input sequence
        double[] vector = {1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0};
        input = new Matrix(vector);
        System.out.print(input);
        System.out.println("Encoding...");
    
        //Encodes the Y0 output stream for the input
        A = new Matrix(input.length, input.length);
        A = A.fillMatrixZero();
        Y0 = A.matrixMultiply(input);
        System.out.print(Y0);

        //Encode the Y1 output stream for the input
        B = new Matrix(input.length, input.length);
        B = B.fillMatrixOne();
        Y1 = B.matrixMultiply(input);
        System.out.print(Y1);

        //Generates a random input stream and encodes both output streams
        input = Matrix.generateInput();
        A = new Matrix(input.length, input.length);
        A = A.fillMatrixZero();
        B = new Matrix(input.length, input.length);
        B = B.fillMatrixOne();
        Y0 = A.matrixMultiply(input);
        Y1 = B.matrixMultiply(input);
        System.out.println("Generating input...");
        System.out.print(input);
        System.out.print("Encoding...");
        System.out.print(Y0);
        System.out.print(Y1);

        //Generates two random output stream and decodes their input streams using
        //Jacobi and Gauss-Seidel iterations. The first random stream is decoded as
        //a Y0 stream and the second as a Y1 stream. Both iterative methods are fun 8 times.
        int maxit = 8;
        System.out.println("Generating random output...");
        output = Matrix.generateOutput();
        A = new Matrix(output.length, output.length);
        A = A.fillMatrixZero();
        B = new Matrix(output.length, output.length);
        B = B.fillMatrixOne();
        System.out.print(output);
        System.out.println("Decoding...");
        input = output.jacobiIteration(A, maxit);
        System.out.print(input);
        input = output.gaussseidelIteration(A, maxit);
        System.out.print(input);

        System.out.println("Generating random output...");
        output = Matrix.generateOutput();
        A = new Matrix(output.length, output.length);
        A = A.fillMatrixZero();
        B = new Matrix(output.length, output.length);
        B = B.fillMatrixOne();
        System.out.print(output);
        System.out.println("Decoding...");
        input = output.jacobiIteration(B, maxit);
        System.out.print(input);
        input = output.gaussseidelIteration(B, maxit);
        System.out.print(input);

        System.out.print("Done");    
    }
}