#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Defining const chars for input file names
const char FILE_NAME_MATRIX[] = "matrix.txt";
const char FILE_NAME_VECTOR[] = "vector.txt";


// Allocating memory for the various array objects
double* AllocateMemory(double* InputTensor)
{
    InputTensor = (double*)calloc(16, sizeof(double));

    return InputTensor;
}

// Reading in a Tensor from an input file
// The file could contain doubles, with '.' or ',' as decimal points.
// Between the doubles any kind of denumerator is allowed, as long
// as they're not the '.', or ',' (decimal points) or the '-' (minus sign) characters.
//
// Also the file should "look like" a Tensor, so every rows and columns
// should have the same amount of numbers. (Imperfect rows/columns are permitted.)
double* ReadInTensor(FILE* Input, double* Tensor, int* NumberOfRows)
{
    // Allocate memory for input tensor
    Tensor = AllocateMemory(Tensor);

    // Temporary storage for numbers in a chararray
    char FileStream;
    char* CurrentNumber = (char*)calloc(16, sizeof(char));

    // Index for elements of the current char array
    int CharIndex = 0;
    // Index for counting elements of output tensor
    int TensorIndex = 0;

    // Index for tensor elements
    double CurrentTensorElement;

    while(1)
    {
        // Read characters one by one from Input file to FileStream char 
        FileStream = getc(Input);

        // If EOF is reached, then append the number in the CurrentNumber buffer to 
        // the Tensor array, then break the loop.
        if(FileStream == EOF)
        {
            // Null terminating string
            CurrentNumber[CharIndex] = '\0';

            // Convert string to double
            CurrentTensorElement = atof(CurrentNumber);

            // Append double to our Tensor
            Tensor[TensorIndex] = CurrentTensorElement;
            
            printf("EOF reached!\n");

            break;
        }

        // Counting rows in the Tensor and incrementing the NumberOfRows
        // integer, everytime a '\n' character is reached
        if(FileStream == '\n')
        {
            *NumberOfRows += 1;
        }

        // If the current char from the input file is a number, or a decimal point/ +- sign, then
        // append it to the CurrentNumber char array buffer
        if((FileStream >= '0' && FileStream <= '9') || FileStream == '.' || FileStream == ',' || FileStream == '-')
        {
            CurrentNumber[CharIndex] = FileStream;
            CharIndex += 1;

            if(CharIndex >= (sizeof(CurrentNumber)/sizeof(CurrentNumber[0])))
            {
                CurrentNumber = (char*)realloc(CurrentNumber, ((sizeof(CurrentNumber)/sizeof(CurrentNumber[0])) * 2));
            }
        }

        // If it isn't, then null terminate the CurrentNumber buffer, convert it to integer,
        // and then append it to the Tensor double* array
        else
        {
            // Null terminating string
            CurrentNumber[CharIndex] = '\0';

            // Convert string to double
            CurrentTensorElement = atof(CurrentNumber);

            // Append double to our Tensor
            Tensor[TensorIndex] = CurrentTensorElement;

            // Relaunch indexing for a new Tensor element
            CharIndex = 0;

            // Temp. counter for elements of the Tensor
            TensorIndex += 1;

            if(TensorIndex >= (sizeof(Tensor)/sizeof(Tensor[0])))
            {
                Tensor = (double*)realloc(Tensor, ((sizeof(Tensor)/sizeof(Tensor[0])) * 2));
            }
        }
    }

    free(CurrentNumber);

    return Tensor;
}


// Checking if the matrix and the vector have correct shapes
//
// To broadcast a matrix and vector together, the vector's lenght should be
// equal to the number of the matrix's columns, if the matrix is multiplied from
// the left side.
// Similarly, if the matrix is multiplied from the right, the vector's lenght
// should be equal to the number of the matrix's rows.
// In the program we multiply our matrix only from the left side
void ConsistencyCheck(int Columns, int RowsVc)
{
    if(Columns == RowsVc)
    {
        printf("Sizeof(Columns): %d\nSizeof(Vector): %d", Columns, RowsVc);
    }

    else
    {
        perror("Matrix and vector are not in a correct shapes!");
        EXIT_FAILURE;
    }
}


// Multiplying given matrix and vector
// If the input matrix have the shape M*N and the input vector have the shape N*1
// Then the output will have shape M*1
// Where M is the number of the rows in the matrix, and N is the number of the columns in the matrix
double* MultiplyThings(int Rows, int Columns, int RowsVc, double* Matrix, double* Vector, double* Output)
{
    ConsistencyCheck(Columns, RowsVc);

    Output = (double*)calloc((Rows * 1), sizeof(double));

    for (int i = 0; i < Rows; i++)
    {
        for(int j = 0; j < Columns; j++)
        {
            Output[i] = (Matrix[Columns * i + j] * Vector[j]);
        }
    }

    return Output;
}


// Printing a matrix with M rows and N columns to the standard output
void PrintOutput(double* Tensor, int Rows, int Columns)
{
    for(int i = 0; i < Rows; i++)
    {
        for(int j = 0; j < Columns; j++)
        {
            printf("%lf ", Tensor[Columns * i  + j]);
        }

        printf("\n");
    }
}


int main()
{

    // ================== Initializing
    // Allocating memory and defining numerators

    // Allocate memory for used array-type variables
    double* Matrix;
    double* Vector;
	double* Output;

    //printf("OK");

    // Index for counting rows of our matrix
    int NumberOfRowsMx = 1;
    int NumberOfRowsVc = 0; // NOT USED, ONLY AUXILIARY

    // ================== Matrix
    // Reading in matrix from a file to the Matrix double* array

     // Opening file - which containing the matrix - in read-only mode
    FILE* Input_Mx = fopen("matrix.txt", "r");

    // Read in matrix from input file
    // The ReadInTensor function returns a double* array
    Matrix = ReadInTensor(Input_Mx, Matrix, &NumberOfRowsMx);

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            printf("%lf ", Matrix[4 * i + j]);
        }
        printf("\n");
    }

    // Close the input file
    fclose(Input_Mx);


    // ================== Vector
    // Reading in vector from a file to the Vector double* array

    // Opening file - which containing the vector - in read-only mode
    FILE* Input_Vc = fopen("vector.txt", "r");

    // Read in vectorfrom input file
    // The ReadInOneDArray function returns a double* array
    Vector = ReadInTensor(Input_Vc, Vector, &NumberOfRowsVc);

    // Close the input file
    fclose(Input_Vc);


    // ================== Remaining calculations
    // Multiply matrix and vector, then print output to stdout

    // Calculate constants for this particular problem
    // Size of the vector
    NumberOfRowsVc = sizeof(Vector) / sizeof(Vector[0]);
    // Size of the full matrix
    int SizeOfMatrix = sizeof(Matrix) / sizeof(Matrix[0]);
    // Column length of the matrix
    int NumberOfColumnsMx = SizeOfMatrix / NumberOfRowsMx;

    // Multiply the matrix and the vector if it's possible
	Output = MultiplyThings(NumberOfRowsMx, NumberOfColumnsMx, NumberOfRowsVc, Matrix, Vector, Output);

    free(Matrix);
    free(Vector);
    
    // Print output
    PrintOutput(Output, NumberOfRowsMx, NumberOfRowsVc);

    free(Output);

    return 0;
}