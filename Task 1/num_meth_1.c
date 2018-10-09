#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Defining const chars for input file names
const char FILE_NAME_MATRIX[] = "../matrix.txt";
const char FILE_NAME_VECTOR[] = "../vector.txt";




// Allocating memory for the various array objects
static double* CallocateMemory()
{
    double* InputTensor;

    InputTensor = (double*)calloc(64, sizeof(double));

    if(InputTensor == NULL)
    {
        perror("Tensor allocation error");
        EXIT_FAILURE;
    }

    return InputTensor;
}

// Allocating memory for char arrays
static char* CallocateMemoryChar()
{
    char* InputCharArray;

    InputCharArray = (char*)calloc(64, sizeof(char));

    if(InputCharArray == NULL)
    {
        perror("Char array allocation error");
        EXIT_FAILURE;
    }

    return InputCharArray;
}


// Reading in a Tensor from an input file
// The file could contain doubles, with '.' or ',' as decimal points.
// Between the doubles any kind of denumerator is allowed, as long
// as they're not the '.', or ',' (decimal points) or the '-' (minus sign) characters.
//
// Also the file should "look like" a Tensor, so every rows and columns
// should have the same amount of numbers. (Imperfect rows/columns are permitted.)
static double* ReadInTensor(FILE* Input, double* Tensor, unsigned int* NumberOfRows, unsigned int* TensorIndex)
{

    // Temporary storage for current character, read from input
    char FileStream;

    // Temporary storage for current number in a char array
    char* CurrentNumber = CallocateMemoryChar();

    printf("%llu\n", sizeof(CurrentNumber));

    // Temporary storage for current tensor element after converting char array (string) to double
    double CurrentTensorElement;

    // Index for elements of the current char array
    unsigned int CharIndex = 0;

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
            Tensor[*TensorIndex] = CurrentTensorElement;

            // Enumerating the final element 
            *TensorIndex += 1;

            printf("\nEOF reached!\n");

            break;
        }

        // Counting rows in the Tensor and incrementing the NumberOfRows
        // integer, everytime a '\n' character is reached
        if(FileStream == '\n')
        {
            (*NumberOfRows) += 1;
        }

        // If the current char from the input file is a number, or a decimal point/ +- sign, then
        // append it to the CurrentNumber char array buffer
        if((FileStream >= '0' && FileStream <= '9') || FileStream == '.' || FileStream == ',' || FileStream == '-')
        {
            CurrentNumber[CharIndex] = FileStream;
            CharIndex += 1;

            if(CharIndex >= (sizeof(CurrentNumber)/sizeof(CurrentNumber[0])))
            {
                CurrentNumber = (char*)realloc(CurrentNumber, (sizeof(CurrentNumber) * 2));

                if(CurrentNumber == NULL)
                {
                    perror("Char reallocation error");
                    EXIT_FAILURE;
                }
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
            Tensor[*TensorIndex] = CurrentTensorElement;

            // Relaunch indexing for a new Tensor element
            CharIndex = 0;

            // Temp. counter for elements of the Tensor
            *TensorIndex += 1;

            if(*TensorIndex >= (sizeof(Tensor)/sizeof(Tensor[0])))
            {
                Tensor = (double*)realloc(Tensor, (sizeof(Tensor) * 2));

                if(Tensor == NULL)
                {
                    perror("Tensor allocation error");
                    EXIT_FAILURE;
                }
            }
        }
    }

    // Close the input file
    fclose(Input);

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
static void ConsistencyCheck(unsigned int Columns, unsigned int RowsVc)
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
static double* MultiplyThings(unsigned int Rows, unsigned int Columns, unsigned int RowsVc, double* Matrix, double* Vector)
{
    // Check before allocating memory for new array
    ConsistencyCheck(Columns, RowsVc);

    // Define and allocate an array for the output
    double* Output;
    Output = (double*)calloc((Rows * 1), sizeof(double));

    // Multiply vector and matrix
    for (unsigned int i = 0; i < Rows; i++)
    {
        for(unsigned int j = 0; j < Columns; j++)
        {
            Output[i] = (Matrix[Columns * i + j] * Vector[j]);
        }
    }

    return Output;
}


// Printing a matrix with M rows and N columns to the standard output
static void PrintOutput(double* Tensor, unsigned int Rows, unsigned int Columns)
{
    for(unsigned int i = 0; i < Rows; i++)
    {
        for(unsigned int j = 0; j < Columns; j++)
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

        // Allocate memory for input tensors
    Matrix = CallocateMemory();
    Vector = CallocateMemory();

    // Index for counting rows of our matrix
    unsigned int NumberOfRowsMx = 1;
    // Index for counting elements of output tensors
    unsigned int TensorIndexMx = 0;
    unsigned int TensorIndexVc = 0;

    unsigned int NumberOfRowsVcAUX = 0; // NOT USED, ONLY AUXILIARY

    // ================== Matrix
    // Reading in matrix from a file to the Matrix double* array

    // Allocate memory for temporary char array
    CurrentNumber = CallocateMemoryChar();

    // Opening file - which containing the matrix - in read-only mode
    FILE* Input_Mx = fopen(FILE_NAME_MATRIX, "r+");

    // Read in matrix from input file
    // The ReadInTensor function returns a double* array
    Matrix = ReadInTensor(Input_Mx, Matrix, CurrentNumber, &NumberOfRowsMx, &TensorIndexMx);

    // Free memory space of current number
    if(CurrentNumber != NULL)
    {
        printf("\nFreeing memory space of CurrentNumber array...\n");
        free(CurrentNumber);
        printf("CurrentNumber array's memory space got freed!\n\n");
    }

    printf("Input Matrix:\n");
    for(unsigned int i = 0; i < 4; i++)
    {
        for(unsigned int j = 0; j < 4; j++)
        {
            printf("%lf ", Matrix[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");


    // ================== Vector
    // Reading in vector from a file to the Vector double* array

    // Opening file - which containing the vector - in read-only mode
    FILE* Input_Vc = fopen(FILE_NAME_VECTOR, "r+");

    // Allocate memory for temporary char array
    CurrentNumber = CallocateMemoryChar();

    // Read in vectorfrom input file
    // The ReadInOneDArray function returns a double* array
    Vector = ReadInTensor(Input_Vc, Vector, CurrentNumber, &NumberOfRowsVcAUX, &TensorIndexVc);

    // Free memory space of current number
    if(CurrentNumber != NULL)
    {
        printf("\nFreeing memory space of CurrentNumber array...\n");
        free(CurrentNumber);
        printf("CurrentNumber array's memory space got freed!\n\n");
    }

    printf("Input vector:\n");
    for(unsigned int k = 0; k < 4; k++)
    {
        printf("%lf ", Vector[k]);
    }
    printf("\n");

    // ================== Remaining calculations
    // Multiply matrix and vector, then print output to stdout

    // Calculate constants for this particular problem
    // Column length of the matrix
    unsigned int NumberOfColumnsMx = TensorIndexMx / NumberOfRowsMx;
    // Size of the vector (renaming variable)
    unsigned int NumberOfRowsVc = TensorIndexVc;


    // Multiply the matrix and the vector if it's possible
	Output = MultiplyThings(NumberOfRowsMx, NumberOfColumnsMx, NumberOfRowsVc, Matrix, Vector);

    free(Matrix);
    free(Vector);

    // Print output
    PrintOutput(Output, NumberOfRowsMx, NumberOfRowsVc);

    free(Output);

    return 0;
}
