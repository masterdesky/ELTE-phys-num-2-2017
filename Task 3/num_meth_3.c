#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Machine epsilon of float32
const double floateps = 1.1920929e-07;



//
// Function prototypes
//
void CheckingIntegrity();
//
void ReadInFileNameSTDIN();
void ReadInIndependentSTDIN();
void ReadInPolyOrderSTDIN();
void CountRowsOfInput();
//
double* CallocateMemory();
//
void LoadFile();
void LoadMeasurementsErrors();
void LoadDependentVariables();
void LoadInitialMatrix();
// Polynomial
void TransposeInitialMatrix();
void MultiplyTransposeAndOriginal();
void MultiplyTransposeAndVector();
void MultiplyGaussedWith_MtxTMulVec();
// Gauss-Jordan
void GaussJordan();
void WriteMtxTMulMtxIntoAnotherMatrix();
void AppendIdentityMatrix();
void SingularMatrixChecker();
void EliminationChecker();
void GaussJordanElimination();
// Output
void PrintFittedParameters();
// Quasi-Main
void FunctionWaveOfGaussJordan();


// Checking file and inputs integrity
void CheckingIntegrity(FILE* InputFile)
{
    if(InputFile == NULL)
    {
        printf("ERROR!\nInput file is missing or corrupted. The program exits.\n\n");
        exit(EXIT_FAILURE);
    }

    else
    {
        printf("Input file is loaded!\n\n");
    }
}


// Scan file name from standard input
void ReadInFileNameSTDIN(char* FILE_NAME)
{
    printf("Please enter the name of the file, which includes the fitting data.   \n"
           "The file name should be something like that: small.dat, or small.txt\n\n"
           "File name: ");
    scanf("%255s", *FILE_NAME);
    printf("\n");
}


void ReadInIndependentSTDIN(int* IndVariables)
{

    int Scanned;

    printf("Please enter the number of the independent variables: ");
    scanf(" %d", &Scanned);
    printf("\n");

    *IndVariables = Scanned;

}


void ReadInPolyOrderSTDIN(int* PolOrder)
{
    int Scanned;

    printf("Order of the polynomial: ");
    scanf(" %d", &Scanned);
    printf("\n");

    *PolOrder = Scanned;
}


// Count the number of rows in the file
void CountRowsOfInput(FILE* InputFile, int* Counter)
{
    char CharacterNow;
    char CharacterPrev = '\0';

    while(1)
    {
        CharacterNow = getc(InputFile);

        if(CharacterNow == EOF)
        {
            Counter++;
        }
        
        if(CharacterNow == '\n' && CharacterPrev != '\n')
        {
            Counter++;
        }
        
        CharacterPrev = CharacterNow;
    }

    //Testing
    printf("Number of lines: %d\n\n", Counter);

    return Counter;
}


// Allocate memory for reading in the full file
double* CallocateMemory(int MemorySize)
{
    // Callocate memory
    double* InputDoubleArray = (double*)calloc(MemorySize, sizeof(double));

    // Checking if callocating was successfull
    if(InputDoubleArray == NULL)
    {
        perror("Allocation error. Program exits");
        EXIT_FAILURE;
    }

    return InputDoubleArray;
}


// Load the full file into a double* array
void LoadFile(FILE* InputFile, double* FullMatrix, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables)
{
    // Load file into matrix with N*(M+2) elements, where
    // N = Number of rows
    // M = Number of independent variables
    // +2 stands for the dependent variables and for the errors columns
    for (unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < NumberOfIndependentVariables+2; j++)
        {
            fscanf(InputFile, "%lf", &FullMatrix[(NumberOfIndependentVariables+2) * i + j]);
        }
    }
}


// Load the errors into a double* array from FullMatrix
void LoadMeasurementsErrors(double* FullMatrix, double* MeasurementsErrors, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables)
{
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = NumberOfIndependentVariables+1; j < NumberOfIndependentVariables+2; j++)
        {
            MeasurementsErrors[i] = FullMatrix[(NumberOfIndependentVariables+2) * i + j];
        }
    }
}


// First load the dependent variables into an array
// Secondly divide them by the corresponding measurement errors
void LoadDependentVariables(double* FullMatrix, double* DependentVariables, double* MeasurementsErrors, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables)
{
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = NumberOfIndependentVariables; j < NumberOfIndependentVariables+1; j++)
        {
            DependentVariables[i] = FullMatrix[(NumberOfIndependentVariables+2) * i + j];
            DependentVariables[i] /= MeasurementsErrors[i];
        }
    }
}


// Load the independent variables into a matrix according to the order of the polynomial fit
// First fill the first column with ones, hence this is the very first (0th) order of the polynomial
// Then fill the remaining spaces with correct powers of the dependent variables
// The matrix should look like this:
//
// | 1 x_1 (x_1)^2 . . (x_1)^n |
// | 1 x_2 (x_2)^2 . . (x_2)^n |
// | 1 x_3 (x_3)^2 . . (x_3)^n |
// | 1 x_4 (x_4)^2 . . (x_4)^n |
// | .  .    .           .     |
// | .  .    .           .     |
//
// Where "n" is the order of the desired polynomial fit
// At last, divide the elements of this matrix with the corresponding measurement errors
void LoadInitialMatrix(double* FullMatrix, double* InitialMatrix, double* DependentVariables, double* MeasurementsErrors, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables, unsigned int DimensionOf_MtxTMulMtx)
{
    // Fill first row with ones
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < 1; j++)
        {
            InitialMatrix[DimensionOf_MtxTMulMtx * i + j] = 1;
        }
    }

    // Load in the datas from the file
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 1; j < NumberOfIndependentVariables; j++)
        {
            InitialMatrix[DimensionOf_MtxTMulMtx * i + j] = pow(DependentVariables[i], j);
            InitialMatrix[DimensionOf_MtxTMulMtx * i + j] /= MeasurementsErrors[i];
        }
    }
}


// Transposing the InitialMatrix
void TransposeInitialMatrix(double* TransposedInitialMatrix, double* InitialMatrix, unsigned int RowsOfInitialMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    //Switch all element M_ij, with the element M_ji
    for(unsigned int j = 0; j < RowsOfInitialMatrix; j++)
    {
        for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
        {
            TransposedInitialMatrix[RowsOfInitialMatrix * i + j] = InitialMatrix[DimensionOf_MtxTMulMtx * j + i];
        }
    }
}


// Multiply transposed matrix and the vector, created from the dependent variables
// Like this: X.T * y
void MultiplyTransposeAndVector(double* MtxTMulVec, double* DependentVariables, double* TransposedInitialMatrix, unsigned int RowsOfInitialMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    //Multiply the vector and the transposed matrix
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = 0; j < RowsOfInitialMatrix; j++)
        {
               MtxTMulVec[i] += (TransposedInitialMatrix[RowsOfInitialMatrix * i + j] * DependentVariables[j]);
        }
    }
}


// Multiply original InitialMatrix and its trasposed variant from the left side
// Like this: X.T * X
void MultiplyTransposeAndOriginal(double* GaussedMatrix, double* MtxTMulMtx, double* InitialMatrix, double* TransposedInitialMatrix, unsigned int RowsOfInitialMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    // Temporary storage for matrix elements
    double Temp = 0;

    // Multiply the two matrices
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxTMulMtx; j++)
        {
            for(unsigned int k = 0; k < RowsOfInitialMatrix; k++)
            {
                Temp += TransposedInitialMatrix[RowsOfInitialMatrix * i + k] * InitialMatrix[DimensionOf_MtxTMulMtx * k + j];
            }

            GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + j] = Temp;
            MtxTMulMtx[DimensionOf_MtxTMulMtx * i + j] = Temp;
            Temp = 0;
        }
    }
}


// Append and identity matrix in the GaussedMatrix's other half
void AppendIdentityMatrix(double* GaussedMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = DimensionOf_MtxTMulMtx; j < 2 * DimensionOf_MtxTMulMtx; j++)
        {
            if(i == (j % DimensionOf_MtxTMulMtx))
            {
                GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + j] = 1;
            }
            else
                GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + j] = 0;
        }
    }
}


// Gauss-Jordan elimination to find inverse of the matrix
// 1. The function first searches for the greatest element in the first column: if it's found in one of the rows, then that row is swapped with the first one.
// 2. In the next steps the function do the same for the other columns too, starting the loop at the second row, then the third, then so on.
// 3. After a step of the loop, the function are normalizing the whole row and performing substitutions to form the required inverse matrix.
// 4. At the end, the diagonal elements will be all just ones, and below the diagonal line will be full of zeros.
void GaussJordan(double* GaussedMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    // Temporary storage for indeces
    // Indicating the index of the row with current greatest element
    unsigned int TempStoreIndex;

    // Temporary storage for a double element of the matrix
    double TempStorage;

    // Temporary index for matrix elements at normalization
    double mat_ij;
    
    // Run column-wise through elements
    // The index (j) here indicates columns
    for(unsigned int j = 0; j < DimensionOf_MtxTMulMtx; j++)
    {
        // Set the jth element of the jth column is initially the greatest element
        TempStoreIndex = j;

        // Search for the greatest element in the jth column, starting from the (j+1)th row
        // It runs through all elements
        for(unsigned int i = (j + 1); j < DimensionOf_MtxTMulMtx; j++)
        {
            if(GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + j] > GaussedMatrix[2 * DimensionOf_MtxTMulMtx * TempStoreIndex + j])
            {
                TempStoreIndex = i;
            }
        }

        // Check if the greatest element is smaller, than the 32bit float machine epsilon
        // If yes, then the matrix is numerically singular, and cannot be inverted
        if(fabs(GaussedMatrix[2 * DimensionOf_MtxTMulMtx * TempStoreIndex + j]) < floateps)
        {
            perror("Given matrix is singular and cannot be inverted! The program exits");
            EXIT_FAILURE;
        }

        // Swapping rows, which has the greatest element in the column
        if(TempStoreIndex != j)
        {
            for(unsigned int k = 0; k < 2 * DimensionOf_MtxTMulMtx; k++)
            {
                TempStorage = GaussedMatrix[2 * DimensionOf_MtxTMulMtx * j + k];
                GaussedMatrix[2 * DimensionOf_MtxTMulMtx * j + k] = GaussedMatrix[2 * DimensionOf_MtxTMulMtx * TempStoreIndex + k];
                GaussedMatrix[2 * DimensionOf_MtxTMulMtx * TempStoreIndex + k] = TempStorage;
            }
        }

        // Perform substitutions and normalization
        for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
        {
            mat_ij = GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + j];

            if(i != j)
            {
                for(unsigned int k = 0; k < 2 * DimensionOf_MtxTMulMtx; k++)
                {
                    GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + k] -= (GaussedMatrix[2 * DimensionOf_MtxTMulMtx * j + k]/GaussedMatrix[2 * DimensionOf_MtxTMulMtx * j + j]) * mat_ij;
                }
            }

            else
            {
                for(unsigned int k = 0; k < 2 * DimensionOf_MtxTMulMtx; k++)
                {
                    GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + k] /= mat_ij ;
                }
            }
        }
    }
}


void EliminationChecker(double* GaussedMatrix, double* MtxTMulMtx, double* IdentityMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    // Temporary storage
    double Temp  = 0;

    // Multiply the two matrices
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxTMulMtx; j++)
        {
            for(unsigned int k = 0; k < DimensionOf_MtxTMulMtx; k++)
            {
                Temp += GaussedMatrix[2 * DimensionOf_MtxTMulMtx * i + k] * MtxTMulMtx[2 * DimensionOf_MtxTMulMtx * k + j];
            }

            if(Temp < floateps)
            {
                IdentityMatrix[2 * DimensionOf_MtxTMulMtx * i + j] = 0;
            }

            else
            {
                IdentityMatrix[2 * DimensionOf_MtxTMulMtx * i + j] = Temp;
            }
            Temp = 0;
        }
    }

    // Testing
    fprintf(stdout, "Should be identity:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxTMulMtx; j++)
        {
            fprintf(stdout, "%g ", IdentityMatrix[2*DimensionOf_MtxTMulMtx * i + j]);
        }
        fprintf(stdout, "\n");
    }
}

void MultiplyGaussedWith_MtxTMulVec(double* FittedParameters, double* MtxTMulVec, double* GaussedMatrix, unsigned int DimensionOf_MtxTMulMtx)
{
    //Multiply them
    for(unsigned int i = 0; i < DimensionOf_MtxTMulMtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxTMulMtx; j++)
        {
               FittedParameters[i] += (GaussedMatrix[2*DimensionOf_MtxTMulMtx * i + j] * MtxTMulVec[j]);
        }
    }
}


/* ------------------------------------------- OUTPUT MANAGEMENT ------------------------------------------- */

void PrintFittedParameters(FILE* OutputFile, double* FittedParameters, unsigned int OrderOfPolynomial, unsigned int NumberOfIndependentVariables)
{
    //Prinf fitting parameters
    fprintf(OutputFile, "Fitted parameters:\n");
    fprintf(OutputFile, "[cnst]: %g\n", FittedParameters[0]);
    for(unsigned int i = 0; i < OrderOfPolynomial; i++)
    {
        for(unsigned int j = 0; j < NumberOfIndependentVariables; j++)
        {
            fprintf(OutputFile, "[x_%d]^%d: %g\n", j+1, i+1, FittedParameters[j + i*NumberOfIndependentVariables + 1]);
        }
    }
}

/* ------------------------------------------- MAIN ------------------------------------------- */

//
//  #   #    #    #  #   #
//  # # #   # #   #  # # #
//  #   #  #   #  #  #   #
//
int main()
{
    // I/O files
    FILE* InputFile;
    FILE* OutputFile;
    
    // Used pointer-type array variables
    double* FullMatrix;                     //Matrix
    double* InitialMatrix;                  //Matrix
    double* DependentVariables;             //Vector
    double* MeasurementsErrors;             //Vector
    double* TransposedInitialMatrix;        //Matrix
    double* MtxTMulMtx;                     //Matrix
    double* MtxTMulVec;                     //Vector
    double* FittedParameters;               //Vector
    double* GaussedMatrix;                  //Matrix
    double* GaussedHalfMatrix;              //Matrix
    double* IdentityMatrix;                 //Matrix

    // Used variables, regarding size of input file and matrix
    // Eg. columns, rows, dependent/independet variables, etc.
    char* FILE_NAME = (char*)calloc(256, sizeof(char));
    unsigned int NumberOfIndependentVariables;
    unsigned int OrderOfPolynomial;
    unsigned int RowsOfInitialMatrix = 0;
    unsigned int DimensionOf_MtxTMulMtx;



    // Input file name
    ReadInFileNameSTDIN(&FILE_NAME);

    // Open files
    InputFile = fopen(FILE_NAME, "r");
    OutputFile = fopen("FittedParameters.txt", "w+");

    // Checking file integrity
    CheckingIntegrity(InputFile);

    // Declare global variables
    ReadInIndependentSTDIN(&NumberOfIndependentVariables);
    ReadInPolyOrderSTDIN(&OrderOfPolynomial);
    CountRowsOfInput(InputFile, &RowsOfInitialMatrix);
    DimensionOf_MtxTMulMtx = OrderOfPolynomial + 1;

    //Test
    printf("\n");
    printf("Number of rows: %d\n", RowsOfInitialMatrix);
    printf("Number of independent variables: %d\n", NumberOfIndependentVariables);
    printf("Dimension of NxN matrix (X_transposed * X): %d\n", DimensionOf_MtxTMulMtx);



    // Allocating memory for double* arrays
    // CallocateMemory function allocates arbitrary sized memory, passed by value to it
    FullMatrix = CallocateMemory(RowsOfInitialMatrix * (NumberOfIndependentVariables + 2));
    InitialMatrix = CallocateMemory(RowsOfInitialMatrix * (OrderOfPolynomial + 1));

    DependentVariables = CallocateMemory(RowsOfInitialMatrix);
    MeasurementsErrors = CallocateMemory(RowsOfInitialMatrix);
    MtxTMulMtx = CallocateMemory(DimensionOf_MtxTMulMtx * DimensionOf_MtxTMulMtx);
    MtxTMulVec = CallocateMemory(DimensionOf_MtxTMulMtx);
    FittedParameters = CallocateMemory(DimensionOf_MtxTMulMtx);

    TransposedInitialMatrix = CallocateMemory(RowsOfInitialMatrix * DimensionOf_MtxTMulMtx);
    IdentityMatrix = CallocateMemory(DimensionOf_MtxTMulMtx * DimensionOf_MtxTMulMtx);

    GaussedMatrix = CallocateMemory(2 * DimensionOf_MtxTMulMtx * DimensionOf_MtxTMulMtx);
    GaussedHalfMatrix = CallocateMemory(DimensionOf_MtxTMulMtx * DimensionOf_MtxTMulMtx);



    // Loading in the designmatrix and the corresponding vectors (Dependent variables and errors)
    LoadFile(InputFile, &FullMatrix, RowsOfInitialMatrix, NumberOfIndependentVariables);
    LoadMeasurementsErrors(FullMatrix, &MeasurementsErrors, RowsOfInitialMatrix, NumberOfIndependentVariables);
    LoadDependentVariables(FullMatrix, &DependentVariables, MeasurementsErrors, RowsOfInitialMatrix, NumberOfIndependentVariables);
    LoadInitialMatrix(FullMatrix, &InitialMatrix, DependentVariables, MeasurementsErrors, RowsOfInitialMatrix, NumberOfIndependentVariables, DimensionOf_MtxTMulMtx);

    free(FullMatrix);

    // Closing input file
    fclose(InputFile);

    // Create the matrix equation
    // X.T * X * a = X.T * y
    TransposeInitialMatrix(&TransposedInitialMatrix, InitialMatrix, RowsOfInitialMatrix, DimensionOf_MtxTMulMtx);
    MultiplyTransposeAndOriginal(&GaussedMatrix, &MtxTMulMtx, InitialMatrix, TransposedInitialMatrix, RowsOfInitialMatrix, DimensionOf_MtxTMulMtx);
    MultiplyTransposeAndVector(&MtxTMulVec, DependentVariables, TransposedInitialMatrix, RowsOfInitialMatrix, DimensionOf_MtxTMulMtx);

    free(TransposedInitialMatrix);

    // Invert MtxTMulMtx with Gauss-Jordan elimination
    AppendIdentityMatrix(&GaussedMatrix, DimensionOf_MtxTMulMtx);
    GaussJordan(&GaussedMatrix, DimensionOf_MtxTMulMtx);
    
    // Calculating the fitting parameters
    MultiplyGaussedWith_MtxTMulVec(FittedParameters, MtxTMulVec, GaussedMatrix, DimensionOf_MtxTMulMtx);

    // Check if Gauss-Jordan elimination was successfull
    EliminationChecker(GaussedMatrix, MtxTMulMtx, IdentityMatrix, DimensionOf_MtxTMulMtx);

    // Output
    PrintFittedParameters(OutputFile, FittedParameters, OrderOfPolynomial, NumberOfIndependentVariables);

    free(InitialMatrix);
    free(DependentVariables);
    free(MeasurementsErrors);
    free(MtxTMulVec);
    free(MtxTMulVec);
    free(FittedParameters);
    free(GaussedMatrix);
    free(GaussedHalfMatrix);
    free(IdentityMatrix);

    return 0;
}
