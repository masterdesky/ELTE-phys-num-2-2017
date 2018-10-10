#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


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
void LoadInitialMatrix();
void LoadDependentVariables();
void LoadMeasurementsErrors();
//
void PrintToStdout();
// Polynomial
void FinishInitialMatrix();
void FinishColumnVector();
void TransposeInitialMatrix();
void MultiplyTransposeAndOriginal();
void MultiplyTransposeAndVector();
void MultiplyGaussedWith_MtxT_Mul_Vec();
// Gauss-Jordan
void WriteMtxTMulMtxIntoAnotherMatrix();
void CreateIdentityMatrix();
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


//Count the number of rows of the matrix
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


// Load the dependent variables into an array
void LoadDependentVariables(double* FullMatrix, double* DependentVariables, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables)
{
    //Temporary storage
    double Temp;

    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = NumberOfIndependentVariables; j < NumberOfIndependentVariables+1; j++)
        {
            DependentVariables[i] = FullMatrix[(NumberOfIndependentVariables+2) * i + j];
        }
    }
}


// Load the independent variables from file into a matrix
void LoadInitialMatrix(double* FullMatrix, double* InitialMatrix, double* DependentVariables, unsigned int RowsOfInitialMatrix, unsigned int NumberOfIndependentVariables, unsigned int DimensionOf_MtxT_Mul_Mtx)
{
    // Fill first row with ones
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < 1; j++)
        {
            InitialMatrix[2 * DimensionOf_MtxT_Mul_Mtx * i + j] = 1;
        }
    }

    // Load in the datas from the file
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 1; j < NumberOfIndependentVariables; j++)
        {
            InitialMatrix[2 * DimensionOf_MtxT_Mul_Mtx * i + j] = pow(DependentVariables[i], j);
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


void PrintToStdout(double* FullMatrix, double* InitialMatrix, double* DependentVariables, double* MeasurementsErrors, double* PolExponents, FILE* OutputFile)
{
    //PRINTING
    fprintf(OutputFile,"InputFile:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < NumberOfIndependentVariables+2; j++)
        {
            fprintf(OutputFile,"%g ",FullMatrix[(NumberOfIndependentVariables+2) * i + j]);
        }
        fprintf(OutputFile,"\n");
    }

    fprintf(OutputFile,"InitialMatrix:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(OutputFile,"%g ",InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(OutputFile,"\n");
    }

    fprintf(OutputFile, "DependentVariables:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        fprintf(OutputFile,"%g\n",DependentVariables[i]);
    }

    fprintf(OutputFile, "MeasurementsErrors:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        fprintf(OutputFile,"%g\n",MeasurementsErrors[i]);
    }

    fprintf(OutputFile, "Exponents:\n");
    for(unsigned int i = 0; i < OrderOfPolynomial; i++)
    {
        fprintf(OutputFile,"%g\n",PolExponents[i]);
    }
}


/* ------------------------------------------- POLYNOMIAL FITTING ------------------------------------------- */

//Build up the final form of the design matrix
void FinishInitialMatrix(double* InitialMatrix, double* MeasurementsErrors, double* PolExponents, FILE* OutputFile)
{
    //Exponentiate matrix elements with the correct exponents
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < OrderOfPolynomial; j++)
        {
            for(unsigned int k = 0; k < NumberOfIndependentVariables; k++)
            {
                InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + k + 1 + j*NumberOfIndependentVariables] = pow(InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + k + 1], PolExponents[j]);
            }
        }
    }

    //Divide the matrix elements by the approximation errors
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] /= MeasurementsErrors[i];
        }
    }
/*
    //TEST
    fprintf(OutputFile, "Design Matrix:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(OutputFile, "%g ",InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(OutputFile, "\n");
    }
    fprintf(OutputFile, "\n\n\n\n");*/
}

void FinishColumnVector(double* DependentVariables, double* MeasurementsErrors)
{
    //Divide the vector elements by the approximation errors
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        DependentVariables[i] /= MeasurementsErrors[i];
    }
}

void TransposeInitialMatrix(double* TransposedMtx, double* InitialMatrix)
{
    //Temporary storage
    double Temp;

    //Switch all element M_ij, with the element M_ji
    for(unsigned int j = 0; j < RowsOfInitialMatrix; j++)
    {
        for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
        {
            Temp = InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + i];
            TransposedMtx[RowsOfInitialMatrix * i + j] = Temp;
        }
    }
}

void MultiplyTransposeAndOriginal(double* MtxT_Mul_Mtx, double* InitialMatrix, double* TransposedMtx, FILE* OutputFile)
{
    //Temporary storage
    double Temp;
/*
    //TEST
    fprintf(OutputFile, "Original Matrix:\n");
    for(unsigned int i = 0; i < RowsOfInitialMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(OutputFile, "%g ", InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(OutputFile, "\n");
    }
    fprintf(OutputFile, "\n\n\n\n");

    //TEST
    fprintf(OutputFile, "Transposed Matrix:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < RowsOfInitialMatrix; j++)
        {
            fprintf(OutputFile, "%g ", TransposedMtx[RowsOfInitialMatrix * i + j]);
        }
        fprintf(OutputFile, "\n");
    }
    fprintf(OutputFile, "\n\n\n\n");
*/

    //Multiply the two matrices
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            for(unsigned int k = 0; k < RowsOfInitialMatrix; k++)
            {
                Temp += TransposedMtx[RowsOfInitialMatrix * i + k] * InitialMatrix[2*DimensionOf_MtxT_Mul_Mtx * k + j];
            }
            MtxT_Mul_Mtx[2*DimensionOf_MtxT_Mul_Mtx * i + j] = Temp;
            Temp = 0;
        }
    }
/*
    //TEST
    fprintf(OutputFile, "MtxT_Mul_Mtx:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(OutputFile, "%g ", MtxT_Mul_Mtx[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(OutputFile, "\n");
    }
    fprintf(OutputFile, "\n\n\n\n");*/
}

void MultiplyTransposeAndVector(double* MtxT_Mul_Vec, double* DependentVariables, double* TransposedMtx, FILE* OutputFile)
{
    //Multiply the vector and the transposed matrix
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < RowsOfInitialMatrix; j++)
        {
               MtxT_Mul_Vec[i] += (TransposedMtx[RowsOfInitialMatrix * i + j] * DependentVariables[j]);
        }
    }
/*
    //TEST
    fprintf(OutputFile, "MtxT_Mul_Vec:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        fprintf(OutputFile, "%g\n", MtxT_Mul_Vec[i]);
    }
    fprintf(OutputFile, "\n\n\n\n");*/
}


void MultiplyGaussedWith_MtxT_Mul_Vec(double* FittedParameters, double* MtxT_Mul_Vec, double* GaussedHalfMatrix)
{
    //Temporary storage
    double Temp;

    //Multiply them
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
               FittedParameters[i] += (GaussedHalfMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] * MtxT_Mul_Vec[j]);
        }
    }
/*
    //TEST
    fprintf(stdout, "Fitted parameters:\n");
    fprintf(stdout, "[cnst]: %g\n", FittedParameters[0]);
    for(unsigned int i = 0; i < OrderOfPolynomial; i++)
    {
        for(unsigned int j = 0; j < NumberOfIndependentVariables; j++)
        {
            fprintf(stdout, "[x_%d]^%d: %g\n", j+1, i+1, FittedParameters[j + i*NumberOfIndependentVariables + 1]);
        }
    }*/
}


/* ------------------------------------------- GAUSS-JORDAN ELIMINATION ------------------------------------------- */

void WriteMtxTMulMtxIntoAnotherMatrix(double* MtxT_Mul_Mtx, double* GaussedMatrix)
{
    //Write into GaussedMatrix
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] = MtxT_Mul_Mtx[2*DimensionOf_MtxT_Mul_Mtx * i + j];
        }
    }
/*
    fprintf(stdout, "Gaussed Input:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }*/
}


//Create and identity matrix in the input matrix's other half
void CreateIdentityMatrix(double* InitialMatrixWithIdentity)
{
/*
    //TEST
    fprintf(stdout, "Matrix before identity\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", InitialMatrixWithIdentity[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }
*/
    //Summon matrix on .5
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = DimensionOf_MtxT_Mul_Mtx; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            InitialMatrixWithIdentity[2*DimensionOf_MtxT_Mul_Mtx * i + j] = 0;
        }
    }
/*
    //TEST
    fprintf(stdout, "Matrix before identity after summon:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout,"%g ", InitialMatrixWithIdentity[(2*DimensionOf_MtxT_Mul_Mtx) * i + j]);
        }
        fprintf(stdout, "\n");
    }
*/
    //Create identity matrix on .5
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = DimensionOf_MtxT_Mul_Mtx; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            if(i == j%DimensionOf_MtxT_Mul_Mtx)
            {
                InitialMatrixWithIdentity[2*DimensionOf_MtxT_Mul_Mtx * i + j] = 1;
            }
            else
                InitialMatrixWithIdentity[2*DimensionOf_MtxT_Mul_Mtx * i + j] = 0;
        }
    }
/*
    fprintf(stdout, "Matrix after identity:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout,"%g ", InitialMatrixWithIdentity[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }*/
}

//Checks if matrix is singular
void SingularMatrixChecker(double* GaussedMatrix, int TempStoreFirst, int j)
{
        if(fabs(GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * TempStoreFirst + j]) < 0.000001)
        {
            printf("Given matrix is singular and cannot be inverted! The program exits.");
            exit(EXIT_FAILURE);
        }
}

void EliminationChecker(double* MtxT_Mul_Mtx, double* GaussedHalfMatrix, double* IdentityMatrix)
{
    //Temporary storage
    double Temp;
/*
    //TEST1
    fprintf(stdout, "MtxT_Mul_Mtx in checking:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", MtxT_Mul_Mtx[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }

    fprintf(stdout, "Gaussed in checking:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", GaussedHalfMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }
*/
    //Multiply the two matrices
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            for(unsigned int k = 0; k < DimensionOf_MtxT_Mul_Mtx; k++)
            {
                Temp += GaussedHalfMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + k] * MtxT_Mul_Mtx[2*DimensionOf_MtxT_Mul_Mtx * k + j];
            }
            if(Temp < 0.000001)
            {
                IdentityMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] = 0;
            }
            else
            {
                IdentityMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] = Temp;
            }
            Temp = 0;
        }
    }

    //TEST2
    fprintf(stdout, "Should be identity:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", IdentityMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }
}

//Execute Gauss-Jordan elimination
void GaussJordanElimination(double* GaussedMatrix, double* GaussedHalfMatrix, FILE* Output)
{
    //Temporary storage for indexes
    int TempNorm;
    int TempStoreFirst;
    double TempStoreSecond;

    //Index for matrix elements
    double max;
    double mat_ji;
/*
    //TEST5
    fprintf(Output, "Matrix with identity:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(Output,"%g ", GaussedMatrix[(2*DimensionOf_MtxT_Mul_Mtx) * i + j]);
        }
        fprintf(Output, "\n");
    }
    fprintf(Output, "\n\n\n\n");
*/
    //Normalize matrix
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        TempNorm = 0;
        for(unsigned int j = 1; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            if(GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] > GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + TempNorm])
            {
                TempNorm = j;
            }
        }

        max = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + TempNorm];

        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] /= max;
        }
    }
/*
    fprintf(Output, "Normalized matrix:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(Output, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(Output, "\n");
    }
    fprintf(Output, "\n\n\n\n");
*/

    //Change the row of the greatest column element with the first row
    for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
    {
        TempStoreFirst = j;

        //Search for the greatest column element in rows
        for(unsigned int i = j + 1; i < DimensionOf_MtxT_Mul_Mtx; i++)
        {
            if(GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j] > GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * TempStoreFirst + j])
            {
                TempStoreFirst = i;
            }
        }

        //Checking for singularity
        SingularMatrixChecker(GaussedMatrix, TempStoreFirst, j);

        //Swapping row which has the greatest column element
        if(TempStoreFirst != j)
        {
            for(unsigned int k = 0; k < 2*DimensionOf_MtxT_Mul_Mtx; k++)
            {
                TempStoreSecond = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + k];
                GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + k] = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * TempStoreFirst + k];
                GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * TempStoreFirst + k] = TempStoreSecond;
            }
        }
    }
/*
    fprintf(Output, "Sorted matrix:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(Output, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(Output, "\n");
    }
    fprintf(Output, "\n\n\n\n");
*/
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        //Perform some magic
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            //Normalize elements with greatest in every row
            if(j == i)
            {
                mat_ji = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + i];
                for(unsigned int k = 0; k < 2*DimensionOf_MtxT_Mul_Mtx; k++)
                {
                    GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + k] /= mat_ji;
                }
            }

            //From every row we substitute the i'th row multiplied with the element M_ji. So in the place M_ji there will be zeros.
            else
            {
                mat_ji = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + i];
                for(unsigned int k = 0; k < 2*DimensionOf_MtxT_Mul_Mtx; k++)
                {
                    GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * j + k] -= (GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + k] / GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + i]) * mat_ji;
                }
            }
        }
    }
/*
    //TESTS
    fprintf(stdout, "Gaussed matrix 0:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }

    fprintf(Output, "Gaussed matrix 1:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(Output, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(Output, "\n");
    }
    fprintf(Output, "\n\n\n\n");

    fprintf(stdout, "Gaussed matrix 2:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = DimensionOf_MtxT_Mul_Mtx; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }
*/
    //Write into GaussedMatrix
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = DimensionOf_MtxT_Mul_Mtx; j < 2*DimensionOf_MtxT_Mul_Mtx; j++)
        {
            GaussedHalfMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j - DimensionOf_MtxT_Mul_Mtx] = GaussedMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j];
        }
    }
/*
    fprintf(stdout, "Gaussed matrix PASSED:\n");
    for(unsigned int i = 0; i < DimensionOf_MtxT_Mul_Mtx; i++)
    {
        for(unsigned int j = 0; j < DimensionOf_MtxT_Mul_Mtx; j++)
        {
            fprintf(stdout, "%g ", GaussedHalfMatrix[2*DimensionOf_MtxT_Mul_Mtx * i + j]);
        }
        fprintf(stdout, "\n");
    }*/
}


/* ------------------------------------------- OUTPUT MANAGEMENT ------------------------------------------- */

void PrintFittedParameters(FILE* OutputFile, double* FittedParameters)
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

//Calling of GJ functions
void FunctionWaveOfGaussJordan(double* MtxT_Mul_Mtx, double* GaussedMatrix, double* IdentityMatrix, double* GaussedHalfMatrix, FILE* OutputFile)
{
    WriteMtxTMulMtxIntoAnotherMatrix(MtxT_Mul_Mtx, GaussedMatrix);
    CreateIdentityMatrix(GaussedMatrix);
    GaussJordanElimination(GaussedMatrix, GaussedHalfMatrix, OutputFile);
    EliminationChecker(MtxT_Mul_Mtx, GaussedHalfMatrix, IdentityMatrix);
}


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
    double* FullMatrix;         //Matrix
    double* InitialMatrix;      //Matrix
    double* DependentVariables; //Vector
    double* MeasurementsErrors; //Vector
    double* TransposedMtx;      //Matrix
    double* MtxT_Mul_Mtx;       //Matrix
    double* MtxT_Mul_Vec;       //Vector
    double* FittedParameters;   //Vector
    double* GaussedMatrix;      //Matrix
    double* GaussedHalfMatrix;  //Matrix
    double* IdentityMatrix;     //Matrix

    // Used variables, regarding size of input file and matrix
    // Eg. columns, rows, dependent/independet variables, etc.
    char* FILE_NAME = (char*)calloc(256, sizeof(char));
    unsigned int NumberOfIndependentVariables;
    unsigned int OrderOfPolynomial;
    unsigned int RowsOfInitialMatrix = 0;
    unsigned int DimensionOf_MtxT_Mul_Mtx;



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
    DimensionOf_MtxT_Mul_Mtx = OrderOfPolynomial + 1;

    //Test
    printf("\n");
    printf("Number of rows: %d\n", RowsOfInitialMatrix);
    printf("Number of independent variables: %d\n", NumberOfIndependentVariables);
    printf("Dimension of NxN matrix (X_transposed * X): %d\n", DimensionOf_MtxT_Mul_Mtx);

    // Allocating memory for double* arrays
    // CallocateMemory function allocates arbitrary sized memory, passed by value to it
    FullMatrix = CallocateMemory(RowsOfInitialMatrix * (NumberOfIndependentVariables + 2));
    InitialMatrix = CallocateMemory(RowsOfInitialMatrix * (OrderOfPolynomial + 1));

    DependentVariables = CallocateMemory(RowsOfInitialMatrix);
    MeasurementsErrors = CallocateMemory(RowsOfInitialMatrix);
    MtxT_Mul_Vec = CallocateMemory(DimensionOf_MtxT_Mul_Mtx);
    FittedParameters = CallocateMemory(DimensionOf_MtxT_Mul_Mtx);

    TransposedMtx = CallocateMemory(2 * RowsOfInitialMatrix * DimensionOf_MtxT_Mul_Mtx);
    MtxT_Mul_Mtx = CallocateMemory(2 * DimensionOf_MtxT_Mul_Mtx * DimensionOf_MtxT_Mul_Mtx);
    GaussedMatrix = CallocateMemory(2 * DimensionOf_MtxT_Mul_Mtx * DimensionOf_MtxT_Mul_Mtx);
    GaussedHalfMatrix = CallocateMemory(2 * DimensionOf_MtxT_Mul_Mtx * DimensionOf_MtxT_Mul_Mtx);
    IdentityMatrix = CallocateMemory(2 * DimensionOf_MtxT_Mul_Mtx * DimensionOf_MtxT_Mul_Mtx);



    //Loading in the matrix and vectors
    LoadFile(InputFile, &FullMatrix, RowsOfInitialMatrix, NumberOfIndependentVariables);
    LoadInitialMatrix(FullMatrix, &InitialMatrix, &DependentVariables, RowsOfInitialMatrix, NumberOfIndependentVariables, DimensionOf_MtxT_Mul_Mtx);
    LoadDependentVariables(FullMatrix, &DependentVariables, RowsOfInitialMatrix, NumberOfIndependentVariables);
    LoadMeasurementsErrors(FullMatrix, &MeasurementsErrors, RowsOfInitialMatrix, NumberOfIndependentVariables);

    //PrintToStdout(FullMatrix, InitialMatrix, DependentVariables, MeasurementsErrors, PolExponents, OutputFile); //TEST

    free(FullMatrix);

    //Closing file
    fclose(InputFile);

    //Performing mathematics and stuff
    FinishInitialMatrix(InitialMatrix, MeasurementsErrors, PolExponents, OutputFile);
    FinishColumnVector(DependentVariables, MeasurementsErrors);
    TransposeInitialMatrix(TransposedMtx, InitialMatrix);
    MultiplyTransposeAndOriginal(MtxT_Mul_Mtx, InitialMatrix, TransposedMtx, OutputFile);
    MultiplyTransposeAndVector(MtxT_Mul_Vec, DependentVariables, TransposedMtx, OutputFile);

    free(TransposedMtx);

    //Invert MtxT_Mul_Mtx with Gauss-Jordan elimination
    FunctionWaveOfGaussJordan(MtxT_Mul_Mtx, GaussedMatrix, IdentityMatrix, GaussedHalfMatrix, OutputFile);

    //Calculating the fitting parameters
    MultiplyGaussedWith_MtxT_Mul_Vec(FittedParameters, MtxT_Mul_Vec, GaussedHalfMatrix);

    //Output
    PrintFittedParameters(OutputFile, FittedParameters);

    free(InitialMatrix);
    free(DependentVariables);
    free(MeasurementsErrors);
    free(MtxT_Mul_Vec);
    free(FittedParameters);
    free(PolExponents);
    free(MtxT_Mul_Mtx);
    free(GaussedMatrix);
    free(GaussedHalfMatrix);
    free(IdentityMatrix);

    return 0;
}
