#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const char FILE_NAME[] = "matrix.txt";

// Machine epsilon of float32
const double floateps = 1.1920929e-07;



// Count the rows of the matrix in the input file
unsigned int CountDimensionOfMatrix(FILE* InputFile)
{
    // Start from 1, thus the last line of the file got also included
    // It is needed, because the while+getc beloew at EOF breaks the loop automatically
    unsigned int Counter = 0;

    // Char variables to store characters of the file
    // getc() reads characters from the file into the CharacterCurrent variable
    // CharacterPrevious stores the previusly read character until the next step of the loop
    char CharacterPrevious = '\0';
    char CharacterCurrent;

    while(1)
    {
        CharacterCurrent = getc(InputFile);

        if(CharacterCurrent == EOF)
        {
            Counter++;
            break;
        }
        

        if(CharacterCurrent == '\n')
        {
            // If there is a double line break ('\n\n'), we don't need to count
            // it as a new row, we just skip it.
            if(CharacterPrevious == '\n')
            {
                continue;
            }

            // Increment the counter by 1, if a new row is found
            else
            {
                Counter++;
            }
        }

        // At the end of the current step of the loop, store the current character for the next round
        CharacterPrevious = CharacterCurrent;
    }

    return Counter;
}


// Returns a double* array with arbitrary size
double* AllocateMemory(unsigned int MemorySize)
{
    double* Matrix = (double*)calloc(MemorySize, sizeof(double));

    return Matrix;
}


// Read in a square matrix from the file if it's possible
double* ReadInMatrix(FILE* InputFile, unsigned int DimensionOfMatrix)
{
    double* InputMatrix = AllocateMemory(DimensionOfMatrix * DimensionOfMatrix);

    // Scan in and print out matrix to the standard output
    printf("Matrix read in:\n");
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j++)
        {
            fscanf(InputFile, "%lf", &InputMatrix[DimensionOfMatrix * i + j]);
            printf("%lf ", InputMatrix[DimensionOfMatrix * i + j]);
        }
        printf("\n");
    }

    printf("\n");

    return InputMatrix;
}


// Create identity matrix with DimensionOfMatrix rows and columns
double* CreateIdentityMatrix(unsigned int DimensionOfMatrix)
{
    double* IdentityMatrix = AllocateMemory(DimensionOfMatrix * DimensionOfMatrix);

    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j++)
        {
            if(i == j)
            {
                IdentityMatrix[DimensionOfMatrix * i + j] = 1;
            }

            else
            {
                IdentityMatrix[DimensionOfMatrix * i + j] = 0;
            }
        }
    }

    return IdentityMatrix;
}


// Create a matrix with M x 2M size, where M are the number of rows
// The righthand M x M part will be initially an identity matrix
// While the lefthand part will be a matrix, read from the input file
double* CreateFullMatrix(FILE* InputFile, unsigned int DimensionOfMatrix)
{

    // Allocate memory for parts of the matrix
    double* FullMatrix = AllocateMemory(2 * DimensionOfMatrix * DimensionOfMatrix);
    double* InputMatrix = AllocateMemory(DimensionOfMatrix * DimensionOfMatrix);
    double* IdentityMatrix = AllocateMemory(DimensionOfMatrix * DimensionOfMatrix);

    // Input matrix is read from a file by the ReadInmatrix function
    InputMatrix = ReadInMatrix(InputFile, DimensionOfMatrix);

    // Close the input file
    fclose(InputFile);

    printf("Matrix cloned:\n");
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j ++)
        {
            printf("%lf ", InputMatrix[DimensionOfMatrix * i + j]);
        }
        printf("\n");
    }

    printf("\n");

    // Identity matrix is arbitrarly created
    IdentityMatrix = CreateIdentityMatrix(DimensionOfMatrix);

    // Include input matrix
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j ++)
        {
            FullMatrix[2 * DimensionOfMatrix * i + j] = InputMatrix[DimensionOfMatrix * i + j];
        }
    }

    // Include identity matrix
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = DimensionOfMatrix; j < 2 * DimensionOfMatrix; j ++)
        {
            FullMatrix[2 * DimensionOfMatrix * i + j] = IdentityMatrix[DimensionOfMatrix * i + j];
        }
    }

    free(InputMatrix);
    free(IdentityMatrix);

    return FullMatrix;
}


// Gauss-Jordan elimination to find inverse of the matrix
// 1. The function first searches for the greatest element in the first column: if it's found in one of the rows, then that row is swapped with the first one.
// 2. In the next steps the function do the same for the other columns too, starting the loop at the second row, then the third, then so on.
// 3. After a step of the loop, the function are normalizing the whole row and performing substitutions to form the required inverse matrix.
// 4. At the end, the diagonal elements will be all just ones, and below the diagonal line will be full of zeros.
double* GaussJordan(FILE* InputFile, unsigned int DimensionOfMatrix)
{
    double* FullGaussMatrix = AllocateMemory(2 * DimensionOfMatrix * DimensionOfMatrix);

    // Create the matrix, which we'll work on
    FullGaussMatrix = CreateFullMatrix(InputFile, DimensionOfMatrix);

    // Temporary storage for indeces
    // Indicating the index of the row with current greatest element
    unsigned int TempStoreIndex;

    // Temporary storage for a double element of the matrix
    double TempStorage;

    // Temporary index for matrix elements at normalization
    double mat_ij;
    
    // Run column-wise through elements
    // The index (j) here indicates columns
    for(unsigned int j = 0; j < DimensionOfMatrix; j++)
    {
        // Set the jth element of the jth column is initially the greatest element
        TempStoreIndex = j;

        // Search for the greatest element in the jth column, starting from the (j+1)th row
        // It runs through all elements
        for(unsigned int i = (j + 1); j < DimensionOfMatrix; j++)
        {
            if(FullGaussMatrix[2 * DimensionOfMatrix * i + j] > FullGaussMatrix[2 * DimensionOfMatrix * TempStoreIndex + j])
            {
                TempStoreIndex = i;
            }
        }

        // Check if the greatest element is smaller, than the 32bit float machine epsilon
        // If yes, then the matrix is numerically singular, and cannot be inverted
        if((FullGaussMatrix[2 * DimensionOfMatrix * TempStoreIndex + j]) < floateps)
        {
            printf("Given matrix is singular and cannot be inverted! The program exits.");
            exit(0);
        }

        // Swapping rows, which has the greatest element in the column
        if(TempStoreIndex != j)
        {
            for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
            {
                TempStorage = FullGaussMatrix[2 * DimensionOfMatrix * j + k];
                FullGaussMatrix[2 * DimensionOfMatrix * j + k] = FullGaussMatrix[2 * DimensionOfMatrix * TempStoreIndex + k];
                FullGaussMatrix[2 * DimensionOfMatrix * TempStoreIndex + k] = TempStorage;
            }
        }

        // Perform substitutions and normalization
        for(unsigned int i = 0; i < DimensionOfMatrix; i++)
        {
            mat_ij = FullGaussMatrix[2 * DimensionOfMatrix * i + j];

            if(i != j)
            {
                for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
                {
                    FullGaussMatrix[2 * DimensionOfMatrix * i + k] -= (FullGaussMatrix[2 * DimensionOfMatrix * j + k]/FullGaussMatrix[2 * DimensionOfMatrix * j + j]) * mat_ij;
                }
            }

            else
            {
                for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
                {
                    FullGaussMatrix[2 * DimensionOfMatrix * i + k] /= mat_ij ;
                }
            }
        }
    }

    return FullGaussMatrix;
}


// Print the input matrix to the standard output
void PrintInverseMatrix(double* FullGaussMatrix, unsigned int DimensionOfMatrix)
{
    printf("\n\nThe invert matrix of the given one:\n");

    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = DimensionOfMatrix; j < 2 * DimensionOfMatrix; j++)
        {
            printf("%lf ", FullGaussMatrix[2 * DimensionOfMatrix * i + j]);
        }
        printf("\n");
    }
}


void GaussEliminationFunctions()
{
    // Variable for containing the size (rank) of the matrix
    unsigned int DimensionOfMatrix;

    // Open input file
    // It should contain a square matrix with double-point elements
    FILE* InputFile = fopen(FILE_NAME, "r");

    // Check if file opened successfully
    if(InputFile == NULL)
    {
        perror("Error during reading \"matrix.txt\"");
        EXIT_FAILURE;
    }

    // Count dimension of matrix
    DimensionOfMatrix = CountDimensionOfMatrix(InputFile);

    // double array, to contain the inverse matrix
    double* FullGaussMatrix = AllocateMemory(2 * DimensionOfMatrix * DimensionOfMatrix);

    // Inverting matrix
    FullGaussMatrix = GaussJordan(InputFile, DimensionOfMatrix);

    // Print the finalized inverse matrix
    PrintInverseMatrix(FullGaussMatrix, DimensionOfMatrix);

    free(FullGaussMatrix);
}


int main()
{
    char Answer;

    for(;;)
    {
        printf("This program will execute Gauss elimination on a given matrix, assuming that\n"
                "it's a square matrix. Please put a .txt file into the same directory with\n"
                "the gauss.c file. The file name must be 'matrix.txt'. If you want to choose\n"
                "an other file name, you can edit it on the very beginning of the .c file.\n"
                "The file should contain only numbers, separated by whitespaces or tabs.\n"
                "If everything's done and OK with your file, write 'Y' then press enter.\n"
                "If not, write 'N' and press enter, and the program exits.\n\n"
                "Are your files ready? (Y/N): ");
        scanf("%c", &Answer);

        switch(Answer)
        {
            case 'y':
            case 'Y':
                GaussEliminationFunctions();
                break;

            case 'n':
            case 'N':
                break;
        }
        break;
    }

    return 0;
}
