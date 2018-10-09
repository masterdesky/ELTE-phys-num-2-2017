#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const char FILE_NAME[] = "matrix.txt";

const float floateps = 1.1920929e-07;


// Count the rows of the matrix in the input file
int CountDimensionOfMatrix(FILE* InputFile)
{
    // Start from 1, thus the last line of the file got also included
    // It is needed, because the while+getc beloew at EOF breaks the loop automatically
    int Counter = 1;

    // Char variables to store characters of the file
    // getc() reads characters from the file into the CharacterCurrent variable
    // CharacterPrevious stores the previusly read character until the next step of the loop
    char CharacterPrevious;
    char CharacterCurrent = getc(InputFile);

    while(CharacterCurrent != EOF)
    {
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
                Counter += 1;
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
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j++)
        {
            fscanf(InputFile, "%lf", &InputMatrix[DimensionOfMatrix * i + j]);
            printf("%lf ", InputMatrix[DimensionOfMatrix * i + j]);
        }
        printf("\n");
    }

    // Close the input file
    fclose(InputFile);

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

    // Identity matrix is arbitrarly created
    IdentityMatrix = CreateIdentityMatrix(DimensionOfMatrix);

    // Include input matrix
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = 0; j < DimensionOfMatrix; j ++)
        {
            FullMatrix[DimensionOfMatrix * i + j] = InputMatrix[DimensionOfMatrix * i + j];
        }
    }

    // Include identity matrix
    for(unsigned int i = 0; i < DimensionOfMatrix; i++)
    {
        for(unsigned int j = DimensionOfMatrix; j < 2 * DimensionOfMatrix; j ++)
        {
            FullMatrix[DimensionOfMatrix * i + j] = IdentityMatrix[DimensionOfMatrix * i + j];
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
double* GaussJordan(FILE* Input, unsigned int DimensionOfMatrix)
{
    double* FullGaussMatrix = AllocateMemory(2 * DimensionOfMatrix * DimensionOfMatrix);

    // Create the matrix, which we'll work on
    FullGaussMatrix = CreateFullMatrix(Input, DimensionOfMatrix);

    // Temporary storage for indeces
    // Indicating rows with current greatest element, and newly found row
    // with even greater element
    int TempStoreFirst;
    int TempStoreSecond;

    // Temporary index for matrix elements at normalization
    int mat_ij;
    
    // Run column-wise through elements
    // The index (j) here indicates columns
    for(unsigned int j = 0; j < DimensionOfMatrix; j++)
    {
        // Set the jth element of the jth column is initially the greatest element
        TempStoreFirst = j;

        // Search for the greatest element in the jth column, starting from the jth row
        // it runs through all elements
        for(unsigned int i = j + 1; j < DimensionOfMatrix; j++)
        {
            if(FullGaussMatrix[DimensionOfMatrix * i + j] > FullGaussMatrix[DimensionOfMatrix * TempStoreFirst + j])
            {
                TempStoreFirst = j;
            }
        }

        // Check if the greatest element is smaller, than the 32bit float machine epsilon
        // If yes, then the matrix is numerically singular, and cannot be inverted
        if((FullGaussMatrix[DimensionOfMatrix * TempStoreFirst + j]) < floateps)
        {
            printf("Given matrix is singular and cannot be inverted! The program exits.");
            exit(0);
        }

        // Swapping rows, which has the greatest element in the column
        if(TempStoreFirst != j)
        {
            for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
            {
                TempStoreSecond = FullGaussMatrix[DimensionOfMatrix * j + k];
                FullGaussMatrix[DimensionOfMatrix * j + k] = FullGaussMatrix[DimensionOfMatrix * TempStoreFirst + k];
                FullGaussMatrix[DimensionOfMatrix * TempStoreFirst + k] = TempStoreSecond;
            }
        }

        // Perform substitutions and normalization
        for(unsigned int i = 0; i < DimensionOfMatrix; i++)
        {
            mat_ij = FullGaussMatrix[DimensionOfMatrix * i + j];

            if(i != j)
            {
                for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
                {
                    FullGaussMatrix[DimensionOfMatrix * i + k] -= (FullGaussMatrix[DimensionOfMatrix * j + k]/FullGaussMatrix[DimensionOfMatrix * j + j]) * mat_ij;
                }
            }

            else
            {
                for(unsigned int k = 0; k < 2 * DimensionOfMatrix; k++)
                {
                    FullGaussMatrix[DimensionOfMatrix * i + k] /= mat_ij ;
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
            printf("%lf ", FullGaussMatrix[DimensionOfMatrix * i + j]);
        }
        printf("\n");
    }
}


void GaussEliminationFunctions()
{
    // Variable for containing the size (rank) of the matrix
    unsigned int DimensionOfMatrix;

    // double array, to contain the inverse matrix
    double* FullGaussMatrix = AllocateMemory(2 * DimensionOfMatrix * DimensionOfMatrix);


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

    // Inverting matrix
    FullGaussMatrix = GaussJordan(InputFile, DimensionOfMatrix);

    // Print the finalized inverse matrix
    PrintInverseMatrix(FullGaussMatrix, DimensionOfMatrix);

    free(FullGaussMatrix);
    
    exit(0);
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
                exit(0);
        }
        printf("\n");
    }

    return 0;
}