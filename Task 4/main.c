#include <stdio.h>
#include <stdlib.h>

//STRUCK FOR GLOBAL DATAS
typedef struct
{
    int RowsOfDataVector;
    double* DataVector;
} DataStruct;

//GLOBAL VARIABLES
const char* FILE_NAME;
double LengthOfSteps = 10;
int NumberOfSteps = 10000000;

//Declare functions
char* InitializingOutput();
double* AllocateMemoryForVariables();
double* AllocateMemoryForDataVector();
void AddDataToStruct();
void Differenciate();
void EulerStep();

//============== ALLOCATING ===============

double* AllocateMemoryForVariables()
{
    double* Variable = (double*)calloc(NumberOfSteps, sizeof(double));

    //Checking for memory
    if(Variable == NULL)
    {
        printf("Variables allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    return Variable;
}

double* AllocateMemoryForDataVector()
{
}

//============== INITIALIZING ===============

char* InitializingOutput()
{
    //Chars for file name
    char FileNameFirst[] = "Results_for_";
    char FileNameLast[] = "_steps";

    //Char*s for file name
    char* StepNumber = (char*)calloc(NumberOfSteps, sizeof(char));
    char* FILE_NAME = (char*)calloc(strlen(FileNameFirst) + strlen(StepNumber) + strlen(FileNameLast), sizeof(char));

    //Write NumberOfSteps integer into StepNumber string
    sprintf(StepNumber, "%d", NumberOfSteps);

    //Build custom file name
    FILE_NAME[0] = '\0';
    strcat(FILE_NAME, FileNameFirst);
    strcat(FILE_NAME, StepNumber);
    strcat(FILE_NAME, FileNameLast);

    //Checking for memory
    if(StepNumber == NULL)
    {
        printf("StepNumber allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    //Checking for memory
    if(FILE_NAME == NULL)
    {
        printf("FILE_NAME allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }
}

DataStruct InitialVariables(DataStruct GivenData)
{
    GivenData.DataVector = (double*)calloc(7, sizeof(double));

    //Checking for memory
    if(GivenData.DataVector == NULL)
    {
        printf("GivenData allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    //Most accurate data i've found on NASA's website
	GivenData.DataVector[0] = 5.9724e24; //Weight of Earth (kg)
	GivenData.DataVector[1] = 7.342e22; //Weight of Moon (kg)
	GivenData.DataVector[2] = 4.055e8; //Apogee (m)
	GivenData.DataVector[3] = 9.64e2; //Min. orbital velocity (m*s^-1)
	GivenData.DataVector[4] = 3.633e8; //Perigee (m)
	GivenData.DataVector[5] = 1.076e3; //Max. orbital velocity (m*s^-1)
	GivenData.DataVector[6] = 6.67384e-11; //Gravitational constant (m^3*kg^-1*s^-2)

	GivenData.RowsOfDataVector = 7;

    return GivenData;
}

DataStruct DependentVariables(DataStruct ProcessedDatam, DataStruct GivenData)
{
    ProcessedData.DataVector = (double*)calloc(4, sizeof(double));

    //Checking for memory
    if(ProcessedData.DataVector == NULL)
    {
        printf("ProcessedData allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    //Declaring independent variables
	ProcessedData.DataVector[0] = GivenData.DataVector[4]; //
	ProcessedData.DataVector[1] = 0; //
	ProcessedData.DataVector[2] = 0; //
	ProcessedData.DataVector[3] = GivenData.DataVector[5]; //

	ProcessedData.RowsOfDataVector = 4;

    return ProcessedData;

}

//============== EULERING ===============

double Differenciate(DataStruct GivenData, DataStruct ProcessedData, double* dy)
{

}

void EulerStep(DataStruct GivenData, DataStruct ProcessedData)
{
    //Running index for loop
    register int i;

    //Euler step
    for(i = 0; i < GivenData.RowsOfDataVector; i++)
    {
        yn[i] += y[i] + LengthOfSteps * Differenciate(GivenData, ProcessedData);
    }
}

int main(void)
{
    //Variables
    FILE* OutputFile;   //File
    double* x;          //Vector
    double* y;          //Vector
    double* yn;         //Vector
    double* dy;         //Vector

    DataStruct GivenData;
    GivenData = InitialVariables(GivenData);

    //OutputFile
    OutputFile = fopen(FILE_NAME, "w+");

    //Allocating
    x = AllocateMemoryForVariables();
    y = AllocateMemoryForVariables();

    //Derivate and Euler
    EulerStep(GivenData, );

    exit(EXIT_SUCCESS);
}
