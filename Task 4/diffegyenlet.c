#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* THE CARTESIAN COORDINATESYSTEM HAS BEEN USED IN THE PROGRAM. THE COORDIANTES ARE INDICATED WITH X AND Y.
   THE X CORDINATE */

//STRUCK FOR GLOBAL DATAS
typedef struct
{
    int RowsOfDataVector;
    double* DataVector;
} DataStruct;

//GLOBAL VARIABLES
const char FILE_NAME[] = "Results.txt";
double LengthOfSteps = 1;
int NumberOfSteps = 24192000;
DataStruct OrbitalValues;

//Declare functions
DataStruct AllocateMemoryForDataVectors();
DataStruct ConstantValues();
DataStruct InitialVariables();
double Differenciate();
void EulerStep();
double EnergyOfMoon();

//============== INITIALIZING ===============

DataStruct AllocateMemoryForDataVectors(int NumberOfVariables)
{
    OrbitalValues.DataVector = (double*)calloc(NumberOfVariables, sizeof(double));

    //Checking for memory
    if(OrbitalValues.DataVector == NULL)
    {
        printf("Allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    return OrbitalValues;
}

//============== INITIALIZING ===============

DataStruct ConstantValues(DataStruct GivenData)
{
    //Most accurate data i've found on NASA's website
	GivenData.DataVector[0] = 5.9724e24; //Weight of Earth (kg)
	GivenData.DataVector[1] = 7.342e22; //Weight of Moon (kg)
	GivenData.DataVector[2] = 4.055e8; //Apoapsis (m)
	GivenData.DataVector[3] = 9.64e2; //Min. orbital velocity (m*s^-1)
	GivenData.DataVector[4] = 3.633e8; //Periapsis (m)
	GivenData.DataVector[5] = 1.076e3; //Max. orbital velocity (m*s^-1)
	GivenData.DataVector[6] = 6.6712819e-11; //Gravitational constant (m^3*kg^-1*s^-2) from http://jlmlasheras.com/en/universal-gravitational-constant-exact-value.html

	GivenData.RowsOfDataVector = 7;

	return GivenData;
}

//WE FOLLOW THE MOTION CONVENTIONALLY STARTING FROM THE PERIAPSIS AT T = 0
DataStruct InitialVariables(DataStruct ProcessedData, DataStruct GivenData)
{
    //DATA AT T = 0, STARTING FROM THE PERIASPSIS
    //THESE DATA WILL EVENTUALLY CHANGING, THE Y COORDINATE WILL HAVE A MAXIMUM AT THE APOAPSIS
    //THE MAXIMUM VELOCITY IN X WILL BE THE MINIMAL ORBITAL VEOLCITY
	ProcessedData.DataVector[0] = GivenData.DataVector[4]; //AT THE PERIAPSIS THE X COORDINATE SHOULD BE THE LENGTH OF THE PERIAPSIS
	ProcessedData.DataVector[1] = 0; //THE Y COORDINATES SHOULD BE 0, MAXIMUM IS THE LENGTH OF APOAPSIS
	ProcessedData.DataVector[2] = 0; //THE VELOCITY OF X SHOULD BE 0 (MAXIMUM IS AT THE APOAPSIS, MINIMAL ORBIT VELOCITY IS THE MAXIMUM HERE)
	ProcessedData.DataVector[3] = GivenData.DataVector[5]; //THE VELOCITY OF Y SHOULD BE THE MAXIMAL ORBIT VELOCITY (AT THE PERIAPSIS)

	ProcessedData.RowsOfDataVector = 4;

	return ProcessedData;
}

//============== EULERING ===============

double Differenciate(DataStruct GivenData, DataStruct ProcessedData, int i)
{
    /* According to Kepler's laws and Newton mechanics, the following equation would
       describe the motion of Moon around the Earth: (d^2)x / (dt^2) = -G*M*x / r^3 and (d^2)y / (dt^2) = -G*M*y / r^3 */

    //r is the distance of the Moon from Earth. It's sqrt(x^2+y^2) in Cartesian coordinates.
    double distance = sqrt(pow(ProcessedData.DataVector[0], 2) + pow(ProcessedData.DataVector[1], 2));
    double dy; //THIS WILL BE THE DERIVATED OBJECT

    //Changes the coordinates
    if(i < 2)
    {
        dy = ProcessedData.DataVector[i+2];
    }

    //Changes the velocities
    if(i >= 2)
    {
        dy = -GivenData.DataVector[6] * GivenData.DataVector[0] * ProcessedData.DataVector[i-2] / pow(distance, 3);
    }

    return dy;
}

void EulerStep(DataStruct GivenData, DataStruct ProcessedData, FILE* OutputFile)
{
    //Running indexes for loops
    register int i;
    register int j;

    //Time
    double Time = 0;

    //ProcessedData.RowsOfDataVector have 4 values, therefore yn will have 4 elements
    double yn[4];

    //Step with the independent variable
    for(j = 0; j <= NumberOfSteps; j++)
    {
        //Fprintf results into the output file for fitting
        fprintf(OutputFile, "%g m\t%g m\t%f m/s\t%f m/s\t%g J\t", ProcessedData.DataVector[0], ProcessedData.DataVector[1], ProcessedData.DataVector[2], ProcessedData.DataVector[3], EnergyOfMoon(ProcessedData, GivenData));
        fprintf(OutputFile, "%g s\n", Time);
        Time += LengthOfSteps;

        //Step with the dependent variable
        //Euler step (stolen from the the website of fiznum)
        for(i = 0; i < ProcessedData.RowsOfDataVector; i++)
        {
            yn[i] = ProcessedData.DataVector[i] + LengthOfSteps * Differenciate(GivenData, ProcessedData, i);
        }

        //This updates the ProcessedData.DataVector with the actual orbital parameters
        for(i = 0; i < ProcessedData.RowsOfDataVector; i++)
        {
            ProcessedData.DataVector[i] = yn[i];
        }
    }
}

double EnergyOfMoon(DataStruct ProcessedData, DataStruct GivenData)
{
    double Energy;

    //E = 1/2 * m * v^2 - m*G*M/r
    Energy = fabs(GivenData.DataVector[2] * sqrt(pow(ProcessedData.DataVector[2], 2) + pow(ProcessedData.DataVector[3], 2)) / 2 - GivenData.DataVector[2] * GivenData.DataVector[6] * GivenData.DataVector[1] / sqrt(pow(ProcessedData.DataVector[0], 2) + pow(ProcessedData.DataVector[1], 2)));

    return Energy;
}

int main(void)
{
    //Variables
    FILE* OutputFile;   //File

    //Declarations
    DataStruct GivenData;
    DataStruct ProcessedData;

    //Allocating
    GivenData = AllocateMemoryForDataVectors(7);
    ProcessedData = AllocateMemoryForDataVectors(4);

    GivenData = ConstantValues(GivenData);
    ProcessedData = InitialVariables(ProcessedData, GivenData);

    //OutputFile
    OutputFile = fopen(FILE_NAME, "w+");

    //Derivate and Euler
    EulerStep(GivenData, ProcessedData, OutputFile);

    //Close file
    fclose(OutputFile);

    //Free memory
    free(GivenData.DataVector);
    free(ProcessedData.DataVector);

    exit(0);
}
