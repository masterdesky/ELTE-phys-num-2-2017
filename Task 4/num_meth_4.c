#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// Struct for orbital datas
typedef struct
{
    int RowsOfDataVector;
    double* DataVector;
} DataStruct;

// Typedef struct's name for calling contained variables
DataStruct OrbitalValues;

// Global control parameters
// 1. Length of Euler step
// 2. Number of total steps
const char FILE_NAME[] = "ResultsSample.txt";
double LengthOfSteps = 100;
int NumberOfSteps = 100000;

// Function prototypes
DataStruct AllocateMemoryForDataVectors();
DataStruct ConstantValues();
DataStruct InitialVariables();
double Differenciate();
void EulerStep();
double EnergyOfMoon();


// Allocating memory
DataStruct AllocateMemoryForDataVectors(int NumberOfVariables)
{
    OrbitalValues.DataVector = (double*)calloc(NumberOfVariables, sizeof(double));

    // Checking if memory callocation was successfull
    if(OrbitalValues.DataVector == NULL)
    {
        printf("Allocation error. Program exits.\n");
        exit(EXIT_FAILURE);
    }

    return OrbitalValues;
}


// Fill a DataStruct with the initial values of orbit parameters, and constants
DataStruct ConstantValues(DataStruct GivenData)
{
    // Most accurate data i've found
	GivenData.DataVector[0] = 5.9724e24;        // Mass of Earth (kg)
	GivenData.DataVector[1] = 7.342e22;         // Mass of Moon (kg)
	GivenData.DataVector[2] = 4.055e8;          // Apoapsis (m)
	GivenData.DataVector[3] = 9.64e2;           // Min. orbital velocity (m*s^-1)
	GivenData.DataVector[4] = 3.633e8;          // Periapsis (m)
	GivenData.DataVector[5] = 1.076e3;          // Max. orbital velocity (m*s^-1)
	GivenData.DataVector[6] = 6.6712819e-11;    // Gravitational constant (m^3*kg^-1*s^-2) from http://jlmlasheras.com/en/universal-gravitational-constant-exact-value.html

	GivenData.RowsOfDataVector = 7;

	return GivenData;
}


// We follow the motion with the ProcessedData.DataVector array
// Data starting at t = 0, when the Moon is the periaspsis
// These data will be updated at every step
DataStruct InitialVariables(DataStruct ProcessedData, DataStruct GivenData)
{
    // 'x'-coordine
    // At the periapsis the 'x' coordinate should be trivially the length of the periapsis
	ProcessedData.DataVector[0] = GivenData.DataVector[4];

    // 'y'-coordinate
    // At the periapsis, the 'y' coordinates should be 0. Its maximum is in the apoapsis, and its trivially the length of apoapsis there
	ProcessedData.DataVector[1] = 0;

    // Velocity, 'x'-component
    // At the periapsis , the velocity's 'x'-component should be trivially 0. Its maximum is in the apoapsis,
    // and that is also the minimum of the orbital velocity.
	ProcessedData.DataVector[2] = 0;

    // Velocity, 'y'-component
    // At the periapsis, the velocity's 'y'-component should be the maximal orbit velocity. Its minimum (0) is in the apoapsis.
	ProcessedData.DataVector[3] = GivenData.DataVector[5];

	ProcessedData.RowsOfDataVector = 4;

	return ProcessedData;
}


// According to Kepler's laws and Newton mechanics, the following equations would
// describe the motion of Moon around the Earth:
// 1. (d^2)x / (dt^2) = -G*M*x / r^3
// and
// 2. (d^2)y / (dt^2) = -G*M*y / r^3
double Differenciate(DataStruct GivenData, DataStruct ProcessedData, register unsigned int i)
{
    // In the equations above, 'r' is the distance of the Moon's and Earth's center of mass
    // It could be easily calculated in Cartesian coordinates as follows:
    // r = sqrt(x^2+y^2)
    double distance = sqrt(pow(ProcessedData.DataVector[0], 2) + pow(ProcessedData.DataVector[1], 2));
    
    // THIS WILL BE THE DERIVATED OBJECT
    double dCoorOrVel;

    // Updating the coordinates
    // Index for coordinates: (i == 0 || i == 1)
    if(i == 0 || i == 1)
    {
        dCoorOrVel = ProcessedData.DataVector[i+2];
    }

    // Updating the velocities
    // Index for velocities: (i == 2 || i == 3)
    if(i == 2 || i == 3)
    {
        dCoorOrVel = -GivenData.DataVector[6] * GivenData.DataVector[0] * ProcessedData.DataVector[i-2] / pow(distance, 3);
    }

    return dCoorOrVel;
}


void EulerStep(DataStruct GivenData, DataStruct ProcessedData, FILE* OutputFile)
{
    // Measuring time
    double Time = 0;

    // ProcessedData.RowsOfDataVector have 4 values, therefore yn will have 4 elements
    double yn[4];

    fprintf(OutputFile, "\'x\'-coordinates\t\t\'y\'-coordinates\t\tVelocity_x\t   Velocity_y\t      Kinetic E\t           Time\n");

    // Step with the independent variable
    for(register unsigned int j = 0; j <= NumberOfSteps; j++)
    {
        //Fprintf results into the output file for fitting
        fprintf(OutputFile, "%f m\t%f m\t%f m/s\t%f m/s\t%f J\t", ProcessedData.DataVector[0], ProcessedData.DataVector[1], ProcessedData.DataVector[2], ProcessedData.DataVector[3], EnergyOfMoon(ProcessedData, GivenData));
        fprintf(OutputFile, "%f s\n", Time);
        Time += LengthOfSteps;

        // Euler step
        for(register unsigned int i = 0; i < ProcessedData.RowsOfDataVector; i++)
        {
            yn[i] = ProcessedData.DataVector[i] + LengthOfSteps * Differenciate(GivenData, ProcessedData, i);
        }

        // This updates the ProcessedData.DataVector with the actual orbital parameters
        for(register unsigned int i = 0; i < ProcessedData.RowsOfDataVector; i++)
        {
            ProcessedData.DataVector[i] = yn[i];
        }
    }
}


// Calculate kinetic energy of the Moon on orbit
double EnergyOfMoon(DataStruct ProcessedData, DataStruct GivenData)
{
    double Energy;

    //E = 1/2 * m * v^2 - m*G*M/r
    Energy = fabs(GivenData.DataVector[2] * sqrt(pow(ProcessedData.DataVector[2], 2) + pow(ProcessedData.DataVector[3], 2)) / 2 - GivenData.DataVector[2] * GivenData.DataVector[6] * GivenData.DataVector[1] / sqrt(pow(ProcessedData.DataVector[0], 2) + pow(ProcessedData.DataVector[1], 2)));

    return Energy;
}


int main()
{
    // Output file to containt orbital data
    FILE* OutputFile;
    OutputFile = fopen(FILE_NAME, "w+");

    // Used structs
    DataStruct GivenData;
    DataStruct ProcessedData;

    // Allocating memory
    GivenData = AllocateMemoryForDataVectors(7);
    ProcessedData = AllocateMemoryForDataVectors(4);

    // Fill structs with required data
    GivenData = ConstantValues(GivenData);
    ProcessedData = InitialVariables(ProcessedData, GivenData);

    // Derivate and Euler
    EulerStep(GivenData, ProcessedData, OutputFile);

    // Close file
    fclose(OutputFile);

    // Free memory
    free(GivenData.DataVector);
    free(ProcessedData.DataVector);

    return 0;
}
