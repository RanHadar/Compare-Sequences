/**
 * @file CompareSequences.c
 * @author  Ran Hadar <ran.hadar1@mail.huji.ac.il>
 * @version 1.0
 * @date 14/11/2018
 *
 * @brief Program that reads sequences from a file and calculates matches.
 *
 * @section DESCRIPTION
 * Input  : txt files.
 * Process: reads sequences from a text file and calculates matches between all of them
 * Output : prints the results to the screen
 */
//************************************  includes ***********************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
//************************************* define *************************************************
#define MAX_LINE 100
#define MAX_SEQ 10
#define MIN_ARGS 2
//*********************************** functions declarations ***********************************
int findMax(int num1, int num2, int num3);
int** setMatrix(const char *str1, const char *str2);
int fillMatrix(int **compMat, const char *seq1, const char *seq2, int match, int mismatch, int gap);
void compSeqs(char *seqArr[MAX_SEQ], int seqCounter, int match, int mismatch, int gap);
char* cleanLine(char fileLine[MAX_LINE]);
void mallocFree(char *seqArr[MAX_SEQ], int seqCounter);
//***********************************************************************************************
/**
 *main function responsible for running the program, gets the args and send to the
 * relevant function for calculation.
 * @param argc - number of args
 * @param argv - txt file, match, mismatch, gap
 * @return
 */
int main(int argc, char *argv[])
{
    char *seqArr[MAX_SEQ] = {};
    char fileLine[MAX_LINE];
    int seqCounter = 0;
    char *strCheck;

    if (argc < MIN_ARGS)
    {
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(argv[1], "r");

    if (fp == NULL)
    {
        exit(EXIT_FAILURE);
    }

    strCheck = fgets(fileLine, MAX_LINE, fp);

    while(strCheck)
    {
        if (fileLine[0] == '>') //seq header was found
        {
            strCheck = fgets(fileLine, MAX_LINE, fp);
            char *str = (char *) malloc(sizeof(fileLine));
            if(str == NULL)
            {
                mallocFree(seqArr, seqCounter);
                exit(EXIT_FAILURE);
            }
            strcat(str, cleanLine(fileLine));
            seqArr[seqCounter] = str;
            strCheck = fgets(fileLine, MAX_LINE, fp);
            while (fileLine[0] != '>' && strCheck) //if its still the same seq
            {
                str = (char*) realloc(str, sizeof(fileLine));
                if(str == NULL)
                {
                    mallocFree(seqArr, seqCounter);
                    exit(EXIT_FAILURE);
                }
                strcat(str, cleanLine(fileLine)); //adds it to the previous one
                seqArr[seqCounter] = str;
                strCheck = fgets(fileLine, MAX_LINE, fp);
            }
            seqCounter++;
        }
        else
        {
            strCheck = fgets(fileLine, MAX_LINE, fp);
        }
    }

    char *temp;
    int match = (int)strtol(argv[2], &temp, 10);
    int mismatch = (int)strtol(argv[3], &temp, 10);
    int gap = (int)strtol(argv[4], &temp, 10);

    compSeqs(seqArr, seqCounter, match, mismatch, gap);
    mallocFree(seqArr, seqCounter); //free all memory after the output

    return 0;
}

/**
 *Responsible for getting a line from the file(seq) and return it clear from \n \r
 * @param fileLine
 * @return clean line
 */
char* cleanLine(char *fileLine)
{
    for(int i = 0; i < (int)strlen(fileLine); i++)
    {
        if(!isalpha(fileLine[i]))
        {
            fileLine[i] = '\0';
        }
    }
    return fileLine;
}

/**
 *Responsible for running matches between the seqs in the file, each couple sends to
 * another function for calculations
 * @param seqArr
 * @param seqCounter
 * @param match
 * @param mismatch
 * @param gap
 */
void compSeqs(char *seqArr[MAX_SEQ], int seqCounter, int match, int mismatch, int gap)
{
    int **compMat;
    int score = 0;
    for(int i = 0; i < seqCounter - 1; i++)
    {
        for(int j = i + 1; j <= seqCounter - 1 ; j++)
        {
            compMat = setMatrix(seqArr[i], seqArr[j]);
            score = fillMatrix(compMat, seqArr[i], seqArr[j], match, mismatch, gap);
            printf("Score for alignment of seq%d to seq%d is %d\n", i + 1, j + 1, score);
            free(compMat);
        }
    }
}

/**
 * Sets a clean multi-array(Matrix) in the memory for calculations
 * @param str1
 * @param str2
 * @return multi_dimensional_array pointer
 */
int** setMatrix(const char *str1, const char *str2)
{
    int **compMat = (int **)malloc(sizeof(int*) * (strlen(str1) + 1));
    if(compMat == NULL)
    {
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < (int)strlen(str1) + 1; i++)
    {
        compMat[i] = (int*)malloc(sizeof(int) * (strlen(str2)) + 1);
    }
    return compMat;
}

/**
 * Responsible for filling the matrix with the relevant calculations for a couple of seqs.
 * @param compMat
 * @param seq1
 * @param seq2
 * @param match
 * @param mismatch
 * @param gap
 * @return the result(right down edge number in the matrix)
 */
int fillMatrix(int **compMat, const char *seq1, const char *seq2, int match, int mismatch, int gap)
{
    int row = 0, colum = 0;
    int tempGup, tempGleft, tempM = 0;

    for(row  = 0; row < (int)strlen(seq1) + 1; row++)
    {
        for(colum = 0; colum  < (int)strlen(seq2) + 1; colum++)
        {
            if(row == 0)
            {
                compMat[row][colum] = gap * colum;
                continue;
            }
            if(colum == 0)
            {
                compMat[row][colum] = gap * row;
                continue;
            }
            else
            {
                if(seq1[row - 1] == seq2[colum - 1])
                {
                    tempM = compMat[row - 1][colum - 1] + match;
                }
                else
                {
                    tempM = compMat[row - 1][colum - 1] + mismatch;
                }
                tempGup =  compMat[row - 1][colum] + gap;
                tempGleft = compMat[row][colum -1] + gap;
                compMat[row][colum] = findMax(tempGleft, tempGup, tempM);
            }
        }
    }
    return compMat[strlen(seq1)][strlen(seq2)];
}

/**
 * Responsible for finding the max num between 3 numbers.
 * @param num1
 * @param num2
 * @param num3
 * @return max num
 */
int findMax(int num1, int num2, int num3)
{

    if(num1 >= num2)
    {
        if(num1 >= num3)
        {
            return num1;
        }
        else
        {
            return num3;
        }
    }
    else
    {
        if(num2 >= num3)
        {
            return num2;
        }
        else
        {
            return num3;
        }
    }
}

/**
 * Responsible for free all memory after use of malloc in other function.
 * @param seqArr
 * @param seqCounter
 */
void mallocFree(char *seqArr[MAX_SEQ], int seqCounter)
{
    for(int i = 0; i < seqCounter; i++)
    {
        free(seqArr[i]);
        seqArr[i] = NULL;
    }
}