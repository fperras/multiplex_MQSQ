#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
static float Pi = 3.1415926535897932384626433;

/*
Multiplex Processing Program is used to process A 3D multiplex Bruker dataset,
organized as t3(1Q) t2(MQphase) t1(MQevolution)
S. Chandra Shekar & Frederic Alain Perras, 01/19/2023
*/
void mPlx(int TD1, int TD2, int TD3, int TD1new, int TD2new, int MQmx);

int main() {
    int TD1, TD2, TD3, MQmx, TD2new, TD1new;

    printf("Please enter the value of TD1\n");
    scanf("%d", &TD1);
    printf("Please enter the value of TD2\n");
    scanf("%d", &TD2);
    printf("Please enter the value of TD3\n");
    scanf("%d", &TD3);

    MQmx  = (TD2-1)/2;
    TD2new= MQmx+1;
    TD1new= TD1*2;

    mPlx(TD1, TD2, TD3, TD1new, TD2new, MQmx);

    system("PAUSE");
    return 0;
} // main

void mPlx(int TD1, int TD2, int TD3, int TD1new, int TD2new, int MQmx) {

    vector<vector<vector<int>>> fDat(TD1,vector<vector<int>>(TD2,vector<int>(TD3,0)));
    vector<vector<vector<float>>> FMQ(TD1new,vector<vector<float>>(TD2new,vector<float>(TD3,0)));

    char filInp[80], filOut[80];
    FILE *fp, *fpLog;
    int TD3half, currntVal, k, k2, kk, l, m,  m2, q;
    float phi, phiMQ, cPhi, sPhi, a0, b0;

    if(TD3%2 != 0)
        printf("WARNING: TD3 not even");

    TD3half= TD3/2;

    printf("Please enter the name of the experimental ser file\n");
    scanf("%s", filInp);
    printf("Please enter the desired name for the processed ser file\n");
    scanf("%s", filOut);

    fpLog= fopen("log.txt", "w");
    fprintf(fpLog,"The original datapoint counts (TD3, TD2, TD1) are (%d, %d, %d) \n", TD3,  TD2,  TD1);
    fprintf(fpLog,"The new datapoint counts (TD3, TD2, TD1) are (%d, %d, %d)\n",  TD3,  TD2new, TD1new);
    fprintf(fpLog,"\nThe following parameters should be set in Topspin:\n");
    fprintf(fpLog,"               F3    F2    F1\n");
    fprintf(fpLog,"    phmod      pk    no    pk            (repeat for 's phmod')\n");
    fprintf(fpLog,"    ftmod     fqc    no    no            (repeat for 's ftmod')\n");
    fprintf(fpLog,"    mc2              QF    echo-antiecho (repeat for 's mc2')\n");
    fprintf(fpLog,"\nThe name of the original ser file is %s\n", filInp);
    fprintf(fpLog,"The name of the processed ser file is %s\n", filOut);
    fclose(fpLog);

    printf("\nThe original datapoint counts (TD3, TD2, TD1) are (%d, %d, %d) \n", TD3,  TD2,  TD1);
    printf("The new datapoint counts (TD3, TD2, TD1) are (%d, %d, %d)\n",  TD3,  TD2new, TD1new);
    printf("\nThe following parameters should be set in Topspin:\n");
    printf("               F3    F2    F1\n");
    printf("    phmod      pk    no    pk            (repeat for 's phmod')\n");
    printf("    ftmod     fqc    no    no            (repeat for 's ftmod')\n");
    printf("    mc2              QF    echo-antiecho (repeat for 's mc2')\n");
    printf("\nThe name of the original ser file is %s\n", filInp);
    printf("The name of the processed ser file is %s\n\n", filOut);

    fp = fopen(filInp,  "rb");
    for(k=0; k<TD1; k++){
        for(l=0; l<TD2; l++){
            for(m=0; m<TD3; m++){
                fread(&currntVal, sizeof(int), 1, fp);
                fDat[k][l][m]= currntVal;
            }
        }
    }
    fclose(fp);

    for(k=0; k<TD1; k++){
        k2=k*2;
        for(q=0; q<MQmx+1; q++){
            for(m=0; m<TD3half; m++){
                m2=m*2;
                FMQ[k2][q][m2]    = 0;
                FMQ[k2][q][m2+1]  = 0;
                FMQ[k2+1][q][m2]  = 0;
                FMQ[k2+1][q][m2+1]= 0;
                for(l=0; l<TD2; l++){
                    phi= l*2*Pi/TD2;
                    phiMQ= -q*phi;
                    cPhi= cos(phiMQ);
                    sPhi= sin(phiMQ);
                    a0= fDat[k][l][m2];
                    b0= fDat[k][l][m2+1];
                    FMQ[k2][q][m2]     += cPhi*a0 +sPhi*b0;
                    FMQ[k2][q][m2+1]   +=-sPhi*a0 +cPhi*b0;
                    FMQ[k2+1][q][m2]   += cPhi*a0 -sPhi*b0;
                    FMQ[k2+1][q][m2+1] += sPhi*a0 +cPhi*b0;
                } // l-loop
            } // m-loop
        } // q-loop
    } // k-loop

    fp= fopen(filOut, "wb"); // output 3d mtrx: (2*TD1)*(MQmx+1)*TD3
    for(k=0; k<TD1; k++){
        k2=k*2;
        for(kk=0; kk<2; kk++) {
            for(l=0; l<TD2new; l++){
                for(m=0; m<TD3half; m++){
                    m2=m*2;
                    currntVal= FMQ[k2+kk][l][m2];
                    fwrite(&currntVal, sizeof(int), 1, fp);
                    currntVal= FMQ[k2+kk][l][m2+1];
                    fwrite(&currntVal, sizeof(int), 1, fp);
                } // m-loop
            } // l-loop
        } // kk-loop
    } // k-loop
    fclose(fp);
} // function mPlx
