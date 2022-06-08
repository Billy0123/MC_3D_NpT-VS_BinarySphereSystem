#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

/*State of optional functions:  (remember to always update text below after changing this file)
 -probDensDistFunMode - INACTIVE
 -FREQUENT UPDATES OF NEIGHBOUR LIST - INACTIVE
 -consfStepByStep - INACTIVE
 -tworzenie IDENTYCZNEGO rozkładu binarnego jak przy gaussie - INACTIVE
 -printActualCycle - INACTIVE
 -checkCorrectNeighbourNumberWhileUpdating - PARTIALLY ACTIVE (in steps: 2), in step 1 big specificMod!
 -obliczanie epsilonów (wyniki) przy GLOBALNYCH średnich Hij, avEij - ACTIVE
 -skalowanie rozmiarow czastek - INACTIVE
*/

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode,neighUpdatingFrequency;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaV=0.1,deltaVND=0.1,
       deltaDiameter,randomStartStep[2],VcpPerParticle,
       neighRadius,neighRadius2,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],dr[4],boxMatrixInnerParameters[11]; //innerParameters: detBoxMatrix,volume,bm1221m1122,bm0122m0221,bm0211m0112,bm1220m1022,bm0022m0220,bm0210m0012,bm1120m1021,bm0021m0120,bm0110m0011;
int alpha,beta;  //parametry sterujące
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferD[20],bufferInitMode[5],bufferFolderIndex[5],bufferAlpha[20],bufferBeta[20],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     loadConfigurationsFileName[200]="Configurations",
     probDensDistFunFileName[200]="ProbDensDistFun",              //probDensDistFunMode
     probDensDistFunResultsFileName[200]="ProbDensDistFunRes",
     loadedJOBID[50]="j-none";
/////////////////  PARTICLE functions{
typedef struct particle {
    double r[3], normR[3];  //x,y,z
    int neighbours[50], neighCounter;
    int type;           //type: 0-black (large), 1-white (small),  tutaj typ w sumie nie byłby potrzebny, ale jak już jest zrobiony, to zostawiam - nie będzie trzeba dodatkowych warunków na rysowanie czarno/białych dysków w mathematice
    double diameter;    //układy binarne: diameter=type?sigma[=1]-deltaDiameter:sigma[=1]+deltaDiameter; w przypadku initMode=7(random) i beta=1 średnice zadane są rozkładem Gaussa - wówczas pole 'type' przyjmuje 1(white) dla diameter<sigma[=1] i 0(black) dla diameter>sigma[=1]
} particle;


void updateBoxMatrixParameters (double box[3][3], double boxInnerParameters[11]) {
    boxInnerParameters[0]=box[0][0]*box[1][1]*box[2][2]+box[0][1]*box[1][2]*box[2][0]+
                          box[0][2]*box[1][0]*box[2][1]-box[0][2]*box[1][1]*box[2][0]-
                          box[0][0]*box[1][2]*box[2][1]-box[0][1]*box[1][0]*box[2][2];
    boxInnerParameters[1]=fabs(boxInnerParameters[0]);
    boxInnerParameters[2]=box[1][2]*box[2][1]-box[1][1]*box[2][2];
    boxInnerParameters[3]=box[0][1]*box[2][2]-box[0][2]*box[2][1];
    boxInnerParameters[4]=box[0][2]*box[1][1]-box[0][1]*box[1][2];
    boxInnerParameters[5]=box[1][2]*box[2][0]-box[1][0]*box[2][2];
    boxInnerParameters[6]=box[0][0]*box[2][2]-box[0][2]*box[2][0];
    boxInnerParameters[7]=box[0][2]*box[1][0]-box[0][0]*box[1][2];
    boxInnerParameters[8]=box[1][1]*box[2][0]-box[1][0]*box[2][1];
    boxInnerParameters[9]=box[0][0]*box[2][1]-box[0][1]*box[2][0];
    boxInnerParameters[10]=box[0][1]*box[1][0]-box[0][0]*box[1][1];
}

void getParticlesDistanceSquared (particle *p1, particle *p2, double boxMatrix[3][3]) {
    double normalizedDR[3],roundedNormalizedDR[3],DR[3];
    for (int i=0;i<3;i++) {
        normalizedDR[i]=p1->normR[i]-p2->normR[i];
        DR[i]=p1->r[i]-p2->r[i];
    }
    for (int i=0;i<3;i++) roundedNormalizedDR[i]=round(normalizedDR[i]);
    for (int j=0;j<3;j++) for (int i=0;i<3;i++) DR[i]-=roundedNormalizedDR[j]*boxMatrix[i][j];
    dr[3]=0; for (int i=0;i<3;i++) {dr[i]=DR[i]; dr[3]+=dr[i]*dr[i];}
}

void adjustNeighRadius (double volume, double mod) {
    neighRadius=mod*cbrt(volume)/cbrt(N);
    neighRadius2=neighRadius*neighRadius;
}

int updateNeighbourList (particle *particles, double boxMatrix[3][3], double volume, double neighModStartValue, bool checkForOnlyCloseNeighbours, int increasingSteps) { //zaczyna od sfery z 'neighModStartValue' (np. 1.4); jeżeli 'checkForOnlyCloseNeighbours==true', to zwiększa ją o czynnik 0.01 'increasingSteps' razy (przy 'increasingSteps==0' lub 'checkForOnlyCloseNeighbours==false', sprawdza TYLKO neighModStartValue)
    bool correctNeighValue=true;
    for (int k=0;k<=increasingSteps;k++) {
        double mod=neighModStartValue+0.01*k;
        adjustNeighRadius(volume,mod);
        for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
        for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
            getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
            if (dr[3]<neighRadius2) {
                particles[i].neighbours[particles[i].neighCounter++]=j;
                particles[j].neighbours[particles[j].neighCounter++]=i;
            }
        }
        int maxFoundNeighCounter=0;
        if (checkForOnlyCloseNeighbours) {
            correctNeighValue=true;
            for (int i=0;i<activeN;i++) {
                if (particles[i].neighCounter>maxFoundNeighCounter) maxFoundNeighCounter=particles[i].neighCounter;
                if (particles[i].neighCounter<12) { //bufor pewien, inną opcją (niż akceptacja więcej niż 12) to interwał mniejszy od 0.01, np. 0.001, ale to znowu będzie mogło, przy wyższych gęstościach, tym samym skutkować. Niech już gdzieś będzie te 13 czy 14 i tyle
                    correctNeighValue=false;
                    printf("FAIL at mod=%.2f: found particles[%d].neighCounter==%d. ",mod,i,particles[i].neighCounter);
                    break;
                }
            }
        }
        if (correctNeighValue) {
            if (checkForOnlyCloseNeighbours) printf("SUCCESS: Updating neighbour list succeeded at neighRadius mod=%.2f (maxFoundNeighCounter=%d)\n",mod,maxFoundNeighCounter);
            break;
        }
    }

    if (correctNeighValue) return 1;
    else {
        printf("ERROR: Updating neighbour list failed (mod: %.2f,0.01,%.2f)\n",neighModStartValue,neighModStartValue+0.01*increasingSteps);
        return 0;
    }
}

void checkSinglePeriodicBoundaryConditions (particle *p, double boxMatrix[3][3]) {
    for (int j=0;j<3;j++) {
        for (int i=0;i<3;i++) p->r[i]-=floor(p->normR[j])*boxMatrix[i][j];
        p->normR[j]=fmod(p->normR[j],1); if (p->normR[j]<0) p->normR[j]++;
    }
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[3][3]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

int createRandomGaps (particle *particles, double boxMatrix[3][3], double volume) {  //TODO: 1) głupie swapowanie - po co przesuwać cząstki w liście, jak wystarczy po prostu 'przerzucać' gapy na (aktualny) koniec listy i zmniejszać jej wielkość (ucinając automatycznie ten ostatni element). 2) nie zamienia wszystkich pól obiektów cząstek (przestarzałe)
    printf("Creating %d random gaps... ",gaps); fflush(stdout);
    adjustNeighRadius(volume,1.4);
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTRandom0to1(randomStartStep)*N);
            for (int j=0;j<i;j++) {
                getParticlesDistanceSquared(&particles[gapsIndexes[i]],&particles[gapsIndexes[j]],boxMatrix);
                if (dr[3]<neighRadius2) {
                    i--;
                    innerAttempt++;
                    break;
                }
                if (j+1==i) innerAttempt=0;
            }
            if (innerAttempt>100000) break;
            if (innerAttempt==0 && i+1==gaps) allReady=1;
        }
    } while (!allReady && attempt<10000000);

    if (attempt>=10000000) {
        printf("ERROR: Couldn't create %d random gaps in %d steps.\n",gaps,attempt);
        return 0;
    } else {
        bool change; do {
            change=false;
            for (int i=gaps-1;i>0;i--) {
                if (gapsIndexes[i-1]>gapsIndexes[i]) {
                    int buffer=gapsIndexes[i];
                    gapsIndexes[i]=gapsIndexes[i-1]; gapsIndexes[i-1]=buffer;
                    change=true;
                }
            }
        } while (change);
        int actualGapIndex=0;
        for (int i=0;i<activeN;i++) {
            if (actualGapIndex<gaps) while (i+actualGapIndex==gapsIndexes[actualGapIndex]) {
                actualGapIndex++;
                if (actualGapIndex==gaps) break;
            }
            for (int j=0;j<3;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
                //TODO: brakujace pola do przepisania
            }
        }
        printf("done\n"); fflush(stdout);
        return 1;
    }
}

double getInitDiameterBinarySystem (int type) {
    double diameter;
    switch (alpha) {
        case 0: diameter=type?1.0-deltaDiameter:1.0+deltaDiameter; break;
        case 1: diameter=type?1.0-deltaDiameter:1.0; break;
    }
    return diameter;
}

double getInitDiameterByDistributionFunction (double (*function)(double), double randomMin, double randomMax) {  //funkcja(double) MUSI być przeskalowana tak, by jej (gęstości prawdopobieństwa) MAX był równy 1
    double diameter;
    do diameter=randomMin+(randomMax-randomMin)*MTRandom0to1(randomStartStep); while (MTRandom0to1(randomStartStep)>function(diameter));
    return diameter;
}

double gaussianFunction (double x) {
    double my=1.0, delta=deltaDiameter;  //my-wartosc srednia rozkladu, delta-odchylenie standardowe
    return /* 1/delta/sqrt(2*M_PI)* */exp(-pow(x-my,2)/(2*delta*delta));  //wykomentowany człon to amplituda funkcji w punkcie MAX - w taki sposób funkcja jest skalowana do maxValue=1
}

double cosinusFunction (double x) {
    return cos((x-1.0)/deltaDiameter);
}

void invertTable (int table[50][50], int size[2]) {
    for (int i=0;i<size[0];i++) for (int j=0;j<size[1];j++)
        if (table[i][j]==0) table[i][j]=1; else table[i][j]=0;
}

int initPositions (particle *particles, double boxMatrix[3][3], double boxParameters[11], double matrixOfParticlesSize[3], int n[3], double matrixCellXYZ[6][6][6][3], double pacFrac) {
    if (generatorStartPoint==0) {
        printf("Setting start position of p-random number generator to actual CPU time (for INIT PROCEDURES)...\n");
        InitRandomMT();
    } else {
        printf("Setting start position of p-random number generator to %ld (for INIT PROCEDURES)...\n",generatorStartPoint);
        InitMT((unsigned int)generatorStartPoint);
    }

    double mod=cbrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]/matrixOfParticlesSize[2]), interval[3][3], actualPosition[3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int layerCounter=0, rowCounter=0, columnCounter=0;

    //w 2D tutaj następowało zdefiniowanie dysków małych/dużych dla struktur uporządkowanych (initMode=0[HD], initMode=1...6, initMode=7[random], initMode=8,...)
    //w przypadku 3D przyjmuję: initMode=0[random], jeżeli chce się uzyskać HD, wystarczy zadać deltaDiameter=0. Jeżeli będzie potrzeba tworzenia jakichś uporządkowanych struktur, należy spojrzeć do programu 2D żeby zobaczyć jak to było robione.
    //Mechanizmy potrzebne do generowania struktur uporządkowanych (odpowiednie zadawanie typu cząstek od tablicy pTC, itd.) JEST tu już zrobione, ale sama tablica pTC NIE jest określana (to powinno znaleźć się w tym miejscu).
    /*int pTC[50][50], typeN[2], horizontalShift;                 //typeN w obu wymiarach [0-x,1-y] musi byc <=50 (wymiar pTC, ew. mozna powiekszyc)
    for (int i=0;i<50;i++) for (int j=0;j<50;j++) pTC[i][j]=0;  //pTC(particleTypeCell) stanowi elementarna komorke typow nadawanych czastkom, horizontalShift - przesuniecie w lewo (wyrazone w liczbie czastek) komorki elementarnej typow
    switch (initMode) {
        case 1: case 2: {typeN[0]=3; typeN[1]=2; horizontalShift=0;    //S1
            pTC[2][0]=1;
            pTC[0][1]=1;} break;
    }
    if ((initMode>0 && initMode<=6 && initMode%2==0) || (initMode>7 && initMode%2==1)) invertTable(pTC,typeN);*/

    for (int i=0;i<N;i++) {
        int cellNumber[3][3]={{columnCounter/n[0],rowCounter/n[1],layerCounter/n[2]},{columnCounter%n[0],rowCounter%n[1],layerCounter%n[2]}}; //cellNumber[0/1][X/Y/Z]: 0-numer komorki, 1-kolumna/rzad W komorce
        for (int j=0;j<3;j++) actualPosition[j]=cellNumber[0][0]*interval[j][0]+cellNumber[0][1]*interval[j][1]+cellNumber[0][2]*interval[j][2]+matrixCellXYZ[cellNumber[1][0]][cellNumber[1][1]][cellNumber[1][2]][j]*cbrt(pacFrac);

        for (int j=0;j<3;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=-(particles[i].r[0]*boxParameters[2]+particles[i].r[1]*boxParameters[3]+particles[i].r[2]*boxParameters[4])/boxParameters[0];
        particles[i].normR[1]=(particles[i].r[0]*boxParameters[5]+particles[i].r[1]*boxParameters[6]+particles[i].r[2]*boxParameters[7])/boxParameters[0];
        particles[i].normR[2]=-(particles[i].r[0]*boxParameters[8]+particles[i].r[1]*boxParameters[9]+particles[i].r[2]*boxParameters[10])/boxParameters[0];

        if (initMode==0) {
            if (beta==0) { //układ binarny dysków BLACK/WHITE w stosunku 50% (za co odpowiada późniejsza część algorytmu)
                //particles[i].type=getInitDiameterByDistributionFunction(gaussianFunction,1.0-6*deltaDiameter,1.0+6*deltaDiameter)<1?1:0;  //tworzenie IDENTYCZNEGO rozkładu binarnego jak przy gaussie (mniejsze/większe od sigma[=1]) jeżeli start generatora będzie ten sam - do tego dochodzi jeszcze poniżej zakomentowanie kodu sprawdzającego udział 50/50
                if (MTRandom0to1(randomStartStep)<0.5) particles[i].type=0; else particles[i].type=1;
                particles[i].diameter=getInitDiameterBinarySystem(particles[i].type);
            } else if (beta==1) { //rozkład Gaussa średnic dysków o wartości średniej \my=sigma[=1] i odchyleniu standardowym \delta=deltaDiameter
                particles[i].diameter=getInitDiameterByDistributionFunction(gaussianFunction,1.0-6*deltaDiameter,1.0+6*deltaDiameter);
                particles[i].type=particles[i].diameter<1?1:0;
            } else if (beta==2) { //rozkład Cossinusowy z MAX w sigma[=1] i okresie 2*Pi*deltaDiameter
                particles[i].diameter=getInitDiameterByDistributionFunction(cosinusFunction,1.0-M_PI/2*deltaDiameter,1.0+M_PI/2*deltaDiameter);
                particles[i].type=particles[i].diameter<1?1:0;
            }
        } /*else {
            particles[i].type=pTC[(columnCounter+rowCounter/typeN[1]*horizontalShift)%typeN[0]][rowCounter%typeN[1]];
            particles[i].diameter=getInitDiameterBinarySystem(particles[i].type);
        }*/

        columnCounter++;
        if (columnCounter*1.000001>=matrixOfParticlesSize[0]*mod) {
            columnCounter=0;
            rowCounter++;
            if (rowCounter*1.000001>=matrixOfParticlesSize[1]*mod) {
                rowCounter=0;
                layerCounter++;
            }
        }
    }

    if (initMode==0) {
        if (beta==0) { //doprowadzenie układu binarnego BLACK/WHITE do udziału po 50% każdego typu
            int blackCounter=0;
            for (int i=0;i<N;i++) if (particles[i].type==0) blackCounter++;
            while (blackCounter>N/2) {
                int randIndex = (int)(MTRandom0to1(randomStartStep)*N);
                if (particles[randIndex].type==0) {particles[randIndex].type=1; particles[randIndex].diameter=getInitDiameterBinarySystem(particles[randIndex].type); blackCounter--;}
            }
            while (blackCounter<N/2) {
                int randIndex = (int)(MTRandom0to1(randomStartStep)*N);
                if (particles[randIndex].type==1) {particles[randIndex].type=0; particles[randIndex].diameter=getInitDiameterBinarySystem(particles[randIndex].type); blackCounter++;}
            }
        } else if (beta==1) { //doprowadzenie do akceptowalnych rozkładów układów o skończonych rozmiarach (N<\infty, przy \infty <\diameter>->1). Testy z publikacji o sferach polidyspersyjnych.
            int attempt=1;
            printf("Tests of Gaussian polydisperse distribution... "); fflush(stdout);
            double diameterSum=0, diameter2Sum=0; for (int i=0;i<N;i++) {diameterSum+=particles[i].diameter; diameter2Sum+=particles[i].diameter*particles[i].diameter;}
            while (fabs(diameterSum/(double)N-1)>0.0000001 ||  //do 0.002 stosowane było 0.00001, po aktualizacji jest 0.0000001 dla 0.001 i o rząd mniej dla 0.0001 itd. (przy delta=1e-6 już bywało bardzo dużo prób)
                   fabs(sqrt(diameter2Sum/(double)N-diameterSum*diameterSum/(double)N/(double)N)/diameterSum*N-deltaDiameter)>0.00000001) { //taktyka 'dummy', ale działa w sensownym czasie NAWET dla deltaDiameter=0.05 przy N=32, więc po co szukać efektywniejszej;  do 0.002 stosowane było 0.000001, po aktualizacji jest 0.00000001 dla 0.001 i o rząd mniej dla 0.0001 itd.
                attempt++; int randIndex = (int)(MTRandom0to1(randomStartStep)*N);
                diameterSum-=particles[randIndex].diameter; diameter2Sum-=particles[randIndex].diameter*particles[randIndex].diameter;
                particles[randIndex].diameter=getInitDiameterByDistributionFunction(gaussianFunction,1.0-6*deltaDiameter,1.0+6*deltaDiameter);
                particles[randIndex].type=particles[randIndex].diameter<1?1:0;
                diameterSum+=particles[randIndex].diameter; diameter2Sum+=particles[randIndex].diameter*particles[randIndex].diameter;
            }
            printf("done (%d attempts)\n",attempt); fflush(stdout);
        }
    }

    //for (int i=0;i<N;i++) particles[i].diameter=1+(particles[i].diameter-1)/pow(10,ASD); //skalowanie rozmiarow czastek - uruchamialem program INITUJAC z tego samego miejsca generatora jak uklady PRZESZLE, przy simulationStage=2; zmienna ASD byla dodana do 'uruchamiacza' 2 na koncu, po kazdym wywolaniu programu byla linijka zmieniajaca nazwe folderu (delta)

    if (gaps>0) return createRandomGaps(particles,boxMatrix,boxParameters[1]);
    else return 1;
}

int getEnergyAll (particle *particles, double boxMatrix[3][3]) {
    int energy=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        double drRoot=sqrt(dr[3]);
        if (drRoot<0.5*(particles[i].diameter+particles[j].diameter)) {
            energy=1;
            j=activeN; i=activeN; break;
        }
    }
    return energy;
}

int getEnergy (particle *particles, particle *dispPart, int index, double boxMatrix[3][3]) {
    int energy=0;
    for (int i=0;i<particles[index].neighCounter;i++) {
        getParticlesDistanceSquared(&particles[particles[index].neighbours[i]],dispPart,boxMatrix);
        double drRoot=sqrt(dr[3]);
        if (drRoot<0.5*(particles[index].diameter+particles[particles[index].neighbours[i]].diameter)) {
            energy=1;
            i=particles[index].neighCounter; break;
        }
    }
    return energy;
}

int attemptToDisplaceAParticle (particle *particles, int index, double boxMatrix[3][3], double boxParameters[11]) {
    int result=1;
    particle displacedParticle;
    for (int i=0;i<3;i++) displacedParticle.r[i]=particles[index].r[i]+(MTRandom0to1(randomStartStep)-0.5)*deltaR;
    displacedParticle.normR[0]=-(displacedParticle.r[0]*boxParameters[2]+displacedParticle.r[1]*boxParameters[3]+displacedParticle.r[2]*boxParameters[4])/boxParameters[0];
    displacedParticle.normR[1]=(displacedParticle.r[0]*boxParameters[5]+displacedParticle.r[1]*boxParameters[6]+displacedParticle.r[2]*boxParameters[7])/boxParameters[0];
    displacedParticle.normR[2]=-(displacedParticle.r[0]*boxParameters[8]+displacedParticle.r[1]*boxParameters[9]+displacedParticle.r[2]*boxParameters[10])/boxParameters[0];

    double newEnPotVicinity=getEnergy(particles,&displacedParticle,index,boxMatrix);

    if (newEnPotVicinity==0) {
        for (int i=0;i<3;i++) {
            particles[index].r[i]=displacedParticle.r[i];
            particles[index].normR[i]=displacedParticle.normR[i];
        }
        checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix);
    } else result=0;
    return result;
}

void cloneParticlesForSpecificBoxMatrix (particle *clonedParticles, particle *particles, double boxMatrix[3][3]) {
    for (int i=0;i<activeN;i++) {
        for (int j=0;j<3;j++) clonedParticles[i].normR[j]=particles[i].normR[j];
        for (int j=0;j<3;j++) clonedParticles[i].r[j]=boxMatrix[j][0]*particles[i].normR[0]+boxMatrix[j][1]*particles[i].normR[1]+boxMatrix[j][2]*particles[i].normR[2];
    }
}

int attemptToChangeVolumeSeparate (particle *particles, double pressure, double boxMatrix[3][3], double boxParameters[11], double deltaV, bool moveType) {
    int result=1;
    double newBoxMatrix[3][3], newBoxMatrixInnerParameters[11], NFactor;
    if (boxParameters[1]/VcpPerParticle/N<pressureRealOfNotFluid) {
        if (moveType) { //logarytmiczny ruch elementów diagonalnych
            NFactor=N+1;
            for (int i=0;i<3;i++) newBoxMatrix[i][i]=exp(log(boxMatrix[i][i])+(MTRandom0to1(randomStartStep)-0.5)*deltaV);
            newBoxMatrix[0][1]=boxMatrix[0][1];
            newBoxMatrix[0][2]=boxMatrix[0][2];
            newBoxMatrix[1][2]=boxMatrix[1][2];
        } else {        //liniowy ruch elementów niediagonalnych
            NFactor=N;
            for (int i=0;i<3;i++) newBoxMatrix[i][i]=boxMatrix[i][i];
            newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
            newBoxMatrix[0][2]=boxMatrix[0][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
            newBoxMatrix[1][2]=boxMatrix[1][2]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        }
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        NFactor=N+1;
        newBoxMatrix[0][0]=exp(log(boxMatrix[0][0])+(MTRandom0to1(randomStartStep)-0.5)*deltaV);
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[2][2]=boxMatrix[2][2]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
        newBoxMatrix[0][2]=boxMatrix[0][2]*modifier;
        newBoxMatrix[1][2]=boxMatrix[1][2]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1]; newBoxMatrix[2][0]=newBoxMatrix[0][2]; newBoxMatrix[2][1]=newBoxMatrix[1][2];
    updateBoxMatrixParameters(newBoxMatrix,newBoxMatrixInnerParameters);

    particle *particlesInNewBox = new particle[activeN];  //dla dużych N (testowane dla N=16384) stack nie wystarcza, trzeba użyć heapu (stack zostaje dla podstawowej tablicy czastek, jest ciut szybszy)
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,newBoxMatrix);
    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            getParticlesDistanceSquared(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]],newBoxMatrix);
            double drRoot=sqrt(dr[3]);
            if (drRoot<0.5*(particles[i].diameter+particles[particles[i].neighbours[j]].diameter)) {
                result=0;
                i=activeN; j=particles[i].neighCounter; break;
            }
        }
    }
    if (result) {
        double arg=-(pressure*(newBoxMatrixInnerParameters[1]-boxParameters[1])-(NFactor*log(newBoxMatrixInnerParameters[1]/boxParameters[1])+log((newBoxMatrix[0][0]+newBoxMatrix[1][1]+newBoxMatrix[2][2])/(boxMatrix[0][0]+boxMatrix[1][1]+boxMatrix[2][2]))));
        if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
        else {
            for (int i=0;i<11;i++) boxParameters[i]=newBoxMatrixInnerParameters[i];
            for (int i=0;i<3;i++) for (int j=0;j<3;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) for (int j=0;j<3;j++) particles[i].r[j]=particlesInNewBox[i].r[j];
        }
    }
    delete [] particlesInNewBox; particlesInNewBox=NULL;
    return result;
}
/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber][1]=strtod(pEnd,&pEnd);
        iterationTable[fileIterateIterationsNumber++][2]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

void addAppendix (char *fileName, char *JOBID, bool jobIdOn) {
    strcpy(buffer,"3D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_D-",5); strncat(buffer,bufferD,20);
    strncat(buffer,"_inM-",6); strncat(buffer,bufferInitMode,5);
    strncat(buffer,"_a-",5); strncat(buffer,bufferAlpha,20);
    strncat(buffer,"_b-",5); strncat(buffer,bufferBeta,20);

    mkdir(buffer,S_IRWXU);
    strncat(buffer,"/",2);
    if (jobIdOn) {
        strncat(buffer,JOBID,50);
        strncat(buffer,"_",2);
    }
    strncat(buffer,fileName,200);
    strcpy(fileName,buffer);
}

void getNextArgument (double prevArg[2], bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            for (int i=0;i<2;i++) prevArg[i]=growing?iterationTable[actIteration][i+1]:iterationTable[fileIterateIterationsNumber-1-actIteration][i+1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg[0]*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>=round(intervalMin[i]*10000) && round(prevArg[0]*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)+round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg[0]/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>round(intervalMin[i]*10000) && round(prevArg[0]*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)-round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)<round(minArg*10000)) growing=-1;
    }
}

bool isLineCorrect(char linia[4096]) {
    sscanf(linia,"%c",linia);
    int actIndex=0, dataIndex=0; while (dataIndex<6) {
        char data[50]="";
        int licznik=0, dotCounter=0;
        while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
        if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;} //test of single dot in a number
        actIndex++; dataIndex++;
        if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10; //test of dot position after first digit
    } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10; //test of max 6 numbers in a row

    return (dataIndex<10);
}

double getAvErrorFromSumEps (double sum, double denominator) {
    return sqrt(sum/denominator);
}

void updateTableAndGetActualMean (double table[100], double & mean, int const & changeIndex, double const & changeValue) {
    mean-=table[changeIndex]*0.01; table[changeIndex]=changeValue; mean+=changeValue*0.01;
}

class doubleSum {//klasa do sumowania zmiennych z wysoka precyzja, na potrzeby koncowych obliczen wartosci srednich
    public:
    double sum, c;
    doubleSum(){sum=0;c=0;}
    doubleSum(const double & a) {sum=a;c=0;}

    doubleSum & operator+=(const double & a) {//Kahan summation algorithm
        double y=a-c, t=sum+y;
        c=(t-sum)-y;
        sum=t;
        return *this;
    }
    doubleSum & operator+=(const doubleSum & a) {//nie wiem czy poprawne, ale ewentualny wpływ i tak powinien byc znikomy
        double y=a.sum-(c+a.c), t=sum+y;
        c=(t-sum)-y;
        sum=t;
        return *this;
    }
    doubleSum & operator*=(const double & a) {
        sum*=a; c*=a;
        return *this;
    }
    doubleSum & operator/=(const double & a) {
        sum/=a; c/=a;
        return *this;
    }

    operator double() const {return sum;}
};



int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[800];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,800,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20]="";
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++]; data[licznik]=' ';
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:initMode=strtol(data,NULL,10);break;
                case 4:deltaDiameter=strtod(data,NULL);break;
                case 5:alpha=strtol(data,NULL,10);break;
                case 6:beta=strtol(data,NULL,10);break;
                case 7:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 8:growing=strtol(data,NULL,10);break;
                case 9:loadedConfiguration=strtol(data,NULL,10);break;
                case 10:loadedArg=strtod(data,NULL);break;
                case 11:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 12:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 13:loadedSetGenerator=strtol(data,NULL,10);break;
                case 14:iterationsNumber=strtol(data,NULL,10);break;
                case 15:intervalSampling=strtol(data,NULL,10);break;
                case 16:intervalOutput=strtol(data,NULL,10);break;
                case 17:saveConfigurations=strtol(data,NULL,10);break;
                case 18:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 19:neighUpdatingFrequency=strtol(data,NULL,10);break;
                case 20:skipFirstIteration=strtol(data,NULL,10);break;
                case 21:useSpecificDirectory=strtol(data,NULL,10);break;
                case 22:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 23:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 24:intervalResults=strtol(data,NULL,10);break;
                case 25:maxDeltaR=strtod(data,NULL);break;
                case 26:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 27:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 28:useFileToIterate=strtol(data,NULL,10);break;
                case 29:startMinPacFrac=strtod(data,NULL);break;
                case 30:startMaxPacFrac=strtod(data,NULL);break;
                case 31:minArg=strtod(data,NULL);break;
                case 32:maxArg=strtod(data,NULL);break;
                case 33:multiplyArgument=strtol(data,NULL,10);break;
                case 34:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-35)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
        }
        fclose(fileConfig);
    } while (testValue!=12345);

    //zlecanie parametrow z poziomu wiersza polecen:
    char JOBID[50]="j-"; int pointNumber=0;
    if (argumentsNumber==1) {
        strncat(JOBID,"none",50);
        if (useFileToIterate) if(createIterationTable()) return 0;
    } else {
        int correctNumberOfArguments=1;
        switch (strtol(arguments[1],NULL,10)) {
            case 0: //ustaw JOBID
                if (argumentsNumber==3) {
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (useFileToIterate) if(createIterationTable()) return 0;
                } else correctNumberOfArguments=0; break;
            case 1: //ustaw JOBID, singleRun dla parametrow zadanych bezposrednio
                if (argumentsNumber==14) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (growing) {
                        startMinPacFrac=strtod(arguments[3],NULL); minArg=strtod(arguments[4],NULL);
                    } else {
                        startMaxPacFrac=strtod(arguments[3],NULL); maxArg=strtod(arguments[4],NULL);
                    }
                    N=strtol(arguments[5],NULL,10);
                    gaps=strtol(arguments[6],NULL,10);
                    deltaDiameter=strtod(arguments[7],NULL);
                    initMode=strtol(arguments[8],NULL,10);
                    growing=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    alpha=strtol(arguments[12],NULL,10);
                    beta=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    deltaDiameter=strtod(arguments[5],NULL);
                    initMode=strtol(arguments[6],NULL,10);
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    pointNumber=strtol(arguments[11],NULL,10);
                    generatorStartPoint=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[9],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=0;
                    alpha=strtol(arguments[12],NULL,10);
                    beta=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, deltaDiameter, initMode, growing, iterationsNumber, useSpecificDirectory, alpha, beta\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, deltaDiameter, initMode, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, alpha, beta\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, deltaDiameter, initMode, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, alpha, beta\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, deltaDiameter, initMode, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, alpha, beta\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, deltaDiameter, initMode, growing, pointNumber, iterationsNumber, useSpecificDirectory, alpha, beta\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[0][i+1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[fileIterateIterationsNumber-1][i+1];
        }
    } else {
        startArg[0]=growing?minArg:maxArg;
        startArg[1]=0;
    }

    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];
    activeN=N-gaps;
    if (fabs(cbrt(N/4.0)-floor(cbrt(N/4.0)))>0.000001) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    }

    //stale wynikajace z zadanych parametrow
    deltaR=maxDeltaR; //jednostka odległości to sigma=1

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps); sprintf(bufferD,"%.3E",deltaDiameter);
    sprintf(bufferInitMode,"%d",initMode); sprintf(bufferAlpha,"%d",alpha); sprintf(bufferBeta,"%d",beta);

    int folderIndex=useSpecificDirectory, checkNext;
    char bufferCheckFolderExisting[200];
    FILE *checkFolderExisting;
    if (!folderIndex) do {
        sprintf(bufferFolderIndex,"%d",++folderIndex);
        strcpy(bufferCheckFolderExisting,resultsFileName);
        addAppendix(bufferCheckFolderExisting,JOBID,false);
        checkFolderExisting = fopen(bufferCheckFolderExisting,"rt");
        if (checkFolderExisting!=NULL) {
            fclose(checkFolderExisting);
            checkNext=1;
        } else checkNext=0;
    } while (checkNext);
    sprintf(bufferFolderIndex,"%d",folderIndex);
    addAppendix(resultsFileName,JOBID,false);
    addAppendix(excelResultsFileName,JOBID,false);
    addAppendix(configurationsFileName,JOBID,true);
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.4E",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(probDensDistFunFileName,JOBID,true); addAppendix(probDensDistFunResultsFileName,JOBID,false);  //probDensDistFunMode

    particle particles[N];
    double args[14];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileAllResults, *fileProbDensDistFun, *fileProbDensDistFunResults;    //probDensDistFunMode
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        fprintf(fileResults,"cycle\tpressure\trho*\tdRho*\tv*\tdV*\tvolume\tdVolume\tH00\tdH00\tH11\tdH11\tH22\tdH22\tH01\tdH01\tH02\tdH02\tH12\tdH12\tB\tdB\tavMy\tdAvMy\tmy1\tdMy1\tmy2\tdMy2\tE\tdE\tnu_100_all\tdNu_100_all\tnu_111_all\tdNu_111_all\tnu_110_1m10\tdNu_110_1m10\tnu_110_001\tdNu_110_001\tS11c\tdS11c\tS12c\tdS12c\tS44c\tdS44c\tC11c\tdC11c\tC12c\tdC12c\tC44c\tdC44c\tS11\tdS11\tS12\tdS12\tS13\tdS13\tS14\tdS14\tS15\tdS15\tS16\tdS16\tS22\tdS22\tS23\tdS23\tS24\tdS24\tS25\tdS25\tS26\tdS26\tS33\tdS33\tS34\tdS34\tS35\tdS35\tS36\tdS36\tS44\tdS44\tS45\tdS45\tS46\tdS46\tS55\tdS55\tS56\tdS56\tS66\tdS66\tB11>0\tB44>0\t-1/2<=B12/B11<=1\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        fprintf(fileExcelResults,"pressure\tv*\tB\tdB\tavMy\tdAvMy\tE\tdE\tnu_100_all\tdNu_100_all\tnu_111_all\tdNu_111_all\tnu_110_1m10\tdNu_110_1m10\tnu_110_001\tdNu_110_001\tS11c\tdS11c\tS12c\tdS12c\tS44c\tdS44c\tC11c\tdC11c\tC12c\tdC12c\tC44c\tdC44c\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[3][3];
    while (growing>=0) {
        double pressureReduced=arg[0], pressure=pressureReduced,            //pressureReduced=\tau*\sigma^3/kT, kT=1[hard] - Dwa punkty widzenia: 1) zmniejszenie sigma ZWIĘKSZA JEDNOSTKE p^*, zatem ten sam STAN FIZYCZNY jak przy sigma=1 bedzie przy mniejszym p^*. pReal to tak naprawde pReducedAdjusted. Inny, równoważny punkt widzenia, to 2) pReal redukuje objetosc, ktora NIE jest wyrazana w jednostkach sigma. Objetosc jest obliczana z boxMatrix, ktory jest inicjowany z czynnikiem *sigma, zatem MA jednostkę, a NIE jest zredukowany. Przeciez gdyby sigma=2, to boxMatrix bylby 2x wiekszy, a 'w jednostkach sigma' (zredukowany) powinien pozostac identyczny
               boxMatrix[3][3],matrixOfParticlesSize[3],unitCellAtCP[3],   //obydwa sprowadzają się do tego, że przy liczeniu prawdopodobieństwa ma być jednostka zredukowana (bezwymiarowa): 1) zakłada, że volume jest zredukowane, więc dostosowuje pReduced do stanu fizycznego; 2) zakłada, że pReduced już jest OK (w końcu jest zredukowane), tylko po prostu objętość NIE jest zredukowana, i trzeba ją zredukować dzieląc przez sigma^3
               matrixCellXYZ[6][6][6][3];
        int n[3]; //n[X/Y/Z], matrixCell[n[XMax]][n[YMax]][n[ZMax]][x/y/z], zatem: n[X/Y/Z](max)=6
        unitCellAtCP[0]=unitCellAtCP[1]=unitCellAtCP[2]=sqrt(2);  //struktura fcc dla kul, jednostka odległości to sigma[=1]
        n[0]=1; n[1]=2; n[2]=2;
        matrixCellXYZ[0][0][0][0]=0; matrixCellXYZ[0][0][0][1]=0; matrixCellXYZ[0][0][0][2]=0;
        matrixCellXYZ[0][1][0][0]=unitCellAtCP[0]/2.0; matrixCellXYZ[0][1][0][1]=unitCellAtCP[1]/2.0; matrixCellXYZ[0][1][0][2]=0;
        matrixCellXYZ[0][0][1][0]=unitCellAtCP[0]/2.0; matrixCellXYZ[0][0][1][1]=0; matrixCellXYZ[0][0][1][2]=unitCellAtCP[2]/2.0;
        matrixCellXYZ[0][1][1][0]=unitCellAtCP[0]; matrixCellXYZ[0][1][1][1]=unitCellAtCP[1]/2.0; matrixCellXYZ[0][1][1][2]=unitCellAtCP[2]/2.0;
        VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]*unitCellAtCP[2]/(double)n[0]/(double)n[1]/(double)n[2];

        matrixOfParticlesSize[0]=1; matrixOfParticlesSize[1]=2; matrixOfParticlesSize[2]=2;
        double NLinearMod = cbrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]/matrixOfParticlesSize[2]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<3;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?cbrt(startMinPacFrac):cbrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0; boxMatrix[2][0]=0.0; boxMatrix[0][2]=0.0; boxMatrix[1][2]=0.0; boxMatrix[2][1]=0.0;
        } else {
            for (int i=0;i<3;i++) for (int j=0;j<3;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        }
        updateBoxMatrixParameters(boxMatrix,boxMatrixInnerParameters);
        double rho=N/boxMatrixInnerParameters[1], pacFrac=1.0/VcpPerParticle/rho;

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E, alpha: %d, beta: %d), deltaDiameter: %.3E\n",N,gaps,growing,startArg[0],rho,pacFrac,alpha,beta,deltaDiameter);
                if (!initPositions(particles,boxMatrix,boxMatrixInnerParameters,matrixOfParticlesSize,n,matrixCellXYZ,pacFrac)) return 0;
                if (!updateNeighbourList(particles,boxMatrix,boxMatrixInnerParameters[1],1.8/*v*=1.45: 1.8-18neigh, 1.95-42neigh*/,false,0)) return 0;  //checkCorrectNeighbourNumberWhileUpdating - if desired type (...,1.3,true,40) else (...,specificMod(i.e.: 1.4),false,0) 1/4
                printf("Checking overlaps in inited configuration... "); fflush(stdout);
                if (getEnergyAll(particles,boxMatrix)==1) {
                    printf("Configuration's init parameters generate overlap(s) [energy=1].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<15) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<14) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==14) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[2][2]=args[7];
                                boxMatrix[1][0]=args[8]; boxMatrix[0][1]=args[8];
                                boxMatrix[2][0]=args[9]; boxMatrix[0][2]=args[9];
                                boxMatrix[1][2]=args[10]; boxMatrix[2][1]=args[10];
                                deltaR=args[11]; deltaV=args[12]; deltaVND=args[13];
                                /*arg[0]=args[3];*/ arg[0]=loadedArg; pressureReduced=arg[0]; pressure=pressureReduced;

                                updateBoxMatrixParameters(boxMatrix,boxMatrixInnerParameters);
                                rho=N/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='b' || character=='w') {pIndex++; particles[pIndex].type=(character=='b')?0:1; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<4;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',' && character!=']') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<3) particles[pIndex].r[i]=strtod(data,NULL);
                                else particles[pIndex].diameter=strtod(data,NULL);
                            }
                            particles[pIndex].normR[0]=-(particles[pIndex].r[0]*boxMatrixInnerParameters[2]+particles[pIndex].r[1]*boxMatrixInnerParameters[3]+particles[pIndex].r[2]*boxMatrixInnerParameters[4])/boxMatrixInnerParameters[0];
                            particles[pIndex].normR[1]=(particles[pIndex].r[0]*boxMatrixInnerParameters[5]+particles[pIndex].r[1]*boxMatrixInnerParameters[6]+particles[pIndex].r[2]*boxMatrixInnerParameters[7])/boxMatrixInnerParameters[0];
                            particles[pIndex].normR[2]=-(particles[pIndex].r[0]*boxMatrixInnerParameters[8]+particles[pIndex].r[1]*boxMatrixInnerParameters[9]+particles[pIndex].r[2]*boxMatrixInnerParameters[10])/boxMatrixInnerParameters[0];
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'b' or 'w'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E, alpha: %d, beta: %d), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaV: %.4E, DeltaVND: %.4E\n",N,gaps,growing,pressureReduced,rho,pacFrac,alpha,beta,args[0],args[1],(long)args[4],deltaR,deltaV,deltaVND);
                //for (int i=0;i<3;i++) for (int j=0;j<3;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d: %.17E,  %.17E,  %.17E,  %.17E,  %d\n",i,particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].diameter,particles[i].type);return 0;
                if (!updateNeighbourList(particles,boxMatrix,boxMatrixInnerParameters[1],1.3,true,40)) return 0;  //checkCorrectNeighbourNumberWhileUpdating - if desired type (...,1.3,true,40) else (...,specificMod(i.e.: 1.4),false,0) 2/4
                printf("Checking overlaps in loaded file... "); fflush(stdout);
                if (getEnergyAll(particles,boxMatrix)==1) {
                    printf("Configuration from loaded file contains overlap(s) [energy=1].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (!onlyMath[0]) {
                if (loadedConfiguration) {
                    if (loadedSetStartGenerator) {
                        printf("Setting start position of p-random number generator to position from file...\n");
                        InitMT((unsigned int)args[0]);
                        randomStartStep[0]=args[0];
                    } else {
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                        if (generatorStartPoint==0) {
                            generatorStartPoint=time(0);
                            printf("Setting start position of p-random number generator to actual CPU time...\n");
                        } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                        InitMT((unsigned int)generatorStartPoint);
                        randomStartStep[0]=generatorStartPoint;
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<args[1];i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {
                    if (generatorStartPoint==0) {
                        generatorStartPoint=time(0);
                        printf("Setting start position of p-random number generator to actual CPU time...\n");
                    } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                    InitMT((unsigned int)generatorStartPoint);
                    randomStartStep[0]=generatorStartPoint;
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)... (%ld cycles)\n",pressureReduced,rho,pacFrac,cyclesOfEquilibration);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E, deltaDiameter: %.3E\n",N,gaps,growing,pressureReduced,deltaDiameter);
            fflush(stdout);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,fullCycleDiagonalMove=activeN+volumeMoveChance/2,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0,
                attemptedNumberVND=0, displacedNumberVND=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0;
            double deltaRTable[100], deltaRMean=deltaR, deltaVTable[100], deltaVMean=deltaV, deltaVNDTable[100], deltaVNDMean=deltaVND;
            for (int i=0;i<100;i++) {deltaRTable[i]=deltaRMean; deltaVTable[i]=deltaVMean; deltaVNDTable[i]=deltaVNDMean;}
            int simulationStage=cyclesOfEquilibration>0?0:cyclesOfMeasurementBuffer>0?1:2;  //0-equilibration, 1-measurement, 2-end
            int cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod;  //TODO: indexScanned w 3D inaczej powinien być ustawiany
            bool volumeMove=false;

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferProbDensDistFun[200],bufferProbDensDistFunResults[200],bufferPressure[100];    //probDensDistFunMode
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferProbDensDistFun,probDensDistFunFileName); strcpy(bufferProbDensDistFunResults,probDensDistFunResultsFileName); //probDensDistFunMode
            sprintf(bufferPressure,"%.4E",pressureReduced);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPressure,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPressure,100); strcpy(bufferSavedConfigurations,bufferConfigurations); strncat(bufferConfigurations,".txt",5); strncat(bufferSavedConfigurations,"_transient.txt",15);
            strncat(bufferProbDensDistFun,"_arg-",6); strncat(bufferProbDensDistFun,bufferPressure,100); strncat(bufferProbDensDistFun,".txt",5); //probDensDistFunMode
            strncat(bufferProbDensDistFunResults,"_arg-",6); strncat(bufferProbDensDistFunResults,bufferPressure,100);    //probDensDistFunMode
            strncat(bufferProbDensDistFunResults,"_N-",4); sprintf(bufferPressure,"%d",N); strncat(bufferProbDensDistFunResults,bufferPressure,100);
            strncat(bufferProbDensDistFunResults,"_inM-",6); sprintf(bufferPressure,"%d",initMode); strncat(bufferProbDensDistFunResults,bufferPressure,100);
            strncat(bufferProbDensDistFunResults,".txt",5);

            fileAllResults = fopen(allResultsFileName,"a");
            if (saveConfigurations) fileSavedConfigurations = fopen(bufferSavedConfigurations,"a");
            //fileProbDensDistFun = fopen(bufferProbDensDistFun,"a");  //probDensDistFunMode - comment if not desired 1/7
            //char confsStepByStepFileName[200]="confsStepByStep.txt"; addAppendix(confsStepByStepFileName,JOBID,false); FILE *fileConfsStepByStep = fopen(confsStepByStepFileName,"a");  //consfStepByStep - comment if not desired 1/3
            if (onlyMath[0]) simulationStage=2;

            timeStart=time(0);
            while (simulationStage<2) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTRandom0to1(randomStartStep)*activeN);
                    volumeMove=false;
                } else randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycle);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(particles,randIndex,boxMatrix,boxMatrixInnerParameters))
                        displacedNumberR++;
                } else {
                    volumeMove=true;
                    if (randIndex<fullCycleDiagonalMove) {
                        attemptedNumberV++;
                        if (attemptToChangeVolumeSeparate(particles,pressure,boxMatrix,boxMatrixInnerParameters,deltaV,true))
                            displacedNumberV++;
                    } else {
                        attemptedNumberVND++;
                        if (attemptToChangeVolumeSeparate(particles,pressure,boxMatrix,boxMatrixInnerParameters,deltaVND,false))
                            displacedNumberVND++;
                    }
                }

                /*bool saveConf=false;    //consfStepByStep - comment if not desired 2/3   -TODO: NIE(dokońca)dostosowane do 3D
                if (cycle<10 || (cycleCounter+1)%fullCycle==0) saveConf=true;
                if (saveConf) {
                    fprintf(fileConfsStepByStep,"multimers[x_,y_,z_,kI_]:={{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},Opacity[If[x==0 && y==0,1,0.3]]",
                            boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
                    for (int i=0;i<activeN;i++) {
                        fprintf(fileConfsStepByStep,",%c[%.12E+x,%.12E+y,%.12E+z,%.17E]",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].diameter);
                    }
                    fprintf(fileConfsStepByStep,"};\nAppendTo[VvsSTEP,{%ld,%.12E}];configurationsList=Append[configurationsList,g[%ld,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                            (long)iStep,fabs(boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1]),(long)iStep,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
                }*/

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++; //printf("cykl: %d, time elapsed: %dsec\n",cycle,time(0)-timeStart);   //printActualCycle

                    // FREQUENT UPDATES OF NEIGHBOUR LIST DURING EQUILIBRATION - na potrzeby silnego ściskania z relatywnie mało-gęstych konfiguracji początkowych
                    //if (cycle%2==0) if (!updateNeighbourList(particles,boxMatrix,boxMatrixInnerParameters[1],1.4,false,0)) return 0;  //checkCorrectNeighbourNumberWhileUpdating - if desired type (...,1.3,true,40) else (...,specificMod(i.e.: 1.4),false,0) 3/4

                    if (cycle%intervalSampling==0) {
                        if (simulationStage==0 && cycle>cyclesOfEquilibration) {
                            simulationStage=1;
                            printf("Equilibration finished after: %ld cycles (%ldsec).\n",cyclesOfEquilibration,time(0)-timeStart);
                            fflush(stdout);
                        }
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV,
                               acceptanceRatioVND = displacedNumberVND/(double)attemptedNumberVND;
                        if (cycle%neighUpdatingFrequency==0) if (!updateNeighbourList(particles,boxMatrix,boxMatrixInnerParameters[1],1.4,false,0)) return 0;  //checkCorrectNeighbourNumberWhileUpdating - if desired type (...,1.3,true,40) else (...,specificMod(i.e.: 1.4),false,0) 4/4

                        if (simulationStage==1) {
                            if (timeEquilibration==0) {
                                timeEquilibration=time(0);

                                printf("Checking overlaps in equilibrated configuration... "); fflush(stdout);
                                if (getEnergyAll(particles,boxMatrix)==1) {
                                    printf("Equilibrated configuration contains overlap(s) [energy=1]. ");
                                    char allResultsErrorFileName[200];
                                    strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                                    if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err).\n");
                                    else printf("Error renaming results file (.err).\n");
                                    return 0;
                                } else  {printf("done\n"); fflush(stdout);}
                            }

                            rho=N/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0) {
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);

                                //Probability Density Distribution Function
                                //probDensDistFunMode - comment if not desired 2/7
                                /*for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                    if (i<particles[i].neighbours[j]) {
                                        getParticlesDistanceSquared(&particles[i],&particles[particles[i].neighbours[j]],boxMatrix);
                                        double drRoot=sqrt(dr[3]);

                                        fprintf(fileProbDensDistFun,"{%d,%d,%.17E},",i,particles[i].neighbours[j],drRoot-0.5*(particles[i].diameter+particles[particles[i].neighbours[j]].diameter));
                                    }
                                }*/
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++) fprintf(fileSavedConfigurations,"%c[%.17E,%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].diameter);
                                fprintf(fileSavedConfigurations,"%c[%.17E,%.17E,%.17E,%.17E]}",particles[activeNMinus1].type==0?'b':'w',particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].r[2],particles[activeNMinus1].diameter);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                printf("Cycle: %ld, ",(cycle+(long)args[4]));
                                printf("simulation time: full-%ldsec, measurement-%ldsec\n",time(0)-timeStart,time(0)-timeEquilibration);
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E, AccRatVND: %.4E, dVND: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV,acceptanceRatioVND,deltaVND);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box22: %.8E\n   box01(10): %.8E, box02(20): %.8E, box12(21): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                fflush(stdout);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2],deltaR,deltaV,deltaVND);
                                for (int i=0;i<activeN;i++) fprintf(fileConfigurations,"%c[%.17E,%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].diameter);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Parallelepiped[{0,0,0},{{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E}}]},{Opacity[0.2],Green,Sphere[{%.12E,%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][0],boxMatrix[1][1],boxMatrix[1][2],boxMatrix[2][0],boxMatrix[2][1],boxMatrix[2][2],particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[2][2]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E, boxMatrix[0][2]=boxMatrix[2][0]=%.12E, boxMatrix[1][2]=boxMatrix[2][1]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
                                fclose(fileConfigurations);
                            }
                        } //else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) {
                                deltaR*=1.05;
                                if (deltaR>maxDeltaR) deltaR=maxDeltaR;  //jednostka odległości to 1; ograniczenie od góry, żeby cząstka nie skakała zbyt mocno po pudle, nie teleportowała się
                            } else {
                                deltaR*=0.95;
                                while (boxMatrix[0][0]+deltaR*0.1==boxMatrix[0][0]) deltaR*=1.05;  //dolna granica deltaR - gdy jest tak małe (w stosunku do boxMatrix[0][0]), że go nie może zmienić ze względu na ograniczoną precyzję (czyli, nie może zmienić pozycji cząstki znajdujące1j się przy granicy pudła), należy ją zwiększać tak długo, aż będzie mogła wprowadzać takie zmiany. Bez uwzględnienia dolnej granicy delt, po jej przekroczeniu, dana delta przestaje być wrażliwa na dostosowywanie (czy się ją delikatnie zwiększy czy zmniejszy, nie zmieni to poziomu akceptacji). Acceptance ruchów będzie wówczas zazwyczaj dość mały (nie 50%, znacznie mniejszy, mniejszy niż 20% - jeżeli byłby duży, to nigdy by do takiej sytuacji nie doszło - bo przecież delta byłaby zwiększana a nie zmniejszana) i będzie dochodziło do dalszego zmniejszania danej delty, aż do -INF. Mnożenie *0.1: 0.5(zakres losowy od -0.5 do 0.5)*0.5(średnia wartość zmiennej losowej od 0 do 1)*0.4(tak dla bezpieczeństwa, żeby nie schodzić do granicznie małych wartości)
                            }
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; //ograniczenie deltaV od góry chyba nie jest potrzebne, w płynie kształt i tak jest blokowany, same ruchy objętościowe nie wprowadzą niczego złego, nawet jak będą bardzo duże (dalej są określane przez entalpię swobodną)
                            else {
                                deltaV*=0.95;
                                while (exp(log(boxMatrix[0][0])+deltaV*0.1)==boxMatrix[0][0]) deltaV*=1.05;     //ograniczenie od dołu jak dla deltaR (ruchy logarytmiczne elementów diagonalnych)
                            }
                            if (acceptanceRatioVND>desiredAcceptanceRatioV) deltaVND*=1.05;
                            else {
                                deltaVND*=0.95;
                                while (boxMatrix[0][1]+deltaVND*0.1==boxMatrix[0][1]) deltaVND*=1.05;   //ograniczenie od dołu jak dla deltaR, ale bazujące na przykładowym elemencie niediagonalnym [tu: 01] (ruchy liniowe elementów niediagonalnych)
                            }

                            int sampleNumberMod100=(cycle/intervalSampling)%100;
                            updateTableAndGetActualMean(deltaRTable,deltaRMean,sampleNumberMod100,deltaR); deltaR=deltaRMean;
                            updateTableAndGetActualMean(deltaVTable,deltaVMean,sampleNumberMod100,deltaV); deltaV=deltaVMean;
                            updateTableAndGetActualMean(deltaVNDTable,deltaVNDMean,sampleNumberMod100,deltaVND); deltaVND=deltaVNDMean;
                        //}
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                        attemptedNumberVND=0; displacedNumberVND=0;
                    }
                    if (simulationStage==1 && cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) simulationStage=2;
                }
            }
            fclose(fileAllResults);
            if (saveConfigurations) fclose(fileSavedConfigurations);
            //fclose(fileProbDensDistFun); //probDensDistFunMode - comment if not desired 3/7
            //fclose(fileConfsStepByStep);    //consfStepByStep - comment if not desired 3/3
            if (timeEquilibration==0) timeEquilibration=time(0);
            printf("Checking overlaps in final configuration... "); fflush(stdout);
            if (!onlyMath[0] && getEnergyAll(particles,boxMatrix)==1) {
                printf("Final configuration contains overlap(s) [energy=1]. ");
                char allResultsErrorFileName[200],configurationErrorFileName[200];
                strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                strcpy(configurationErrorFileName,bufferConfigurations); strncat(configurationErrorFileName,".err",5);
                if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err). ");
                else printf("Error renaming results file (.err). ");
                if (rename(bufferConfigurations,configurationErrorFileName)==0) printf("Configuration file successfully renamed (.err).\n");
                else printf("Error renaming configuration file (.err).\n");
                return 0;
            } else {printf("done\n"); fflush(stdout);}
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie liczby linii danych (potrzebne do podziału na zespoły i obliczenia średnich błędów)
            printf("Calculation of data lines... "); fflush(stdout);
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[4096]; double dataLicznik=0; int faultyLines=0, onlyMathLinesBuffer=onlyMath[1];
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) {
                if (fgets(linia,300,fileAllResults)!=NULL && !isLineCorrect(linia)) onlyMathLinesBuffer++;
            }
            while(fgets(linia,300,fileAllResults)!=NULL) {
                if (isLineCorrect(linia)) dataLicznik++; else faultyLines++;
            }
            fclose(fileAllResults);
            printf("done (Found %ld data lines [%d faulty lines occurred].",(long)dataLicznik,faultyLines);
            if ((long)dataLicznik%10>0) printf(" Last %ld lines won't be considered, due to calculations of averages in 10 sets.)\n",(long)dataLicznik%10); else printf(")\n");
            dataLicznik-=(long)dataLicznik%10;

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);
            doubleSum avVolumeSet[10], avBoxMatrixSet[10][6], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            /*for (int i=0;i<10;i++) { //doubleSum w konstruktorze inicjalizuje zmienne (0)
                avVolumeSet[i]=0; for (int j=0;j<6;j++) avBoxMatrixSet[i][j]=0;
                avRhoSet[i]=0; avPacFracSet[i]=0;
            }*/
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[6]; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10;
                    else dataD[dataIndex++]=strtod(data,NULL);
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<6;i++) avBoxMatrixSet[setIndex][i]+=dataD[i];
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            doubleSum avVolume/*=0*/, avBoxMatrix[6]/*={0,0,0,0,0,0}*/, avRho/*=0*/, avPacFrac/*=0*/; //doubleSum w konstruktorze inicjalizuje zmienne (0)
            for (int i=0;i<10;i++) {
                for (int j=0;j<6;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
                avVolumeSet[i]=fabs(avBoxMatrixSet[i][0]*avBoxMatrixSet[i][1]*avBoxMatrixSet[i][2]+avBoxMatrixSet[i][3]*avBoxMatrixSet[i][5]*avBoxMatrixSet[i][4]+
                                    avBoxMatrixSet[i][4]*avBoxMatrixSet[i][3]*avBoxMatrixSet[i][5]-avBoxMatrixSet[i][4]*avBoxMatrixSet[i][1]*avBoxMatrixSet[i][4]-
                                    avBoxMatrixSet[i][0]*avBoxMatrixSet[i][5]*avBoxMatrixSet[i][5]-avBoxMatrixSet[i][3]*avBoxMatrixSet[i][3]*avBoxMatrixSet[i][2]);
                avVolume+=avVolumeSet[i];
                avRhoSet[i]=N/avVolumeSet[i]; avRho+=avRhoSet[i];
                avPacFracSet[i]=1.0/VcpPerParticle/avRhoSet[i]; avPacFrac+=avPacFracSet[i];
            }
            avVolume*=0.1; avRho*=0.1; avPacFrac*=0.1; for (int i=0;i<6;i++) avBoxMatrix[i]*=0.1;
            //obliczenie bledow mierzonych wielkosci
            doubleSum dAvVolume/*=0*/, dAvBoxMatrix[6]/*={0,0,0,0,0,0}*/, dAvRho/*=0*/, dAvPacFrac/*=0*/; //doubleSum w konstruktorze inicjalizuje zmienne (0)
            for (int i=0;i<10;i++) {double epsilon=avVolume-avVolumeSet[i]; dAvVolume+=epsilon*epsilon;} dAvVolume=getAvErrorFromSumEps(dAvVolume,90.0); //10*9 (n(n-1))
            for (int j=0;j<6;j++) {for (int i=0;i<10;i++) {double epsilon=avBoxMatrix[j]-avBoxMatrixSet[i][j]; dAvBoxMatrix[j]+=epsilon*epsilon;} dAvBoxMatrix[j]=getAvErrorFromSumEps(dAvBoxMatrix[j],90.0);}
            for (int i=0;i<10;i++) {double epsilon=avRho-avRhoSet[i]; dAvRho+=epsilon*epsilon;} dAvRho=getAvErrorFromSumEps(dAvRho,90.0);
            for (int i=0;i<10;i++) {double epsilon=avPacFrac-avPacFracSet[i]; dAvPacFrac+=epsilon*epsilon;} dAvPacFrac=getAvErrorFromSumEps(dAvPacFrac,90.0);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            doubleSum e0000Set[10], e0001Set[10], e0002Set[10], e0011Set[10], e0012Set[10], e0022Set[10],
                   e0101Set[10], e0102Set[10], e0111Set[10], e0112Set[10], e0122Set[10],
                   e0202Set[10], e0211Set[10], e0212Set[10], e0222Set[10],
                   e1111Set[10], e1112Set[10], e1122Set[10],
                   e1212Set[10], e1222Set[10],
                   e2222Set[10],
                   H00,H11,H22,H01,H02,H12,
                   H00_2,H11_2,H22_2,H01_2,H01_3,H01_4,H02_2,H02_3,H02_4,H12_2,H12_3,H12_4,
                   H12_2mH11tH22,H11tH12,H00tH12,H00tH11,H01tH12,H01tH02,H02tH11,H02tH12,H11tH22,H00tH02,
                   denominatorComponent,denominatorInv;
            /*for (int i=0;i<10;i++) { //doubleSum w konstruktorze inicjalizuje zmienne (0)
                e0000Set[i]=0; e0001Set[i]=0; e0002Set[i]=0; e0011Set[i]=0; e0012Set[i]=0; e0022Set[i]=0;
                e0101Set[i]=0; e0102Set[i]=0; e0111Set[i]=0; e0112Set[i]=0; e0122Set[i]=0;
                e0202Set[i]=0; e0211Set[i]=0; e0212Set[i]=0; e0222Set[i]=0;
                e1111Set[i]=0; e1112Set[i]=0; e1122Set[i]=0;
                e1212Set[i]=0; e1222Set[i]=0;
                e2222Set[i]=0;
            }*/
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;

            //obliczanie epsilonów (wyniki) przy GLOBALNYCH średnich Hij, avEij - comment if NOT desired 1/2
            H00=avBoxMatrix[0]; H11=avBoxMatrix[1]; H22=avBoxMatrix[2];
            H01=avBoxMatrix[3]; H02=avBoxMatrix[4]; H12=avBoxMatrix[5];
            H00_2=H00*H00; H11_2=H11*H11; H22_2=H22*H22; H01_2=H01*H01; H01_3=H01*H01_2; H01_4=H01_2*H01_2; H02_2=H02*H02; H02_3=H02*H02_2; H02_4=H02_2*H02_2; H12_2=H12*H12; H12_3=H12*H12_2; H12_4=H12_2*H12_2;
            H12_2mH11tH22=H12_2-H11*H22; H11tH12=H11*H12; H00tH12=H00*H12; H00tH11=H00*H11; H01tH12=H01*H12; H01tH02=H01*H02; H02tH11=H02*H11; H02tH12=H02*H12; H11tH22=H11*H22; H00tH02=H00*H02;
            denominatorComponent=H02_2*H11-2*H01tH02*H12+H01_2*H22+H00*H12_2mH11tH22; denominatorInv=1/(2*denominatorComponent*denominatorComponent);
            //

            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    //obliczanie epsilonów (wyniki) przy GLOBALNYCH średnich Hij, avEij - uncomment if NOT desired 2/2
                    /*H00=avBoxMatrixSet[setIndex][0]; H11=avBoxMatrixSet[setIndex][1]; H22=avBoxMatrixSet[setIndex][2];
                    H01=avBoxMatrixSet[setIndex][3]; H02=avBoxMatrixSet[setIndex][4]; H12=avBoxMatrixSet[setIndex][5];
                    H00_2=H00*H00; H11_2=H11*H11; H22_2=H22*H22; H01_2=H01*H01; H01_3=H01*H01_2; H01_4=H01_2*H01_2; H02_2=H02*H02; H02_3=H02*H02_2; H02_4=H02_2*H02_2; H12_2=H12*H12; H12_3=H12*H12_2; H12_4=H12_2*H12_2;
                    H12_2mH11tH22=H12_2-H11*H22; H11tH12=H11*H12; H00tH12=H00*H12; H00tH11=H00*H11; H01tH12=H01*H12; H01tH02=H01*H02; H02tH11=H02*H11; H02tH12=H02*H12; H11tH22=H11*H22; H00tH02=H00*H02;
                    denominatorComponent=H02_2*H11-2*H01tH02*H12+H01_2*H22+H00*H12_2mH11tH22; denominatorInv=1/(2*denominatorComponent*denominatorComponent);*/
                    //
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                double h[6],h00,h11,h22,h01,h02,h12;
                int dataIndex=0; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10) {
                        if ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.')) dataIndex=10;
                        else h[dataIndex++]=strtod(data,NULL);
                    }
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    h00=h[0]; h11=h[1]; h22=h[2]; h01=h[3]; h02=h[4]; h12=h[5];
                    double h00_2=h00*h00,h11_2=h11*h11,h22_2=h22*h22,h01_2=h01*h01,h02_2=h02*h02,h12_2=h12*h12,
                           h12th22=h12*h22,h01th02=h01*h02,h00th01=h00*h01,h11th12=h11*h12,h01th12=h01*h12,h02th12=h02*h12,h00th02=h00*h02,h01th11=h01*h11,
                           h00ph11=h00+h11,h00ph22=h00+h22,h11ph22=h11+h22,

                           e00=(-H02_4*H11_2+4*H01*H02_3*H11tH12+H12_2*(-2*h01*H01*h12*H12+(h00_2-H00_2+h01_2)*H12_2+H01_2*(h12_2+h22_2))-2*H12*((h00_2-H00_2+h01_2)*H11tH12-h01*H01*(H11*h12+(h00ph11)*H12)+H01_2*(H00tH12+h12*(h11ph22)))*H22+(-H01_4-2*h01*H01*(h00ph11)*H11+(h00-H00)*(h00+H00)*H11_2+h01_2*(H01_2+H11_2)+H01_2*(h11_2+2*H00tH11+h12_2))*H22_2+H02_2*((h01_2-4*H01_2+h11_2+h12_2)*H12_2+H11_2*(h12_2+h22_2+2*H00*H22)-2*H11*(H12*(H00tH12+h12*(h11ph22))+H01_2*H22))+2*H02*(H01tH12*(-H11*(h12_2+h22_2)+H12*(2*H00tH12+h12*(h11ph22)))+2*H01_3*H12*H22-H01*(-h11*H11*h12+h11_2*H12+(h01_2+2*H00tH11+h12_2)*H12-H11*h12th22)*H22+h01*(H11*h12-(h00ph11)*H12)*H12_2mH11tH22)+h02_2*(H02_2*H11_2-2*H01tH02*H11tH12+H01_2*H12_2+H12_2mH11tH22*H12_2mH11tH22)+2*h02*(-h01*(H02tH11-H01tH12)*(H02tH12-H01*H22)-H12_2mH11tH22*(-h00*H02tH11+h00*H01tH12+H02*h12*H12-H02tH11*h22+H01tH12*h22-H01*h12*H22)))*denominatorInv,
                           e11=(h01_2*H02_4+H02_4*h11_2-H02_4*H11_2+H02_4*h12_2+2*H00*h01th02*H02_2*H12-2*h00th01*H02_3*H12-2*h01*H02_3*h11*H12-2*h02*H02_3*h12*H12+2*H00*H02_2*h11th12*H12+H00_2*h02_2*H12_2-2*h00*H00*h02*H02*H12_2+h00_2*H02_2*H12_2+h01_2*H02_2*H12_2+h02_2*H02_2*H12_2-2*H00*H02_2*H11*H12_2-2*H00*h01*H02*h12*H12_2+H00_2*h12_2*H12_2-H00_2*H12_4+2*H00*H02_2*h12*H12*h22-2*H00*h02*H02*H12_2*h22+H00_2*H12_2*h22_2+4*H01_3*H02tH12*H22-2*H00*(h01_2*H02_2+H02_2*(h11_2-H11_2+h12_2)+h01*(H00*h02-H02*(h00ph11))*H12-h02*H02*h12*H12+H00tH12*(-H11tH12+h12*(h11ph22)))*H22-H01_4*H22_2+H00_2*(h01_2+h11_2-H11_2+h12_2)*H22_2+H01_2*(H02_2*(h02_2+h12_2-4*H12_2+h22_2)-2*(H00*H12_2+H02*(H02tH11+h01th12+h02*(h00ph22)))*H22+(h00_2+h01_2+h02_2+2*H00tH11)*H22_2)+2*H01*(-H02*(h01*H02*(h02*H02-h12*H12)-h02*H02tH12*(h00ph22)+H00tH12*(h02_2+h12_2-2*H12_2+h22_2)+H02_2*(-2*H11tH12+h12*(h11ph22)))+(H02*(h01*H02*(h00ph11)+h02*H02*h12-(h00_2+h01_2+h02_2)*H12)+H00*(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*(-2*H11tH12+h12*(h11ph22))))*H22-H00*(h01*(h00ph11)+h02th12)*H22_2))*denominatorInv,
                           e22=(H02_2*(h00_2+h01_2+h02_2-H02_2)*H11_2+2*H00*H02tH11*(-h00th02*H11-h01*H11*h12+h00th01*H12+h01th11*H12+h02th12*H12-H02*H12_2-h02*H11*h22+H02tH11*H22)-2*H01_3*(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*h12*(h11ph22)-2*H02tH12*H22)+H01_4*(h02_2+h12_2+h22_2-H22_2)+H01_2*(-2*H00*h02_2*H11-2*H00tH11*h12_2+2*H00*h11th12*H12+h00_2*H12_2+h02_2*H12_2+2*h01*(H02tH11*h12+H00*h02*H12+H02*(h00ph11)*H12)+h01_2*(H02_2+H12_2)+2*H00*h12*H12*h22-2*H00tH11*h22_2+2*h02*H02*(h12*H12+H11*(h00ph22))-2*H00*H12_2*H22+2*H00tH11*H22_2+H02_2*(h11_2+h12_2-4*H12_2-2*H11tH22))+H00_2*(h02_2*H11_2-2*h01th02*H11tH12+H12_2*(h01_2+h11_2+h12_2-H12_2)+2*H11tH12*(-h12*(h11ph22)+H12*H22)+H11_2*(h12_2+h22_2-H22_2))-2*H01*(H02tH11*(h01*H02*(h00ph11)+h02*H02*h12+(h00_2+h01_2+h02_2-2*H02_2)*H12)+H00*(h01_2*H02tH12+h01*(-h02*H02tH11+H12*(-H11*h12+(h00ph11)*H12))+h02*H12*(h12*H12-H11*(h00ph22))+H02*(-h11*H11*h12+h11_2*H12+h12_2*H12-2*H12_3-H11*h12th22+2*H11tH12*H22))))*denominatorInv,
                           e01=(h01th02*H02_3*H11+H02_3*h11*H11*h12-h01_2*H02_3*H12-H02_3*h11_2*H12+H00*h02_2*H02tH11*H12-h00th02*H02_2*H11tH12-h01*H02_2*H11*h12*H12-H02_3*h12_2*H12+H00*H02tH11*h12_2*H12-H00*h01th02*H02*H12_2+2*h00th01*H02_2*H12_2+2*h01*H02_2*h11*H12_2+2*h02*H02_2*h12*H12_2-H00tH02*h11th12*H12_2+h00*H00*h02*H12_3-h00_2*H02*H12_3-h01_2*H02*H12_3-h02_2*H02*H12_3+H00*h01th12*H12_3+H02_3*H11*h12th22-h02*H02_2*H11tH12*h22-H00tH02*h12*H12_2*h22+H00*h02*H12_3*h22+H00*H02tH11*H12*h22_2+(H02tH11*(-H02*(h01*(h00ph11)+h02th12)+(h00_2+h01_2+h02_2)*H12)+H00*(h01_2*H02tH12-h01*(h02*H02tH11+H12*(H11*h12+(h00ph11)*H12))+H02*(-h11*H11*h12+h11_2*H12+h12_2*H12-H11*h12th22)-h02*H12*(h12*H12+H11*(h00ph22))))*H22+H00tH11*(h01*(h00ph11)+h02th12)*H22_2+H01_2*(H02tH12*(h02_2+h12_2+h22_2)-(h01*(h02*H02+h12*H12)+h02*H12*(h00ph22)+H02*h12*(h11ph22))*H22+(h01*(h00ph11)+h02th12)*H22_2)+H01*(-(H02_2*H11+H00*H12_2)*(h02_2+h12_2+h22_2)+(H02_2*(h11_2+h12_2)+2*h01*(H02tH11*h12+H00*h02*H12-H02*(h00ph11)*H12)+h01_2*(H02_2+H12_2)+2*h02*H02*(-h12*H12+H11*(h00ph22))+H12*((h00_2+h02_2)*H12+2*H00*h12*(h11ph22)))*H22-((h00_2+h01_2+h02_2)*H11+H00*(h01_2+h11_2+h12_2))*H22_2))*denominatorInv,
                           e02=(-H00*(h02_2*H02*H11_2+h01*H12_2*(H11*h12-(h00ph11)*H12)+H02*((h01_2+h11_2+h12_2)*H12_2-2*H11*h12*H12*(h11ph22)+H11_2*(h12_2+h22_2))+h02*H12*(-2*h01*H02tH11+H12*(-h12*H12+H11*(h00ph22))))+H00tH11*(h00th02*H11+h01*H11*h12-h00th01*H12-h01th11*H12-h02th12*H12+h02*H11*h22)*H22+H02tH11*(h00th02*H02tH11+h01*H02tH11*h12-h00th01*H02tH12-h01*H02*h11*H12-h02*H02*h12*H12+h00_2*H12_2+h01_2*H12_2+h02_2*H12_2+h02*H02tH11*h22-(h00_2+h01_2+h02_2)*H11tH22)+H01_2*(h02_2*H02tH11+2*h01th12*H12_2+2*h02*H12_2*(h00ph22)+H02tH11*(h12_2+h22_2)-H02*(h01_2+h11_2+h12_2)*H22-h01*(H11*h12+(h00ph11)*H12)*H22-h02*(h12*H12+H11*(h00ph22))*H22)+H01_3*(-H12*(h02_2+h12_2+h22_2)+(h01th02+h12*(h11ph22))*H22)+H01*(H02_2*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22))-H12*((h00_2+h02_2)*H12_2-H00*(h02_2*H11-h12*H12*(h11ph22)+H11*(h12_2+h22_2)))+((h00_2+h02_2)*H11tH12+H00*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22)))*H22+h01_2*H12*(H02_2-H12_2+(H00+H11)*H22)-2*h02*H02tH11*(H12*(h00ph22)-h12*H22)-h01*(2*H02tH11*(h12*H12-(h00ph11)*H22)+h02*(H02_2*H11+H00*(H12_2+H11tH22)))))*denominatorInv,
                           e12=(H02_2*H11*(h01*H02*(h00ph11)+h02*H02*h12-(h00_2+h01_2+h02_2)*H12)+H01_3*(-H02*(h02_2+h12_2+h22_2)+(h01th12+h02*(h00ph22))*H22)+H01_2*(2*h01th02*H02_2+2*H02_2*h11th12+H00*h02_2*H12+H00*h12_2*H12+2*H02_2*h12th22+H00tH12*h22_2-(h01*H02*(h00ph11)+h02*H02*h12+(h00_2+h01_2+h02_2)*H12+H00*(h01th02+h12*(h11ph22)))*H22)+H00tH02*(h01_2*H02tH12+h02*H12*(-h12*H12+2*H11*(h00ph22))+H02*(-h11*H11*h12+h11_2*H12+h12*(h12*H12-H11*h22))-h02*H11*h12*H22-h01*(h02*H02tH11-2*H11*h12*H12+(h00ph11)*H12_2+(h00ph11)*H11tH22))+H00_2*(-h02_2*H11tH12-H11tH12*(h12_2+h22_2)+H11*h12*(h11ph22)*H22+h01th02*(H12_2+H11tH22)+H12*(h12*H12*(h11ph22)-(h01_2+h11_2+h12_2)*H22))+H01*(-H02_3*(h11_2+h12_2)-h02*H02_2*H11*(h00ph22)+h01_2*H02*(-H02_2+H12_2+(H00+H11)*H22)-H00*h02*(H12_2*(h00ph22)+(-2*h12*H12+H11*(h00ph22))*H22)+H02*((h00_2+h02_2)*(H12_2+H11tH22)+H00*(h02_2*H11-2*h12*H12*(h11ph22)+H11*(h12_2+h22_2)+(h11_2+h12_2)*H22))-h01*(H02_2*H11*h12+2*H00*h02*H02tH12+H00*(-2*(h00ph11)*H12*H22+h12*(H12_2+H11tH22)))))*denominatorInv;

                    e0000Set[setIndex]+=e00*e00; e0001Set[setIndex]+=e00*e01; e0002Set[setIndex]+=e00*e02; e0011Set[setIndex]+=e00*e11; e0012Set[setIndex]+=e00*e12; e0022Set[setIndex]+=e00*e22;
                    e0101Set[setIndex]+=e01*e01; e0102Set[setIndex]+=e01*e02; e0111Set[setIndex]+=e01*e11; e0112Set[setIndex]+=e01*e12; e0122Set[setIndex]+=e01*e22;
                    e0202Set[setIndex]+=e02*e02; e0211Set[setIndex]+=e02*e11; e0212Set[setIndex]+=e02*e12; e0222Set[setIndex]+=e02*e22;
                    e1111Set[setIndex]+=e11*e11; e1112Set[setIndex]+=e11*e12; e1122Set[setIndex]+=e11*e22;
                    e1212Set[setIndex]+=e12*e12; e1222Set[setIndex]+=e12*e22;
                    e2222Set[setIndex]+=e22*e22;
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            doubleSum e0000/*=0*/, e0001/*=0*/, e0002/*=0*/, e0011/*=0*/, e0012/*=0*/, e0022/*=0*/,
                   e0101/*=0*/, e0102/*=0*/, e0111/*=0*/, e0112/*=0*/, e0122/*=0*/,
                   e0202/*=0*/, e0211/*=0*/, e0212/*=0*/, e0222/*=0*/,
                   e1111/*=0*/, e1112/*=0*/, e1122/*=0*/,
                   e1212/*=0*/, e1222/*=0*/,
                   e2222/*=0*/; //doubleSum w konstruktorze inicjalizuje zmienne (0)
            for (int i=0;i<10;i++) {
                e0000Set[i]/=dataLicznik*0.1; e0000+=e0000Set[i]; e0001Set[i]/=dataLicznik*0.1; e0001+=e0001Set[i]; e0002Set[i]/=dataLicznik*0.1; e0002+=e0002Set[i];
                e0011Set[i]/=dataLicznik*0.1; e0011+=e0011Set[i]; e0012Set[i]/=dataLicznik*0.1; e0012+=e0012Set[i]; e0022Set[i]/=dataLicznik*0.1; e0022+=e0022Set[i];
                e0101Set[i]/=dataLicznik*0.1; e0101+=e0101Set[i]; e0102Set[i]/=dataLicznik*0.1; e0102+=e0102Set[i]; e0111Set[i]/=dataLicznik*0.1; e0111+=e0111Set[i];
                e0112Set[i]/=dataLicznik*0.1; e0112+=e0112Set[i]; e0122Set[i]/=dataLicznik*0.1; e0122+=e0122Set[i];
                e0202Set[i]/=dataLicznik*0.1; e0202+=e0202Set[i]; e0211Set[i]/=dataLicznik*0.1; e0211+=e0211Set[i]; e0212Set[i]/=dataLicznik*0.1; e0212+=e0212Set[i];
                e0222Set[i]/=dataLicznik*0.1; e0222+=e0222Set[i];
                e1111Set[i]/=dataLicznik*0.1; e1111+=e1111Set[i]; e1112Set[i]/=dataLicznik*0.1; e1112+=e1112Set[i]; e1122Set[i]/=dataLicznik*0.1; e1122+=e1122Set[i];
                e1212Set[i]/=dataLicznik*0.1; e1212+=e1212Set[i]; e1222Set[i]/=dataLicznik*0.1; e1222+=e1222Set[i];
                e2222Set[i]/=dataLicznik*0.1; e2222+=e2222Set[i];
            }
            e0000*=0.1; e0001*=0.1; e0002*=0.1; e0011*=0.1; e0012*=0.1; e0022*=0.1;
            e0101*=0.1; e0102*=0.1; e0111*=0.1; e0112*=0.1; e0122*=0.1;
            e0202*=0.1; e0211*=0.1; e0212*=0.1; e0222*=0.1;
            e1111*=0.1; e1112*=0.1; e1122*=0.1;
            e1212*=0.1; e1222*=0.1;
            e2222*=0.1;
            //obliczenie bledow iloczynow elementow tensora odkształceń
            doubleSum dE0000/*=0*/, dE0001/*=0*/, dE0002/*=0*/, dE0011/*=0*/, dE0012/*=0*/, dE0022/*=0*/,
                   dE0101/*=0*/, dE0102/*=0*/, dE0111/*=0*/, dE0112/*=0*/, dE0122/*=0*/,
                   dE0202/*=0*/, dE0211/*=0*/, dE0212/*=0*/, dE0222/*=0*/,
                   dE1111/*=0*/, dE1112/*=0*/, dE1122/*=0*/,
                   dE1212/*=0*/, dE1222/*=0*/,
                   dE2222/*=0*/; //doubleSum w konstruktorze inicjalizuje zmienne (0)
            for (int i=0;i<10;i++) {double epsilon=e0000-e0000Set[i]; dE0000+=epsilon*epsilon;} dE0000=getAvErrorFromSumEps(dE0000,90.0); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {double epsilon=e0001-e0001Set[i]; dE0001+=epsilon*epsilon;} dE0001=getAvErrorFromSumEps(dE0001,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0002-e0002Set[i]; dE0002+=epsilon*epsilon;} dE0002=getAvErrorFromSumEps(dE0002,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0011-e0011Set[i]; dE0011+=epsilon*epsilon;} dE0011=getAvErrorFromSumEps(dE0011,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0012-e0012Set[i]; dE0012+=epsilon*epsilon;} dE0012=getAvErrorFromSumEps(dE0012,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0022-e0022Set[i]; dE0022+=epsilon*epsilon;} dE0022=getAvErrorFromSumEps(dE0022,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0101-e0101Set[i]; dE0101+=epsilon*epsilon;} dE0101=getAvErrorFromSumEps(dE0101,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0102-e0102Set[i]; dE0102+=epsilon*epsilon;} dE0102=getAvErrorFromSumEps(dE0102,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0111-e0111Set[i]; dE0111+=epsilon*epsilon;} dE0111=getAvErrorFromSumEps(dE0111,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0112-e0112Set[i]; dE0112+=epsilon*epsilon;} dE0112=getAvErrorFromSumEps(dE0112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0122-e0122Set[i]; dE0122+=epsilon*epsilon;} dE0122=getAvErrorFromSumEps(dE0122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0202-e0202Set[i]; dE0202+=epsilon*epsilon;} dE0202=getAvErrorFromSumEps(dE0202,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0211-e0211Set[i]; dE0211+=epsilon*epsilon;} dE0211=getAvErrorFromSumEps(dE0211,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0212-e0212Set[i]; dE0212+=epsilon*epsilon;} dE0212=getAvErrorFromSumEps(dE0212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e0222-e0222Set[i]; dE0222+=epsilon*epsilon;} dE0222=getAvErrorFromSumEps(dE0222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1111-e1111Set[i]; dE1111+=epsilon*epsilon;} dE1111=getAvErrorFromSumEps(dE1111,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1112-e1112Set[i]; dE1112+=epsilon*epsilon;} dE1112=getAvErrorFromSumEps(dE1112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1122-e1122Set[i]; dE1122+=epsilon*epsilon;} dE1122=getAvErrorFromSumEps(dE1122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1212-e1212Set[i]; dE1212+=epsilon*epsilon;} dE1212=getAvErrorFromSumEps(dE1212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1222-e1222Set[i]; dE1222+=epsilon*epsilon;} dE1222=getAvErrorFromSumEps(dE1222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e2222-e2222Set[i]; dE2222+=epsilon*epsilon;} dE2222=getAvErrorFromSumEps(dE2222,90.0);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            //Eijkl - bezwymiarowe [strain], volume - [sigma^3], kT=1[hard]
            double s0000=e0000*avVolume, dS0000=fabs(e0000*dAvVolume)+fabs(dE0000*avVolume),
                   s0001=e0001*avVolume, dS0001=fabs(e0001*dAvVolume)+fabs(dE0001*avVolume),
                   s0002=e0002*avVolume, dS0002=fabs(e0002*dAvVolume)+fabs(dE0002*avVolume),
                   s0011=e0011*avVolume, dS0011=fabs(e0011*dAvVolume)+fabs(dE0011*avVolume),
                   s0012=e0012*avVolume, dS0012=fabs(e0012*dAvVolume)+fabs(dE0012*avVolume),
                   s0022=e0022*avVolume, dS0022=fabs(e0022*dAvVolume)+fabs(dE0022*avVolume),
                   s0101=e0101*avVolume, dS0101=fabs(e0101*dAvVolume)+fabs(dE0101*avVolume),
                   s0102=e0102*avVolume, dS0102=fabs(e0102*dAvVolume)+fabs(dE0102*avVolume),
                   s0111=e0111*avVolume, dS0111=fabs(e0111*dAvVolume)+fabs(dE0111*avVolume),
                   s0112=e0112*avVolume, dS0112=fabs(e0112*dAvVolume)+fabs(dE0112*avVolume),
                   s0122=e0122*avVolume, dS0122=fabs(e0122*dAvVolume)+fabs(dE0122*avVolume),
                   s0202=e0202*avVolume, dS0202=fabs(e0202*dAvVolume)+fabs(dE0202*avVolume),
                   s0211=e0211*avVolume, dS0211=fabs(e0211*dAvVolume)+fabs(dE0211*avVolume),
                   s0212=e0212*avVolume, dS0212=fabs(e0212*dAvVolume)+fabs(dE0212*avVolume),
                   s0222=e0222*avVolume, dS0222=fabs(e0222*dAvVolume)+fabs(dE0222*avVolume),
                   s1111=e1111*avVolume, dS1111=fabs(e1111*dAvVolume)+fabs(dE1111*avVolume),
                   s1112=e1112*avVolume, dS1112=fabs(e1112*dAvVolume)+fabs(dE1112*avVolume),
                   s1122=e1122*avVolume, dS1122=fabs(e1122*dAvVolume)+fabs(dE1122*avVolume),
                   s1212=e1212*avVolume, dS1212=fabs(e1212*dAvVolume)+fabs(dE1212*avVolume),
                   s1222=e1222*avVolume, dS1222=fabs(e1222*dAvVolume)+fabs(dE1222*avVolume),
                   s2222=e2222*avVolume, dS2222=fabs(e2222*dAvVolume)+fabs(dE2222*avVolume),
                   //Voigt notation
                   S11=s0000, dS11=dS0000, S12=s0011, dS12=dS0011, S13=s0022, dS13=dS0022, S14=2*s0012, dS14=2*dS0012, S15=2*s0002, dS15=2*dS0002, S16=2*s0001, dS16=2*dS0001,
                   S22=s1111, dS22=dS1111, S23=s1122, dS23=dS1122, S24=2*s1112, dS24=2*dS1112, S25=2*s0211, dS25=2*dS0211, S26=2*s0111, dS26=2*dS0111,
                   S33=s2222, dS33=dS2222, S34=2*s1222, dS34=2*dS1222, S35=2*s0222, dS35=2*dS0222, S36=2*s0122, dS36=2*dS0122,
                   S44=4*s1212, dS44=4*dS1212, S45=4*s0212, dS45=4*dS0212, S46=4*s0112, dS46=4*dS0112,
                   S55=4*s0202, dS55=4*dS0202, S56=4*s0102, dS56=4*dS0102,
                   S66=4*s0101, dS66=4*dS0101,
                   //cubic average elements - compliances
                   S11c=(S11+S22+S33)/3.0, dS11c=(dS11+dS22+dS33)/3.0,
                   S12c=(S12+S13+S23)/3.0, dS12c=(dS12+dS13+dS23)/3.0,
                   S44c=(S44+S55+S66)/3.0, dS44c=(dS44+dS55+dS66)/3.0,
                   //cubic average elements - stiffness
                   C11c=pressure+(S11c+S12c)/(S11c*S11c+S11c*S12c-2*S12c*S12c), dC11c=(2*dS12c*fabs(S12c)*fabs(2*S11c+S12c)+dS11c*fabs(S11c*S11c+2*S11c*S12c+3*S12c*S12c))/pow(fabs((S11c-S12c)*(S11c+2*S12c)),2),
                   C12c=-pressure-S12c/(S11c*S11c+S11c*S12c-2*S12c*S12c), dC12c=(dS11c*fabs(S12c)*fabs(2*S11c+S12c)+dS12c*fabs(S11c*S11c+2*S12c*S12c))/pow(fabs((S11c-S12c)*(S11c+2*S12c)),2),
                   C44c=pressure+1/S44c, dC44c=fabs(dS44c/S44c/S44c),
                   //isotropic elastic moduli
                   B=1/(3*S11c+6*S12c), dB=(fabs(dS11c)+2*fabs(dS12c))/(3*pow(S11c+2*S12c,2)),
                   my1=1/(2*S11c-2*S12c), dMy1=(fabs(dS11c)+fabs(dS12c))/(2*pow(S11c-S12c,2)),
                   my2=1/S44c, dMy2=fabs(dS44c/S44c/S44c),
                   avMy=(my1+my2)/2.0, dAvMy=(dMy1+dMy2)/2.0,
                   E=(9*B*avMy)/(3*B+avMy), dE=(9*(3*fabs(B*B*dAvMy)+fabs(dB*avMy*avMy)))/pow(3*B+avMy,2),
                   //Poisson ratio in key directions: [100], [111], [110][1m10], [110][001]
                   nu_100_all=-S12c/S11c, dNu_100_all=fabs(dS12c/S11c)+fabs((dS11c*S12c)/S11c/S11c),
                   nu_111_all=(-2*S11c-4*S12c+S44c)/(2*(S11c+2*S12c+S44c)), dNu_111_all=(3*(fabs(dS44c*(S11c+2*S12c))+(fabs(dS11c)+2*fabs(dS12c))*fabs(S44c)))/(2*pow(S11c+2*S12c+S44c,2)),
                   nu_110_1m10=1-(4*(S11c+S12c))/(2*(S11c+S12c)+S44c), dNu_110_1m10=(4*fabs(dS44c)*fabs(S11c+S12c)+4*(fabs(dS11c)+fabs(dS12c))*fabs(S44c))/pow(2*(S11c+S12c)+S44c,2),
                   nu_110_001=-((4*S12c)/(2*(S11c+S12c)+S44c)), dNu_110_001=(4*(2*fabs(dS11c)*fabs(S12c)+fabs(dS44c*S12c)+fabs(dS12c*(2*S11c+S44c))))/pow(2*(S11c+S12c)+S44c,2),
                   //Bij for stability conditions
                   B11=(S11c+S12c)/(S11c*S11c+S11c*S12c-2*S12c*S12c),   //>0
                   B44=1/S44c,                                          //>0
                   B12dB11=-S12c/(S11c+S12c);                           //>=-1/2 && <=1
            printf("done\n");

            //wyznaczenie Probability Density Distribution Function
            /*printf("Determination of probability density distribution function... "); fflush(stdout);     //probDensDistFunMode - comment if not desired 4/7
            fileProbDensDistFun=fopen(bufferProbDensDistFun,"rt");
            double **avEdgeDistance=new double*[activeN]; for (int i=0;i<activeN;i++) avEdgeDistance[i]=new double[activeN];  //tablica tworzona dynamicznie zapisywana jest w HEAP a nie STACK (stack może zostać przepełniony i jest stackoverflow/segmentation fault   - pamiętać o delete[] tabel
            double licznik=0;
            for (int i=0;i<activeN;i++) for (int j=0;j<activeN;j++) avEdgeDistance[i][j]=0;
            char data[3][50]={"","",""}; int actIndex=0,character,dataType;
            while ((character=fgetc(fileProbDensDistFun))!=EOF) {
                if (character=='{') dataType=0;
                else if (character!=',' && character!='}') data[dataType][actIndex++]=character;
                else {
                    data[dataType++][actIndex++]=' '; actIndex=0;

                    if (dataType==3) {
                        int index[2]={(int)strtol(data[0],NULL,10),(int)strtol(data[1],NULL,10)};
                        if (index[0]==0 && index[1]==1) licznik++;  //i=0, j=1  identyfikuja pelen 'cykl' par (to pierwsza para jaka powinna zawsze byc zapisana wg metody updateNeighbours())
                        avEdgeDistance[index[0]][index[1]]+=strtod(data[2],NULL);
                    }
                    if (character=='}') fgetc(fileProbDensDistFun);
                }
            }
            fclose(fileProbDensDistFun);
            for (int i=0;i<activeN;i++) for (int j=0;j<activeN;j++) avEdgeDistance[i][j]/=licznik;
            printf("done\n");*/

            long timeEndMath=time(0);





/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) fileConfigurations = fopen(bufferConfigurations,"w");
            //fileProbDensDistFunResults = fopen(bufferProbDensDistFunResults,"w");     //probDensDistFunMode - comment if not desired 5/7

            fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,(double)avRho,(double)dAvRho,(double)avPacFrac,(double)dAvPacFrac,(double)avVolume,(double)dAvVolume,(double)avBoxMatrix[0],(double)dAvBoxMatrix[0],(double)avBoxMatrix[1],(double)dAvBoxMatrix[1],(double)avBoxMatrix[2],(double)dAvBoxMatrix[2],(double)avBoxMatrix[3],(double)dAvBoxMatrix[3],(double)avBoxMatrix[4],(double)dAvBoxMatrix[4],(double)avBoxMatrix[5],(double)dAvBoxMatrix[5],B,dB,avMy,dAvMy,my1,dMy1,my2,dMy2,E,dE,nu_100_all,dNu_100_all,nu_111_all,dNu_111_all,nu_110_1m10,dNu_110_1m10,nu_110_001,dNu_110_001,S11c,dS11c,S12c,dS12c,S44c,dS44c,C11c,dC11c,C12c,dC12c,C44c,dC44c,S11,dS11,S12,dS12,S13,dS13,S14,dS14,S15,dS15,S16,dS16,S22,dS22,S23,dS23,S24,dS24,S25,dS25,S26,dS26,S33,dS33,S34,dS34,S35,dS35,S36,dS36,S44,dS44,S45,dS45,S46,dS46,S55,dS55,S56,dS56,S66,dS66,B11,B44,B12dB11);
            fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,(double)avPacFrac,B,dB,avMy,dAvMy,E,dE,nu_100_all,dNu_100_all,nu_111_all,dNu_111_all,nu_110_1m10,dNu_110_1m10,nu_110_001,dNu_110_001,S11c,dS11c,S12c,dS12c,S44c,dS44c,C11c,dC11c,C12c,dC12c,C44c,dC44c);

            if (!onlyMath[0]) {
                rho=N/boxMatrixInnerParameters[1]; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2],deltaR,deltaV,deltaVND);
                for (int i=0;i<activeN;i++) fprintf(fileConfigurations,"%c[%.17E,%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].r[2],particles[i].diameter);
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Parallelepiped[{0,0,0},{{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E},{%.17E,%.17E,%.17E}}]},{Opacity[0.2],Green,Sphere[{%.12E,%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][0],boxMatrix[1][1],boxMatrix[1][2],boxMatrix[2][0],boxMatrix[2][1],boxMatrix[2][2],particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].r[2],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[2][2]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E, boxMatrix[0][2]=boxMatrix[2][0]=%.12E, boxMatrix[1][2]=boxMatrix[2][1]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[2][2],boxMatrix[0][1],boxMatrix[0][2],boxMatrix[1][2]);
            }

            /*for (int i=0;i<activeN;i++) {  //probDensDistFunMode - comment if not desired 6/7
                for (int j=0;j<activeN;j++) if (avEdgeDistance[i][j]!=0) fprintf(fileProbDensDistFunResults,"%d,%d,%.17E\n",i,j,avEdgeDistance[i][j]);
                delete [] avEdgeDistance[i];
            } delete [] avEdgeDistance;*/

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) fclose(fileConfigurations);
            //fclose(fileProbDensDistFunResults);  //probDensDistFunMode - comment if not desired 7/7
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec, time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n"); fflush(stdout);
}
