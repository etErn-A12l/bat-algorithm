#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>

#define POP_NO 40         // 44-43       // number of bats
#define DIM 71            // number of dimensions
#define ITERATIONS 200000 // number of iterations

const char *fName = "ftv70.txt";
const int opt_fit = 1950;
int **matrix;
int ROW;

int tst = 0;
typedef struct
{
    float pulse, vel[DIM], loud, freq;
    long fit;
    int pos[DIM];
} BAT;

void Initiate(BAT[]);
void calFitness(BAT[]);
int bestFitness(BAT[]);
void adjustFreq(BAT[]);
void updtVel(BAT[], BAT[], int);
void updtPos(BAT *, BAT *, int);
void genLclSol(BAT[], int);
void updateLoudPulse(BAT[], int);
int **read_matrix(const char *, int *);
void posShuffle(int[]);
void updateGnome(BAT[], BAT[]);
int fitness(int[]);
void formatTime(double);
void swap(int *, int, int);

void InverseMutation(int *);
void Cyclic_Crossover(int *, int *);
void CancyCrushUpdate(int *);
void VelocityMutation(int *, float *);
void Random_Mutation(int *);
int hasDuplicate(BAT[]);
float avg_vel(float[]);
void lin_kernighan(int *);
void simulated_annealing(int *, double, double);
void rgibnnm(int *);
int distance(int, int);
void reverse(int *, int, int);

int main(int argc, char const *argv[])
{
    // system("cls");
    clock_t start_time, end_time;
    start_time = clock(); // Record the starting time
    matrix = read_matrix(fName, &ROW);

    // for (int i = 0; i < ROW; i++)
    // {
    //     for (int j = 0; j < ROW; j++)
    //         printf("%d\t", matrix[i][j]);
    //     printf("\n\n");
    // }

    BAT x[POP_NO];

    srand(time(NULL)); // seedind with the value of current time

    // Generate random frequencies and loud
    float A = (float)rand() / (float)RAND_MAX;

    int Itr = 0, bestGen = 0;

    Initiate(x);

    calFitness(x);

    int bst_IDX = bestFitness(x); // Best fit in solutions
    BAT bestBat = x[bst_IDX];

    // =================================

    int test_count = 0;
    // printf("\n\n");
    // while (Itr++ < ITERATIONS)
    do
    {
        // system("cls");

        hasDuplicate(x);
        BAT y[POP_NO]; // New Solutions

        // Prepare New Solution
        for (int i = 0; i < POP_NO; i++)
            memcpy(y[i].pos, x[i].pos, ROW * sizeof(int));

        // printf("\n\n\n\t\t======= FITNESS OF %d ITERATION =======\n\n\n", Itr);
        // for (int i = 0; i < POP_NO; i++)
        // {
        //     printf("\n  BAT %d : \t%d  ", i + 1, x[i].fit);
        //     // printf("\n\n================\n\n");
        //     // for (int k = 0; k < ROW; k++)
        //     //     printf(" %d", x[i].pos[k]);
        // }

        printf("\r  Iteration : %6d | Best Fitness : %6ld | From BAT : %4d", Itr, x[bst_IDX].fit, bst_IDX + 1);

        adjustFreq(y);
        updtVel(x, y, bst_IDX);
        updtPos(x, y, bst_IDX);

        hasDuplicate(y);

        calFitness(y);

        for (int i = 0; i < POP_NO; i++)
        {
            float r = (float)rand() / (float)(RAND_MAX);

            // if (x[i].pulse > r)
            // {
            //     // Generating new local solutions
            //     genLclSol(x, i);
            //     test_count++;
            // }

            if (r < x[i].loud && y[i].fit < x[i].fit)
                // reduce loud and increase puse rate
                updateLoudPulse(x, Itr);
        }

        calFitness(x);

        // ==== RANKING THE NEW SOLUTIONS ====

        for (int i = 0; i < POP_NO; i++)
        {
            // printf("\n%d < %d ?\n", y[i].fit, x[i].fit);
            if (y[i].fit < x[i].fit) // Conditionally accept New Solution
            {
                // x[i].pulse = y[i].pulse;
                // x[i].loud = y[i].loud;
                // x[i].freq = y[i].freq;
                x[i].fit = y[i].fit;
                memcpy(x[i].pos, y[i].pos, ROW * sizeof(int));
                memcpy(x[i].vel, y[i].vel, ROW * sizeof(int));
            }
        }

        bst_IDX = bestFitness(x);

        if (x[bst_IDX].fit < bestBat.fit)
        {
            bestGen = Itr;
            // bestBat.pulse = x[bst_IDX].pulse;
            // bestBat.loud = x[bst_IDX].loud;
            bestBat.freq = x[bst_IDX].freq;
            bestBat.fit = x[bst_IDX].fit;
            memcpy(bestBat.pos, x[bst_IDX].pos, ROW * sizeof(int));
        }
        Itr++;
    } while (bestBat.fit > opt_fit);

    printf("\n\n  ");
    for (int c = 0; c < 100; c++)
        printf("-");
    printf("\n\n  Best Fitness at %d ITERATION with %ld FITNESS !\n", bestGen + 1, bestBat.fit);
    printf("\n  Need to Rdeuce More %ld Fitness :(\n", bestBat.fit - opt_fit);
    printf("\n  The GENOME WAS :\n\n ");
    for (int i = 0; i < ROW; i++)
        printf(" %i", bestBat.pos[i]);
    printf("\n\n  ");
    for (int c = 0; c < 100; c++)
        printf("-");
    printf("\n\n");
    end_time = clock(); // Record the ending time
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    formatTime(execution_time);
    // printf("\n\n  Test Count: %d", tst);
    printf("\n\n");

    return 0;
}

// =================================================================
//                       Initiator Functions
// =================================================================

int **read_matrix(const char *fName, int *rows)
{
    FILE *file = fopen(fName, "r");
    char line[10000];
    int row_count = 0;

    if (file == NULL)
    {
        printf("Error opening file\n");
        return NULL;
    }

    // Count rows and columns
    // while (fgets(line, sizeof(line), file))
    // {
    //     if (strstr(line, "9999999") != NULL)
    //         row_count++;
    // }
    // row_count++;
    row_count = DIM;
    printf("\n  %d x %d Cost Matrix\n\n", row_count, row_count);

    // Allocate memory for the 2D array
    int **matrix = (int **)malloc(row_count * sizeof(int *));
    for (int i = 0; i < row_count; i++)
        matrix[i] = (int *)malloc(row_count * sizeof(int));

    // Reset file pointer to beginning
    fseek(file, 0L, SEEK_SET);

    // Read the matrix values
    for (int i = 0; i < row_count; i++)
    {
        for (int j = 0; j < row_count; j++)
            fscanf(file, "%d", &matrix[i][j]);
    }

    fclose(file);

    *rows = row_count;
    return matrix;
}

void Initiate(BAT bats[])
{
    // Assigning random values
    for (int i = 0; i < POP_NO; i++)
    {
        int count = 0;

        for (int j = 0; j < ROW; j++) // Initializing 0 - ROW
        {
            bats[i].pos[j] = count++;
            bats[i].vel[j] = (float)rand() / (float)(RAND_MAX);
        }
        posShuffle(bats[i].pos);

        bats[i].pulse = ((float)rand() / (float)(RAND_MAX)) * (1 - 0) + 0;
        bats[i].loud = ((float)rand() / (float)(RAND_MAX)) * (2 - 1) + 1;
        bats[i].freq = ((float)rand() / (float)(RAND_MAX)) * (2 - 0) + 0;
    }
}

void calFitness(BAT bats[])
{
    for (int i = 0; i < POP_NO; i++)
        bats[i].fit = fitness(bats[i].pos);
}

// =================================================================
//                       Updater Functions
// =================================================================

void updtPos(BAT *bats, BAT *nBats, int B_ind)
{
    // printf("\n");
    for (int i = 0; i < POP_NO; i++)
    {
        InverseMutation(nBats[i].pos);
        if (bats[i].fit < fitness(nBats[i].pos))
        {
            simulated_annealing(nBats[i].pos, 100, 0.29);

            // lin_kernighan(nBats[i].pos);
            if (bats[i].fit < fitness(nBats[i].pos))
                // printf(" %d", i);
                rgibnnm(nBats[i].pos);
        }
    }
    // printf("\n");
}

void adjustFreq(BAT nBats[])
{
    for (int i = 0; i < POP_NO; i++)
        nBats[i].freq = ((float)rand() / (float)(RAND_MAX)) * (2 - 0) + 0;
}

void updtVel(BAT bats[], BAT nBats[], int Index)
{
    for (int i = 0; i < POP_NO; i++)
        for (int j = 0; j < DIM; j++)
            nBats[i].vel[j] = bats[i].vel[j] + (bats[i].pos[j] - bats[Index].pos[j]) * bats[i].freq;
}

void genLclSol(BAT bats[], int ind)
{
    /*
    float avgLoud, sum = 0;
    for (int i = 0; i < POP_NO; i++)
        sum += bats[i].loud;

    avgLoud = sum / POP_NO;



    float Eta = ((float)rand() / (float)(RAND_MAX)) * 2.0f - 1.0f; // Random number between 1 and -1
    // float Eta = 0.6294;

    for (int j = 0; j < DIM; j++)
    {
        bats[ind].pos[j] = bats[ind].pos[j] + Eta * avgLoud;

        // Boundary check
        if (bats[ind].pos[j] < 0)
            bats[ind].pos[j] = 0;
        else if (bats[ind].pos[j] > 2)
            bats[ind].pos[j] = 2;
    }

    */
    for (int i = 0; i < POP_NO; i++)
        lin_kernighan(bats[i].pos);
}

void updateLoudPulse(BAT bats[], int Iteration)
{
    float ALPHA = 0.9;
    float t = Iteration, GAMMA = 0.9, r0;

    for (int i = 0; i < POP_NO; i++)
    {
        // Reducing Loudness
        bats[i].loud = ALPHA * bats[i].loud;

        // Increasing Pulse Rate
        r0 = bats[i].pulse * 0.001 / 100; // 0.001 % of initial pulse rate
        bats[i].pulse = r0 * (1 - exp(-GAMMA * t));
    }
}

void updateGnome(BAT bats[], BAT nBats[])
{
    for (int i = 0; i < POP_NO; i++)
    {
        for (int j = 0; j < ROW; j++)
            nBats[i].pos[j] = bats[i].pos[j];

        posShuffle(nBats[i].pos);
    }
}

// apply the 2-opt move to the tour
// void two_opt(BAT bat, int **matrix)
// {
//     printf("\nEntered Two Opt");
//     int improved = 1;
//     while (improved)
//     {
//         improved = 0;
//         for (int i = 0; i < ROW - 2; i++)
//         {
//             for (int j = i + 2; j < ROW; j++)
//             {
//                 int delta = matrix[bat.pos[i]][bat.pos[j]] + matrix[bat.pos[i + 1]][bat.pos[j + 1]] - matrix[bat.pos[i]][bat.pos[i + 1]] - matrix[bat.pos[j]][bat.pos[j + 1]];
//                 if (delta < 0)
//                 {
//                     reverse(i + 1, j, bat.pos);
//                     improved = 1;
//                 }
//             }
//         }
//     }
// }

// apply the 2-opt move to the tour
void two_opt(BAT bat, int **matrix)
{
    // printf("\nEntered Two Opt");
    int mid = ROW / 2;
    int max = ROW - 1;
    int lower = (rand() % (mid - 0 + 1)) + 0;
    int upper = (rand() % (max - (mid + 1) + 1)) + mid;
    // printf("\nLower: %d\nUpper: %d", lower, upper);

    // printf("\n  OLD Fitness : %d", fitness(bat.pos, matrix));
    // reverse(lower, upper, bat.pos);
    while (lower < upper)
    {
        int tmp = bat.pos[lower];
        bat.pos[lower++] = bat.pos[upper];
        bat.pos[upper--] = tmp;
    }

    // printf("\n  New Fitness : %d", bat.fit);
}

void Random_Mutation(int *pos)
{
    int n = (ROW / 2) - 2;
    if (n % 2 == 1)
        n += 1; // Making Even number
    int indices[n];
    for (int k = 0; k < n; k++)
        indices[k] = rand() % ROW;

    for (int k = 0; k < n; k += 2)
        swap(pos, indices[k], indices[k + 1]);
}

void CancyCrushUpdate(int *pos)
{
    int min_fit = fitness(pos);
    int index = 0, temp;

    for (int q = 0; q < ROW / 2; q++)
    {
        swap(pos, q, ROW - 1 - q);
        int anew_fit = fitness(pos);
        if (anew_fit < min_fit)
        {
            index = q;
            min_fit = anew_fit;
        }
        // Revert
        swap(pos, q, ROW - 1 - q);
    }

    if (min_fit != fitness(pos))
        swap(pos, index, ROW - 1 - index);
}

void InverseMutation(int *pos)
{
    int lower = rand() % ROW, upper = rand() % ROW;

    while (lower == upper)
        upper = rand() % ROW;

    if (lower < upper)
        reverse(pos, lower, upper);

    else
    {
        int prev_low = lower;
        while (prev_low != upper)
        {
            swap(pos, lower, upper);
            lower = (lower + 1) % ROW;
            if (upper - 1 == -1)
                upper = ROW - 1;
            else
                --upper;
        }
    }
}

void VelocityMutation(int *posi, float *velo)
{
    float vel = avg_vel(velo);
    int jd = (int)(ROW * vel);
    if (jd % 2 != 0)
        jd -= 1;

    int pos[jd];
    for (int f = 0; f < jd; f++)
        pos[f] = rand() % ROW;

    for (int d = 0; d < jd; d += 2)
        swap(posi, pos[d], pos[d + 1]);
}

void Cyclic_Crossover(int *pos, int *tour)
{
    int child[ROW];
    for (short e = 0; e < ROW; e++)
        child[e] = -1;
    child[0] = pos[0];
    short z = 0, j, a = tour[z];
    while (a != pos[0])
    {
        a = tour[z];
        for (j = 0; j < ROW; j++)
        {
            if (pos[j] == a)
            {
                z = j;
                child[z] = pos[z];
            }
        }
    }
    for (int k = 1; k < ROW; k++)
    {
        if (-1 == child[k])
            child[k] = tour[k];
    }
    for (int k = 0; k < ROW; k++)
        pos[k] = child[k];
}

void simulated_annealing(int *tour, double initial_temp, double cooling_rate)
{
    int n = ROW;
    int *new_tour = malloc(n * sizeof(int));
    memcpy(new_tour, tour, n * sizeof(int));

    double temp = initial_temp;
    while (temp > 1)
    {
        int i = rand() % n;
        int j = rand() % n;
        while (i == j)
        {
            j = rand() % n;
        }
        if (i > j)
        {
            int temp = i;
            i = j;
            j = temp;
        }

        reverse(new_tour, i, j);

        long current_cost = 0;
        for (int k = 0; k < n; k++)
            current_cost += matrix[tour[k]][tour[(k + 1) % n]];

        long new_cost = 0;
        for (int k = 0; k < n; k++)
            new_cost += matrix[new_tour[k]][new_tour[(k + 1) % n]];

        if (new_cost < current_cost)
            memcpy(tour, new_tour, n * sizeof(int));

        else
        {
            double p = exp((current_cost - new_cost) / temp);
            if ((double)rand() / RAND_MAX < p)
                memcpy(tour, new_tour, n * sizeof(int));

            else
                memcpy(new_tour, tour, n * sizeof(int));
        }

        temp *= cooling_rate;
    }

    free(new_tour);
}

void lin_kernighan(int *p)
{
    int improved = 1;
    int n = 20;
    while (improved)
    {
        improved = 0;
        for (int i = 1; i < n - 2; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                int delta = matrix[p[i - 1]][p[i]] + matrix[p[j]][p[(j + 1) % n]] - matrix[p[i - 1]][p[j]] - matrix[p[i]][p[(j + 1) % n]];
                if (delta > 0)
                {
                    reverse(p, i, j);
                    improved = 1;
                }
            }
        }
    }
}

// void rgibnnm(int *tour)
// {
//     // printf("Entered");
//     int n = ROW;
//     int x = 0, improved = 0;
//     int s_tour[ROW];
//     memcpy(s_tour, tour, n * sizeof(int));
//     while (!improved && x < n)
//     {
//         // Select a random gene
//         int gene = x;
//         long old_fit = fitness(tour);
//         // Find the nearest neighbor of the selected gene
//         int nearest = -1;
//         int min_dist = INT_MAX;
//         for (int i = 0; i < n; i++)
//         {
//             if (i != gene)
//             {
//                 int dist = distance(tour[gene], tour[i]);
//                 if (dist < min_dist)
//                 {
//                     min_dist = dist;
//                     nearest = i;
//                 }
//             }
//         }

//         // Insert the selected gene beside its nearest neighbor
//         if (nearest != -1)
//         {
//             int temp = tour[gene];
//             if (gene < nearest)
//             {
//                 for (int i = gene; i < nearest; i++)
//                     tour[i] = tour[i + 1];

//                 tour[nearest] = temp;
//             }
//             else
//             {
//                 for (int i = gene; i > nearest + 1; i--)
//                     tour[i] = tour[i - 1];

//                 tour[nearest + 1] = temp;
//             }
//             // if (nearest != ROW - 1)
//             // swap(tour, gene, (nearest + 1) % ROW);
//             // else
//             //     swap(tour, gene, nearest - 1);
//         }
//         long new_fit = fitness(tour);
//         if (new_fit < old_fit)
//             improved = 1;
//         else
//         {
//             memcpy(tour, s_tour, n * sizeof(int));
//             x++;
//         }
//     }
//     // printf("...Exited");
// }

void rgibnnm(int *tour)
{
    // printf("Entered");
    int n = ROW;
    int x = rand() % n;

    // Select a random gene
    int gene = x;
    // Find the nearest neighbor of the selected gene
    int nearest = -1;
    int min_dist = INT_MAX;
    for (int i = 0; i < n; i++)
    {
        if (i != gene)
        {
            int dist = distance(tour[gene], tour[i]);
            if (dist < min_dist)
            {
                min_dist = dist;
                nearest = i;
            }
        }
    }

    // Insert the selected gene beside its nearest neighbor
    if (nearest != -1)
    {
        int temp = tour[gene];
        if (gene < nearest)
        {
            for (int i = gene; i < nearest; i++)
                tour[i] = tour[i + 1];

            tour[nearest] = temp;
        }
        else
        {
            for (int i = gene; i > nearest + 1; i--)
                tour[i] = tour[i - 1];

            tour[nearest + 1] = temp;
        }
        // if (nearest != ROW - 1)
        // swap(tour, gene, (nearest + 1) % ROW);
        // else
        //     swap(tour, gene, nearest - 1);
    }
    // printf("...Exited");
}

void reverse(int *tour, int i, int j)
{
    while (i < j)
    {
        swap(tour, i, j);
        i++;
        j--;
    }
}

void swap(int *pos, int i, int j)
{
    int temp = pos[i];
    pos[i] = pos[j];
    pos[j] = temp;
}

// =================================================================
//                       Helper Functions
// =================================================================

int bestFitness(BAT bats[])
{
    // Finding the best fit
    int Index = 0;
    for (int i = 0; i < POP_NO; i++)
        Index = (bats[i].fit < bats[Index].fit) ? i : Index;

    return Index;
}

void posShuffle(int pos[])
{
    for (int j = ROW - 1; j > 0; j--)
    {
        int u = rand() % (j + 1);
        int temp = pos[j];
        pos[j] = pos[u];
        pos[u] = temp;
    }
}

int distance(int x, int y)
{
    return matrix[x][y];
}

int fitness(int pos[])
{
    int fit = 0;
    for (int i = 0; i < ROW - 1; i++)
    {
        if (pos[i] != pos[i + 1])
            fit += matrix[pos[i]][pos[i + 1]];
    }
    fit += matrix[pos[ROW - 1]][pos[0]];
    return fit;
}

int hasDuplicate(BAT bats[])
{
    for (int h = 0; h < POP_NO; h++)
    {
        int i, j;
        for (i = 0; i < ROW; i++)
        {
            for (j = i + 1; j < ROW; j++)
            {
                if (bats[h].pos[i] == bats[h].pos[j])
                    return 1;
            }
        }
    }
    // If no duplicates are found, return false
    return 0;
}

float avg_vel(float x[])
{
    float s = 0;
    for (int g = 0; g < ROW; g++)
        s += x[g];
    return s / ROW;
}

void formatTime(double seconds)
{
    int hours = (int)seconds / 3600;
    int minutes = ((int)seconds % 3600) / 60;
    int remainingSeconds = ((int)seconds % 3600) % 60;
    printf("  Execution Time : ");
    if (hours > 0)
        printf("%d hour(s), ", hours);
    if (minutes > 0 || hours > 0)
        printf("%d minute(s), ", minutes);
    printf("%d second(s)\n", remainingSeconds);
}