#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <pthread.h>
#include <assert.h>

#define POP_NO 42         // 44-43       // number of bats
#define DIM 56           // number of dimensions
#define ITERATIONS 200000 // number of iterations
#define MAX_THREAD 6      // Number of threads
#define COOL_RATE 0.99
#define TEMP 100
#define NUM_THREADS 6
#define PAIRS_PER_THREAD 7

const char *data_file = "/home/xron/Drives/Docs_And_Media/Documents/Everything/Code-More/Application Development/Programming Languages/C/Others/HK mam TSP Prob 2/Datasets/ftv55.txt";
const char *sol_file = "/home/xron/Drives/Docs_And_Media/Documents/Everything/Code-More/Application Development/Programming Languages/C/Others/HK mam TSP Prob 2/Solutions/ftv55_sol.txt";

const int opt_fit = 1608;
int **matrix;
int ROW;

int tst = 0;
typedef struct
{
    float pulse, vel[DIM], loud, freq;
    long fit;
    int pos[DIM];
} BAT;

typedef struct
{
    BAT *bats;
    BAT *nBats;
    int start_index;
    int end_index;
} ThreadData;

void Initiate(BAT[]);
void calFitness(BAT[]);
int bestFitness(BAT[]);
void adjustFreq(BAT[]);
void updtVel(BAT[], BAT[], int);
void updtPos(BAT *, BAT *, int);
void genLclSol(BAT[], int, int);
void updateLoudPulse(BAT[], int);
int **read_matrix(const char *, int *);
void posShuffle(int[]);
void updateGnome(BAT[], BAT[]);
int fitness(int[]);
void formatTime(double);
void swap(int *, int, int);

void *processPairs(void *);
void InverseMutation(int *);
void simulated_annealing(int *, double, double);
void Cyclic_Crossover(int *, int *);
void rgibnnm(int *);
int distance(int, int);
void reverse(int *, int, int);
float max_velocity(float *);
int reverseImprovesTour(int *, int, int, int);
void threeOpt(int *);
void twoOpt(int *);

void Simul_Annel_3_Opt(int *, double, double);
void apply_3_opt(int *, int, int, int);
void Swap_Operator(int *tour);

int main(int argc, char const *argv[])
{
    system("clear");
    clock_t start_time, end_time;
    start_time = clock(); // Record the starting time
    matrix = read_matrix(data_file, &ROW);

    BAT x[POP_NO];

    srand(time(NULL)); // seedind with the value of current time

    // Generate random frequencies and loud
    float A = (float)rand() / (float)RAND_MAX;

    long Itr = 0, bestGen = 0;

    Initiate(x);

    calFitness(x);

    int bst_IDX = bestFitness(x); // Best fit in solutions
    BAT bestBat = x[bst_IDX];
    double execution_time;
    // =================================

    int test_count = 0;
    // printf("\n\n");
    // while (Itr++ < ITERATIONS)
    do
    {
        // system("clear");

        // hasDuplicate(x);
        BAT y[POP_NO]; // New Solutions

        // Prepare New Solution
        for (int i = 0; i < POP_NO; i++)
            memcpy(y[i].pos, x[i].pos, ROW * sizeof(int));

        // printf("\n\n\n\t\t======= FITNESS OF %ld ITERATION =======\n\n\n", Itr);
        // for (int i = 0; i < POP_NO; i++)
        // {
        //     printf("\n  BAT %d : \t%ld  ", i + 1, x[i].fit);
        //     // printf("\n\n================\n\n");
        //     // for (int k = 0; k < ROW; k++)
        //     //     printf(" %d", x[i].pos[k]);
        // }
        end_time = clock();
        double exectime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        printf("\r  Iteration : %6ld | Best Fitness : %6ld | Optimum : %4d", Itr, x[bst_IDX].fit, opt_fit);
        // formatTime(exectime);

        adjustFreq(y);
        updtVel(x, y, bst_IDX);
        updtPos(x, y, bst_IDX);

        // hasDuplicate(y);

        calFitness(y);

        for (int i = 0; i < POP_NO; i++)
        {
            float r = (float)rand() / (float)(RAND_MAX);

            if (x[i].pulse > r)
            {
                // Generating new local solutions
                genLclSol(x, i, bst_IDX);
                test_count++;
            }

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

                // Write Best Found Solution into Disk
                FILE *solf = fopen(sol_file, "w");
                fprintf(solf, "\n\n\t\t\t\t----* %d x %d Matrix *----\n\n\n", ROW, ROW);
                fprintf(solf, "\t  Iteration : %6ld | Best Fitness : %6ld | Optimum : %4d \n\n\nPATH:\n\n", Itr, x[bst_IDX].fit, opt_fit);

                for (int h = 0; h < ROW; h++)
                    fprintf(solf, "%d  ", x[bst_IDX].pos[h]);

                fprintf(solf, "\n\n\n================================\n\n\n");
                for (int f = 0; f < POP_NO; f++)
                    fprintf(solf, "\n  BAT %d : \t%ld  ", f + 1, x[f].fit);

                fclose(solf);
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
    printf("\n\n  Best Fitness at %ld ITERATION with %ld FITNESS !\n", bestGen + 1, bestBat.fit);
    printf("\n  Need to Rdeuce More %ld Fitness :(\n", bestBat.fit - opt_fit);
    printf("\n  The GENOME WAS :\n\n ");
    for (int i = 0; i < ROW; i++)
        printf(" %i", bestBat.pos[i]);
    printf("\n\n  ");
    for (int c = 0; c < 100; c++)
        printf("-");
    printf("\n\n");
    end_time = clock(); // Record the ending time
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    formatTime(execution_time);
    // printf("\n\n  Test Count: %d", tst);
    printf("\n\n");

    return 0;
}

// =================================================================
//                       Initiator Functions
// =================================================================

int **read_matrix(const char *data_file, int *rows)
{
    FILE *file = fopen(data_file, "r");
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

// void updtPos(BAT *bats, BAT *nBats, int B_ind)
// {
//     for (int i = 0; i < POP_NO; i++)
//     {
//         float cool = max_velocity(bats[i].vel);
//         simulated_annealing(nBats[i].pos, TEMP, cool);
//         if (bats[i].fit < fitness(nBats[i].pos))
//         {
//             Simul_Annel_3_Opt(nBats[i].pos, TEMP, cool);

//             if (bats[i].fit < fitness(nBats[i].pos))
//             {
//                 InverseMutation(nBats[i].pos);
//                 if (bats[i].fit < fitness(nBats[i].pos))
//                 {
//                     rgibnnm(nBats[i].pos);
//                     // twoOpt(nBats[i].pos);
//                     // threeOpt(nBats[i].pos);
//                     Swap_Operator(nBats[i].pos);
//                     //     Simul_Annel_3_Opt(nBats[i].pos, TEMP, cool);
//                     //     if (bats[i].fit < fitness(nBats[i].pos))
//                     //     {
//                     //         threeOpt(nBats[i].pos);
//                     //         if (bats[i].fit < fitness(nBats[i].pos))
//                     //         {
//                     //             twoOpt(nBats[i].pos);
//                     //         }
//                     //     }
//                 }
//             }
//         }
//     }
// }

void updtPos(BAT *bats, BAT *nBats, int B_ind)
{
    pthread_t threads[NUM_THREADS];
    ThreadData threadData[NUM_THREADS];

    int pairs_per_thread = PAIRS_PER_THREAD;
    int remaining_pairs = POP_NO % NUM_THREADS;

    int pair_count = 0;
    for (int i = 0; i < NUM_THREADS; i++)
    {
        threadData[i].bats = bats;
        threadData[i].nBats = nBats;

        int pairs_to_process = pairs_per_thread;
        if (remaining_pairs > 0)
        {
            pairs_to_process++;
            remaining_pairs--;
        }

        threadData[i].start_index = pair_count;
        threadData[i].end_index = pair_count + pairs_to_process - 1;

        pthread_create(&threads[i], NULL, processPairs, (void *)&threadData[i]);

        pair_count += pairs_to_process;
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }
}

void *processPairs(void *data)
{
    ThreadData *threadData = (ThreadData *)data;
    BAT *bats = threadData->bats;
    BAT *nBats = threadData->nBats;
    int start_index = threadData->start_index;
    int end_index = threadData->end_index;

    for (int i = start_index; i <= end_index; i++)
    {
        float cool = max_velocity(bats[i].vel);

        InverseMutation(nBats[i].pos);
        // twoOpt(nBats[i].pos);
        // Simul_Annel_3_Opt(nBats[i].pos, TEMP, cool);

        if (bats[i].fit < fitness(nBats[i].pos))
        {
            Simul_Annel_3_Opt(nBats[i].pos, TEMP, cool);

            // lin_kernighan(nBats[i].pos);
        }
        if (bats[i].fit < fitness(nBats[i].pos))
        {
            simulated_annealing(nBats[i].pos, TEMP, cool);

            // memcpy(nBats[i].pos, bats[i].pos, ROW * sizeof(int));
        }
        if (bats[i].fit < fitness(nBats[i].pos))
        {
            rgibnnm(nBats[i].pos);
            Swap_Operator(nBats[i].pos);
        }
    }

    // int i = rand() % (ROW - 1);
    // int j = rand() % (ROW - 1);
    // int k = rand() % (ROW - 1);

    // while (i == j || j == k || i == k)
    // {
    //     j = rand() % (ROW - 1);
    //     k = rand() % (ROW - 1);
    // }

    // int small = (i <= j && i <= k) ? i : ((j <= i && j <= k) ? j : k);
    // int big = (i >= j && i >= k) ? i : ((j >= i && j >= k) ? j : k);
    // int middle = (i != small && i != big) ? i : ((j != small && j != big) ? j : k);

    // apply_3_opt(nBats[i].pos, small, middle, big);
    // InverseMutation(nBats[i].pos);
    // Simul_Annel_3_Opt(nBats[i].pos, TEMP, cool);

    pthread_exit(NULL);
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
        {
            // float vel = bats[i].vel[j] + (bats[i].pos[j] - bats[Index].pos[j]) * bats[i].freq;
            // while (vel > 1)
            //     vel /= 10;
            float vel = bats[i].vel[j] + (distance(bats[i].pos[j], bats[Index].pos[j]) * bats[i].freq);
            int i_vel = (int)vel;
            nBats[i].vel[j] = vel - i_vel;
        }
}

void genLclSol(BAT bats[], int ind, int bst)
{
    Cyclic_Crossover(bats[ind].pos, bats[bst].pos);
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

void InverseMutation(int *pos)
{
    int lower = rand() % ROW, upper = rand() % ROW;

    while (lower == upper)
        upper = rand() % ROW;

    reverse(pos, lower, upper);
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
            j = rand() % n;

        // if (i > j)
        // {
        //     int temp = i;
        //     i = j;
        //     j = temp;
        // }

        // reverse(new_tour, i, j);

        /***** Modified *****/

        reverse(new_tour, i, j);

        /***** Modified *****/

        long current_cost = fitness(tour);
        long new_cost = fitness(new_tour);

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

void rgibnnm(int *tour)
{
    // printf("Entered");
    int n = ROW;
    int x = 0, improved = 0;
    int s_tour[ROW];
    memcpy(s_tour, tour, n * sizeof(int));
    while (!improved && x < n)
    {
        // Select a random gene
        int gene = x;
        long old_fit = fitness(tour);
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
        long new_fit = fitness(tour);
        if (new_fit < old_fit)
            improved = 1;
        else
        {
            memcpy(tour, s_tour, n * sizeof(int));
            x++;
        }
    }
    // printf("...Exited");
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

void twoOpt(int *tour)
{
    // printf("\nEntered two opt");
    bool improved = false;
    int tmp_tour[ROW];
    // while (improved)
    // {
    // improved = false;
    for (int i = 0; i < ROW - 1; i++)
    {
        for (int j = i + 1; j < ROW; j++)
        {
            memcpy(tmp_tour, tour, ROW * sizeof(int));
            reverse(tmp_tour, i, j);
            long Before = fitness(tour);
            long After = fitness(tmp_tour);

            if (After < Before)
            {
                reverse(tour, i, j);
                improved = true;
                // break;
                continue;
            }
            else
            {
                memcpy(tmp_tour, tour, ROW * sizeof(int));
                reverse(tmp_tour, j, i);
                After = fitness(tmp_tour);
                if (After < Before)
                {
                    reverse(tour, j, i);
                    improved = true;
                }
            }
        }
        // if (improved)
        //     break;
    }
    // }
}

void Swap_Operator(int *tour)
{
    bool improved = false;
    int tmp_tour[ROW];
    memcpy(tmp_tour, tour, ROW * sizeof(int));
    for (int i = 0; i < ROW - 1; i++)
    {
        for (int j = i + 1; j < ROW; j++)
        {
            swap(tmp_tour, i, j);
            long Before = fitness(tour);
            long After = fitness(tmp_tour);

            if (After < Before)
            {
                // printf("\nSwap Worked !");
                swap(tour, i, j);
                improved = true;
                break;
            }
            else
                swap(tmp_tour, i, j);
        }
        if (improved)
            break;
    }
}

void threeOpt(int *tour)
{
    // printf("\nEntered three opt");
    bool improved = false;
    // while (improved)
    // {
    //     improved = false;
    for (int i = 0; i < ROW - 2; i++)
    {
        for (int j = i + 1; j < ROW - 1; j++)
        {
            for (int k = j + 1; k < ROW; k++)
            {
                if (reverseImprovesTour(tour, i, j, k) == 1)
                {
                    reverse(tour, i, j);
                    reverse(tour, j + 1, k);
                    improved = true;
                    // printf("\nFound Improve: %d", fitness(tour));
                    // break;
                }
            }
            // if (improved)
            //     break;
        }
        // if (improved)
        //     break;
    }
    // }
    // printf("\nDid not Improve");
}

void apply_3_opt(int *path, int i, int j, int k)
{
    // Apply the 3-opt move to the solution
    int tmp_path[ROW];
    memcpy(tmp_path, path, ROW * sizeof(int));

    reverse(tmp_path, i + 1, j);
    reverse(tmp_path, j + 1, k);

    memcpy(path, tmp_path, ROW * sizeof(int));
}

void Simul_Annel_3_Opt(int *path, double temperature, double cooling_rate)
{

    int tour[ROW];
    memcpy(tour, path, ROW * sizeof(int));
    double temp = temperature;
    while (temp > 1)
    {
        // Generate a new solution by applying a 3-opt move
        int i = rand() % (ROW - 1);
        int j = rand() % (ROW - 1);
        int k = rand() % (ROW - 1);

        while (i == j || j == k || i == k)
        {
            j = rand() % (ROW - 1);
            k = rand() % (ROW - 1);
        }

        int small = (i <= j && i <= k) ? i : ((j <= i && j <= k) ? j : k);
        int big = (i >= j && i >= k) ? i : ((j >= i && j >= k) ? j : k);
        int middle = (i != small && i != big) ? i : ((j != small && j != big) ? j : k);

        // printf("\nEEntering 3-opt: i = %d, j = %d, k = %d", i, j, k);

        apply_3_opt(tour, small, middle, big);

        long current_cost = fitness(path);
        long new_cost = fitness(tour);
        // printf("\nAAAAAA");
        if (new_cost < current_cost)
            memcpy(path, tour, ROW * sizeof(int));

        else
        {
            double p = exp((current_cost - new_cost) / temp);
            if ((double)rand() / RAND_MAX < p)
                memcpy(path, tour, ROW * sizeof(int));

            // else
            //     memcpy(tour, path, ROW * sizeof(int));
        }

        temp *= cooling_rate;
    }
}

void reverse(int *tour, int i, int j)
{
    if (i < j)
    {
        while (i < j)
        {
            swap(tour, i, j);
            i++;
            j--;
        }
    }
    else
    {
        int prev_low = i;
        while (prev_low != j)
        {
            swap(tour, i, j);
            i = (i + 1) % ROW;
            if (j - 1 == -1)
                j = ROW - 1;
            else
                --j;
        }
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
        assert(pos[i] != pos[i + 1]);
        if (pos[i] != pos[i + 1])
            fit += distance(pos[i], pos[i + 1]);
    }
    fit += distance(pos[ROW - 1], pos[0]);
    return fit;
}

int reverseImprovesTour(int *tour, int i, int j, int k)
{
    int tmp_tour[ROW];
    memcpy(tmp_tour, tour, ROW * sizeof(int));

    long distBefore = fitness(tour);
    // Reverse the segment from i+1 to j
    reverse(tmp_tour, i + 1, j);
    // Reverse the segment from j+1 to k
    reverse(tmp_tour, j + 1, k);

    long distAfter = fitness(tmp_tour);

    // Compare total distances
    if (distAfter < distBefore)
    {
        // printf("\nReturned 1");
        memcpy(tour, tmp_tour, ROW * sizeof(int));
        return 1; // Improvement
    }
    else
    {
        return 0; // No improvement
    }
}

void formatTime(double seconds)
{
    int hours = (int)seconds / 3600;
    int minutes = ((int)seconds % 3600) / 60;
    int remainingSeconds = ((int)seconds % 3600) % 60;
    printf("  Exec Time : ");
    if (hours > 0)
        printf("%d hr, ", hours);
    if (minutes > 0 || hours > 0)
        printf("%d min, ", minutes);
    printf("%d sec", remainingSeconds);
}

float max_velocity(float *vel)
{
    float max = 0;
    for (int i = 0; i < ROW; i++)
        max = (vel[i] > max) ? vel[i] : max;

    return max;
}