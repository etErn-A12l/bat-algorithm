#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int **read_matrix(const char *filename, int *rows);
int row_number(const char *filename);

int main()
{
    const char *filename = "F:/Documents/TEMP/HK Mam/matrix.txt";
    int rows;
    int **matrix = read_matrix(filename, &rows);
    int cols = rows;

    // Use the matrix as needed
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n\n");
    }

    // Free allocated memory
    for (int i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}

int **read_matrix(const char *filename, int *rows)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error opening file\n");
        return NULL;
    }

    // Count rows and columns
    int row_count = row_number(filename);

    printf("\nRows: %d\nCols: %d\n\n\n", row_count, row_count);

    // Allocate memory for the 2D array
    int **matrix = (int **)malloc(row_count * sizeof(int *));
    for (int i = 0; i < row_count; i++)
    {
        matrix[i] = (int *)malloc(row_count * sizeof(int));
    }

    // Reset file pointer to beginning
    fseek(file, 0L, SEEK_SET);

    // Read the matrix values
    for (int i = 0; i < row_count; i++)
    {
        for (int j = 0; j < row_count; j++)
        {
            fscanf(file, "%d", &matrix[i][j]);
        }
    }

    fclose(file);

    *rows = row_count;
    return matrix;
}

int row_number(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error opening file\n");
        return 0;
    }

    // Count rows and columns
    int row_count = 1;
    int ch;
    while ((ch = fgetc(file)) != EOF)
    {
        if (ch == '\n')
            row_count++;
    }

    printf("\nRows: %d\n\n\n", row_count);
    fclose(file);
    return row_count;
}