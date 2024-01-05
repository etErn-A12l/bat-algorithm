#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *file;
    int count = 0;
    int num;

    // Open the file
    file = fopen("F:/Documents/TEMP/Yo/Yo/ft70.txt", "r");
    if (file == NULL) {
        printf("Failed to open the file.\n");
        return 1;
    }

    // Count the integers in the file
    while (fscanf(file, "%d", &num) == 1) {
        count++;
    }

    // Close the file
    fclose(file);

    // Print the count
    printf("Number of integers in the file: %d\n", count);

    return 0;
}