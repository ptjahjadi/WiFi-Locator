/*Wi-Fi Locator. The University of Melbourne.
Written in C by Patrick Tjahjadi. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_WAPS 99
#define MAX_COORDS 99
#define CONST_FSPL 32.45
#define TWENTY_FSPL 20
#define GRID_SIZE 78
#define MAP_HEIGHT 39
#define STRENGTH_THRESHOLD -45
#define CELL_CENTRE 0.5
#define PLUS 0
#define ONE_START -5
#define ONE_LIMIT -15
#define THREE_START -25
#define THREE_LIMIT -35

typedef struct {
    double power[MAX_WAPS];
    double freq[MAX_WAPS];
    double x_coord[MAX_WAPS];
    double y_coord[MAX_WAPS];
} wap_t;

typedef struct {
    double x_coord[MAX_COORDS];
    double y_coord[MAX_COORDS];
    double signal_strength[MAX_COORDS][MAX_WAPS];
} coordinates_t;

void read_waps (int* num_waps, double* signal_strength, wap_t* sequence);
void calc_max_signal_strength(int num_waps, double* signal_strength);
void assess_points (wap_t sequence, int* obs_points, int num_waps);
void find_all_points(wap_t sequence, int num_waps);
void map_plotting(wap_t sequence, int num_waps);

int main(int argc, char* argv[]) {
    int num_waps = 0;
    int obs_points = 0;
    double signal_strength[MAX_WAPS];
    wap_t sequence;
    read_waps(&num_waps, signal_strength, &sequence);
    calc_max_signal_strength(num_waps, signal_strength);
    assess_points(sequence, &obs_points, num_waps);
    find_all_points(sequence, num_waps);
    map_plotting(sequence, num_waps);
    return 0;
}

void read_waps (int* num_waps, double* signal_strength, wap_t* sequence) {
    double fspl = 0;
    double euclidean_distance = 0;
    char letter;
    while (scanf("%c", &letter)) {
        if (letter == 'P') {
            scanf("%lf %lf %lf %lf\n", &sequence->power[*num_waps], &sequence->freq[*num_waps], &sequence->x_coord[*num_waps], &sequence->y_coord[*num_waps]);
            euclidean_distance = sqrt(pow(sequence->x_coord[*num_waps], 2) + pow(sequence->y_coord[*num_waps], 2));
            fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence->freq[*num_waps]) + CONST_FSPL;
            signal_strength[*num_waps] = sequence->power[*num_waps] - fspl;
            *num_waps += 1;
        }
        else {
            break;
        }
    }
}

void calc_max_signal_strength(int num_waps, double* signal_strength) {
    int i;
    double max_strength = signal_strength[0];
    for (i = 0; i < num_waps; i++) {
        if (signal_strength[i] > max_strength) {
            max_strength = signal_strength[i];
        }
    }
    printf("Stage 1\n");
    printf("==========\n");
    printf("Number of WAPs: %02d\n", num_waps);
    printf("Maximum signal strength at (00.0, 00.0): %.2lf dBm\n", max_strength);
}

void assess_points(wap_t sequence, int* obs_points, int num_waps) {
    char letter;
    int i, j;
    double max_strength, fspl;
    double euclidean_distance = 0;
    coordinates_t observations;
    do {
        scanf("%lf %lf\n", &observations.x_coord[*obs_points], &observations.y_coord[*obs_points]);
        *obs_points += 1;
        scanf("%c", &letter);
    } while (letter == 'O');
    printf("\nStage 2\n");
    printf("==========\n");
    for (i = 0; i < *obs_points; i++) {
        for(j = 0; j < num_waps; j++) {
            euclidean_distance = sqrt(pow((sequence.x_coord[j] - observations.x_coord[i]), 2) + pow((sequence.y_coord[j] - observations.y_coord[i]), 2));
            fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence.freq[j]) + CONST_FSPL;
            observations.signal_strength[i][j] = sequence.power[j] - fspl;
        }
        max_strength = observations.signal_strength[i][0];
        for (j = 1; j < num_waps; j++) {
            if (observations.signal_strength[i][j] > max_strength) {
                max_strength = observations.signal_strength[i][j];
            }
        }
        printf("Maximum signal strength at (%.1lf, %.1lf): %.2lf dBm\n", observations.x_coord[i], observations.y_coord[i], max_strength);
    }
}

void find_all_points(wap_t sequence, int num_waps) {
    double point_signal_strengths[GRID_SIZE-1][GRID_SIZE-1][num_waps];
    int grid_area = (GRID_SIZE-1) * (GRID_SIZE-1);
    double max_strength, fspl, euclidean_distance;
    int points_below = 0;
    double percentage_below;
    int x, y, i;
    for (x = 0; x < GRID_SIZE - 1; x++) {
        for (y = 0; y < GRID_SIZE - 1; y++) {
            for (i = 0; i < num_waps; i++) {
                euclidean_distance = sqrt(pow((sequence.x_coord[i] - (x + 1)), 2) + pow((sequence.y_coord[i] - (y + 1)), 2));
                fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence.freq[i]) + CONST_FSPL;
                point_signal_strengths[x][y][i] = sequence.power[i] - fspl;
            }
            max_strength = point_signal_strengths[x][y][0];
            for (i = 0; i < num_waps; i++) {
                if (point_signal_strengths[x][y][i] > max_strength) {
                    max_strength = point_signal_strengths[x][y][i];
                }
            }
            if (max_strength <= STRENGTH_THRESHOLD) {
                points_below += 1;
            }
        }
    }
    percentage_below = 100 * (double)points_below / (double)grid_area;
    printf("\nStage 3\n");
    printf("==========\n");
    printf("%d points sampled\n", grid_area);
    printf("%04d points (%.2lf%%) with maximum signal strength <= %d dBm\n", points_below, percentage_below, STRENGTH_THRESHOLD);
}

void map_plotting(wap_t sequence, int num_waps) {
    double cell_signal_strengths[GRID_SIZE][MAP_HEIGHT][num_waps];
    double max_strength, fspl, euclidean_distance;
    int x, y, i;
    int counter;
    printf("\nStage 4\n");
    printf("==========\n");
    counter = 0;
    for (y = MAP_HEIGHT*2; y > 0; y -= 2) {
        for (x = 0; x < GRID_SIZE; x++) {
            for (i = 0; i < num_waps; i++) {
                euclidean_distance = sqrt(pow((sequence.x_coord[i] - (x + CELL_CENTRE)), 2) + pow((sequence.y_coord[i] - (y - 1)), 2));
                fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence.freq[i]) + CONST_FSPL;
                cell_signal_strengths[x][counter][i] = sequence.power[i] - fspl;
            }
            max_strength = cell_signal_strengths[x][counter][0];
            for (i = 0; i < num_waps; i++) {
                if (cell_signal_strengths[x][counter][i] > max_strength) {
                    max_strength = cell_signal_strengths[x][counter][i];
                }   
            }
            if (max_strength > PLUS) {
                printf("+");
            }
            else if (max_strength < ONE_START && max_strength > ONE_LIMIT) {
                printf("1");
            }
            else if (max_strength < THREE_START && max_strength > THREE_LIMIT) {
                printf("3");
            }
            else if (max_strength <= STRENGTH_THRESHOLD) {
                printf("-");
            }
            else {
                printf(" ");    
            }

        }
        counter += 1;
        printf("\n");
    }
}