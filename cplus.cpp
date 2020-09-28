/*Wi-Fi Locator. The University of Melbourne.
Written in C++ by Patrick Tjahjadi. */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

// Hash define
#define MAX_WAPS 99 // Define maximum wireless access points in a grid as 99
#define TWENTY_FSPL 20 // Define transmission power of home wi-fi routers as 20
#define CONST_FSPL 32.45 // Define constant free-space path loss as 32.45
#define MAX_COORDS 99 // Define maximum number of observation points as 99
#define GRID_SIZE 77 // Define maximum length of grid as 77
#define MAP_HEIGHT 39 // Define maximum heigh of grid as 39
#define THRESHOLD -45 // Define threshold of transmission as -45dBm
#define PERCENTAGE 100 // Define percentage as out of 100
#define CELL_CENTRE 0.5 // Define centre of a space as 0.5 (between 0 and 1)
#define PLUS 0 // Define "+" as 0 for mapping
#define ONE_START -5 // Define "1" upper limit as -5 for mapping
#define ONE_LIMIT -15 // Define "1" lower limit as -15 for mapping
#define THREE_START -25 // Define "3" upper limit as -25 for mapping
#define THREE_LIMIT -35 // Define "3" lower limit as -35 for mapping

// Struct for wireless access points, containing power, frequency, x and y coordinates
struct wap_t {
    double power[MAX_WAPS];
    double freq[MAX_WAPS];
    double x_coord[MAX_WAPS];
    double y_coord[MAX_WAPS];
}; 

// Struct for observation points, containing x and y coordinates and signal strength
struct coordinates_t{
    double x_coord[MAX_COORDS];
    double y_coord[MAX_COORDS];
    double signal_strength[MAX_COORDS][MAX_WAPS];
};

// ---------------------------------Function Definitions---------------------------------
void read_waps(int* num_waps, double* signal_strength, wap_t* sequence);
void find_max_signal_strength(int num_waps, double* signal_strength);
void assess_points (wap_t sequence, int* obs_points, int num_waps);
void find_all_points(int num_waps, wap_t sequence);
double calculate_signal_strength_at_zero(wap_t* sequence, int* num_waps);
double calculate_signal_strength(wap_t sequence, coordinates_t observations, int i, int j);
void map_plotting(wap_t sequence, int num_waps);
// --------------------------------------------------------------------------------------

int main (int argc, char *argv[]) {
    int num_waps = 0;
    int obs_points = 0;
    double signal_strength[MAX_WAPS];
    wap_t sequence;
    read_waps(&num_waps, signal_strength, &sequence);
    find_max_signal_strength(num_waps, signal_strength);
    assess_points(sequence, &obs_points, num_waps);
    find_all_points(num_waps, sequence);
    map_plotting(sequence, num_waps);
    return 0;
}

// Function to read properties of wireless access points
void read_waps (int* num_waps, double* signal_strength, wap_t* sequence) {
    char letter;
    while(cin>>letter) {
        if (letter == 'P') {
            cin>>sequence->power[*num_waps]>>sequence->freq[*num_waps]
            >>sequence->x_coord[*num_waps]>>sequence->y_coord[*num_waps];

            signal_strength[*num_waps] = calculate_signal_strength_at_zero(sequence, num_waps);
            *num_waps += 1;
        }
        else {
            break;
        }
    }
}

// Function to find the maximum signal strength at (0, 0)
void find_max_signal_strength(int num_waps, double* signal_strength) {
    cout<<"Stage 1"<<endl<<"=========="<<endl;
    cout<<"Number of WAPs: "<<num_waps<<endl;
    cout<<fixed<<setprecision(2)<<"Maximum signal strength at (00.0, 00.0): "
    <<*max_element(signal_strength, signal_strength+num_waps)<<" dBm"<<endl;
}

// Function to calculate the signal strength at (0, 0)
double calculate_signal_strength_at_zero(wap_t* sequence, int* num_waps) {
    double euclidean_distance = 0;
    double fspl = 0;
    double signal_strength;
    euclidean_distance = sqrt(pow(sequence->x_coord[*num_waps], 2) + pow(sequence->y_coord[*num_waps], 2));
    fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence->freq[*num_waps]) + CONST_FSPL;
    signal_strength = sequence->power[*num_waps] - fspl;
    return signal_strength;
}

// Function to obtain data of observation points and its signal strengths
void assess_points (wap_t sequence, int* obs_points, int num_waps) {
    coordinates_t observations;
    int i = 0;
    int j = 0;
    char letter;
    do {
        cin>>observations.x_coord[i]>>observations.y_coord[i];
        *obs_points += 1;
        cin>>letter;
        i += 1;
    } while (letter == 'O');
    cout<<endl<<"Stage 2"<<endl;
    cout<<"=========="<<endl;
    for (i = 0; i < *obs_points; i++) {
        for (j = 0; j < num_waps; j++) {
            observations.signal_strength[i][j] = calculate_signal_strength(sequence, observations, i, j);
        }
        cout<<fixed<<setprecision(2)<<"Maximum signal strength at ("<<observations.x_coord[i]
        <<", "<<observations.y_coord[i]<<"): "<<*max_element(observations.signal_strength[i],
        observations.signal_strength[i]+num_waps)<<"dBm"<<endl;
    }
}

// Function to calculate the signal strength at a given point
double calculate_signal_strength(wap_t sequence, coordinates_t observations, int i, int j) {
    double euclidean_distance, fspl, signal_strength;

    euclidean_distance = sqrt(pow(sequence.x_coord[j] - observations.x_coord[i], 2) + 
    pow(sequence.y_coord[j] - observations.y_coord[i], 2));

    fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL * log10(sequence.freq[j]) + CONST_FSPL;
    signal_strength = sequence.power[j] - fspl;
    return signal_strength;
}

// Function to calculate the signal strengths of all points in the grid
void find_all_points(int num_waps, wap_t sequence) {
    int i, j, k;
    int num_points = GRID_SIZE * GRID_SIZE;
    int threshold_points = 0;
    double euclidean_distance, fspl;
    double* signal_strengths;
    double** max_signal_strengths = (double**)malloc(GRID_SIZE * sizeof(double*));
    // Determine the number of points in the grid
    for (i = 0; i < GRID_SIZE; i++) {
        max_signal_strengths[i] = (double*)malloc(GRID_SIZE * sizeof(double));
    }
    cout<<endl<<"Stage 3"<<endl;
    cout<<"=========="<<endl;
    cout<<num_points<<" points sampled"<<endl;
    // Determine the signal strengths of each point
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            double* signal_strengths = (double*)malloc(num_waps * sizeof(double));
            for (k = 0; k < num_waps; k++) {
                euclidean_distance = sqrt(pow(sequence.x_coord[k] - (i + 1), 2) + 
                pow(sequence.y_coord[k] - (j + 1), 2));

                fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL *
                log10(sequence.freq[k]) + CONST_FSPL;

                signal_strengths[k] = sequence.power[k] - fspl;
            }
            max_signal_strengths[i][j] = signal_strengths[0];
            for (k = 0; k < num_waps; k++) {
                if(signal_strengths[k] > max_signal_strengths[i][j]) {
                    max_signal_strengths[i][j] = signal_strengths[k];
                }
            }
            if (max_signal_strengths[i][j] <= THRESHOLD) {
                threshold_points += 1;
            }
        }
    }
    cout<<threshold_points<<" points (";
    cout<<fixed<<setprecision(2)<<(double(threshold_points) * PERCENTAGE /(GRID_SIZE * GRID_SIZE))
    <<"%) with maximum signal strength <="<<THRESHOLD<<" dBm"<<endl;

}

// Function to visualise a map that shows signal strenghts at various locations
void map_plotting(wap_t sequence, int num_waps){
    double max_strength, fspl, euclidean_distance;
    double cell_signal_strength[GRID_SIZE+1][MAP_HEIGHT][num_waps];
    int x, y, i, counter;
    cout<<endl<<"Stage 4"<<endl;
    cout<<"=========="<<endl;
    for (y = (MAP_HEIGHT * 2); y > 0; y -= 2) {
        counter = 0;
        for (x = 0; x <= GRID_SIZE; x += 1) {
            for (i = 0; i < num_waps; i += 1) {
                euclidean_distance = sqrt(pow((sequence.x_coord[i] - (x + CELL_CENTRE)), 2) + 
                pow((sequence.y_coord[i] - (y - 1)) , 2));

                fspl = TWENTY_FSPL * log10(euclidean_distance) + TWENTY_FSPL *
                log10(sequence.freq[i]) + CONST_FSPL;

                cell_signal_strength[x][counter][i] = sequence.power[i] - fspl;
            }
            max_strength = cell_signal_strength[x][counter][0];
            for (i = 0; i < num_waps; i += 1) {
                if (max_strength < cell_signal_strength[x][counter][i]) {
                    max_strength = cell_signal_strength[x][counter][i];
                }
            }
            // Output "+" if signal strength is higher than 0
            if (max_strength > PLUS) {
                cout<<"+";
            }
            // Output "1" if signal strength is between -5 and -15
            else if (max_strength < ONE_START && max_strength > ONE_LIMIT) {
                cout<<"1";
            }
            // Output "3" if signal strength is between -25 and -35
            else if (max_strength < THREE_START && max_strength > THREE_LIMIT) {
                cout<<"3";
            }
            // Output "-" if signal strenght is lower than -45
            else if (max_strength <= THRESHOLD) {
                cout<<"-";
            }
            else {
                cout<<" ";    
            }
            counter += 1;    
        }
        cout<<endl;
    }
}
