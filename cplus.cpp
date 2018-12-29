#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

#define MAX_WAPS 99
#define TWENTY_FSPL 20
#define CONST_FSPL 32.45

struct wap_t {
    double power[MAX_WAPS];
    double freq[MAX_WAPS];
    double x_coord[MAX_WAPS];
    double y_coord[MAX_WAPS];
};

void read_waps (int* num_waps, double* signal_strength, wap_t* sequence);
void calc_max_signal_strength(int num_waps, double* signal_strength);


int main (int argc, char *argv[]) {
    int num_waps = 0;
    double signal_strength[MAX_WAPS];
    wap_t sequence;
    read_waps(&num_waps, signal_strength, &sequence);
    calc_max_signal_strength(num_waps, signal_strength);
    return 0;
}

void read_waps (int* num_waps, double* signal_strength, wap_t* sequence) {
    double fspl = 0;
    double euclidean_distance = 0;
    char letter;
    while(cin>>letter) {
        if (letter == 'P') {
            cin>>sequence->power[*num_waps]>>sequence->freq[*num_waps]>>sequence->x_coord[*num_waps]>>sequence->y_coord[*num_waps];
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
    cout<<"Stage 1"<<endl<<"=========="<<endl;
    cout<<"Number of WAPs: "<<num_waps<<endl;
    cout<<fixed<<setprecision(2)<<"Maximum signal strength at (00.0, 00.0): "<<max_strength<<" dBm"<<endl;
}

