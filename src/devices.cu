#include <iomanip>
#include <iostream>
#include <math.h>
#include "devices.h"
using namespace std;

size_t num_digits(int number) {
    if ( number == 0 ) {
        return 1;
    }
    int sign = number < 0 ? 1 : 0;
    return floor(log10(abs(number))) + 1 + sign;
}

enum columns{Id,Name,Version,MPCs,Memory};

void show_devices() {
    const char *pipe = " \u2502", *cross = "\u253c", *bar = "\u2500";
    int devices;
    cudaGetDeviceCount(&devices);
    printf("DEVICES\n");
    const int num_cols = 5;
    const char *headers[] = {"Id", "Name", "Version", "MPCs", "Memory"};
    size_t lens[num_cols];
    for( int i = 0; i < num_cols; ++i ) {
        lens[i] = strlen(headers[i]);
    }
    lens[Id] = max(lens[Id], num_digits(devices));
    for( int i = 0; i < devices; ++i ) {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        lens[Name] = max(lens[Name], strlen(props.name));
        lens[Version] = max(lens[Version],
                num_digits(props.major) + num_digits(props.minor) + 2);
        lens[MPCs] = max(lens[MPCs], num_digits(props.multiProcessorCount));
        lens[Memory] = max(lens[Memory], num_digits(props.totalGlobalMem >> 20));
    }
    for(int i=0;i<num_cols;++i){
        ++lens[i];
    }
    cout << setw(lens[Id]) << headers[Id];
    for( int i = 1; i < num_cols; ++i ) {
        cout << pipe << setw(lens[i]) << headers[i];
    }
    cout << endl;
    for( int i = 0; i < num_cols; ++i ) {
        if ( i > 0 ) {
            cout << cross;
        }
        for ( int j = 0; j < lens[i] + 1; ++j ) {
            cout << bar;
        }
    }
    cout << endl;

    char *version = (char*)malloc((lens[1]+1) * sizeof(char));
    for( int i = 0; i < devices; ++i ) {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);

        cout << setw(lens[Id]) << i << pipe;
        cout << setw(lens[Name]) << props.name << pipe;
        sprintf(version, "v%d.%d", props.major, props.minor);
        cout << setw(lens[Version]) << version << pipe;
        cout << setw(lens[MPCs]) << props.multiProcessorCount << pipe;
        cout << setw(lens[Memory]) << (props.totalGlobalMem >> 20) << endl;
    }
    free(version);
}
