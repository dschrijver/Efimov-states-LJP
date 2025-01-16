#define _POSIX_SOURCE
#include <stdio.h> // printf(), fflush(), stdout
#include <unistd.h> // fork(), pipe(), close()
#include <signal.h> // kill(), SIGKILL
#include <stdlib.h> // srand()
#include <time.h> // time(), clock(), clock_t, CLOCKS_PER_SEC, sleep(), time_t, struct tm, localtime
#include <math.h> // sqrt()

#include "src/child_processes.h" // child_process()
#include "src/update_network.h" // update_prep_network()
#include "src/matrix_mat.h" // matuniform(), vecuniform(), printmat(), printrowvec()
#include "main.h"


int main(void) {
    srand(time(NULL));

    // Use the following if the calculation starts with a random initialization
    matuniform((double)-1/sqrt((double)N_in), (double)1/sqrt((double)N_in), N_in, N_hid, w1);
    vecuniform((double)-1/sqrt((double)N_hid), (double)1/sqrt((double)N_hid), N_hid, w2);

    int pids[CHILD_PROCESSES];

    // Create pipes
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        pipe(pipes[2*i]); // Sends to child process
        pipe(pipes[2*i + 1]); // Receives from child process
    }

    // Create child processes
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        fflush(stdout); // To help printing
        pids[i] = fork();
        if (pids[i] == 0) {
            child_process(i);
            // The child process will not return here
        }
    }

    // Close pipes that are not needed
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        // Main process sends to pipes[2*i][1]
        close(pipes[2*i][0]);
        // Main process reads from pipes[2*i + 1][0]
        close(pipes[2*i + 1][1]);
    }

    sleep(1); // Wait for all processes to be started, not actually necessary

    clock_t begin, end;
    double time_spent;
    time_t now;
    struct tm *now_tm;
    int hour, minute, second;

    // Preparation
    begin = clock();
    now = time(NULL);
    now_tm = localtime(&now);
    hour = now_tm->tm_hour;
    minute = now_tm->tm_min;
    second = now_tm->tm_sec;
    printf("\nStarted preparation... (%02d:%02d:%02d)\n\n", hour, minute, second);
    fflush(stdout);
    double h = prep_network();
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nPreparation took %.1f minutes\n\n", time_spent/(double)60.0);
    fflush(stdout);

    // Minimization
    begin = clock();
    now = time(NULL);
    now_tm = localtime(&now);
    hour = now_tm->tm_hour;
    minute = now_tm->tm_min;
    second = now_tm->tm_sec;
    printf("\nStarted minimization... (%02d:%02d:%02d)\n\n", hour, minute, second);
    fflush(stdout);
    h = minimize_energy(h);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nMinimization took %.1f minutes\n\n", time_spent/(double)60.0);
    fflush(stdout);

    // Sweep
    begin = clock();
    now = time(NULL);
    now_tm = localtime(&now);
    hour = now_tm->tm_hour;
    minute = now_tm->tm_min;
    second = now_tm->tm_sec;
    printf("\nStarted sweep... (%02d:%02d:%02d)\n\n", hour, minute, second);
    fflush(stdout);
    (void)sweep(h);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nSweep took %.1f minutes\n\n", time_spent/(double)60.0);
    fflush(stdout);

    // Print final network parameters
    printf("w1 = ");
    printmat(N_in, N_hid, w1);
    printf("b1 = ");
    printrowvec(N_hid, b1);
    printf("w2 = ");
    printrowvec(N_hid, w2);
    fflush(stdout);
    
    // Kill child processes
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        kill(pids[i], SIGKILL);
    }

    return 0;
}