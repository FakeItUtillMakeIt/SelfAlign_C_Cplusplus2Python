
#ifndef ALIGN_H
#define ALIGN_H


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void mexFunction(int arg_count, double* arg_alignment_trace, double* arg_bad_levels, double* arg_cumulate_score,
    double* arg_alignment_matrix, double* arg_score_matrix, long int arg_num_rows, int arg_num_cols,
    double* arg_step_type_penalties, double* arg_modes, long int arg_num_modes, long int arg_lookback,
    char arg_alignment_type, double* arg_slip_location, bool arg_use_periodic_boundaries);


long int mod(long int a, long int b);

void align(
    double* score_matrix,
    long int num_rows,
    long int num_cols,
    double* step_penalties,
    double* bad_penalties,
    double* slip_penalties,
    double* modes,
    long int num_modes,
    long int lookback,
    double* slip_locations,
    int use_periodic_boundaries,
    double* alignment_trace,
    double* cumulate_score,
    double* level_is_bad,
    double* alignment_matrix,
    int num_step_penalties_per_row);

void selfAlign(
    double* score_matrix,
    long int num_rows,
    long int num_cols,
    double* step_penalties,
    double SA_penalty,
    double* alignment_trace,
    double* cumulate_score,
    double* alignment_matrix,
    int num_step_penalties_per_row);
    
int ret();

#endif