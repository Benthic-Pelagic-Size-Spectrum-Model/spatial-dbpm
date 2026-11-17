// memory_outputs.h and memory_outputs.c contain the functions that return model
// output to R as arrays held in memory.

#ifndef MEMORY_OUTPUTS_H
#define MEMORY_OUTPUTS_H

#include "common.h"

void setup_output_arrays(RUN *run, GRID *grid, OUTPUT_ARRAYS *out);
void capture_output(RUN *run, GRID *grid, COMMUNITY *community, int time);

#endif
