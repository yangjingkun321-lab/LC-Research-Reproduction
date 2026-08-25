#ifndef RL_SERVER_H
#define RL_SERVER_H

#include "state.h"

// Run one interactive RL episode.
//
// The GlobalState is kept alive for the whole process.
// Commands are read from stdin.
//
// Supported commands:
//
// STATE
// STEP <loop_id>
// QUIT
//
// Protocol responses are prefixed by [RL] so that Python can
// distinguish them from normal LoopyCuts profiler output.
int run_rl_server(GlobalState &state);

#endif // RL_SERVER_H