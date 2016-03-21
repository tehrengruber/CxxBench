#define CXX_BENCH_TIMER_START(state)    \
    state.start();                      \
    state.timer.start();

#define CXX_BENCH_TIMER_STOP(state)     \
    state.timer.stop();                 \
    state.stop();