/* Shared measurement driver for the standalone PolyBench-Tiramisu wrappers.
 *
 * Included at the end of each generated *_wrapper.cpp, which must define
 *
 *     static void tpb_single_measurement(int dump_arrays);
 *
 * performing exactly one PolyBench-style measurement: aligned allocation,
 * deterministic init_array initialization (from the *_init.h sidecar),
 * polybench_timer_start (cache flush + gettimeofday), one kernel
 * invocation, polybench_timer_stop, and printing the elapsed time in
 * MILLISECONDS to stdout.
 *
 * Protocol (mirrors TiraLib's PolybenchHarness wrapper):
 *   - each measurement runs in a freshly exec'ed process, like
 *     PolyBench's time_benchmark.sh re-running the benchmark binary;
 *   - number of measurements: argv[1] if given, else $MIN_RUNS, else 1;
 *   - one time (ms) per measurement on stdout, space-separated;
 *   - if $TIRALIB_DUMP_ARRAYS is set, the first measurement dumps the
 *     live-out arrays to stderr in PolyBench's POLYBENCH_DUMP_ARRAYS
 *     format (for correctness checking against PolyBench/C).
 */
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

int main(int argc, char **argv) {
    /* Child mode: one measurement in a fresh process, then exit. */
    if (argc > 1 && std::strcmp(argv[1], "--single-run") == 0) {
        int dump = (argc > 2 && std::strcmp(argv[2], "--dump") == 0);
        tpb_single_measurement(dump);
        return 0;
    }

    int n_runs = 1;
    const char *min_runs_env = std::getenv("MIN_RUNS");
    if (min_runs_env != NULL)
        n_runs = std::atoi(min_runs_env);
    if (argc > 1)
        n_runs = std::atoi(argv[1]);
    if (n_runs <= 0) {
        std::fprintf(stderr, "usage: %s [n_measurements]\n", argv[0]);
        return 2;
    }
    int dump_requested = (std::getenv("TIRALIB_DUMP_ARRAYS") != NULL);

    for (int i = 0; i < n_runs; ++i) {
        pid_t pid = fork();
        if (pid < 0) {
            std::perror("polybench-tiramisu wrapper: fork");
            return 1;
        }
        if (pid == 0) {
            if (dump_requested && i == 0) {
                execl(argv[0], argv[0], "--single-run", "--dump",
                      (char *)NULL);
                execl("/proc/self/exe", argv[0], "--single-run", "--dump",
                      (char *)NULL);
            } else {
                execl(argv[0], argv[0], "--single-run", (char *)NULL);
                execl("/proc/self/exe", argv[0], "--single-run",
                      (char *)NULL);
            }
            std::perror("polybench-tiramisu wrapper: execl");
            _exit(127);
        }
        int status = 0;
        while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {
        }
        if (WIFSIGNALED(status)) {
            std::fprintf(stderr,
                         "polybench-tiramisu wrapper: measurement killed by "
                         "signal %d\n",
                         WTERMSIG(status));
            return 1;
        }
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0)
            return WEXITSTATUS(status);
    }
    std::printf("\n");
    return 0;
}
