#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include <util.h>
#include <barvinok.h>

static void time_diff(struct tms *before, struct tms *after)
{
    long ticks = sysconf(_SC_CLK_TCK);
    printf("User: %g; Sys: %g\n", 
	    (0.0 + after->tms_utime - before->tms_utime) / ticks,
	    (0.0 + after->tms_stime - before->tms_stime) / ticks);
}

int main()
{
    Value cm, cb;
    struct tms tms_before, tms_between, tms_after;
    Polyhedron *A;
    Matrix *M;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, 600);
    Matrix_Free(M);
    value_init(cm);
    value_init(cb);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    times(&tms_before);
    manual_count(A, &cm);
    times(&tms_between);
    barvinok_count(A, &cb);
    times(&tms_after);
    printf("manual: ");
    value_print(stdout, P_VALUE_FMT, cm);
    puts("");
    time_diff(&tms_before, &tms_between);
    printf("Barvinok: ");
    value_print(stdout, P_VALUE_FMT, cb);
    puts("");
    time_diff(&tms_between, &tms_after);
    value_clear(cm);
    value_clear(cb);
    Polyhedron_Free(A);
    return 0;
}
