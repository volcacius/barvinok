#ifndef BARVINOK_OPTIONS_H
#define BARVINOK_OPTIONS_H

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options {
    /* PolyLib options */
    unsigned	MaxRays;

    /* lexmin options */
		/* Check for integer point in domain
		 * 0: no
		 * 1: yes
		 */
    int		emptiness_check;
};

struct barvinok_options *barvinok_options_new_with_defaults();

#if defined(__cplusplus)
}
#endif

#endif
