noinst_PROGRAMS = @bv_omega_programs@
EXTRA_PROGRAMS = occ

occ_SOURCES = \
	Exit.cc \
	lex.yy.cc \
	y.tab.cc
LDADD = -lcode_gen -lomega
LDFLAGS = @OMEGA_LDFLAGS@
CPPFLAGS = @OMEGA_CPPFLAGS@

YACC_FLAGS = -d

y.tab.hh y.tab.cc: parser.y
	$(YACC) -o y.tab.cc -v -d $<

lex.yy.cc: parser.l y.tab.hh
	$(LEX) -olex.yy.cc -i $<

