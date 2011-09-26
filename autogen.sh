#!/bin/sh
autoreconf -i --no-recursive
if test -f polylib/autogen.sh; then
	(cd polylib; ./autogen.sh)
fi
if test -f isl/autogen.sh; then
	(cd isl; ./autogen.sh)
fi
if test -f pet/autogen.sh; then
	(cd pet; ./autogen.sh)
fi
if test -f isl-polylib/autogen.sh; then
	(cd isl-polylib; ./autogen.sh)
fi
if test -f cloog/autogen.sh; then
	(cd cloog; ./autogen.sh)
fi
