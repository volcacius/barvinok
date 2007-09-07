#!/bin/sh
libtoolize -c --force
aclocal -I m4
autoheader
automake -a -c --foreign
autoconf
if test -f polylib/autogen.sh; then
	(cd polylib; ./autogen.sh)
fi
if test -f piplib/autogen.sh; then
	(cd piplib; ./autogen.sh)
fi
(cd bernstein; ./autogen.sh)
