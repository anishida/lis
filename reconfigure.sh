#!/bin/sh

# Shell script for the GNU Build System
# To remake ./configure, edit ./configure.in and run the following command:
#
# > ./reconfigure.sh
#

# Check that the utility is present.

type autoreconf >/dev/null 2>&1 || { echo >&2 "You need autoreconf to remake configure script. Aborting."; exit 1; }

# Run `autoconf' (and `autoheader', `aclocal', `automake', `autopoint', 
# and `libtoolize' where appropriate) repeatedly to remake 
# the GNU Build System files.

autoreconf
