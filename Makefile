# SYSTEM COMMANDS

CHMOD = chmod ug+rwx
ECHO = echo
RM = rm
RSYNC = rsync -auvz

# VARIABLES

BINDIR=../../bin
CONFFILES = *.conf
HTMLFILES = *.html
INIFILES = *.ini
PNGFILES = *.png
TXTFILES=*.txt
PYPROGS = *.py
RPROGS = *.R
SHPROGS = *.sh
CPROGS =
CPPPROGS =
HPLUSPLUSFILES = *.h
OBJS=*.o
INTERMEDIATES = *.pyc

# Compiling C programs

CC=gcc
CDEBUGFLAGS=-g -Wall
COPTFLAGS=-O3
CLDFLAGS=
CLIBINC=
CZLIB=-lz
CMATHS=-lm
CPP=g++
CPPFLAGS=-O3

# COMMANDS

dummy:
	${ECHO} 'Please specify which command should be executed, e.g., make install'

permissions:
	cp ../admin/chmod_files.sh ../../bin ; ../../bin/chmod_files.sh

htmlfiles: ${HTMLFILES}
	${RSYNC} ${HTMLFILES} ${BINDIR} >/dev/null

inifiles: ${INIFILES}
	${RSYNC} ${INIFILES} ${BINDIR} >/dev/null

txtfiles: ${TXTFILES}
	${RSYNC} ${TXTFILES} ${BINDIR} >/dev/null

pngfiles: ${PNGFILES}
	${RSYNC} ${PNGFILES} ${BINDIR} >/dev/null

pyprogs: ${PYPROGS}
	${RSYNC} ${PYPROGS} ${BINDIR} >/dev/null

rprogs: ${RPROGS}
	${RSYNC} ${RPROGS} ${BINDIR} >/dev/null

shprogs: ${SHPROGS}
	${RSYNC} ${SHPROGS} ${BINDIR} >/dev/null

cprogs: ${CPROGS}
	${RSYNC} ${CPROGS} ${BINDIR} >/dev/null

cppprogs: $(CPPPROGS)
	${RSYNC} ${CPPPROGS} ${BINDIR} >/dev/null

install: pyprogs rprogs

clean:
	${RM} ${INTERMEDIATES} ${CPROGS} ${CPPPROGS}
