# Process this file with autoconf to produce a configure script.
AC_INIT(WashUDAQ, 1.1)
AC_CONFIG_SRCDIR(digittes/digiTES.cpp)
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])

AC_EXEEXT
AC_OBJEXT

# we use libtool, thus we need to include it
LT_INIT

#  This should track the version in
#  AM_INIT_AUTOMAKE it will be used to
#  define the version information in shared
#  images.



# Checks for programs.


AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AM_PROG_LEX
AC_PROG_YACC
AC_PROG_GCC_TRADITIONAL
AC_PROG_LIBTOOL

AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADER(time.h)
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h memory.h stddef.h stdlib.h strings.h string.h sys/file.h sys/ioctl.h sys/socket.h sys/time.h unistd.h values.h windows.h machine/hal_sysinfo.h sys/proc.h sys/mtio.h sys/mman.h sys/ipc.h sys/shm.h machine/param.h sys/un.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STAT
AC_HEADER_TIME
AC_STRUCT_TM
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SIGNAL

# Checks for library functions.
AC_FUNC_MMAP
AC_REPLACE_FUNCS([atexit dup2 ftruncate getpagesize memset putenv])
AC_REPLACE_FUNCS([strchr strcspn strerror strrchr bzero mkstemp stat]) 
AC_REPLACE_FUNCS([malloc memcmp strftime alloca vfork strtok_r])
AC_CHECK_FUNCS([select socket strspn fork drand48 seed48]) 
AC_CHECK_FUNCS([pipe sbrk shm_open shmget fcntl])
AC_CHECK_FUNCS([floor getcwd gethostname gettimeofday mbrlen])
AC_CHECK_FUNCS(time)
AC_CHECK_LIB(m, pow)
AC_CHECK_LIB(m, sqrt)
AC_CHECK_LIB(m, logb)
AC_CHECK_FUNCS([strdup strstr strtol tzset mktime strtod timezone])






AM_PATH_CPPUNIT
AM_CONDITIONAL([CPPUNIT], [test -n "${CPPUNIT_LIBS}"] )
AM_CONDITIONAL([CPPUNITDONOTBUILD], [test -n ""])
AM_CONDITIONAL([THIS_IS_NOT_DEFINED], [test -n ""])

#
#   We'd like to compile with -pedantic if supported
#   Intel looks like gcc so we check explicitly for a gcc
#   named compiler
#


AC_SUBST(CFLAGS)

#  Build and install iniparser

(cd iniparser; make clean default)
(cd iniparser; prefix=${prefix} ./install)




# create output files.


AC_OUTPUT([Makefile
	include/Makefile
	lib/Makefile
	digittes/Makefile
	readout/Makefile
	 ])
