AC_DEFUN([FLA_REQUIRE_FC_VENDOR],
[
	AC_REQUIRE([FLA_REQUIRE_FC])

	dnl Ascertain the compiler "vendor".
	FC_VENDOR=$(echo "$FC" | egrep -o 'gfortran|f77|g77|xlf|frt|pgf77|cf77|fort77|fl32|af77|xlf90|f90|pgf90|pghpf|epcf90|g95|xlf95|f95|fort|ifort|ifc|efc|pgfortran|pgf95|lf95|ftn|nagfor' | { read first rest ; echo $first ; })

	if test "${FC_VENDOR}" = "" ; then
		AC_MSG_ERROR([configure can't determine the compiler vendor for $FC. Please submit a bug report to the FLAME developers.])
	else
		AC_MSG_NOTICE([[setting the FC_VENDOR environment variable to ${FC_VENDOR}.]])
	fi

	dnl Substitute the user-defined FFLAGS into the autoconf output files.
	AC_SUBST(FC_VENDOR)
])
