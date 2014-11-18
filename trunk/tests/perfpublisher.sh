#!/bin/bash
# Script to format tests results into a single xml file.
# See https://wiki.jenkins-ci.org/display/JENKINS/PerfPublisher+Plugin
# -----
# 2014/11/17 - Written by AB <Alexis.Breust@imag.fr>

XMLFILE=$1
tests=$2

#========#
# Header #
#========#

echo '<?xml version="1.0" encoding="UTF-8"?>' >> $XMLFILE
echo '<report name="tests-report" categ="tests">' >> $XMLFILE

#=======#
# Start #
#=======#

echo '<start>' >> $XMLFILE
echo '<date format="YYYYMMDD" val="'$(date +%Y%m%d)'" />' >> $XMLFILE
echo '<time format="HHMMSS" val="'$(date +%H%M%S)'" />' >> $XMLFILE
echo '</start>' >> $XMLFILE

#=======#
# Tests #
#=======#

for test in $tests
do
	echo '[Compiling]' $test

	COMPILESTART=$(date +%s%3N)
	COMPILELOG=$(make $test 2>&1; echo 'Returned state: '$?)
	COMPILEEND=$(date +%s%3N)
	COMPILETIME=$(($COMPILEEND - $COMPILESTART))
	COMPILECHECK=$(echo $COMPILELOG | grep -o '[^ ]*$')
	
	if [[ $COMPILECHECK -ne 0 ]]
	then
		#Compilation failure
		EXECUTED='no'
		PASSED='no'
		STATE='0'
		EXECUTIONLOG=''
		EXECUTIONTIME='0.0'
		COMPILETIMERELEVANT='false'
		EXECUTIONTIMERELEVANT='false'
		ERRORLOG='Does not compile.'
		echo '-> Does not compile.'
	else
		#Compilation success
		EXECUTED='yes'
		COMPILETIMERELEVANT='true'

		echo '[Executing]' $test

		EXECUTIONSTART=$(date +%s%3N)
		EXECUTIONLOG=$(./$test  2>&1; echo 'Returned state: '$?)
		EXECUTIONEND=$(date +%s%3N)
		EXECUTIONTIME=$(($EXECUTIONEND - $EXECUTIONSTART))
		EXECUTIONCHECK=$(echo $EXECUTIONLOG | grep -o '[^ ]*$')
		if [[ $EXECUTIONCHECK -ne 0 ]]
		then
			#Execution failure
			PASSED='no'
			STATE='0'
			EXECUTIONTIMERELEVANT='false'
			ERRORLOG='Execution failure.'
			echo '-> Execution failure.'
		else
			#Execution success
			PASSED='yes'
			STATE='100'
			EXECUTIONTIMERELEVANT='true'
			ERRORLOG=''
		fi
	fi

	echo '<test name="'$test'" executed="'$EXECUTED'">' >> $XMLFILE
	echo '<targets><target>TEST</target></targets>' >> $XMLFILE
	echo '<result>' >> $XMLFILE
	
	# Logs
	echo '<success passed="'$PASSED'" state="'$STATE'" />' >> $XMLFILE
	echo '<errorlog><![CDATA['$ERRORLOG']]></errorlog>' >> $XMLFILE
	echo '<log name="Compile output"><![CDATA['"$COMPILELOG"']]></log>' >> $XMLFILE
	echo '<log name="Execution output"><![CDATA['"$EXECUTIONLOG"']]></log>' >> $XMLFILE
	
	# Times
	echo '<compiletime unit="ms" mesure="'$COMPILETIME'" isRelevant="'$COMPILETIMERELEVANT'" />' >> $XMLFILE
	echo '<executiontime unit="ms" mesure="'$EXECUTIONTIME'" isRelevant="'$EXECUTIONTIMERELEVANT'" />' >> $XMLFILE
	
	echo '</result>' >> $XMLFILE
	echo '</test>' >> $XMLFILE
done

#========#
# Footer #
#========#

echo '</report>' >> $XMLFILE

#==========#
# Epilogue #
#==========#

echo 'Results correctly exported to' $XMLFILE

