#!/bin/bash
# Script to format benchmarks results into a single xml file.
# See https://wiki.jenkins-ci.org/display/JENKINS/PerfPublisher+Plugin
# -----
# 2014/11/17 - Written by AB <Alexis.Breust@imag.fr>

XMLFILE=$1
benchmarks=$2

#========#
# Header #
#========#

echo '<?xml version="1.0" encoding="UTF-8"?>' >> $XMLFILE
echo '<report name="benchmarks-report" categ="benchmarks">' >> $XMLFILE

#=======#
# Start #
#=======#

echo '<start>' >> $XMLFILE
echo '<date format="YYYYMMDD" val="'$(date +%Y%m%d)'" />' >> $XMLFILE
echo '<time format="HHMMSS" val="'$(date +%H%M%S)'" />' >> $XMLFILE
echo '</start>' >> $XMLFILE

#============#
# Benchmarks #
#============#

for benchmark in $benchmarks
do
	echo '[Compiling]' $benchmark

	COMPILESTART=$(date +%s%3N)
	COMPILELOG=$(make $benchmark 2>&1; echo 'Returned state: '$?)
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
		PERFORMANCEFLOPS='0.0'
		COMPILETIMERELEVANT='false'
		EXECUTIONTIMERELEVANT='false'
		PERFORMANCEFLOPSRELEVANT='false'
		ERRORLOG='Does not compile.'
		echo '-> Does not compile.'
	else
		#Compilation success
		EXECUTED='yes'
		COMPILETIMERELEVANT='true'

		echo '[Executing]' $benchmark

		EXECUTIONLOG=$(./$benchmark 2>&1)
		if [[ $EXECUTIONLOG != "Time:"* ]]
		then
			#Execution failure
			PASSED='no'
			STATE='0'
			EXECUTIONTIME='0.0'
			PERFORMANCEFLOPS='0.0'
			EXECUTIONTIMERELEVANT='false'
			PERFORMANCEFLOPSRELEVANT='false'
			ERRORLOG='Unexpected output.'
			echo '-> Unexpected output.'
		else
			#Execution success
			PASSED='yes'
			STATE='100'
			EXECUTIONTIME=$(echo $EXECUTIONLOG | cut -d' ' -f2)
			PERFORMANCEFLOPS=$(echo $EXECUTIONLOG | cut -d' ' -f4)
			EXECUTIONTIMERELEVANT='true'
			if [[ $PERFORMANCEFLOPS != "Irrelevant" ]]
			then
				PERFORMANCEFLOPSRELEVANT='true'
			else
				PERFORMANCEFLOPSRELEVANT='false'
				PERFORMANCEFLOPS='0.0'
			fi
			ERRORLOG=''
		fi
	fi

	echo '<test name="'$benchmark'" executed="'$EXECUTED'">' >> $XMLFILE
	echo '<targets><target>BENCHMARK</target></targets>' >> $XMLFILE
	echo '<result>' >> $XMLFILE
	
	# Logs
	echo '<success passed="'$PASSED'" state="'$STATE'" />' >> $XMLFILE
	echo '<errorlog><![CDATA['$ERRORLOG']]></errorlog>' >> $XMLFILE
	echo '<log name="Compile output"><![CDATA['"$COMPILELOG"']]></log>' >> $XMLFILE
	echo '<log name="Execution output"><![CDATA['"$EXECUTIONLOG"']]></log>' >> $XMLFILE
	
	# Times
	echo '<compiletime unit="ms" mesure="'$COMPILETIME'" isRelevant="'$COMPILETIMERELEVANT'" />' >> $XMLFILE
	echo '<executiontime unit="s" mesure="'$EXECUTIONTIME'" isRelevant="'$EXECUTIONTIMERELEVANT'" />' >> $XMLFILE
	echo '<performance unit="GFLOPS" mesure="'$PERFORMANCEFLOPS'" isRelevant="'$PERFORMANCEFLOPSRELEVANT'" />' >> $XMLFILE
	
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

