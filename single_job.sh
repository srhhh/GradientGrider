#!/bin/bash

echo ""
echo ""
echo "Now running the background job"
echo ""

###############################################################################################################################################
###############################################################################################################################################

#Enter the name of the original source code files

#The parameters file
oldPARAMETERS="PARAMETERS"

#The analysis file
oldANALYSIS="ANALYSIS"

#The physics file
oldPHYSICS="PHYSICS"

#The variables file
oldVARIABLES="VARIABLES"

#The makefile for grid creation
oldMAKEGRID="make_makeGridwithNewTrajectories2"

#The makefile for grid analysis
oldMAKEANALYSIS="make_checkNewTrajectorieswithMultipleGrids2"

###############################################################################################################################################
###############################################################################################################################################

#Enter the name of the source code files to be made in the new directory
#These must be different from the original source code file names

#The parameters file
newPARAMETERS=PARAMETERS_new

#The analysis file
newANALYSIS=ANALYSIS_new

#The makefile for grid creation
newMAKEGRID=make_makeGridwithNewTrajectories2_new

#The makefile for grid analysis
newMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids2_new

###############################################################################################################################################
###############################################################################################################################################

#Change the below variables for whatever grid you want
#If you want to change other variables (ex. spacing1) you would
#have to make more sed statements below to change that

#How many children generations do you want?
Norder_max=1

#The ratio of child-level to parent-level cell spacings with respect to var1
scaling1_0=004
scaling1_1=004

#The ratio of child-level to parent-level cell spacings with respect to var2
scaling2_0=004
scaling2_1=004

#The number of frames a parent-level cell accepts before it is subdivided
overcrowd0=10050
overcrowd1=10000
overcrowd2=01010

#The number of trajectories simulated and added to a new grid
Ntraj_max=0050

#The number of grids to add to a new library
Ngrid_max=1

#Whether to add duplicate copies or use a labelling scheme
force_Duplicates=.false.

#Whether to use the labelling scheme or drop it altogether
force_NoLabels=.true.

#The default flags to be used for analyses
#Of course, you don't want all analyses to be the same so go down to each analysis and change
#what you want each individual one to do
heatmap_flag=.false.
trueSA_flag=.false.
trueED_flag=.false.
testtraj_flag=.true.
testtrajRMSD_flag=.false.
percentthreshold_flag=.true.
testtrajSA_flag=.false.
testtrajSAheatmap_flag=.false.
Nsort=1
outer_threshold=.100000d0
outer_threshold1=.100000d0
outer_threshold2=.100000d0
outer_threshold3=.100000d0
outer_threshold4=.100000d0
outer_threshold5=.050000d0
inner_threshold="0.00000d0"
inner_threshold1="0.00000d0"
inner_threshold2="0.00000d0"
inner_threshold3="0.00000d0"
inner_threshold4="0.00000d0"
inner_threshold5="0.00000d0"
alpha_ratio="1.0d0"
alpha_ratio1="1.0d-2"
alpha_ratio2="1.0d0"
alpha_ratio3="1.0d-2"
alpha_ratio4="1.0d0"
alpha_ratio5="1.0d-2"
R1_threshold="3.0000d-6"
R2_threshold="1.0000d-3"
force_Permutations=.true.
interpolation_flag=.true.
gather_interpolation_flag=.true.
reject_flag=.true.
accept_first=.false.
accept_worst=.false.
grid_addition=1
Ngrid_cap=1
Norder_cap=1
#Ngrid_cap=${Ngrid_max}
Ntrajectories=5
Naccept_max=05
Nthreads=1

#Names of the experiments
#exp1name=RMSD100_7trajtest_DEC02_empirical
exp1name=exp001
exp2name=exp101
exp3name=exp312
exp4name=exp013
exp5name=exp014

#This flag states whether we are continuing an old experiment
continue_analysis=.true.

#If we want to use a fixed set of initial conditions,
#specify which experiment they come from here
useoldinitialbonddata_flag=.false.
#initialbondfolder="001reject.10000"
initialbondfolder=exp051/

#If we are using frames from a trajectory then list this
#true and read from this file
readtrajectory_flag=.true.
readtrajectoryfile="readtrajectories_allPT.txt"
readtrajectoryfolder="processed_traj"

#readtrajectoryfile1="readtrajectories_single.txt"
#readtrajectoryfolder1="traj_single"
#readtrajectoryfile1="readtrajectories_fewlongPT.txt"
readtrajectoryfile1="readtrajectories_somePT.txt"
readtrajectoryfolder1="some_processed_traj"
readtrajectoryfile2="readtrajectories2.txt"
readtrajectoryfolder2="traj2"
readtrajectoryfile3="readtrajectories3.txt"
readtrajectoryfolder3="traj3"
readtrajectoryfile4="readtrajectories4.txt"
readtrajectoryfolder4="traj4"
readtrajectoryfile5="readtrajectories5.txt"
readtrajectoryfolder6="traj5"

#If you have special set of parameters you want to compare, list them here
#These will be compared at each compilation

#These are the options for comparison:
#ScatteringAngle (rad)
#AbsoluteEnergyChange (eV)
#RelativeEnergyChange (eV)
#RotationalEnergyChange (eV)

#If the comparison lower and upper limits are the same, the program will
#use whatever the minimum and maximum is of the data (bad if outliers exist)
#comparison_flag=notScatteringAngle
#comparison_lowerlimit="0.0d0"
#comparison_upperlimit="0.2500d0"
#comparison_upperlimit="0.0180d0"

#declare -a prefixes
#prefixes[0]="001accept.15000"
#prefixes[1]="002accept.15000"
#prefixes[2]="004accept.15000"
#prefixes[3]="008accept.15000"
#prefixes[4]="001reject.15000"
#prefixes[0]="004accept.05000"
#prefixes[1]="004accept.01000"
#prefixes[2]="004accept.00500"
#prefixes[0]="001omegaA.50005"
#prefixes[1]="001omegaA.37505"
#prefixes[2]="001omegaA.25005"
#prefixes[3]="001omegaA.12505"
#prefixes[4]="001reject.05000"
#prefixes[0]="exp002"
#prefixes[1]="exp003"
#prefixes[2]="exp003"

###############################################################################################################################################
###############################################################################################################################################

#Enter what path the original source code is on and what you would like
#the library housing all the grids (a folder) to be called

#The name of the new library (folder)
#newGRID=HH_${scaling1_0}_${scaling2_0}_${overcrowd0}_${Ntraj_max}_1
#newGRID="HBrCO2_Sep25"
#newGRID="HBrCO2_Jul30"
newGRID="HBrCO2_Apr1"

#If you want to make a new grid, set this to 1; otherwise, set it to zero
newGRID_flag=0

#How often you want to check the progress of the new grid's creation
#(has an intrinsic minimum of Ntraj_max/10)
#Just set this to a very large number if no progress checks are wanted
newGRID_check_min=10000

#The number of post-grid analyses you would like done
#These are separate from the comparison and the post-grid-making analysis
Nanalyses=1

#The path that has the original source code
currentPATH=$(pwd)/src

#The file that keeps the times of the code
bashout="timing.dat"

#DO NOT TOUCH THESE
gridPATH=$(pwd)/$newGRID
newSOURCE=SOURCE
newPATH=$gridPATH/$newSOURCE

###############################################################################################################################################
###############################################################################################################################################

#This job accepts may accept one argument.
#This argument must be the name of:
# 1) a comparison analysis file
# 2) an interpolation analysis file

#Temporary file for additonal information about the experiments
comparisonsfile=comparison.txt

if [ "$1" != "" ]; then
	#The first line differentiates between 1 and 2
	comparison_flag=$(head -n 1 "$1")

	#This is the option for interpolation:
	#InterpolationTDD (natural numbers, starts at 0)
	#InterpolationCED (real numbers, centered at 1)

        #These are the options for comparison:
        #ScatteringAngle (rad)
        #AbsoluteEnergyChange (eV)
        #RelativeEnergyChange (eV)
        #RotationalEnergyChange (eV)

#       if [ $comparison_flag == "ScatteringAngle" ] || [ $comparison_flag == "AbsoluteEnergyChange" ] || [ $comparison_flag == "RelativeEnergyChange" ] || [ $comparison_flag == "RotationalEnergyChange" ]
#then
#	gather_interpolation_flag=.false.
#elif [ $comparison_flag == "InterpolationTDD" ] || [ $comparison_flag == "InterpolationRED" ] || [ $comparison_flag == "InterpolationAED" ] || [ $comparison_flag == "InterpolationIED" ] || [ $comparison_flag == "InterpolationR1D" ] || [ $comparison_flag == "InterpolationRSV1D" ] || [ $comparison_flag == "InterpolationRSV2D" ]
#       then
#	gather_interpolation_flag=.true.
#else
#	echo "" 
#	echo "Argument supplied does not have a first-line option that is offered!" 
#	echo "" 
#	exit
#fi

	read -r comparison_lowerlimit comparison_upperlimit <<<$(sed -n '2p' "$1")
	read -r Ntrajectories <<<$(sed -n '3p' "$1")

#decimal_match="^[\$0-9][\$0-9]*\.[\$0-9][\$0-9]*[eEdD]-[\$0-9]*$"
	decimal_match='^[0-9][0-9]*\.[0-9][0-9]*[eEdD][-+]?[0-9]*$'
#integer_match="^[\$0-9][\$0-9]*$"
	integer_match='^[\$0-9][\$0-9]*$'

	if [[ $comparison_lowerlimit =~ $decimal_match ]] ||  [[ $comparison_upperlimit =~ $decimal_match ]]
	then
		:
	else
		echo "" 
		echo "Argument supplied does not have a valid second-line option!" 
		echo "   Must be two positive decimals separated by whitespace" 
		echo "" 
		exit
	fi

	if [[ $Ntrajectories =~ $integer_match ]]
	then
		:
	else
		echo "" 
		echo "Argument supplied does not have a valid third-line option!" 
		echo "   Must be a whole number" 
		echo "" 
		exit
	fi

	comparison_lowerlimit=$comparison_lowerlimit
	comparison_upperlimit=$comparison_upperlimit

	rm -f $gridPATH/$comparisonsfile
        declare -a prefixes
	while read -r aline
        do
		read -r expname extra <<<$aline
		if [ ! -d $gridPATH/$expname ]; then
        		echo "" 
        		echo "Argument supplied does not have a valid experiment for comparison!"
        		echo "   The folder \"$expname\" was not found" 
        		echo "" 
        		exit
		fi
		prefixes[${#prefixes[@]}]="$expname"/
		echo "$extra" >>$gridPATH/$comparisonsfile
	done<<<"$(tail -n +4 $1)"

	if [ ${#prefixes[@]} == 0 ]; then
		echo "" 
		echo "Argument supplied does not have any experiments for comparison!"
		echo "" 
		exit
        fi
else
	comparison_flag=""
	comparison_lowerlimit="0.0d0"
	comparison_upperlimit="0.0d0"
fi

###############################################################################################################################################
###############################################################################################################################################
#		GRID CREATION AND SUMMARY ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#If you want to create a new grid, do that first
if [ $newGRID_flag -eq 1 ]
then

#If there is another folder of the same name delete that folder first
rm -r $gridPATH/

#Make the directories, copy all the original source code, etc.
mkdir $gridPATH/
mkdir $newPATH/
cp $currentPATH/*.f90 $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/

#temporarily do this
cp -r $(pwd)/dropoffs/ $gridPATH/

mkdir $gridPATH/startup/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile $gridPATH/startup/readtrajectories.txt
mkdir $gridPATH/startup/traj/
cp $readtrajectoryfolder/* $gridPATH/startup/traj/

fi

#Make changes to the parameters file as specified in the variables above
#Unless you want to change MORE variables, don't touch this
sed "s|Ntraj_max = [0-9]*|Ntraj_max = $Ntraj_max| 
     s|Norder_max = [0-9]*|Norder_max = $Norder_max|
     s|Ngrid_max = [0-9]*|Ngrid_max = $Ngrid_max|
     s|trajectory_text_length = [0-9]*|trajectory_text_length = ${#Ntraj_max}|
     s|gridpath_length = .*|gridpath_length = $((${#gridPATH}+1))|
     s|Ngrid_check_min = .*|Ngrid_check_min = $newGRID_check_min|
     s|force_Duplicates = .*|force_Duplicates = $force_Duplicates|
     s|force_NoLabels = .*|force_NoLabels = $force_NoLabels|
     s|gridpath0 = .*|gridpath0 = \\&\\n\"$gridPATH/\"|
     s|character([0-9]*),parameter :: parametersfile = .*|character($((${#newPARAMETERS}+4))),parameter :: parametersfile = \"$newPARAMETERS.f90\"|" <$currentPATH/$oldPARAMETERS.f90 >$newPATH/$newPARAMETERS.f90

longtext="/default_scaling/{N;N;N;N;s|\\(default_scaling\\s=\\s&\\n\\s*reshape(\\s&\\)"
longtext="$longtext""\\n.*\\n.*\\n"
longtext="$longtext""|\\1\\n"
longtext="$longtext""                (/ $scaling1_0, $scaling2_0, 4, 4, \\&\\n"
longtext="$longtext""                   $scaling1_1, $scaling2_1, 10, 10, \\&\\n"
longtext="$longtext""|}"

sed -i "$longtext" $newPATH/$newPARAMETERS.f90

longtext="/default_overcrowd/{N;s|\\(default_overcrowd\\s=\\s&\\)"
longtext="$longtext""\\n.*"
longtext="$longtext""|\\1\\n"
longtext="$longtext""        (/ $overcrowd0, $overcrowd1, $overcrowd2, 50 /)"
longtext="$longtext""|}"

sed -i "$longtext" $newPATH/$newPARAMETERS.f90

#Make changes to the analysis file as specified in the variables above
#This analysis file controls any analysis done during grid creation:
#particularly things like rejecting frames and not checking too many cells.
#
#It also controls the first post-grid-making analysis for the trajectories
#in the grid that was just made (heatmap, scattering angle, etc)
sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s|Norder_cap = [0-9]*|Norder_cap = $Norder_cap|
     s/heatmap_flag = .*/heatmap_flag = .true./
     s/trueSA_flag = .*/trueSA_flag = .true./
     s/trueED_flag = .*/trueED_flag = .true./
     s/testtraj_flag = .*/testtraj_flag = .false./
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/comparison_flag = .*/comparison_flag = .false./
     s/percentthreshold_flag = .*/percentthreshold_flag = .false./
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = .true./
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = .false./
     s/testtrajSA_flag = .*/testtrajSA_flag = .false./
     s|expfolder_length = .*|expfolder_length = 8|
     s|expfolder = .*|expfolder = \"startup/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/startup/$newANALYSIS.f90

#DO NOT TOUCH THIS
sed "s/$oldPARAMETERS\\.o/$newPARAMETERS.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/startup/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEGRID >$newPATH/$newMAKEGRID

#DO NOT TOUCH THIS
sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/startup/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS



#Now, all we need to do is go into the new folder and make the output files
cd $newPATH

make -f $newPATH/$newMAKEGRID

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "Grid Creation" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "GRIDMAKING %E  %U  %S  %P  %O" ./a.out
##valgrind --leak-check=yes ./a.out

exit

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "Post Analysis" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS00 %E  %U  %S  %P  %O" ./a.out

fi

###############################################################################################################################################
###############################################################################################################################################


###############################################################################################################################################
###############################################################################################################################################
#		COMPARISON (USUALLY DONE AFTER ANALYSES)
###############################################################################################################################################
###############################################################################################################################################

#This prevents a new parameters file from being copied into the library
#While at the same time, updating all of the .f90 files
#Remember, the parameters file should never change after creation!
shopt -s extglob
cp $currentPATH/!($oldPARAMETERS|$newPARAMETERS|$oldPHYSICS|$oldVARIABLES|$oldANALYSIS|$newANALYSIS)+(.f90) $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/
shopt -s extglob

#Whatever trajectories we want to check, we need to accumulate and pass
#that to the fortran code somehow
allprefixes=""
alllengths=()
for prefix in "${prefixes[@]}"
do
      allprefixes=$allprefixes$prefix
      alllengths+=("${#prefix}")
done
alllengths_statement=$(IFS=, ; echo "${alllengths[*]}")

#if [ $comparison_flag == "ScatteringAngle" ] || [ $comparison_flag == "AbsoluteEnergyChange" ] || [ $comparison_flag == "RelativeEnergyChange" ] || [ $comparison_flag == "RotationalEnergyChange" ] || [ $comparison_flag == "InterpolationTDD" ] || [ $comparison_flag == "InterpolationRED" ] || [ $comparison_flag == "InterpolationAED" ] || [ $comparison_flag == "InterpolationIED" ] || [ $comparison_flag == "InterpolationR1D" ] || [ $comparison_flag == "InterpolationRSV1D" ] || [ $comparison_flag == "InterpolationRSV2D" ]
#then

if [ "$comparison_flag" != "" ]
then

rm -r $gridPATH/comparison/
mkdir $gridPATH/comparison/

cp $1 $gridPATH/comparison/

#Change the comparison analysis as specified in the variables above
sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = .false./
     s/trueED_flag = .*/trueED_flag = .false./
     s/testtraj_flag = .*/testtraj_flag = .false./
     s/continue_analysis = .*/continue_analysis = .false./
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = .false./
     s/percentthreshold_flag = .*/percentthreshold_flag = .false./
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio/
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .true./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/comparison_SATRVname = .*/comparison_SATRVname = \"$comparison_flag\"/
     s/comparison_number = .*/comparison_number = ${#prefixes[@]}/
     s|character([^)]*),parameter :: allprefixes = .*|character(${#allprefixes}),parameter :: allprefixes = \"$allprefixes\"|
     s|alllengths = .*|alllengths = (/$alllengths_statement/)|
     s/testheatmapSA_flag = .*/testheatmapSA_flag = .false./
     s/testtrajSA_flag = .*/testtrajSA_flag = .false./
     s|expfolder_length = .*|expfolder_length = 11|
     s|expfolder = .*|expfolder = \"comparison/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/comparison/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/comparison/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newPATH

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "Comparison" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "COMPARISON %E  %U  %S  %P  %O" ./a.out

rm $gridPATH/$comparisonsfile

fi

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 1 ]
then
exit
fi

#Now, we come to the part of the code which may run in parallel
#We have to check if the user input a correct number of threads

#Is it a positive integer? (made entirely of digits?)

if [[ "$Nthreads" =~ ^[0-9]+$ ]]
then
:
else
echo ""
echo "WARNING"
echo "Number of threads ($Nthreads) is not a positive integer."
echo "Exiting program before compiling analysis."
echo ""
exit
fi

#Is it greater than zero and less than or equal to the total number of cores?

Ncores=$(nproc --all)
if [[ $Nthreads -gt 0 && $Nthreads -le $Ncores ]]
then
:
else
echo ""
echo "WARNING"
echo "Number of threads ($Nthreads) may not be compatible on this system with $Ncores cores."
echo "Exiting program before compiling analysis."
echo ""
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FIRST ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#This prevents a new parameters file from being copied into the library
#While at the same time, updating all of the .f90 files
#Remember, the parameters file should never change after creation!
shopt -s extglob
cp $currentPATH/!($oldPARAMETERS|$newPARAMETERS|$oldPHYSICS|$oldVARIABLES|$oldANALYSIS|$newANALYSIS)+(.f90) $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/
shopt -s extglob

#rm -r $gridPATH/$exp1name/
mkdir -p $gridPATH/$exp1name/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile1 $gridPATH/$exp1name/readtrajectories.txt
mkdir $gridPATH/$exp1name/traj/
cp $readtrajectoryfolder1/* $gridPATH/$exp1name/traj/

fi

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/Nthreads = [0-9]*/Nthreads = $Nthreads/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold1/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold1/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio1/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .false./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/
     s|expfolder_length = .*|expfolder_length = $((${#exp1name}+1))|
     s|expfolder = .*|expfolder = \"$exp1name/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/$exp1name/$newANALYSIS.f90


sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/$exp1name/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newPATH

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "$Ntrajectories Trajectories Check ($exp1name)" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS01 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 2 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		SECOND ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#rm -r $gridPATH/$exp2name/
mkdir -p $gridPATH/$exp2name/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile2 $gridPATH/$exp2name/readtrajectories.txt
mkdir $gridPATH/$exp2name/traj/
cp $readtrajectoryfolder2/* $gridPATH/$exp2name/traj/

fi

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/Nthreads = [0-9]*/Nthreads = $Nthreads/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold2/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold2/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio2/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .false./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/
     s|expfolder_length = .*|expfolder_length = $((${#exp2name}+1))|
     s|expfolder = .*|expfolder = \"$exp2name/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/$exp2name/$newANALYSIS.f90


sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/$exp2name/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "$Ntrajectories Trajectories Check ($exp2name)" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS02 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 3 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		THIRD ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#rm -r $gridPATH/$exp3name/
mkdir -p $gridPATH/$exp3name/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile3 $gridPATH/$exp3name/readtrajectories.txt
mkdir $gridPATH/$exp3name/traj/
cp $readtrajectoryfolder3/* $gridPATH/$exp3name/traj/

fi

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/Nthreads = [0-9]*/Nthreads = $Nthreads/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold3/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold3/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio3/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .false./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/
     s|expfolder_length = .*|expfolder_length = $((${#exp3name}+1))|
     s|expfolder = .*|expfolder = \"$exp3name/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/$exp3name/$newANALYSIS.f90


sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/$exp3name/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "$Ntrajectories Trajectories Check ($exp3name)" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS03 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 4 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FOURTH ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#rm -r $gridPATH/$exp4name/
mkdir -p $gridPATH/$exp4name/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile4 $gridPATH/$exp4name/readtrajectories.txt
mkdir $gridPATH/$exp4name/traj/
cp $readtrajectoryfolder4/* $gridPATH/$exp4name/traj/

fi

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/Nthreads = [0-9]*/Nthreads = $Nthreads/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold4/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold4/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio4/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .false./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/
     s|expfolder_length = .*|expfolder_length = $((${#exp4name}+1))|
     s|expfolder = .*|expfolder = \"$exp4name/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/$exp4name/$newANALYSIS.f90


sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/$exp4name/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "$Ntrajectories Trajectories Check ($exp4name)" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS04 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 5 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FIFTH ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#rm -r $gridPATH/$exp5name/
mkdir -p $gridPATH/$exp5name/

if [ $readtrajectory_flag == ".true." ]
then

cp $readtrajectoryfile5 $gridPATH/$exp5name/readtrajectories.txt
mkdir $gridPATH/$exp5name/traj/
cp $readtrajectoryfolder5/* $gridPATH/$exp5name/traj/

fi

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/continue_analysis = .*/continue_analysis = $continue_analysis/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondfolder_length = .*|initialbondfolder_length = $((${#initialbondfolder}))|
     s|initialbondfolder = .*|initialbondfolder = \"$initialbondfolder\"|
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/Naccept_max = .*/Naccept_max = $Naccept_max/
     s/Nthreads = [0-9]*/Nthreads = $Nthreads/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/outer_threshold_SI = .*/outer_threshold_SI = $outer_threshold5/
     s/inner_threshold_SI = .*/inner_threshold_SI = $inner_threshold5/
     s/R1_threshold_SI = .*/R1_threshold_SI = $R1_threshold/
     s/R2_threshold_SI = .*/R2_threshold_SI = $R2_threshold/
     s/Nsort = .*/Nsort = $Nsort/
     s/alpha_ratio = .*/alpha_ratio = $alpha_ratio5/
     s|force_Permutations = .*|force_Permutations = $force_Permutations|
     s/interpolation_flag = .*/interpolation_flag = $interpolation_flag/
     s/gather_interpolation_flag = .*/gather_interpolation_flag = $gather_interpolation_flag/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s|readtrajectory_flag = .*|readtrajectory_flag = $readtrajectory_flag|
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/grid_addition = .*/grid_addition = $grid_addition/
     s/comparison_flag = .*/comparison_flag = .false./
     s|character([^)]*),parameter :: comparison_file = .*|character(${#comparisonsfile}),parameter :: comparison_file = \"$comparisonsfile\"|
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/
     s|expfolder_length = .*|expfolder_length = $((${#exp5name}+1))|
     s|expfolder = .*|expfolder = \"$exp5name/\"|
     s|character([0-9]*),parameter :: analysisfile = .*|character($((${#newANALYSIS}+4))),parameter :: analysisfile = \"$newANALYSIS.f90\"|" <$currentPATH/$oldANALYSIS.f90 >$gridPATH/$exp5name/$newANALYSIS.f90


sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s|EXPSOURCE = .*|EXPSOURCE = $gridPATH/$exp5name/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo "" >> $newPATH/$bashout
echo $(date) >> $newPATH/$bashout
echo "$Ntrajectories Trajectories Check ($exp5name)" >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS05 %E  %U  %S  %P  %O" ./a.out
