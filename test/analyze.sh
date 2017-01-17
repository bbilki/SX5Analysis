#!/bin/sh

scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/admin/RunList.txt Files/
wait
scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results/PlotList.txt Files/
wait

if [ -f Files/RunDiff.txt ]
then
	rm Files/RunDiff.txt
fi

if [ ! -d Data ]
then
	mkdir Data
fi

if [ ! -d Histos ]
then
	mkdir Histos
fi

if [ ! -d Plots ]
then
	mkdir Plots
fi

if [ ! -d NTuples ]
then
	mkdir NTuples
fi

runfilename="Files/RunList.txt"
plotfilename="Files/PlotList.txt"
rt=( 0 1 2 1 2 1 2 2 2 1 2 3 3 3 2 3 3 2 )
while read -r -a l
do
	if [ ${l[0]} -gt 1214 ]
	then
		if [ ! -f Data/SX5_${l[0]}.root ]
		then
			scp daq@cmshcal35.cern.ch:/data/spool/SX5_${l[0]}.root Data/
			wait
			r=0
			if [ ${l[2]} -ne -1 ]
			then
				r=${rt[${l[2]}]}
				if [ ! -d Plots/${l[1]} ]
				then
					mkdir Plots/${l[1]}
					wait
					mkdir Histos/${l[1]}
					wait
				fi
			fi
			if [ ${l[4]} -ne -1 ]
			then
				r=${rt[${l[4]}]}
				if [ ! -d Plots/${l[3]} ]
				then
					mkdir Plots/${l[3]}
					wait
					mkdir Histos/${l[3]}
					wait
				fi
			fi
			if [ ${l[6]} -ne -1 ]
			then
				r=${rt[${l[6]}]}
				if [ ! -d Plots/${l[5]} ]
				then
					mkdir Plots/${l[5]}
					wait
					mkdir Histos/${l[5]}
					wait
				fi
			fi
			if [ $r -gt 0 ]
			then
				echo "cmsRun sx5analyzer_cfg.py ${l[0]} $r"
				cmsRun sx5analyzer_cfg.py ${l[0]} $r
				wait
	# 			mv H_${l[0]}.root Histos/
				mv N_${l[0]}.root NTuples/
				wait
				echo ${l[0]} >> Files/RunDiff.txt
			fi
		fi
	fi
done < "$runfilename"

cd Files
./plotter 1
wait
cd ..
# 
# # while read -r -a l
# # do
# # 	if [ ${l[1]} -ne 17 ]
# # 	then
# # 		scp -r -P 53222 Plots/${l[0]} hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results
# # 		wait
# # 	fi
# # done < "$plotfilename"
# 
# 
# # scp -r -P 53222 Plots/* hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results
# # wait
# 
cd Plots
tar -cf plots.tar *
wait
gzip plots.tar
wait
scp -P 53222 plots.tar.gz hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results
wait
ssh -p 53222 hfSX5@feynman.physics.uiowa.edu "cd /var/www/html/HFPhase1SX5Ops/results/ ; tar -zxf /var/www/html/HFPhase1SX5Ops/results/plots.tar.gz"
wait
rm plots.tar.gz
cd -




