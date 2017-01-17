#!/bin/sh

if [ $# -eq 0 ]
then
	scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/admin/RunList.txt Files/
	wait
	scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results/PlotList.txt Files/
	wait
	
	cd Plots
	tar -cvf plots.tar *
	wait
	gzip plots.tar
	wait
	scp -P 53222 plots.tar.gz hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results
	wait
	ssh -p 53222 hfSX5@feynman.physics.uiowa.edu "cd /var/www/html/HFPhase1SX5Ops/results/ ; tar -zxvf /var/www/html/HFPhase1SX5Ops/results/plots.tar.gz"
	wait
	rm plots.tar.gz
	cd -
else
	if [ $1 -eq 1 ]
	then
		scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/admin/RunList.txt Files/
		wait
		scp -P 53222 hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFPhase1SX5Ops/results/PlotList.txt Files/
		wait
	elif [ $1 -eq 2 ]
	then
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
	fi
fi