#!/bin/bash
# Script written by Anne-Morwenn Pastier, ISTerre, Grenoble

# Compiles REEF code
./makefile

# Parametrization
sl='waelbroeck2002'		# Sea-level curve = filename without the extension
u=-0.30				# Vertical rate * 10
r=07				# Reef growth rate
v=50				# Eroded volume
p=03				# Initial slope
zr=50				# Maximum depth for reef growth

# Absolute path to save model outputs
path=./Storage/umoy

# Simulations loop 

#while read sl
#do

#for u0 in {-27..-30}											# Vertical rate loop # ! Input vertical rate has to be integer between -99 and 99
#do
#	u=$(echo "scale=2;$u0/10" | bc)									# Computes decimal vertical rate for model input

#	for r in {05..09}										# Maximum reef growth rate loop	do
#		while read v										# Erosion rate loop (reading from file V0s)
#		do
#		      for p in {02..03}									# Initial slope loop (reading from file pentes)
#			do
				if [ ! -d $path$u/Rgr$r/V0$v/Pente$p ]					# Checks if the simulation is already done
				then

					# Writing input parameters for simulation in param0.dat
					echo $sl.dat  '! RSL variations file' > param0.dat
					echo $p '! Slope (%)' >> param0.dat
					echo $r '! Maximum reef growth rate (mm/yr)' >> param0.dat
					echo $zr '! Maximum reef growth depth (m)' >> param0.dat
					echo '1' '! Optimal depth for vertical gradient (m)' >> param0.dat
					echo '2' '! Optimal depth for horizontal gradient (m)' >> param0.dat 
					echo '3' '! Maximum depth for wave erosion (m)' >> param0.dat
					echo $u '! Uplift rate (mm/yr)' >> param0.dat
					echo $v '! Eroded volume (mm3/yr)' >> param0.dat
					echo '1000' '! Temporal increment (yr)' >> param0.dat
					echo '.true.' '! Reef construction or not ?' >> param0.dat
					echo '.true.' '! Wave erosion or not ?' >> param0.dat
				
				
					# Simulation
					./terrasse
	
					# Final figure with GMT
					./REEF-morpho-series.sh
					# Saving final figure
					mv morpho.ps umoy$u-rgr$r-V$v-P$p-$sl.ps

					# Compressing results
					tar zcf umoy$u-rgr$r-V$v-P$p-$sl.tar.gz out* sedim* rsl* topo* box* param* umoy*.ps  
					# Creating directory for saving
					if [ ! -d $path$u ] 
					then
						mkdir $path$u
					fi
					if [ ! -d $path$u/Rgr$r ] 
					then
						mkdir $path$u/Rgr$r
					fi
					if [ ! -d $path$u/Rgr$r/V0$v ] 
					then
						mkdir $path$u/Rgr$r/V0$v
					fi
					if [ ! -d $path$u/Rgr$r/V0$v/Pente$p ] 
					then
						mkdir $path$u/Rgr$r/V0$v/Pente$p
					fi
					
					# Storing outputs
					mv umoy$u-rgr$r-V$v-P$p-$sl.tar.gz umoy$u-rgr$r-V$v-P$p*.ps $path$u/Rgr$r/V0$v/Pente$p #vrelRSL* modern* vol* crois* superf* flat*

					# Cleaning from preceding simulation
					rm out* sedim* topo* box* rsl0* umoy*.ps modern* volreef* crois* superf* flat* vrelRSL* #abras* const*

				else
					echo "u$u r$r v$v p$p already done"
			fi
#			done 
#		done
#	done 
#done
#done < sls
