#!/bin/sh

# Script plotting final figure for REEF model
# Written for GMT 5.1

# Variabes initialisation
PROJ=-JX18/10	
read para < param
#mv volvitesse* volvitesse$para
#mv croissvitesse* croissvitesse$para

xmin=`awk '{if (NR==1) print $1}' box$para`
xmax=`awk '{if (NR==1) print $2}' box$para`
ymin=`awk '{if (NR==1) print $3}' box$para`
ymax=`awk '{if (NR==1) print $4}' box$para`
BOX=-R$xmin/$xmax/$ymin/$ymax

# Increment for axes scales
dx=`awk '{if(NR==2) print $1}' box$para`
dy=`awk '{if(NR==2) print $2}' box$para`
g=g$dy
dy=$dy$g

# Temporal increment
tmin=0
tmax=22
tstep=1

# Input files names
ext0="out"
ext2="sedim"
# Output file name
out=morpho.ps

# Recovering simulation parameters
u=`echo $para | cut -b1-5`
r=`echo $para | cut -b6-7`
v=`echo $para | cut -b8-9`
p=`echo $para | cut -b10-11`


# Finale figure plot

gmt5 gmtset MAP_GRID_PEN_PRIMARY grey

gmt5 psbasemap $PROJ $BOX -Bpxa$dx+l"Distance (m)" -Bpya$g+l"Elevation (m)" -BWS+glightgrey -V -K -P > $out

for i in `ls $ext0* | sort -r`
do
	MIS=`echo $i | cut -b4-5`
	mis=`echo $MIS | awk '{sub(/^0*/,"");}1'`
	AGE=`echo $i | cut -b6-9`
	age=`echo $AGE | awk '{sub(/^0*/,"");}1'`
	
	# Computes homemade color scale: blue for glacial periods, red for interglacial periods, the oldest the palest, from MIS 11 to present (older=palest blue)
#	if [ $(( $mis % 2 )) -eq 0 ]
#	then
#		hue=240
#		sat=`echo "scale=4;1-($mis-2)/22" | bc`	
#	else
#		hue=355
#		sat=`echo "scale=4;1-($mis-1)/21" | bc`	
#	fi

	# New colorscale
      	if [ $AGE -eq "0000" ]
	      then
		col=255/128/111
	else

		if [ $age -le 1500 ] && [ $age -ge 425 ]
		then
			col=220/220/220	
		fi
		if [ $age -ge 359 -a $age -le 425 ]		# MIS 11
		then
			col=217/131/251	
		fi
		if [ $age -ge 338 ] && [ $age -le 358 ]		# MIS10
		then
			col=161/233/237	
		fi
		if [ $age -ge 299 ] && [ $age -le 337 ]		# MIS 9
		then
			col=239/173/222
		fi
		if [ $age -ge 280 ] && [ $age -le 298 ]		# MIS 8e
		then
			col=98/125/148
		fi
		if [ $age -ge 247 ] && [ $age -le 279 ]		# MIS 8c
		then
			col=132/168/197
		fi
		if [ $age -ge 227 ] && [ $age -le 246 ]		# MIS 7e
		then
			col=252/62/67
		fi
		if [ $age -ge 202 ] && [ $age -le 226 ]		# MIS 7c
		then
			col=252/134/137
		fi
		if [ $age -ge 189 ] && [ $age -le 201 ]		# MIS 7a
		then
			col=239/173/175
		fi
		if [ $age -ge 134 ] && [ $age -le 188 ]		# MIS 6
		then
			col=130/181/161
		fi
		if [ $age -ge 110 ] && [ $age -le 133 ]		# MIS 5e
		then
			col=238/107/46
		fi
		if [ $age -ge 89 ] && [ $age -le 109 ]		# MIS 5c
		then
			col=242/165/129
		fi
		if [ $age -ge 74 ] && [ $age -le 88 ]		# MIS 5a
		then
			col=255/210/189
		fi
		if [ $age -ge 64 ] && [ $age -le 73 ]		# MIS 4
		then
			col=90/127/199
		fi
		if [ $age -ge 31 ] && [ $age -le 63 ]		# MIS 3
		then
			col=207/189/255
		fi
		if [ $age -ge 13 ] && [ $age -le 30 ]		# MIS 2
		then
			col=131/231/255
		fi
		if [ $age -ge 0 ] && [ $age -le 12 ]		# MIS 1
		then
			col=255/128/111
		fi
	fi


	gmt5 psxy $i -G$col -L $PROJ $BOX -O -K >> $out
	gmt5 psxy $ext2$MIS$AGE* -G200-0-0.5 -t95 $PROJ $BOX -O -K >> $out
done
	
# Plots final topograhy (hides eroded aerial areas)
gmt5 psxy topo$para -G255 -L $PROJ $BOX -W -O -K >> $out 
gmt5 psxy rsl0$para $PROJ $BOX -Wblue -O -K >> $out
gmt5 pstext parame.dat -R0/50/-2/110 $PROJ -F+f8p+jBL -O -K >> $out

# colorscale
gmt5 psscale -D9/15/18/0.7h -G0/450 -CMIS.cpt -O -K >> $out

ps2pdf $out
