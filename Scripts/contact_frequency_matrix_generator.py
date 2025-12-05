set VAR1 [open "int.dat" w]

puts $VAR1 "frame\ index\ resid"

set nf [molinfo top get numframes]

for {set i 0} {$i < $nf} {incr i} {

	set A1 [atomselect top "(alpha and pfrag 6 to 11) and same residue as within 4 of pfrag 0 to 5" frame $i]
		
	set B1 [$A1 get resid]
	
	set L1 [llength $B1]
	
  		for {set j1 0} {$j1 < $L1} {incr j1} {
   		 set C1 [lindex $B1 $j1]
    		puts $VAR1 "$i $j1 $C1"
		}

}

close $VAR1
