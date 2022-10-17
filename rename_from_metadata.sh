#!/bin/bash
#pipe in the metadata with e.g. < metadata_XXXX.csv 
#Metadata should be a csv file, with the Library name including a trailing 
#underscore to capture different lenghts
#first field, the metadata second - e.g. XXXX_A_ NL-BZ17-24A

#Main

#read the variables lib and meta
while read lib meta 
#get the files in i
	do for i in *.gz 
#compare the from to the filename in i and if it matches then
 		do if [ "${i%_all*}" = "$lib" ] 
 		then
 #do the rename
  			mv $i ${meta}_${i} ; 
 #end the if
 		fi 
 # end the 2nd do
	done
# end the 1st do
done

