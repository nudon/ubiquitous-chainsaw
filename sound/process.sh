#!/bin/sh
space=" "
empty=""
hyphen="-"

convertAll() {
    ext=$1
    match=./*."$ext"
    for a in $match; do
	#for some reason file name found is the reg-exp search term
	if  [ "$a" != "$match" ]; then
	    z="${a//[[:punct:]]/$empty}"
	    b="${z[@]/%$ext/.wav}"
	    c="${b//$space/$empty}"
	    d="${c//[[:digit:]]/}"
	    echo "converting $a to a wav and deleting original file" 
	    ffmpeg -i "$a" -f wav "$d"
	    rm "$a"
	fi
    done

}

convertAll "mp3"
convertAll "flac"
convertAll "m4a"
