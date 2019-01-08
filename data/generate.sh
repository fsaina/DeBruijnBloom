#!/bin/bash

ALPHABET='ATGC'

if [ $# -eq 0 ]
  then
    echo "No arguments supplied, please provide the length for the generating sequence. E.g. ./generate 10"
    exit 1
fi

if [ "$1" -lt "1" ]
  then 
    echo "ERROR: Length of the sequence to generate has to be greater than 0"
    exit 1

fi

if [[ "$OSTYPE" == "linux-gnu" ]]; then
        # Linux
        cat /dev/urandom | tr -dc $ALPHABET | fold -w $1 | head -n 1 | sed -e 's/^/>r1|/' | tr '|' '\n'

elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
        cat /dev/urandom | env LC_CTYPE=C tr -dc $ALPHABET | fold -w $1 | head -n 1 | sed -e 's/^/>r1|/' | tr '|' '\n'
fi
