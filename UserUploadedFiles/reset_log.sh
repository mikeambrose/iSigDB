#!bin/bash
# scheduled once every 24 hours
d="./logs.txt"
rm $d
echo "IP	nickname	size
" > $d
