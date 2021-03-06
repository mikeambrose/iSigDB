#!bin/bash
# checks that the directory doesn't get too large
max_size=3000000;
dir_size=$(du -sb './files' | cut -f1)
if [[ $dir_size -ge $max_size ]]; then
    echo tooooo big
fi

# checks that there is a certain amount of free space
min_space=100000000;
free_space=$(df './files' | awk '{print $4}' | tail -1)
if [[ $free_space -le $min_space ]]; then
    echo not enough room
fi
