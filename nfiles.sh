### find and count the number of files in each directory (including subdirectories) starting from the current directory and 
### store the results in a file named nfiles.out, sorted by the number of files in descending order.

find . -type d -print0 | while read -d '' -r dir; do
    files=("$dir"/*)
    printf "%5d files in directory %s\n" "${#files[@]}" "$dir"
done > nfilestmp
sort -nrk1,1 nfilestmp > nfiles.out
rm nfilestmp
