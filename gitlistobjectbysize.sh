#!/bin/bash -e
# work over each commit and append all files in tree to $tempFile
tempFile=$(mktemp)
IFS=$'\n'
for commitSHA1 in $(git rev-list --all); do
	git ls-tree -r --long "$commitSHA1" >>"$tempFile"
done

# sort files by SHA1, de-dupe list and finally re-sort by filesize
sort --key 3 "$tempFile" | \
	uniq | \
	sort --key 4 --numeric-sort --reverse | \
	awk '!visited[$5]++' | \
	awk '{ split( "B KiB MiB GiB" , v ); s=1; while( $4>1024 ){ $4/=1024; s++ } print int($4) " " v[s] "\t\t" $5 }' | \
	tac

# remove temp file
rm "$tempFile"